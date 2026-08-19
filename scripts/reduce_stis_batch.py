import os, sys, glob, json, argparse, shutil
import numpy as np

# our-own stis extractor, one SN at a time. runs in wsl (surf_uv).
# usage: python reduce_stis_batch.py <base_dir> <sn_name> [--outroot output] [--z Z] [--expl-mjd MJD] [--targname STR]
# writes the nested tree: output/<sn>/STIS/<CCD|MAMA|ECHELLE>/<date>_day<phase>/<grating>/{...}
os.environ['CRDS_SERVER_URL'] = 'https://hst-crds.stsci.edu'
os.environ['CRDS_PATH'] = '/mnt/c/Users/eluru/CALTECH_SURF_2026/WORK/crds_cache'
os.environ['oref'] = '/mnt/c/Users/eluru/CALTECH_SURF_2026/WORK/crds_cache/references/hst/stis/'
# hstcal's calstis executables (cs6.e etc.) live in the env bin; ensure it's on PATH even when this is
# spawned via the raw env python (conda activate not run), else stistools can't find cs6.e and EVERY
# x1d fails with FileNotFoundError - whole gratings vanish silently (this was the real 2022MOX defect).
os.environ['PATH'] = os.path.dirname(sys.executable) + os.pathsep + os.environ.get('PATH', '')

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import paths
import crds
import stistools.x1d
import stistools.defringe
from astropy.io import fits
from astropy.time import Time
from stis_extract import find_trace, adaptive_bg, adaptive_extrsize, draw_extraction, extrsize_for, is_echelle
import coadd as co

p = argparse.ArgumentParser()
p.add_argument('base')
p.add_argument('sn')
p.add_argument('--outroot', default=paths.OUT)
p.add_argument('--z', type=float, default=0.0)
p.add_argument('--expl-mjd', type=float, default=None)   # explosion MJD for the day-phase label
p.add_argument('--targname', default=None)
p.add_argument('--skip-proposid', default=None, type=int,  # skip exposures from a specific program
               help='skip exposures whose PROPOSID equals this integer (e.g. 17205 for the early saturated 2023IXF visits)')
p.add_argument('--sync-only', action='store_true',
               help='gather exposures + CRDS assign_bestrefs, then exit (serial cache warm-up before the parallel run)')
a = p.parse_args()

CAL = ('LAMP', 'FLAT', 'BIAS', 'DARK', 'WAVE')


def detector_of(h):
    return 'CCD' if str(h.get('DETECTOR', '')).upper() == 'CCD' else 'MAMA'


def is_cal(h):
    it = str(h.get('IMAGETYP', '')).upper(); tn = str(h.get('TARGNAME', '')).upper()
    return any(x in it for x in CAL) or any(x in tn for x in ('FLAT', 'LAMP', 'BIAS', 'DARK', 'WAVE'))


def epoch_dir(mjd):
    date = Time(mjd, format='mjd').iso[:10] if mjd else 'unknown'
    if a.expl_mjd and mjd:
        return f'{date}_day{int(round(mjd - a.expl_mjd))}'
    return date


def find_fringe_flat(sci_root):
    visit = sci_root[:6]
    for raw in glob.glob(f'{a.base}/mastDownload/HST/{visit}*/*_raw.fits'):
        h = fits.getheader(raw, 0)
        if 'FLAT' in str(h.get('TARGNAME', '')).upper() and str(h.get('OPT_ELEM', '')).upper() == 'G750L':
            return raw
    return None


def defringe_g750l(sci_input, sci_root, flat_raw, work):
    d = os.path.dirname(sci_input)
    wav = f'{d}/{sci_root}_wav.fits'
    if not os.path.exists(wav):
        return None
    flat_root = os.path.basename(flat_raw).replace('_raw.fits', '')
    os.makedirs(work, exist_ok=True)
    for src in (sci_input, wav, flat_raw):
        shutil.copy(src, work)
    cwd = os.getcwd(); os.chdir(work)
    try:
        sci_local = os.path.basename(sci_input)
        for f in (f'{flat_root}_nsp.fits', f'{flat_root}_frr.fits', f'{sci_root}_drj.fits'):
            if os.path.exists(f): os.remove(f)
        stistools.defringe.normspflat(f'{flat_root}_raw.fits', f'{flat_root}_nsp.fits',
                                      do_cal=True, wavecal=f'{sci_root}_wav.fits')
        stistools.defringe.mkfringeflat(sci_local, f'{flat_root}_nsp.fits', f'{flat_root}_frr.fits',
                                        beg_shift=-2.0, end_shift=1.0, shift_step=0.1)
        stistools.defringe.defringe(sci_local, f'{flat_root}_frr.fits', overwrite=True)
        drj = os.path.abspath(f'{sci_root}_drj.fits')
    except Exception as e:
        print('defringe error', sci_root, repr(e)[:150]); drj = None
    finally:
        os.chdir(cwd)
    return drj if drj and os.path.exists(drj) else None


def x1d_run(inp, out, **kw):
    if os.path.exists(out):
        os.remove(out)
    try:
        stistools.x1d.x1d(inp, output=out, verbose=False, **kw)
    except Exception as e:
        print('x1d error', repr(e)[:120])
    return os.path.exists(out)


def write_1d(path, hdr, dat):
    # per spectral order: rest_wvl obs_wvl flux err net dq mask
    with open(path, 'w') as f:
        f.write(hdr + '\n')
        f.write('# rest_wvl obs_wvl flux error net dq mask\n')
        for r in dat:
            for w, fl, er, nt, dq in zip(r['WAVELENGTH'], r['FLUX'], r['ERROR'], r['NET'], r['DQ']):
                f.write(f'{w/(1+a.z):.4f} {w:.4f} {fl:.6e} {er:.6e} {nt:.6e} {int(dq)} 0\n')


# ---- gather science exposures (crj preferred over flt) ----
seen = {}
for f in sorted(glob.glob(f'{a.base}/mastDownload/HST/*/*_crj.fits') +
                glob.glob(f'{a.base}/mastDownload/HST/*/*_flt.fits')):
    root = os.path.basename(f).replace('_crj.fits', '').replace('_flt.fits', '')
    if root in seen and seen[root]['inp'].endswith('_crj.fits'):
        continue
    h = fits.getheader(f, 0)
    if is_cal(h):
        continue
    if a.targname and a.targname.upper() not in str(h.get('TARGNAME', '')).upper():
        continue
    if a.skip_proposid and int(h.get('PROPOSID', 0) or 0) == a.skip_proposid:
        continue
    seen[root] = {'root': root, 'inp': f, 'grat': str(h.get('OPT_ELEM', '?')).upper(),
                  'det': detector_of(h), 'mjd': float(h.get('TEXPSTRT', 0) or h.get('EXPSTART', 0) or 0),
                  'exptime': float(h.get('TEXPTIME', 0) or 0)}
exps = list(seen.values())
# drop failed guide-star-acq / aborted exposures (texptime < 5s): they are pure noise (median good
# exposure is ~850s). record them in the manifest so nothing is silently lost. an epoch whose
# exposures ALL fail then just produces no coadd (as the paper dropped the 2023ixf day-50 of4304).
FAIL_EXPTIME = 5.0
failed = [e for e in exps if e['exptime'] < FAIL_EXPTIME]
exps = [e for e in exps if e['exptime'] >= FAIL_EXPTIME]
print(f'{len(exps)} STIS science exposures for {a.sn} ({len(failed)} dropped: texptime<{FAIL_EXPTIME}s)')
if not exps:
    print('BATCHDONE'); sys.exit(0)

# crds bestrefs + reference sync. run the full catalog at --workers 1 (serial): a concurrent
# sync_references=True writes a ref into the shared oref while another worker's x1d reads it
# (FileNotFoundError, silently dropping gratings = the 2022MOX defect). serial = one writer = no race.
crds.assign_bestrefs([e['inp'] for e in exps], sync_references=True)
if a.sync_only:
    print('CRDS-SYNCED', len(exps), 'exposures for', a.sn); sys.exit(0)

# ---- defringe G750L up front (needs the drj for both the finder and x1d) ----
for e in exps:
    e['defr'] = None
    if e['grat'] == 'G750L':
        flat = find_fringe_flat(e['root'])
        if flat:
            drj = defringe_g750l(e['inp'], e['root'], flat, f"{a.base}/defringe_work/{e['root']}")
            if drj:
                e['inp'] = drj; e['defr'] = 'defringed'
            else:
                e['defr'] = 'no_defringe(chain_failed)'
        else:
            e['defr'] = 'no_defringe(no_flat)'

# ---- pass 1: finder on every exposure; pick a per-source extrsize per ccd grating ----
best = {}   # grating -> (detect, extrsize)
for e in exps:
    if is_echelle(e['grat']):
        e['res'] = None; e['es'] = extrsize_for(e['grat'], e['det']); e['es_mode'] = 'echelle'; continue
    sci2d = fits.getdata(e['inp'], 1)
    res = find_trace(sci2d, e['det'], grating=e['grat'])
    es, es_mode = adaptive_extrsize(sci2d, res['a2center'], e['grat'], e['det'], flags=res['flags'])
    e['res'] = res; e['es'] = es; e['es_mode'] = es_mode
    if e['det'] == 'CCD' and es_mode.startswith('adaptive'):
        if e['grat'] not in best or res['detect'] > best[e['grat']][0]:
            best[e['grat']] = (res['detect'], es)
# lock the per-source adaptive value (from the best epoch) across all epochs of that grating
for e in exps:
    if e['grat'] in best:
        e['es'] = best[e['grat']][1]

# assign epoch labels by clustering within each detector on a rolling ~1-day window, so a single
# visit that straddles ut midnight (two calendar dates a few hours apart) stays one epoch.
for det in set(e['det'] for e in exps):
    dexps = sorted((e for e in exps if e['det'] == det), key=lambda e: e['mjd'])
    start = None; prev = None
    for e in dexps:
        if start is None or e['mjd'] - prev > 0.6:   # 0.6d: merge same-visit midnight-straddle, keep daily epochs apart
            start = e['mjd']
        e['epoch'] = epoch_dir(start)
        prev = e['mjd']

# ---- pass 2: extract, write the tree ----
manifest = {'sn': a.sn, 'z': a.z, 'instrument': 'STIS', 'epochs': {},
            'failed_exposures': [{'root': e['root'], 'grating': e['grat'], 'detector': e['det'],
                                  'exptime': round(e['exptime'], 1)} for e in failed]}
scaling = []
groups = {}   # (det, epoch) -> {grating: [spec, ...]}
gflags = {}   # (det, epoch, grating) -> set of qc flags, for the plot titles
for e in exps:
    det = e['det'] if not is_echelle(e['grat']) else 'ECHELLE'
    ep = e['epoch']
    outdir = f"{a.outroot}/{a.sn}/STIS/{det}/{ep}/{e['grat']}"
    os.makedirs(outdir, exist_ok=True)
    es = e['es']
    flags = list(e['res']['flags']) if e['res'] else ['echelle']
    if e['defr']: flags.append(e['defr'])

    x1d_out = f"{outdir}/{e['root']}_x1d.fits"
    if is_echelle(e['grat']):
        ok = x1d_run(e['inp'], x1d_out)   # calstis order extraction, defaults
        a2_used = np.nan
    else:
        a2 = e['res']['a2center']
        b1, b2 = adaptive_bg(fits.getdata(e['inp'], 1), int(round(a2)), es)
        o1, o2 = b1 - int(round(a2)), b2 - int(round(a2))
        kw = dict(extrsize=es, bk1offst=o1, bk2offst=o2, bk1size=10, bk2size=10)
        # default-first: trust calstis' own trace search; only fall back to our finder center
        # (maxsrch=0) when calstis returns nothing ("cannot extract"), the faint-trace case the
        # finder was built for. the finder + flags still drive the aperture, bg and qc; this only
        # sets which a2center x1d anchors on.
        ok = x1d_run(e['inp'], x1d_out, **kw)                               # calstis auto-search
        if not ok:
            print('RETRY-FINDER', e['root'], f'a2center={a2}')
            ok = x1d_run(e['inp'], x1d_out, a2center=a2, maxsrch=0, **kw)    # our finder center
    if not ok:
        print('FAIL no x1d', e['root'], e['grat']); continue

    dat = fits.getdata(x1d_out, 1)
    if dat is None or len(dat) == 0:      # x1d wrote an empty table (no order extracted)
        print('FAIL empty x1d', e['root'], e['grat']); continue
    if not is_echelle(e['grat']):
        a2_used = float(dat['A2CENTER'][0]) if 'A2CENTER' in dat.columns.names and len(dat) else np.nan
        if e['res'] and not np.isnan(a2_used) and abs(a2_used - e['res']['a2center']) > 5:
            flags.append(f"center_mismatch({a2_used:.0f}vs{e['res']['a2center']:.0f})")

    hdr = (f"# {a.sn} {e['root']} {e['grat']} {det} z={a.z} a2center={a2_used:.1f} "
           f"extrsize={es} exptime={e['exptime']:.1f} mode={e['es_mode']} flags={','.join(flags) or 'none'}")
    write_1d(f"{outdir}/{a.sn}_{e['root']}_{e['grat']}_1d.txt", hdr, dat)

    if e['res'] is not None:
        png = f"{outdir}/{a.sn}_{e['root']}_{e['grat']}_2d.png"
        draw_extraction(fits.getdata(e['inp'], 1), e['res'], e['det'], es,
                        title=f"{a.sn} {e['root']} {e['grat']} {e['exptime']:.0f}s", save=png)

    w = np.concatenate([r['WAVELENGTH'] for r in dat])
    fl = np.concatenate([r['FLUX'] for r in dat])
    er = np.concatenate([r['ERROR'] for r in dat])
    dq = np.concatenate([r['DQ'] for r in dat])
    # store OBSERVED-frame wavelengths throughout; redshift is applied once, at plot time.
    groups.setdefault((det, ep), {}).setdefault(e['grat'], []).append((np.asarray(w, float), fl, er, dq))
    gflags.setdefault((det, ep, e['grat']), set()).update(
        f for f in flags if f in ('finder_fail', 'extended_host', 'low_detect')
        or f.startswith('center_mismatch') or f.startswith('no_defringe'))
    print('X1D', e['root'], e['grat'], 'es', es, 'flags', flags or 'none')

# ---- per-filter + per-epoch coadds ----
import matplotlib
matplotlib.use('Agg'); import matplotlib.pyplot as plt

def _ylim(f):
    # robust linear y-limits: clip to the 1-99 pct so noisy tails don't flatten the plot.
    g = f[np.isfinite(f)]
    if len(g) < 5:
        return None
    lo, hi = np.nanpercentile(g, [1, 99])
    if not np.isfinite(hi - lo) or hi <= lo:
        return None
    pad = 0.08 * (hi - lo)
    return (min(lo, 0.0) - pad, hi + pad)

def _save1d(base, w, f, e, title, z):
    # write the 1d ascii (obs-frame w,f,e) + a linear rest-frame plot for one tier (native/resel).
    m = np.isfinite(f)
    if m.sum() < 2:
        return
    np.savetxt(f'{base}.txt', np.column_stack([w[m], f[m], e[m]]), header='obs_wvl flux error', comments='# ')
    fig, axp = plt.subplots(figsize=(11, 4))
    axp.plot(w[m] / (1 + z), f[m], lw=0.6, color='k')
    yl = _ylim(f[m])
    if yl:
        axp.set_ylim(*yl)
    axp.set_xlabel('rest wavelength (A)'); axp.set_ylabel('flux (uncorrected)'); axp.set_title(title)
    fig.tight_layout(); fig.savefig(f'{base}.png', dpi=110); plt.close(fig)

for (det, ep), grats in groups.items():
  try:
    epdir = f'{a.outroot}/{a.sn}/STIS/{det}/{ep}'
    nat_legs, res_legs = [], []   # (grating, w, f, e) per tier for the cross-grating merge
    for g, specs in grats.items():
        gf = gflags.get((det, ep, g), set())
        tag = f'  [{",".join(sorted(gf))}]' if gf else ''
        gp = f'{epdir}/{g}'; os.makedirs(gp, exist_ok=True)
        nw, nf, ne = co.coadd_native(specs, g)          # raw/native 1px (the priority product)
        if len(nw):
            _save1d(f'{gp}/{a.sn}_{ep}_{g}_native', nw, nf, ne, f'{a.sn} {ep} {g} native{tag}', a.z)
            nat_legs.append((g, nw, nf, ne))
        rw, rf, re = co.coadd_resel(specs, g)           # resel (2 native px, near-Nyquist)
        if len(rw):
            _save1d(f'{gp}/{a.sn}_{ep}_{g}_resel', rw, rf, re, f'{a.sn} {ep} {g} resel{tag}', a.z)
            res_legs.append((g, rw, rf, re))

    # cross-grating epoch coadds (native + resel): each on its own union grid, NO inter-grating scaling
    # (flux cal agrees), inverse-variance combine so a grating's drooped low-throughput edge self-down-weights.
    epf = set().union(*[gflags.get((det, ep, g), set()) for g in grats]) if grats else set()
    neg = False
    for tier, legs in (('native', nat_legs), ('resel', res_legs)):
        if not legs:
            continue
        gw, gmf, gme, recs = co.align_and_merge(legs, prefer=('G230LB',))
        if not len(gw):
            continue
        if tier == 'native':
            for r in recs:
                scaling.append({'epoch': ep, 'detector': det, 'grating': r['leg'], 'scale': r['scale']})
            fin = gmf[np.isfinite(gmf)]
            neg = len(fin) > 20 and float(np.nanmedian(fin)) < 0   # caught nothing real -> non-detection
        flg = sorted(epf | ({'NEG_CONTINUUM'} if neg else set()))
        tg = f'  [{",".join(flg)}]' if flg else ''
        _save1d(f'{epdir}/{a.sn}_{ep}_epochcoadd_{tier}', gw, gmf, gme, f'{a.sn} {ep} cross-grating {tier}{tg}', a.z)
    manifest['epochs'][f'{det}/{ep}'] = {'gratings': list(grats), 'detector': det, 'neg_continuum': bool(neg)}
  except Exception as e:
    print('epoch coadd error', det, ep, repr(e)[:150])

# ---- write the stis manifest fragment + scaling ----
sndir = f'{a.outroot}/{a.sn}'
os.makedirs(sndir, exist_ok=True)
with open(f'{sndir}/{a.sn}_stis_manifest.json', 'w') as f:
    json.dump(manifest, f, indent=1)
with open(f'{sndir}/{a.sn}_stis_scaling.json', 'w') as f:
    json.dump(scaling, f, indent=1)

# ---- cleanup heavy intermediates (keep raw downloads + the x1d/products) ----
for d in glob.glob(f'{a.base}/defringe_work/*'):
    shutil.rmtree(d, ignore_errors=True)

print('BATCHDONE')
