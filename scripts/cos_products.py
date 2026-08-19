import os, sys, glob, json, argparse
import numpy as np
from astropy.io import fits
from astropy.time import Time
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# build the cos data tree from the already-run default calcos products in data/cos_catalog/reduced/.
# per SN: group x1dsum by detector (FUV/NUV) + date + grating, apply the artifact masks, write
# 1d ascii + per-filter coadd + a 2D cross-dispersion viz + 1d plot, plus a manifest.
# usage: python cos_products.py <SN|ALL> [--reduced data/cos_catalog/reduced] [--outroot output] [--z ...] [--expl-mjd ...]

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, HERE)
import paths

NUVC_ART = (2150.0, 2350.0)   # known nuvc artifact stripe
# per-grating useful wavelength range (COS Instrument Handbook). data outside is detector-edge or
# second-order garbage: e.g. the G230L NUVA blue stripe (<1650, spurious ~1e-9 spikes) and the G140L
# noisy red edge (>1950, extends to ~2500). replaces the old blue-only COS_FUV_WMIN with a full range.
# VERIFY / extend against the handbook when a new grating/cenwave appears.
COS_RANGE = {'G130M': (1130.0, 1470.0), 'G160M': (1380.0, 1780.0), 'G140L': (1130.0, 1950.0),
             'G140M': (1130.0, 1470.0), 'G230L': (1650.0, 3200.0), 'G185M': (1650.0, 2150.0),
             'G225M': (2050.0, 2550.0), 'G285M': (2450.0, 3250.0)}
# geocoronal airglow windows (OBSERVED frame): H I Lya 1215.67 + the O I 1302/1305/1306 triplet
# (STScI COS airglow page; He II quasars Syphers 2012 - contamination is "primarily H I Lya + O I 1302").
# in G140L the Lya / N I 1200 blend is broad and low-res, so mask the whole region.
AIRGLOW = [(1213.5, 1218.0), (1300.5, 1307.5)]          # medium-res FUV: narrow Lya + O I 1302
AIRGLOW_G140L = [(1194.0, 1226.0), (1295.0, 1312.0)]    # G140L low-res: broad Lya/N I + broad O I


def airglow_windows(grating):
    return AIRGLOW_G140L if str(grating).upper() == 'G140L' else AIRGLOW


def mask_row(segment, wvl, grating=None):
    # per-pixel mask CODE (0 keep; 1 NUVC artifact; 2 geocoronal airglow; 4 out-of-grating-range).
    # the coadd drops any nonzero code; the 1d.txt keeps the code so every mask is recoverable.
    seg = str(segment).upper(); grat = str(grating).upper()
    wvl = np.asarray(wvl, float)
    m = np.zeros(len(wvl), int)
    rng = COS_RANGE.get(grat)
    if rng is not None:
        m[(wvl < rng[0]) | (wvl > rng[1])] = 4
    if seg == 'NUVC':
        if grat == 'G230L':
            m[:] = 1                       # NUVC stripe is pure noise in G230L (validated 2/2 NUV sources)
        else:
            m[(wvl >= NUVC_ART[0]) & (wvl <= NUVC_ART[1])] = 1
    if seg in ('FUVA', 'FUVB'):
        for lo, hi in airglow_windows(grat):   # Lya + O I airglow (grating-dependent width)
            m[(wvl >= lo) & (wvl <= hi)] = 2
    return m


def epoch_dir(mjd, expl):
    date = Time(mjd, format='mjd').iso[:10] if mjd else 'unknown'
    return f'{date}_day{int(round(mjd - expl))}' if (expl and mjd) else date


def counts_2d(reddir, root, outpng, grat):
    # cross-dispersion profile from the counts image(s), so the tree has a 2D-style viz for cos too.
    cf = sorted(glob.glob(f'{reddir}/*_counts_a.fits') + glob.glob(f'{reddir}/*_counts.fits'))
    if not cf:
        return False
    try:
        img = fits.getdata(cf[0])
    except (IndexError, OSError, TypeError):
        return False
    if img is None:
        return False
    prof = np.nansum(img, axis=1)
    fig, ax = plt.subplots(1, 2, figsize=(13, 4), gridspec_kw={'width_ratios': [3, 1]})
    vlo, vhi = np.nanpercentile(img, [40, 99.5])
    ax[0].imshow(img, origin='lower', aspect='auto', vmin=vlo, vmax=vhi, cmap='gray')
    ax[0].set_xlabel('dispersion (px)'); ax[0].set_ylabel('cross-disp (px)')
    ax[0].set_title(f'{grat} counts ({os.path.basename(cf[0])})', fontsize=9)
    ax[1].plot(prof, np.arange(len(prof)), 'k', lw=0.6)
    ax[1].set_xlabel('counts'); ax[1].set_title('cross-disp profile', fontsize=9)
    fig.tight_layout(); fig.savefig(outpng, dpi=100); plt.close(fig)
    return True


def prods(d):
    # x1dsum is the fp-pos-combined product; fall back to per-exposure x1d when calcos ran per-corrtag
    # and never wrote an x1dsum (e.g. SN2010AL). coadd_resel then fp-combines the per-exposure x1ds.
    xs = sorted(glob.glob(f'{d}/*_x1dsum.fits'))
    return xs if xs else sorted(glob.glob(f'{d}/*_x1d.fits'))


def process_sn(sn, reduced, outroot, z, expl):
    import coadd as co
    dirs = sorted(glob.glob(f'{reduced}/{sn}_*'))
    if not dirs:
        print(f'no reduced cos products for {sn}'); return 0
    # drop failed/aborted associations (exptime < 5s), the same gate as stis. record them.
    good, failed_cos = [], []
    for d in dirs:
        xs = prods(d)
        if not xs:
            continue
        try:
            et = float(fits.getheader(xs[0], 1).get('EXPTIME', 0)
                       or fits.getheader(xs[0], 0).get('EXPTIME', 0) or 0)
        except Exception:
            et = 0.0
        if et < 5.0:
            h0 = fits.getheader(xs[0], 0)
            failed_cos.append({'root': os.path.basename(d).split('_')[-2],
                               'grating': str(h0.get('OPT_ELEM', '?')).upper(),
                               'detector': str(h0.get('DETECTOR', '')).upper(), 'exptime': round(et, 1)})
        else:
            good.append(d)
    dirs = good
    manifest = {'sn': sn, 'z': z, 'instrument': 'COS', 'epochs': {}, 'failed_exposures': failed_cos}
    groups = {}   # (det, ep, grat) -> {'specs': [...], 'counts_dir': d, 'root': root}
    n = 0
    # cluster associations per detector on a rolling ~1-day window so a visit that straddles ut
    # midnight (two calendar dates a few hours apart) stays one epoch.
    byd = {}
    for d in dirs:
        xs = prods(d)
        if not xs:
            continue
        h0 = fits.getheader(xs[0], 0)
        det = str(h0.get('DETECTOR', 'FUV')).upper()
        mjd = float(h0.get('EXPSTART', 0) or fits.getheader(xs[0], 1).get('EXPSTART', 0) or 0)
        byd.setdefault(det, set()).add(mjd)
    ep_of = {}
    for det, mjds in byd.items():
        start = None; prev = None
        for mjd in sorted(mjds):
            if start is None or mjd - prev > 0.6:   # 0.6d: merge same-visit midnight-straddle, keep daily epochs apart
                start = mjd
            ep_of[(det, mjd)] = epoch_dir(start, expl)
            prev = mjd
    for d in dirs:
        for x1dsum in prods(d):
            root = os.path.basename(x1dsum).split('_')[0]
            kind = 'x1dsum' if 'x1dsum' in os.path.basename(x1dsum) else 'x1d'
            h0 = fits.getheader(x1dsum, 0)
            det = str(h0.get('DETECTOR', 'FUV')).upper()
            grat = str(h0.get('OPT_ELEM', '?')).upper()
            mjd = float(h0.get('EXPSTART', 0) or fits.getheader(x1dsum, 1).get('EXPSTART', 0) or 0)
            ep = ep_of.get((det, mjd), epoch_dir(mjd, expl))
            outdir = f'{outroot}/{sn}/COS/{det}/{ep}/{grat}'
            os.makedirs(outdir, exist_ok=True)
            try:
                dat = fits.getdata(x1dsum, 1)
            except (IndexError, OSError, TypeError):
                print(f'skip empty {kind} {os.path.basename(x1dsum)}'); continue
            if dat is None or len(dat) == 0:
                continue
            if not np.any(np.asarray(dat['FLUX'], float) != 0):   # uncalibrated x1d (all-zero flux): skip
                print(f'skip zero-flux {kind} {os.path.basename(x1dsum)}'); continue
            txt = f'{outdir}/{sn}_{ep}_{grat}_{root}_{kind}_1d.txt'
            w_all, f_all, e_all = [], [], []
            with open(txt, 'w') as f:
                f.write(f'# {sn} {root} {grat} {det} z={z}  cos default calcos {kind}\n')
                f.write('# segment rest_wvl obs_wvl flux error gross dq mask\n')
                for r in dat:
                    seg = str(r['SEGMENT'])
                    w = np.asarray(r['WAVELENGTH'], float); fl = np.asarray(r['FLUX'], float)
                    er = np.asarray(r['ERROR'], float)
                    gr = np.asarray(r['GROSS'], float) if 'GROSS' in dat.columns.names else np.zeros_like(w)
                    dq = np.asarray(r['DQ'], int) if 'DQ' in dat.columns.names else np.zeros(len(w), int)
                    mk = mask_row(seg, w, grat)
                    for i in range(len(w)):
                        f.write(f'{seg} {w[i]/(1+z):.4f} {w[i]:.4f} {fl[i]:.6e} {er[i]:.6e} '
                                f'{gr[i]:.4e} {int(dq[i])} {mk[i]}\n')
                    keep = (mk == 0)
                    # store OBSERVED-frame wvl; redshift applied once at plot time (1d.txt keeps both cols)
                    w_all.append(w[keep]); f_all.append(fl[keep]); e_all.append(er[keep])
            # guard against empty arrays (all-masked segment, e.g. E0102, zero-flux stripes)
            if not any(len(x) > 0 for x in w_all):
                continue
            w_all = [x for x in w_all if len(x) > 0]
            f_all = [x for x in f_all if len(x) > 0]
            e_all = [x for x in e_all if len(x) > 0]
            # per-stripe legs (NOT concatenated): each cos stripe is a distinct wvl segment. concatenating
            # them made the FCR bridge the inter-stripe gaps with bogus flux (the day214 Mg II 2x-low bug).
            g = groups.setdefault((det, ep, grat), {'specs': [], 'counts_dir': d, 'root': root, 'nexp': 0})
            for wk, fk, ek in zip(w_all, f_all, e_all):
                g['specs'].append((wk, fk, ek, np.zeros(len(wk), int)))
            g['nexp'] += 1
            n += 1

    # per (det, ep, grat): coadd the associations (fp-pos/epoch) at native + resel, write both + 2d viz
    def mask_cos(ax, fx, grat, det):
        # null out-of-range + airglow windows, then cut 100x-local-median junk (geocoronal Lya spikes)
        gm = np.zeros(len(ax), bool)
        _r = COS_RANGE.get(grat)
        if _r is not None:
            gm |= (ax < _r[0]) | (ax > _r[1])
        if det == 'FUV':
            for lo, hi in airglow_windows(grat):
                gm |= (ax >= lo) & (ax <= hi)
        fx = np.where(gm, np.nan, fx)
        if np.any(np.isfinite(fx)):
            from scipy.ndimage import generic_filter
            med = generic_filter(np.where(np.isfinite(fx), fx, 0.0), np.nanmedian, size=51)
            scale = np.maximum(np.abs(med), np.nanpercentile(np.abs(fx[np.isfinite(fx)]), 25) or 1e-30)
            fx = np.where(np.abs(fx) > 100 * scale, np.nan, fx)
        return fx

    def ylim_cos(f):
        gg = f[np.isfinite(f)]
        if len(gg) < 5:
            return None
        lo, hi = np.nanpercentile(gg, [1, 99])
        if not np.isfinite(hi - lo) or hi <= lo:
            return None
        pad = 0.08 * (hi - lo)
        return (min(lo, 0.0) - pad, hi + pad)

    def save_cos(base, w, f, e, title):
        m = np.isfinite(f)
        if m.sum() < 2:
            return
        np.savetxt(f'{base}.txt', np.column_stack([w[m], f[m], e[m]]), header='obs_wvl flux error', comments='# ')
        fig, axp = plt.subplots(figsize=(11, 4))
        axp.plot(w[m] / (1 + z), f[m], lw=0.5, color='k')
        yl = ylim_cos(f[m])
        if yl:
            axp.set_ylim(*yl)
        axp.set_xlabel('rest wavelength (A)'); axp.set_ylabel('flux (uncorrected)'); axp.set_title(title)
        fig.tight_layout(); fig.savefig(f'{base}.png', dpi=110); plt.close(fig)

    epoch_nat, epoch_res = {}, {}   # (det, ep) -> [(grat, w, f, e)] per tier for the cross-grating merge
    for (det, ep, grat), g in groups.items():
        outdir = f'{outroot}/{sn}/COS/{det}/{ep}/{grat}'; os.makedirs(outdir, exist_ok=True)
        for tier, (ax, fx, ex), store in (('native', co.coadd_native(g['specs'], grat), epoch_nat),
                                          ('resel', co.coadd_resel(g['specs'], grat, cos=True), epoch_res)):
            if len(ax) == 0:
                continue
            fx = mask_cos(ax, fx, grat, det)
            m = np.isfinite(fx)
            if m.sum() < 2:
                continue
            save_cos(f'{outdir}/{sn}_{ep}_{grat}_{tier}', ax[m], fx[m], ex[m],
                     f'{sn} {ep} COS {grat} {tier} ({g["nexp"]} exp)')
            store.setdefault((det, ep), []).append((grat, ax[m], fx[m], ex[m]))
        counts_2d(g['counts_dir'], g['root'], f'{outdir}/{sn}_{ep}_{grat}_2d.png', grat)
        manifest['epochs'].setdefault(f'{det}/{ep}', {'gratings': [], 'detector': det})
        if grat not in manifest['epochs'][f'{det}/{ep}']['gratings']:
            manifest['epochs'][f'{det}/{ep}']['gratings'].append(grat)
        print('COS', sn, ep, det, grat, f'({g["nexp"]} exp)')

    # per (det, ep): cross-grating epoch coadds (native + resel), union grid, no scaling, ivar combine.
    # one cos epoch is one product/line (fixes the two-line-per-epoch case, e.g. 2009IP G130M + G160M).
    for tier, store in (('native', epoch_nat), ('resel', epoch_res)):
        for (det, ep), gl in store.items():
            epdir = f'{outroot}/{sn}/COS/{det}/{ep}'; os.makedirs(epdir, exist_ok=True)
            gw, gmf, gme, recs = co.align_and_merge([(gg, w, f, e) for gg, w, f, e in gl])
            if not len(gw):
                continue
            fin = gmf[np.isfinite(gmf)]
            neg = len(fin) > 20 and float(np.nanmedian(fin)) < 0
            if tier == 'native':
                manifest['epochs'].setdefault(f'{det}/{ep}', {'gratings': [], 'detector': det})
                manifest['epochs'][f'{det}/{ep}']['neg_continuum'] = bool(neg)
            tg = '  [NEG_CONTINUUM]' if neg else ''
            save_cos(f'{epdir}/{sn}_{ep}_epochcoadd_{tier}', gw, gmf, gme,
                     f'{sn} {ep} COS {det} cross-grating {tier}{tg}')

    sndir = f'{outroot}/{sn}'
    os.makedirs(sndir, exist_ok=True)
    with open(f'{sndir}/{sn}_cos_manifest.json', 'w') as f:
        json.dump(manifest, f, indent=1)
    return n


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('sn')
    ap.add_argument('--reduced', default=os.path.join(ROOT, 'data', 'cos_catalog', 'reduced'))
    ap.add_argument('--outroot', default=paths.OUT)
    ap.add_argument('--z', type=float, default=0.0)
    ap.add_argument('--expl-mjd', type=float, default=None)
    a = ap.parse_args()
    if a.sn == 'ALL':
        sns = sorted({os.path.basename(d).split('_')[0] for d in glob.glob(f'{a.reduced}/*')})
        sns = [s for s in sns if s not in ('N132D-KNOT', 'E0102-HOTSPOT')]   # omit the snrs
        for s in sns:
            process_sn(s, a.reduced, a.outroot, a.z, a.expl_mjd)
    else:
        process_sn(a.sn, a.reduced, a.outroot, a.z, a.expl_mjd)
    print('COSDONE')


if __name__ == '__main__':
    main()
