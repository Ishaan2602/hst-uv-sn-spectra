import os, sys, glob, json, argparse, re
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import cm, colors
import matplotlib as mpl
from astropy.time import Time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import paths

# the pre-reduced early 1d files are ordinal-named (1st..5th); their phases are known days
EARLY_PHASE = {'1st': 3.25, '2nd': 4.25, '3rd': 5.25, '4th': 8.25, '5th': 11.25}

# SN-level product assembly: read every per-epoch coadd in the tree and build the time-series
# waterfall, a merged manifest, and a flattened scaling.csv.
# --early-skip-dir: dir of pre-reduced 1D ascii files (e.g. data/earlytime_2023ixf/) whose
#   spectra should be included in the timeseries without being re-reduced.
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import coadd as co

DATE_RE = re.compile(r'(\d{4}-\d{2}-\d{2})(?:_day(-?\d+))?')


def load_txt(path):
    try:
        a = np.loadtxt(path, comments='#')
        if a.ndim < 2 or a.shape[0] < 2:
            return None, None, None
        e = a[:, 2] if a.shape[1] > 2 else None
        return a[:, 0], a[:, 1], e
    except Exception:
        return None, None, None


def phase_num(phase_str):
    try:
        return float(phase_str)
    except (TypeError, ValueError):
        return None


def gather_epochs(sndir, tier, early_dir=''):
    # collect per-(detector,epoch) <tier> cross-grating coadds (stis + cos) + the early pre-reduced
    # 1Ds (native tier only). each carries an instrument tag (STIS CCD/MAMA, COS FUV/NUV) + date/phase.
    specs = []
    for inst_dir in ('STIS', 'COS'):
        for p in glob.glob(f'{sndir}/{inst_dir}/*/*/*_epochcoadd_{tier}.txt'):
            ep = os.path.basename(os.path.dirname(p))
            det = os.path.basename(os.path.dirname(os.path.dirname(p)))
            m = DATE_RE.search(ep)
            w, f, e = load_txt(p)
            if w is None:
                continue
            ph = m.group(2) if m and m.group(2) else ''
            specs.append({'date': m.group(1) if m else '', 'phase': ph, 'phase_num': phase_num(ph),
                          'inst': f'{inst_dir} {det}', 'w': w, 'f': f, 'e': e, 'early': False})
    # pre-reduced early-time spectra (native-resolution 1Ds, e.g. days 3-11 for SN2023IXF): native tier
    # only. files have a text header line (not #-commented), format: WAVELENGTH FLUX ERROR ...
    if tier == 'native' and early_dir and os.path.isdir(early_dir):
        for p in sorted(glob.glob(f'{early_dir}/*.txt')):
            try:
                rows = []
                with open(p) as fh:
                    for line in fh:
                        parts = line.split()
                        if len(parts) >= 2:
                            try:
                                rows.append((float(parts[0]), float(parts[1])))
                            except ValueError:
                                pass
                if len(rows) < 5:
                    continue
                arr = np.array(rows)
                w, f = arr[:, 0], arr[:, 1]
                fname = os.path.basename(p)
                ordn = next((k for k in EARLY_PHASE if k in fname), None)
                ph = EARLY_PHASE.get(ordn)
                specs.append({'date': '', 'phase': f'{ph:.1f}' if ph else '', 'phase_num': ph,
                              'inst': 'STIS CCD', 'w': w, 'f': f, 'e': None, 'early': True})
            except Exception:
                pass
    return specs


def cluster_epochs(specs, sndir, tier, expl_mjd=None):
    # group the per-detector epoch spectra into global epochs on a rolling ~1-day window, across
    # detectors AND instruments: a date carrying both a ccd g750m stub and a mama uv epoch, or a
    # contemporaneous stis + cos visit, collapses to ONE epoch line (like the gold standard).
    # cross-instrument merge each cluster onto its OWN union grid (no coarsening, no scaling - the
    # reference trims + median-combines) and save epochs/<label>_<tier>.txt. also stash a COMMON_AXIS
    # realign per epoch for the display-only time-series waterfall. returns epochs earliest-first.
    for s in specs:
        if s['date']:
            try:
                s['mjd'] = float(Time(s['date']).mjd)
            except Exception:
                s['mjd'] = None
        elif s['phase_num'] is not None and expl_mjd:
            s['mjd'] = expl_mjd + s['phase_num']
            s['date'] = Time(s['mjd'], format='mjd').iso[:10]   # fill the date for the label
        else:
            s['mjd'] = None
    dated = sorted([s for s in specs if s['mjd'] is not None], key=lambda s: s['mjd'])
    clusters, cur, start, prev = [], [], None, None
    for s in dated:
        # 0.6-day gap (not 1.0): merges a single visit's exposures/gratings that straddle ut
        # midnight (hours apart) but keeps genuine daily-cadence epochs separate (e.g. the 2023ixf
        # early days 3/4/5, which are one day apart and are distinct epochs in the gold standard).
        if start is None or s['mjd'] - prev > 0.6:
            if cur:
                clusters.append(cur)
            cur, start = [s], s['mjd']
        else:
            cur.append(s)
        prev = s['mjd']
    if cur:
        clusters.append(cur)
    clusters += [[s] for s in specs if s['mjd'] is None]   # undated (no expl date) each their own

    epdir = f'{sndir}/epochs'
    os.makedirs(epdir, exist_ok=True)
    out = []
    for cl in clusters:
        # cross-instrument merge on the union grid, no scaling (flux cal agrees; the stis ccd leg is
        # only tagged as anchor). ivar where every leg has an error, else nanmedian (early 1Ds have none).
        gw, merged, merr, _ = co.align_and_merge([(s['inst'], s['w'], s['f'], s.get('e')) for s in cl],
                                                  prefer=('STIS CCD', 'CCD'))
        rep = min(cl, key=lambda s: s['phase_num'] if s['phase_num'] is not None else 1e9)
        insts = sorted({s['inst'] for s in cl})
        if not len(gw):
            continue
        mm = np.isfinite(merged)
        if mm.any() and rep['date']:
            tag = f"_day{rep['phase']}" if rep['phase'] else ''
            np.savetxt(f"{epdir}/{rep['date']}{tag}_{tier}.txt",
                       np.column_stack([gw[mm], merged[mm], merr[mm]]),
                       header=f'obs_wvl flux error  |  merged epoch ({tier}): {",".join(insts)}', comments='# ')
        # the waterfall (time series) is the one place the coarse COMMON_AXIS still lives: realign there.
        fcommon = co.resample_fcr(gw[mm], merged[mm], co.COMMON_AXIS)[0] if mm.sum() > 2 else np.full_like(co.COMMON_AXIS, np.nan)
        out.append({'date': rep['date'], 'phase': rep['phase'], 'phase_num': rep['phase_num'],
                    'insts': insts, 'merged_common': fcommon})
    out.sort(key=lambda s: (0, s['phase_num']) if s['phase_num'] is not None else (1, s['date']))
    return out


def waterfall(sn, epochs, out, z=0.0):
    # one line per global epoch (merged across detectors/instruments), earliest on top. log10(f/med)
    # like the gold standard; nan coverage gaps break the line (no connectors). the vertical offset
    # is DYNAMIC: each epoch gets room proportional to its own robust spread, so noisy late epochs
    # stop bleeding into the next one (the fixed 1.6 gap was too small for the noisy uv legs).
    # products are stored observed-frame; convert the axis to rest here (the single z-application).
    if not epochs:
        return
    n = len(epochs)
    have_phase = [e['phase_num'] for e in epochs if e['phase_num'] is not None]
    use_phase = len(have_phase) == n and n > 1
    if use_phase:
        cnorm = colors.LogNorm(vmin=max(1.0, min(have_phase)), vmax=max(have_phase) * 1.1)

    axc = co.COMMON_AXIS / (1 + z)
    ys, extents = [], []
    for e in epochs:
        f = np.asarray(e['merged_common'], float)
        with np.errstate(invalid='ignore', divide='ignore'):
            norm = np.nanmedian(f[np.isfinite(f) & (f > 0)]) or 1.0
            y = np.log10(f / norm)
        ys.append(y)
        fin = np.isfinite(y)
        ext = (np.nanpercentile(y[fin], 98) - np.nanpercentile(y[fin], 2)) if fin.any() else 1.0
        extents.append(float(np.clip(ext, 0.6, 2.5)))
    pad = 0.7
    base = [0.0] * n
    for i in range(1, n):
        base[i] = base[i - 1] - (extents[i - 1] / 2 + pad + extents[i] / 2)

    fig, ax = plt.subplots(figsize=(12, max(4, 0.5 * (sum(extents) + pad * n) + 1.5)))
    ally = []
    for i, e in enumerate(epochs):
        y = ys[i] + base[i]
        c = cm.turbo(cnorm(e['phase_num'])) if use_phase else cm.turbo(i / max(1, n - 1))
        ax.plot(axc, y, color=c, lw=0.5)          # nan in y breaks the line at coverage gaps
        fin = np.isfinite(y)
        if not fin.any():
            continue
        lab = e['date'] or (f"day {e['phase']}" if e['phase'] else '?')
        if e['phase'] and e['date']:
            lab += f" (d{e['phase']})"
        lab += f" [{'+'.join(e['insts'])}]"
        idx = np.where(fin)[0][-1]
        ax.text(axc[idx] * 1.01, y[idx], f'  {lab}', fontsize=7, va='center', color=c)
        ally.append(y[fin])

    if ally:
        ally = np.concatenate(ally)
        ax.set_ylim(np.nanmin(ally) - 0.6, np.nanmax(ally) + 0.6)
    ax.set_xlabel('rest wavelength (A)')
    ax.set_ylabel('log(normalized flux) - offset')
    ax.set_title(f'{sn} UV time series ({n} epochs)')
    if use_phase:
        sm = mpl.cm.ScalarMappable(cmap='turbo', norm=cnorm); sm.set_array([])
        cb = fig.colorbar(sm, ax=ax, pad=0.01, fraction=0.02)
        cb.set_label('phase (days)', fontsize=9)
    fig.tight_layout(); fig.savefig(out, dpi=130, bbox_inches='tight'); plt.close(fig)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('sn')
    ap.add_argument('--outroot', default=paths.OUT)
    ap.add_argument('--early-skip-dir', default='', dest='early_skip_dir')
    ap.add_argument('--expl-mjd', type=float, default=None, dest='expl_mjd')
    ap.add_argument('--z', type=float, default=0.0)
    a = ap.parse_args()
    sndir = f'{a.outroot}/{a.sn}'
    if not os.path.isdir(sndir):
        print(f'no products dir for {a.sn}'); return
    specs_n = gather_epochs(sndir, 'native', a.early_skip_dir)
    epochs = cluster_epochs(specs_n, sndir, 'native', a.expl_mjd)   # native drives the time-series waterfall
    cluster_epochs(gather_epochs(sndir, 'resel', a.early_skip_dir), sndir, 'resel', a.expl_mjd)
    waterfall(a.sn, epochs, f'{sndir}/{a.sn}_timeseries.png', a.z)
    # master_coadd removed: coadding epochs of an evolving SN is not physically meaningful.

    # merge the stis + cos manifests
    manifest = {'sn': a.sn, 'instruments': {}, 'epochs': {}}
    for inst in ('stis', 'cos'):
        mp = f'{sndir}/{a.sn}_{inst}_manifest.json'
        if os.path.exists(mp):
            m = json.load(open(mp))
            manifest['instruments'][inst.upper()] = {'z': m.get('z')}
            manifest['epochs'].update(m.get('epochs', {}))
    manifest['n_epochs'] = len(epochs)
    json.dump(manifest, open(f'{sndir}/{a.sn}_manifest.json', 'w'), indent=1)

    print(f'built products for {a.sn}: {len(epochs)} epochs')


if __name__ == '__main__':
    main()
