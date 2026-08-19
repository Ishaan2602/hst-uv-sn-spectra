import glob, os
import numpy as np
from astropy.io import fits
from collections import defaultdict

# fast scan: count stis x1ds per grating from the folder names (no fits reads), then read ONE
# representative file per grating for the native per-px dispersion (it's an instrument constant).
# shows which modes get downsampled onto the coarse COMMON_AXIS (echelle + medium-res) vs ~1:1.
roots = ['output', 'output2']
files = []
for r in roots:
    files += glob.glob(os.path.join(r, '*', 'STIS', '*', '*', '*', '*_x1d.fits'))

by_grat = defaultdict(lambda: {'files': [], 'src': set()})
for f in files:
    p = f.split(os.sep)
    grat = p[-2].upper(); src = p[-6]
    by_grat[grat]['files'].append(f)
    by_grat[grat]['src'].add(src)


def common_step(wc):
    if wc < 1650: return 1.0
    if wc < 3050: return 1.4
    if wc < 5600: return 2.7
    return 4.9


print(f'{"grat":8s} {"nfile":>5s} {"nsrc":>4s} {"disp_A/px":>9s} {"nord":>4s} '
      f'{"wrange":>15s} {"cmn_step":>8s} {"downsamp":>8s}')
rows = []
for g, d in by_grat.items():
    disp = np.nan; nord = 0; wlo = whi = np.nan
    for f in d['files'][:8]:   # try a few in case the first is empty
        try:
            dat = fits.getdata(f, 1)
        except Exception:
            continue
        if dat is None or len(dat) == 0 or 'WAVELENGTH' not in dat.columns.names:
            continue
        wl = np.asarray(dat['WAVELENGTH'], float)
        if wl.ndim == 1:
            wl = wl[None, :]
        ds = [np.nanmedian(np.abs(np.diff(row))) for row in wl if row.size > 1]
        if ds:
            disp = float(np.nanmedian(ds)); nord = int(wl.shape[0])
            wlo = float(np.nanmin(wl)); whi = float(np.nanmax(wl))
            break
    cs = common_step((wlo + whi) / 2) if np.isfinite(wlo) else np.nan
    rows.append((g, len(d['files']), len(d['src']), disp, nord, wlo, whi, cs,
                 cs / disp if disp and np.isfinite(disp) else np.nan))

for g, nf, ns, disp, nord, wlo, whi, cs, ds in sorted(rows, key=lambda r: -(r[8] if np.isfinite(r[8]) else 0)):
    print(f'{g:8s} {nf:5d} {ns:4d} {disp:9.4f} {nord:4d} {wlo:6.0f}-{whi:5.0f} {cs:8.2f} {ds:7.1f}x')
