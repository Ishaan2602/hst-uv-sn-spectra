# discover -> pick a stis ccd epoch -> download -> reduce (wsl) -> coadd.
# heavy reduction (x1d/defringe) goes through reduce_epoch.py in the surf_uv wsl env.
# run as a script for a demo, or import + call run().
import os, subprocess, glob
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.coordinates import SkyCoord
import astropy.units as u
from astroquery.mast import Observations

WSL_BASE = '/mnt/c/Users/eluru/CALTECH_SURF_2026/WORK'
STIS_CCD_GRATINGS = ('G230LB', 'G430L', 'G750L')


def find_uv_sne(limit=None):
    # discovery front-end: query by proposer classification, not name (catches disguised names like TESS-SN)
    tabs = []
    for inst in ('STIS/CCD', 'STIS/NUV-MAMA', 'STIS/FUV-MAMA', 'COS/NUV', 'COS/FUV'):
        t = Observations.query_criteria(obs_collection='HST', instrument_name=inst,
                                        target_classification='*SUPERNOVA*')
        if len(t):
            tabs.append(t['target_name', 's_ra', 's_dec', 'instrument_name', 'filters'])
    from astropy.table import vstack
    allobs = vstack(tabs)
    # dedupe to unique SNe by coordinate (0.5 arcmin)
    uniq, seen = [], []
    for row in allobs:
        c = (float(row['s_ra']), float(row['s_dec']))
        if all((c[0]-s[0])**2 + (c[1]-s[1])**2 > (0.5/60)**2 for s in seen):
            seen.append(c); uniq.append(row)
    return uniq[:limit] if limit else uniq


def pick_stis_ccd_epoch(target, ra=None, dec=None):
    # find one visit with all 3 ccd gratings + the contemporaneous g750l ccdflat.
    # use coordinate query if ra/dec given (name resolver fails for recent SNe).
    if ra is not None and dec is not None:
        coord = SkyCoord(ra, dec, unit='deg')
        t = Observations.query_criteria(coordinates=coord, radius=5*u.arcsec,
                                        obs_collection='HST', instrument_name='STIS/CCD')
    else:
        t = Observations.query_criteria(objectname=target, obs_collection='HST',
                                        instrument_name='STIS/CCD', radius='0.05 deg')
    visits = {}
    for r in t:
        oid, grat, tname = str(r['obs_id']), str(r['filters']), str(r['target_name']).upper()
        v = oid[:6]
        visits.setdefault(v, {})
        if 'FLAT' in tname and grat == 'G750L':
            visits[v]['flat'] = oid
        elif grat in STIS_CCD_GRATINGS:
            visits[v].setdefault(grat.lower(), oid)
    for v, d in visits.items():
        if all(g.lower() in d for g in STIS_CCD_GRATINGS):
            return {k: d[k] for k in ('g230lb', 'g430l', 'g750l') if k in d} | ({'flat': d['flat']} if 'flat' in d else {})
    raise RuntimeError(f"no single STIS/CCD visit with all 3 gratings for {target}")


def fetch(ids, target):
    roots = list(ids.values())
    pl = Observations.get_product_list(Observations.query_criteria(obs_id=roots))
    keep = [f.endswith(('_flt.fits', '_raw.fits', '_wav.fits')) for f in pl['productFilename']]
    Observations.download_products(pl[keep], download_dir=f'../data/{target}', mrp_only=False)


def reduce(ids, target):
    base = f'{WSL_BASE}/data/{target}/mastDownload/HST'
    args = ' '.join(f'{k}={v}' for k, v in ids.items())
    cmd = (f"source ~/miniforge3/etc/profile.d/conda.sh && conda activate surf_uv && "
           f"python {WSL_BASE}/scripts/reduce_epoch.py {base} {args}")
    subprocess.run(['wsl.exe', '-d', 'Debian', '--', 'bash', '-lc', cmd], check=True)


def _clean(w, f):
    o = np.argsort(w); w, f = w[o], f[o]
    keep = np.concatenate([[True], np.diff(w) > 0])
    return w[keep], f[keep]


def coadd(ids, target, save=True):
    hst = f'../data/{target}/mastDownload/HST'
    newax = np.concatenate([np.arange(1650, 3050, 1.4), np.arange(3050, 5600, 2.7), np.arange(5600, 10260, 4.9)])
    stack, fig_ax = [], plt.subplots(figsize=(13, 5))
    fig, ax = fig_ax
    for g in ('g230lb', 'g430l', 'g750l'):
        if g not in ids:
            continue
        r = ids[g]
        t = fits.getdata(f'{hst}/{r}/{r}_x1d.fits', 1)[0]
        good = (t['DQ'] & 16 != 16) & (t['DQ'] & 512 != 512) & np.isfinite(t['FLUX']) & (t['FLUX'] != 0)
        w, f = _clean(np.asarray(t['WAVELENGTH'][good], float), np.asarray(t['FLUX'][good], float))
        o = np.interp(newax, w, f, left=np.nan, right=np.nan)
        stack.append(o); ax.plot(newax, o, lw=0.5, alpha=0.5, label=g.upper()) # ISSUE HERE
    comb = np.nanmedian(np.vstack(stack), axis=0)
    ax.plot(newax, comb, lw=0.8, color='k', label='combined')
    ax.set_xlabel('wavelength (A)'); ax.set_ylabel('flux (erg/s/cm2/A)')
    ax.set_title(f'{target}  STIS UV-optical coadd (pipeline)'); ax.legend(ncol=4, fontsize=8)
    plt.tight_layout()
    if save:
        os.makedirs('../output', exist_ok=True)
        np.savetxt(f'../output/{target}_coadd.csv', np.column_stack([newax, comb]),
                   delimiter=',', header='wavelength_A,flux', comments='')
        fig.savefig(f'../output/{target}_coadd.png', dpi=120)
        print(f'saved ../output/{target}_coadd.csv + .png')
    return newax, comb


def run(target, ids=None):
    ids = ids or pick_stis_ccd_epoch(target)
    print(f'{target}: {ids}')
    fetch(ids, target)
    reduce(ids, target)
    return coadd(ids, target)


if __name__ == '__main__':
    # first-pass demo on SN2023ixf (data may already be cached locally)
    run('SN2023ixf', ids={'g230lb': 'oezt01040', 'g430l': 'oezt010h0',
                          'g750l': 'oezt010e0', 'flat': 'oezt010d0'})
