import os, sys, argparse, warnings
import numpy as np
import pandas as pd
from astropy.coordinates import SkyCoord
import astropy.units as u
from astroquery.mast import Observations
warnings.filterwarnings('ignore')

# download all STIS science data for one SN (or ALL) by coordinate.
# coord query (name resolver fails on recent SNe), verify each obs target within 5" of the catalog
# coord (kills cone-search false positives), pull raw+wav (+ any same-field ccdflat for g750l defringe).
# resumable (mast skips existing). runs on native windows python. cos is already downloaded+reduced.
HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
from paths import CATALOG as CLEAN
DATA = os.path.join(ROOT, 'data')


def download_sn(name, ra, dec):
    coord = SkyCoord(ra, dec, unit='deg')
    try:
        obs = Observations.query_criteria(coordinates=coord, radius=8 * u.arcsec,
                                          obs_collection='HST',
                                          instrument_name=['STIS/CCD', 'STIS/NUV-MAMA', 'STIS/FUV-MAMA'])
    except Exception as e:
        print(f'{name}: query err {repr(e)[:100]}'); return 0
    if obs is None or len(obs) == 0:
        print(f'{name}: no STIS obs'); return 0
    # keep spectroscopy on-target (within 5"); ccdflats share the field so they survive the cut
    sep = coord.separation(SkyCoord(obs['s_ra'], obs['s_dec'], unit='deg'))
    obs = obs[(sep < 5 * u.arcsec) & (obs['dataproduct_type'] == 'spectrum')]
    if len(obs) == 0:
        print(f'{name}: no on-target STIS spectra'); return 0
    pl = Observations.get_product_list(obs)
    # flt/crj are the x1d inputs; raw+wav are needed for the g750l ccdflat defringe chain
    keep = np.array([str(f).endswith(('_flt.fits', '_crj.fits', '_wav.fits', '_raw.fits'))
                     for f in pl['productFilename']])
    pl = pl[keep]
    if len(pl) == 0:
        print(f'{name}: no raw/wav products'); return 0
    Observations.download_products(pl, download_dir=os.path.join(DATA, name), mrp_only=False)
    print(f'{name}: downloaded {len(pl)} products ({len(obs)} obs)')
    return len(pl)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('sn')                       # a single SN name (as in the clean catalog) or ALL
    ap.add_argument('--limit', type=int, default=None)
    a = ap.parse_args()
    df = pd.read_csv(CLEAN)
    if a.sn != 'ALL':
        df = df[df['name'] == a.sn]
        if len(df) == 0:
            print(f'{a.sn} not in clean catalog'); return
    if a.limit:
        df = df.head(a.limit)
    for _, r in df.iterrows():
        # skip cos-only targets (no stis grating listed)
        if 'STIS' not in str(r.get('instr', '')):
            print(f"{r['name']}: cos-only, skip stis download"); continue
        download_sn(str(r['name']), float(r['ra']), float(r['dec']))
    print('DOWNLOADDONE')


if __name__ == '__main__':
    main()
