import os, sys, glob
# generic stis 1d extractor for the late 2023ixf epochs (MAMA G140L/G230L + CCD G750M).
# no defringe / no ocrreject - MAMA is photon-counting (x1d off the flt); G750M at Ha is below the
# fringing regime so we x1d off the mast crj. defringe stays in reduce_epoch_ts.py for G750L.
# usage: python reduce_stis_generic.py <visit_base_dir> [<visit_base_dir> ...]
#   each base has mastDownload/HST/<obs>/<obs>_{raw,flt,crj,wav}.fits
os.environ['CRDS_SERVER_URL'] = 'https://hst-crds.stsci.edu'
os.environ['CRDS_PATH'] = '/mnt/c/Users/eluru/CALTECH_SURF_2026/WORK/crds_cache'
os.environ['oref'] = '/mnt/c/Users/eluru/CALTECH_SURF_2026/WORK/crds_cache/references/hst/stis/'

import crds
import stistools.x1d
from astropy.io import fits

sci = []
for base in sys.argv[1:]:
    for raw in sorted(glob.glob(f'{base}/mastDownload/HST/*/*_raw.fits')):
        d = os.path.dirname(raw)
        root = os.path.basename(raw).replace('_raw.fits', '')
        if 'IXF' not in str(fits.getheader(raw, 0).get('TARGNAME', '')):
            continue  # skip lamp flats / non-science
        crj = f'{d}/{root}_crj.fits'
        inp = crj if os.path.exists(crj) else f'{d}/{root}_flt.fits'   # crj (ccd) else flt (mama)
        sci.append((root, inp))

print(f'{len(sci)} science exposures to extract')
crds.assign_bestrefs([inp for _, inp in sci], sync_references=True)

for root, inp in sci:
    out = f'{os.path.dirname(inp)}/{root}_x1d.fits'
    if os.path.exists(out):
        os.remove(out)
    try:
        stistools.x1d.x1d(inp, output=out, verbose=False)
        print('X1D', root, os.path.basename(inp))
    except Exception as e:
        print('FAIL', root, repr(e))
print('ALLDONE')
