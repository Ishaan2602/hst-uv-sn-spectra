import os, sys, glob
# generic stis x1d for late 2023ixf epochs: mama (g140l/g230l) + ccd g750m.
# no defringe/ocrreject. mama is photon-counting, g750m is below the fringing regime.
# usage: reduce_stis_generic.py <base_dir> [<base_dir>...] [--targname=SN2023IXF]
#   --targname: filter by substring of TARGNAME (case-insensitive). if omitted, extracts all science.
os.environ['CRDS_SERVER_URL'] = 'https://hst-crds.stsci.edu'
os.environ['CRDS_PATH'] = '/mnt/c/Users/eluru/CALTECH_SURF_2026/WORK/crds_cache'
os.environ['oref'] = '/mnt/c/Users/eluru/CALTECH_SURF_2026/WORK/crds_cache/references/hst/stis/'

import crds
import stistools.x1d
from astropy.io import fits

# parse optional --targname=... arg
targ_filter = None
dirs = []
for a in sys.argv[1:]:
    if a.startswith('--targname='):
        targ_filter = a.split('=', 1)[1].upper()
    else:
        dirs.append(a)

sci = []
for base in dirs:
    for raw in sorted(glob.glob(f'{base}/mastDownload/HST/*/*_raw.fits')):
        d = os.path.dirname(raw)
        root = os.path.basename(raw).replace('_raw.fits', '')
        h = fits.getheader(raw, 0)
        targname = str(h.get('TARGNAME', '')).upper()
        imagetyp = str(h.get('IMAGETYP', '')).upper()
        # skip lamp/calibration exposures (IMAGETYP contains LAMP, FLAT, BIAS, DARK, etc.)
        if any(x in imagetyp for x in ('LAMP', 'FLAT', 'BIAS', 'DARK', 'SFLAT', 'WAVECAL')):
            continue
        # optional target name filter
        if targ_filter and targ_filter not in targname:
            continue
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
