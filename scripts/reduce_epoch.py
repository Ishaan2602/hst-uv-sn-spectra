import os, sys, shutil
# generalized stis epoch reducer, runs in wsl (surf_uv). called by pipeline.py.
# usage: python reduce_epoch.py <base_dir> g230lb=<id> g430l=<id> g750l=<id> flat=<id>
# any grating may be omitted. single CRSPLIT=1 exposures extract straight off the flt;
# g750l is defringed first if a contemporaneous flat is given.
os.environ['CRDS_SERVER_URL'] = 'https://hst-crds.stsci.edu'
os.environ['CRDS_PATH'] = '/mnt/c/Users/eluru/CALTECH_SURF_2026/WORK/crds_cache'
os.environ['oref'] = '/mnt/c/Users/eluru/CALTECH_SURF_2026/WORK/crds_cache/references/hst/stis/'

import crds
import stistools.x1d, stistools.defringe

base = sys.argv[1]
ids = dict(a.split('=') for a in sys.argv[2:])   # e.g. {'g230lb':'oezt01040', 'flat':'oezt010d0'}

# sync refs for everything we'll touch
files = [f"{base}/{r}/{r}_flt.fits" for k, r in ids.items() if k != 'flat']
if 'flat' in ids:
    files.append(f"{base}/{ids['flat']}/{ids['flat']}_raw.fits")
crds.assign_bestrefs([f for f in files if os.path.exists(f)], sync_references=True)


def x1d_flt(root):
    out = f"{base}/{root}/{root}_x1d.fits"
    if os.path.exists(out):
        os.remove(out)
    stistools.x1d.x1d(f"{base}/{root}/{root}_flt.fits", output=out, verbose=False)
    print("X1D", root)


for g in ('g230lb', 'g430l'):
    if g in ids:
        x1d_flt(ids[g])

if 'g750l' in ids:
    r = ids['g750l']
    d = f"{base}/{r}"
    if 'flat' in ids:                              # defringe with the contemporaneous ccdflat
        flat = ids['flat']
        shutil.copy(f"{base}/{flat}/{flat}_raw.fits", f"{d}/{flat}_raw.fits")
        os.chdir(d)
        for f in (f"{flat}_nsp.fits", f"{flat}_frr.fits", f"{r}_drj.fits", f"{r}_x1d.fits"):
            if os.path.exists(f):
                os.remove(f)
        stistools.defringe.normspflat(f"{flat}_raw.fits", f"{flat}_nsp.fits", do_cal=True, wavecal=f"{r}_wav.fits")
        stistools.defringe.mkfringeflat(f"{r}_flt.fits", f"{flat}_nsp.fits", f"{flat}_frr.fits",
                                        beg_shift=-2.0, end_shift=1.0, shift_step=0.1)
        stistools.defringe.defringe(f"{r}_flt.fits", f"{flat}_frr.fits", overwrite=True)
        stistools.x1d.x1d(f"{r}_drj.fits", output=f"{r}_x1d.fits", verbose=False)
        print("X1D g750l defringed", r)
    else:
        x1d_flt(r)
print("EPOCHDONE")
