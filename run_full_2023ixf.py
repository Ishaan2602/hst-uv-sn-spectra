import os, sys, shutil
# runs in wsl (surf_uv). full self-reduction of one stis epoch of 2023ixf: all 3 ccd gratings,
# g750l defringed with the contemporaneous ccdflat. syncs crds refs itself then does the work.
os.environ['CRDS_SERVER_URL'] = 'https://hst-crds.stsci.edu'
os.environ['CRDS_PATH'] = '/mnt/c/Users/eluru/CALTECH_SURF_2026/WORK/crds_cache'
os.environ['oref'] = '/mnt/c/Users/eluru/CALTECH_SURF_2026/WORK/crds_cache/references/hst/stis/'

import crds
import stistools.x1d, stistools.defringe

base = '/mnt/c/Users/eluru/CALTECH_SURF_2026/WORK/Data/2023ixf/mastDownload/HST'
sci = {'g230lb': 'oezt01040', 'g430l': 'oezt010h0', 'g750l': 'oezt010e0'}
flat = 'oezt010d0'

# make sure refs for this epoch are in the cache + headers point at them
files = [f"{base}/{r}/{r}_flt.fits" for r in sci.values()] + [f"{base}/{flat}/{flat}_raw.fits"]
crds.assign_bestrefs([f for f in files if os.path.exists(f)], sync_references=True)


def x1d_flt(root):
    out = f"{base}/{root}/{root}_x1d.fits"
    if os.path.exists(out):
        os.remove(out)
    stistools.x1d.x1d(f"{base}/{root}/{root}_flt.fits", output=out, verbose=False)
    print("X1D", root)


# blue gratings: extract straight off the flt (single exposure, nothing to CR-combine)
x1d_flt(sci['g230lb'])
x1d_flt(sci['g430l'])

# g750l: defringe the flt with the ccdflat, then extract the defringed frame
r = sci['g750l']
d = f"{base}/{r}"
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
print("FULLDONE")
