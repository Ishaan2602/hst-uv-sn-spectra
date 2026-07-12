import os, sys, shutil
# crj-aware stis epoch reducer for the 2023ixf time series, runs in wsl (surf_uv).
# usage: python reduce_epoch_ts.py <base_dir> g230lb=<id> g430l=<id> g750l=<id> flat=<id>
# the main epochs (of43) are CRSPLIT>1, so mast gives a cr-combined crj per grating - we
# extract off that crj (same as the paper: pipeline cr-rejection + our own x1d). single
# CRSPLIT=1 exposures (early oezt01) have no crj so we fall back to the flt. g750l is
# defringed with the contemporaneous ccdflat before extraction.
os.environ['CRDS_SERVER_URL'] = 'https://hst-crds.stsci.edu'
os.environ['CRDS_PATH'] = '/mnt/c/Users/eluru/CALTECH_SURF_2026/WORK/crds_cache'
os.environ['oref'] = '/mnt/c/Users/eluru/CALTECH_SURF_2026/WORK/crds_cache/references/hst/stis/'

import crds
import stistools.x1d, stistools.defringe

base = sys.argv[1]
ids = dict(a.split('=') for a in sys.argv[2:])   # e.g. {'g230lb':'of4301010', 'flat':'of4301040'}


def sci_input(root):
    # prefer the cr-combined crj (CRSPLIT>1 main epochs); fall back to the flt (single exposure)
    crj = f"{base}/{root}/{root}_crj.fits"
    return crj if os.path.exists(crj) else f"{base}/{root}/{root}_flt.fits"


# sync refs for everything we'll touch
files = [sci_input(r) for k, r in ids.items() if k != 'flat']
if 'flat' in ids:
    files.append(f"{base}/{ids['flat']}/{ids['flat']}_raw.fits")
crds.assign_bestrefs([f for f in files if os.path.exists(f)], sync_references=True)


def x1d_from(root, inp):
    out = f"{base}/{root}/{root}_x1d.fits"
    if os.path.exists(out):
        os.remove(out)
    stistools.x1d.x1d(inp, output=out, verbose=False)
    if not os.path.exists(out):
        # calstis "Cannot extract" (no exception, just no output): the auto trace-find bails when the
        # blue flux is very low - the day-66 g230lb hits this, same as the paper. retry with the fixed
        # E1 trace center (893.5, from the successful epochs) and no search, i.e. a manual extraction.
        print("RETRY", root, "auto extract failed, forcing a2center=893.5 maxsrch=0")
        stistools.x1d.x1d(inp, output=out, a2center=893.5, maxsrch=0, extrsize=7, verbose=False)
    print("X1D", root, os.path.basename(inp))


for g in ('g230lb', 'g430l'):
    if g in ids:
        x1d_from(ids[g], sci_input(ids[g]))

if 'g750l' in ids:
    r = ids['g750l']
    d = f"{base}/{r}"
    inp = os.path.basename(sci_input(r))           # crj if present, else flt
    if 'flat' in ids:                              # defringe with the contemporaneous ccdflat
        flat = ids['flat']
        shutil.copy(f"{base}/{flat}/{flat}_raw.fits", f"{d}/{flat}_raw.fits")
        os.chdir(d)
        for f in (f"{flat}_nsp.fits", f"{flat}_frr.fits", f"{r}_drj.fits", f"{r}_x1d.fits"):
            if os.path.exists(f):
                os.remove(f)
        stistools.defringe.normspflat(f"{flat}_raw.fits", f"{flat}_nsp.fits", do_cal=True, wavecal=f"{r}_wav.fits")
        stistools.defringe.mkfringeflat(inp, f"{flat}_nsp.fits", f"{flat}_frr.fits",
                                        beg_shift=-2.0, end_shift=1.0, shift_step=0.1)
        stistools.defringe.defringe(inp, f"{flat}_frr.fits", overwrite=True)
        stistools.x1d.x1d(f"{r}_drj.fits", output=f"{r}_x1d.fits", verbose=False)
        print("X1D g750l defringed", r)
    else:
        x1d_from(r, sci_input(r))
print("EPOCHDONE")
