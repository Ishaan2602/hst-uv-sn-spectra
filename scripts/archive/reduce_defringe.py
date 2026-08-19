import os, sys, shutil
import stistools.defringe, stistools.x1d

# runs in wsl (surf_uv), oref exported. G750L defringing for SN2024iss.
# copies the needed files into a flat work dir (the tools write lots of siblings), then:
# normspflat (normalize the contemporaneous fringe flat) -> mkfringeflat (match it to the science)
# -> defringe (divide out) -> x1d on both the fringed (crj) and defringed (drj).
base, work = sys.argv[1], sys.argv[2]
# sci/flat rootnames default to the 2024iss pair, override for other targets
sci  = sys.argv[3] if len(sys.argv) > 3 else 'of8b02040'
flat = sys.argv[4] if len(sys.argv) > 4 else 'of8b02050'

os.makedirs(work, exist_ok=True)
shutil.copy(f"{base}/{sci}/{sci}_crj.fits", f"{work}/{sci}_crj.fits")
shutil.copy(f"{base}/{sci}/{sci}_wav.fits", f"{work}/{sci}_wav.fits")
shutil.copy(f"{base}/{flat}/{flat}_raw.fits", f"{work}/{flat}_raw.fits")
os.chdir(work)

# clear anything from a previous run (the pipeline wont overwrite)
for f in [f"{flat}_nsp.fits", f"{flat}_frr.fits", f"{flat}_crj.fits", f"{sci}_drj.fits",
          f"{sci}_x1d_fr.fits", f"{sci}_x1d_df.fits"]:
    if os.path.exists(f):
        os.remove(f)

stistools.defringe.normspflat(f"{flat}_raw.fits", f"{flat}_nsp.fits", do_cal=True, wavecal=f"{sci}_wav.fits")
stistools.defringe.mkfringeflat(f"{sci}_crj.fits", f"{flat}_nsp.fits", f"{flat}_frr.fits",
                                beg_shift=-2.0, end_shift=1.0, shift_step=0.1)
stistools.defringe.defringe(f"{sci}_crj.fits", f"{flat}_frr.fits", overwrite=True)

stistools.x1d.x1d(f"{sci}_crj.fits", output=f"{sci}_x1d_fr.fits", verbose=False)   # fringed
stistools.x1d.x1d(f"{sci}_drj.fits", output=f"{sci}_x1d_df.fits", verbose=False)   # defringed
print("DEFRINGEDONE")
