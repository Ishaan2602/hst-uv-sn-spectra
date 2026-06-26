import os, argparse
import stistools.ocrreject, stistools.x1d

# runs in wsl (surf_uv env) since ocrreject/x1d need the hstcal binaries.
# parametrized so we can play with extraction knobs (extrsize, bg offsets/sizes, manual center).
p = argparse.ArgumentParser()
p.add_argument("obsdir")                              # wsl path to the obsid folder
p.add_argument("root")                                # rootname, e.g. of8b02010
p.add_argument("--out", default="230_x1d_ours.fits")  # output x1d name (lands in obsdir)
for k in ("extrsize", "a2center", "maxsrch", "bk1offst", "bk2offst", "bk1size", "bk2size"):
    p.add_argument(f"--{k}", type=float, default=None)
a = p.parse_args()

flt = f"{a.obsdir}/{a.root}_flt.fits"
crj_ours = f"{a.obsdir}/230_crj_ours.fits"
x1d_out = f"{a.obsdir}/{a.out}"

# ocrreject only needs doing once, reuse the crj if it's already sitting there
if not os.path.exists(crj_ours):
    stistools.ocrreject.ocrreject(flt, crj_ours, verbose=False)

# pipeline wont overwrite, so clear this output first
if os.path.exists(x1d_out):
    os.remove(x1d_out)

# only pass the knobs we actually set, let stistools default the rest
kw = {k: v for k, v in vars(a).items()
      if k in ("extrsize", "a2center", "maxsrch", "bk1offst", "bk2offst", "bk1size", "bk2size") and v is not None}
stistools.x1d.x1d(crj_ours, output=x1d_out, verbose=False, **kw)
print("wrote:", os.path.basename(x1d_out), "with", kw)
