import os, sys
import stistools.ocrreject, stistools.x1d

# param sweep driver for the stis extraction play. runs in wsl (surf_uv env).
# literal args only - shell var/loop expansion gets eaten crossing git-bash -> wsl, so we loop here in python.
obsdir, root = sys.argv[1], sys.argv[2]
flt = f"{obsdir}/{root}_flt.fits"
crj = f"{obsdir}/230_crj_ours.fits"
if not os.path.exists(crj):
    stistools.ocrreject.ocrreject(flt, crj, verbose=False)

def run(out, **kw):
    p = f"{obsdir}/{out}"
    if os.path.exists(p):
        os.remove(p)
    stistools.x1d.x1d(crj, output=p, verbose=False, **kw)
    print("wrote", out)

bg = dict(a2center=894, maxsrch=0, bk1offst=-14, bk2offst=14, bk1size=10, bk2size=10)
run("230_x1d_closebg.fits", extrsize=5, **bg)          # bg straddling trace, for the extraction-region viz
for e in (1, 3, 5, 7, 9, 11):                           # aperture width sweep
    run(f"230_x1d_ext{e}.fits", extrsize=e, **bg)
for c in (890, 894, 898):                               # miscentering test (true center ~894)
    run(f"230_x1d_cen{c}.fits", extrsize=7, **{**bg, "a2center": c})
print("ALLDONE")
