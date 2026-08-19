import os, sys, shutil
from astropy.io import fits
import calcos

# runs in wsl (surf_uv env), lref exported before calling.
# FUV height sweep: edit the XTRACTAB source-box HEIGHT for FUVA+FUVB, point copies of both
# rawtags at it (boxcar), run calcos per height.
# args: rawtag_a, xtractab_src, workdir, h1 h2 ...
rawtag_a, xt_src, workdir = sys.argv[1], sys.argv[2], sys.argv[3]
heights = [int(h) for h in sys.argv[4:]]
rawtag_b = rawtag_a.replace('_rawtag_a', '_rawtag_b')

h0 = fits.getheader(rawtag_a)
os.makedirs(workdir, exist_ok=True)

for hgt in heights:
    sub = f"{workdir}/h{hgt}"
    if os.path.exists(sub):
        shutil.rmtree(sub)
    os.makedirs(sub)

    new_xt = f"{sub}/edit_1dx.fits"
    with fits.open(xt_src) as f:
        d = f[1].data
        sel = (d['OPT_ELEM'] == h0['OPT_ELEM']) & (d['CENWAVE'] == h0['CENWAVE']) & (d['APERTURE'] == 'PSA')
        d['HEIGHT'][sel] = hgt
        f.writeto(new_xt, overwrite=True)

    for rt in (rawtag_a, rawtag_b):
        dst = f"{sub}/{os.path.basename(rt)}"
        shutil.copy(rt, dst)
        fits.setval(dst, 'XTRACTAB', value=new_xt, ext=0)
        fits.setval(dst, 'XTRCTALG', value='BOXCAR', ext=0)

    calcos.calcos(f"{sub}/{os.path.basename(rawtag_a)}", verbosity=0, outdir=f"{sub}/out")
    print("wrote h", hgt)
print("SWEEPDONE")
