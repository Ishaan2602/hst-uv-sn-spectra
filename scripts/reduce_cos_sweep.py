import os, sys, shutil
from astropy.io import fits
import calcos

# runs in wsl (surf_uv env), lref exported before calling.
# sweeps the COS source-box HEIGHT: for each height, edit the XTRACTAB, point a corrtag copy at it, run calcos.
# args: corrtag_src, xtractab_src, workdir, h1 h2 h3 ...
corrtag_src, xt_src, workdir = sys.argv[1], sys.argv[2], sys.argv[3]
heights = [int(h) for h in sys.argv[4:]]

h0 = fits.getheader(corrtag_src)
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

    ct = f"{sub}/{os.path.basename(corrtag_src)}"
    shutil.copy(corrtag_src, ct)
    fits.setval(ct, 'XTRACTAB', value=new_xt, ext=0)

    calcos.calcos(ct, verbosity=0, outdir=f"{sub}/out")
    print("wrote h", hgt)
print("SWEEPDONE")
