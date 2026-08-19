import os, sys, shutil
from astropy.io import fits
import calcos

# runs in wsl (surf_uv env), lref exported before calling.
# custom COS extraction: edit the source-box HEIGHT in the XTRACTAB, point a copy of the
# corrtag at the edited table, re-run calcos.
# args: corrtag_src, outdir, new_height, xtractab_src, workdir
corrtag_src, outdir = sys.argv[1], sys.argv[2]
new_height = float(sys.argv[3])
xt_src, workdir = sys.argv[4], sys.argv[5]

os.makedirs(workdir, exist_ok=True)
h0 = fits.getheader(corrtag_src)

# edit the xtractab rows matching this mode (all 3 NUV stripes)
new_xt = os.path.join(workdir, "edit_1dx.fits")
with fits.open(xt_src) as f:
    d = f[1].data
    sel = (d['OPT_ELEM'] == h0['OPT_ELEM']) & (d['CENWAVE'] == h0['CENWAVE']) & (d['APERTURE'] == 'PSA')
    print("orig heights:", list(d['HEIGHT'][sel]), "-> new", new_height)
    d['HEIGHT'][sel] = new_height
    f.writeto(new_xt, overwrite=True)

# copy the corrtag and aim its header at our edited table
ct = os.path.join(workdir, os.path.basename(corrtag_src))
shutil.copy(corrtag_src, ct)
fits.setval(ct, 'XTRACTAB', value=new_xt, ext=0)

if os.path.exists(outdir):
    shutil.rmtree(outdir)
os.makedirs(outdir, exist_ok=True)
calcos.calcos(ct, verbosity=0, outdir=outdir)
print("custom calcos done, height", new_height)
