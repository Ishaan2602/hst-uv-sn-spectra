import os, sys, glob, shutil
import numpy as np
from astropy.io import fits
import calcos

# HEIGHT sweep for one cos association. edits xtractab, restages, re-runs calcos per height step.
# usage: python reduce_cos_sweep.py <base_dir> <asn_root> <h1,h2,...>

base = sys.argv[1]
root = sys.argv[2]
heights = [int(h) for h in sys.argv[3].split(',')]

mast = os.path.join(base, 'mastDownload', 'HST')
asn = os.path.join(mast, root, f'{root}_asn.fits')
work_base = os.path.join(base, 'work')
sweep_out = os.path.join(base, 'sweep', root)
os.makedirs(sweep_out, exist_ok=True)

# get grating/cenwave/XTRACTAB from the first rawtag member
mem = fits.getdata(asn)
members = [m['MEMNAME'].strip().lower() for m in mem if m['MEMTYPE'].strip().upper().startswith('EXP')]
rt0 = None
for mn in members:
    for suf in ['rawtag_a.fits', 'rawtag.fits']:
        p = os.path.join(mast, mn, f'{mn}_{suf}')
        if os.path.exists(p): rt0 = p; break
    if rt0: break
assert rt0, f'no rawtag found for {root}'
h0 = fits.getheader(rt0, ext=0)
opt_elem = h0['OPT_ELEM']; cenwave = int(h0['CENWAVE'])
xt_ref = h0['XTRACTAB'].split('$')[-1]
xt_src = os.path.join('crds_cache/references/hst/cos', xt_ref)
print(f'{root}: {opt_elem} cenwave={cenwave}, XTRACTAB={xt_ref}')

def find_file(mn, suf):
    p = os.path.join(mast, mn, f'{mn}_{suf}')
    return p if os.path.exists(p) else None

logf = open(os.path.join(sweep_out, 'sweep.log'), 'a')
def log(*a):
    msg=' '.join(str(x) for x in a); print(msg); logf.write(msg+'\n'); logf.flush()

for h in heights:
    if h % 2 == 0: h += 1   # HEIGHT must be odd
    out = os.path.join(sweep_out, f'h{h:03d}')
    if glob.glob(os.path.join(out, '*_x1dsum*.fits')):
        log(f'skip h={h} (x1dsum exists)'); continue

    # edit the XTRACTAB: copy + change HEIGHT for PSA rows matching this mode
    new_xt = os.path.join(sweep_out, f'xtractab_h{h:03d}.fits')
    with fits.open(xt_src) as f:
        d = f[1].data
        m = np.array([str(a).strip()=='PSA' and str(e).strip()==opt_elem
                      for a,e in zip(d['APERTURE'], d['OPT_ELEM'])])
        m &= (d['CENWAVE'] == cenwave)
        orig = list(d['HEIGHT'][m])
        d['HEIGHT'][m] = h
        f.writeto(new_xt, overwrite=True)
    log(f'h={h}: edited HEIGHT {orig} -> {h}')

    # stage work dir: asn + rawtags + spt
    work = os.path.join(work_base, f'{root}_sweep_h{h:03d}')
    if os.path.exists(work): shutil.rmtree(work)
    os.makedirs(work)
    shutil.copy(asn, work)
    for mn in members:
        for suf in ['rawtag.fits','rawtag_a.fits','rawtag_b.fits','spt.fits']:
            f = find_file(mn, suf)
            if f: shutil.copy(f, work)
    # point each rawtag in workdir at our edited XTRACTAB (absolute path)
    new_xt_abs = os.path.abspath(new_xt)
    for rt in glob.glob(os.path.join(work, '*rawtag*.fits')):
        fits.setval(rt, 'XTRACTAB', value=new_xt_abs, ext=0)

    os.makedirs(out, exist_ok=True)
    asn_in = os.path.join(work, os.path.basename(asn))
    try:
        calcos.calcos(asn_in, verbosity=0, outdir=out)
        log(f'h={h}: ok -> {out}')
    except Exception as e:
        log(f'h={h}: FAIL {repr(e)[:200]}')
    finally:
        shutil.rmtree(work, ignore_errors=True)

logf.close()
print('sweep done')
