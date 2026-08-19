import os, sys, glob, shutil, subprocess, time
import pandas as pd
from astropy.io import fits
import calcos

# per-association default calcos batch for the whole cos catalog.
# reads asn_manifest.csv, stages each asn + member rawtags, runs bestrefs + calcos -> x1dsum.
# resumable: skips an asn if its x1dsum already exists. optional 2nd arg = limit for smoke tests.
# run from wsl surf_uv with lref + CRDS env exported.

base = os.path.abspath(sys.argv[1]) if len(sys.argv) > 1 else 'data/cos_catalog'
limit = int(sys.argv[2]) if len(sys.argv) > 2 else None   # smoke-test: process only N associations
mast = os.path.join(base, 'mastDownload', 'HST')
mani = pd.read_csv(os.path.join(base, 'asn_manifest.csv'))
if limit: mani = mani.head(limit)
outroot = os.path.join(base, 'reduced')
workroot = os.path.join(base, 'work')
os.makedirs(outroot, exist_ok=True)
os.makedirs(workroot, exist_ok=True)

logf = open(os.path.join(base, 'reduce_cos_batch.log'), 'a')
def log(*a):
    msg = ' '.join(str(x) for x in a)
    print(msg); logf.write(msg + '\n'); logf.flush()

def find_file(root, suffix):
    p = os.path.join(mast, root, f'{root}_{suffix}')
    return p if os.path.exists(p) else None

log(f'\n=== batch start {time.ctime()}   {len(mani)} associations ===')
done = fail = skip = 0
for _, r in mani.iterrows():
    root = r['root']; sn = r['sn']; grat = r['grating']
    asn = os.path.join(mast, root, f'{root}_asn.fits')   # asn lives in its own rootname folder
    outdir = os.path.join(outroot, f'{sn}_{root}_{grat}')

    if not os.path.exists(asn):
        log(f'FAIL {sn} {root} - asn not found at {asn}'); fail += 1; continue
    if glob.glob(os.path.join(outdir, '*_x1dsum*.fits')):
        log(f'skip {sn} {root} ({grat}) - x1dsum exists'); skip += 1; continue

    # stage asn + members into a clean workdir
    work = os.path.join(workroot, root)
    if os.path.exists(work): shutil.rmtree(work)
    os.makedirs(work)
    shutil.copy(asn, work)
    mem = fits.getdata(asn)
    members = [m['MEMNAME'].strip().lower() for m in mem if m['MEMTYPE'].strip().upper().startswith('EXP')]
    staged = []
    for mn in members:
        got = False
        for suf in ['rawtag.fits', 'rawtag_a.fits', 'rawtag_b.fits', 'spt.fits']:
            f = find_file(mn, suf)
            if f:
                shutil.copy(f, work); got = got or 'rawtag' in suf
        staged.append(got)
    if not any(staged):
        log(f'FAIL {sn} {root} - no member rawtags found'); fail += 1; continue

    # bestrefs: update the rawtag headers to the best ref files and sync them into the crds cache
    raws = glob.glob(os.path.join(work, '*_rawtag*.fits'))
    try:
        subprocess.run(['crds', 'bestrefs', '--files', *raws,
                        '--sync-references=1', '--update-bestrefs'],
                       check=True, capture_output=True, text=True)
    except subprocess.CalledProcessError as e:
        log(f'FAIL {sn} {root} - bestrefs: {e.stderr[-300:]}'); fail += 1; continue

    # calcos wants an empty outdir
    if os.path.exists(outdir): shutil.rmtree(outdir)
    os.makedirs(outdir)
    asn_in = os.path.join(work, os.path.basename(asn))
    try:
        calcos.calcos(asn_in, verbosity=0, outdir=outdir)
        log(f'OK   {sn} {root} ({grat})'); done += 1
    except Exception as e:
        log(f'FAIL {sn} {root} - calcos: {repr(e)[:300]}'); fail += 1
    finally:
        shutil.rmtree(work, ignore_errors=True)

log(f'=== batch done {time.ctime()}  ok={done} fail={fail} skip={skip} ===')
logf.close()
