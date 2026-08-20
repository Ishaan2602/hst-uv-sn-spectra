import os, sys, glob, json, argparse, subprocess, time, shutil
import pandas as pd
from concurrent.futures import ProcessPoolExecutor, as_completed

# top-level driver for the full-catalog run. per SN: download STIS -> reduce STIS ->
# build COS products -> assemble SN-level products. resumable (MAST skips existing files),
# per-SN output override (old STIS tree deleted before each re-run so nothing stale persists).
# runs in the surf_uv wsl env (has astroquery + hstcal + astropy).
#
# parallelization: each SN is fully independent (separate data/ and output/ trees), so we
# run up to --workers processes in parallel. each worker writes to a tmp per-SN log file;
# the main process appends those to run_full_catalog.log serially as each worker finishes,
# so the log stays coherent and the existing 'grep built products' progress check still works.
# default 4 workers keeps us under MAST's concurrent-request threshold.
HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
from paths import CATALOG as CLEAN, OUT
DATA = os.path.join(ROOT, 'data')
COSRED = os.path.join(DATA, 'cos_catalog', 'reduced')
PY = sys.executable

# explosion-epoch override (mjd) for the day-phase labels. 2023IXF has a well-determined explosion
# epoch (better than discovery); every other SN falls back to its TNS discovery date (tns_disc_mjd
# from the cleaned catalog), i.e. the labels are "days since TNS discovery".
EXPL_MJD = {'SN2023IXF': 60082.79}

# SNe where early-time data is pre-reduced and should not be re-downloaded or reduced.
EARLY_SKIP = {'SN2023IXF': os.path.join(DATA, 'earlytime_2023ixf')}


def sh(args, log):
    log.write('$ ' + ' '.join(args) + '\n'); log.flush()
    r = subprocess.run(args, cwd=ROOT, capture_output=True, text=True)
    log.write(r.stdout[-4000:] + '\n' + r.stderr[-2000:] + '\n'); log.flush()
    return r.returncode


def script(name):
    return os.path.join(HERE, name)


def cos_dirs_for(name):
    # read-only glob; safe to call from multiple workers simultaneously
    return glob.glob(os.path.join(COSRED, f'{name}_*'))


def run_sn_worker(row_dict, redo=False, out=OUT, cos_only=False):
    """Process one SN. Called by each pool worker. Returns (name, logpath, elapsed_s).

    Writes ALL output to a per-SN tmp log (OUT/_run_sn_<name>.log) so concurrent
    workers never interleave writes to the shared main log.  The main process appends
    the tmp log to run_full_catalog.log after this function returns.
    """
    t0 = time.time()
    r = row_dict
    name = str(r['name']); z = float(r['z']); expl = EXPL_MJD.get(name)
    if expl is None:                                   # fall back to the TNS discovery epoch
        dm = r.get('tns_disc_mjd')
        try:
            expl = float(dm) if dm is not None and str(dm).strip() not in ('', 'nan', 'None') else None
        except (TypeError, ValueError):
            expl = None
    has_stis = 'STIS' in str(r.get('instr', ''))
    has_cos = len(cos_dirs_for(name)) > 0
    base = os.path.join(DATA, name)

    logpath = os.path.join(out, f'_run_sn_{name}.log')
    with open(logpath, 'w') as log:
        log.write(f'\n===== {name}  z={z}  STIS={has_stis} COS={has_cos} =====\n')

        if has_stis and not cos_only:
            sh([PY, script('download_all.py'), name], log)
            stis_out = os.path.join(out, name, 'STIS')
            if os.path.isdir(stis_out):
                shutil.rmtree(stis_out)
            if os.path.isdir(base):
                args = [PY, script('reduce_stis_batch.py'), base, name, '--outroot', out, '--z', str(z)]
                if expl: args += ['--expl-mjd', str(expl)]
                if name == 'SN2023IXF':
                    args += ['--skip-proposid', '17205']
                sh(args, log)

        if has_cos:
            cos_out = os.path.join(out, name, 'COS')
            if os.path.isdir(cos_out): shutil.rmtree(cos_out)
            args = [PY, script('cos_products.py'), name, '--outroot', out, '--z', str(z)]
            if expl: args += ['--expl-mjd', str(expl)]
            sh(args, log)

        if os.path.isdir(os.path.join(out, name)):
            args = [PY, script('build_products.py'), name, '--outroot', out, '--z', str(z),
                    '--early-skip-dir', EARLY_SKIP.get(name, '')]
            if expl:
                args += ['--expl-mjd', str(expl)]
            sh(args, log)

    return name, logpath, time.time() - t0


def presync_crds(rows, logf):
    # warm the CRDS cache ONCE, up front, so the parallel workers never sync concurrently - parallel
    # assign_bestrefs(sync_references=True) races on the shared cache and silently drops gratings (very
    # likely the 2022MOX missing-legs defect). download only the sources not already local (a local
    # check, no MAST query), then do ONE batched bestrefs over every downloaded STIS science frame.
    import glob as _glob
    os.environ.setdefault('CRDS_SERVER_URL', 'https://hst-crds.stsci.edu')
    os.environ.setdefault('CRDS_PATH', os.path.join(ROOT, 'crds_cache'))
    names = [str(r['name']) for r in rows if 'STIS' in str(r.get('instr', ''))]

    def frames(name):
        hst = os.path.join(DATA, name, 'mastDownload', 'HST')
        return (_glob.glob(os.path.join(hst, '*', '*_crj.fits')) +
                _glob.glob(os.path.join(hst, '*', '*_flt.fits')))

    for name in names:                       # download only what's missing locally (skip the MAST query)
        if not frames(name):
            sh([PY, script('download_all.py'), name], logf)
    # one file per exposure root (crj+flt share refs) so bestrefs isn't computed twice per exposure
    files, seen = [], set()
    for name in names:
        for f in sorted(frames(name)):
            root = os.path.basename(f).replace('_crj.fits', '').replace('_flt.fits', '')
            if root not in seen:
                seen.add(root); files.append(f)
    logf.write(f'\n#### crds pre-sync: {len(files)} STIS exposures, {len(names)} SNe (batched) ####\n'); logf.flush()
    if files:
        import crds
        try:
            crds.assign_bestrefs(files, sync_references=True)
        except Exception as e:
            logf.write(f'crds pre-sync error: {repr(e)[:300]}\n')
    logf.write(f'crds pre-sync done {time.ctime()}\n'); logf.flush()
    print(f'[presync] warmed CRDS cache: {len(files)} exposures across {len(names)} SNe', flush=True)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--limit', type=int, default=None)
    ap.add_argument('--only', default=None, help='comma-separated SN names')
    ap.add_argument('--redo', action='store_true')
    ap.add_argument('--workers', type=int, default=4,
                    help='parallel worker processes (default 4; cap ~6 for MAST rate limits)')
    ap.add_argument('--outroot', default=OUT,
                    help='per-SN product root (pass a candidate output dir for a new reduction run; catalog stays in output)')
    ap.add_argument('--presync', action='store_true',
                    help='optional one-time serial CRDS cache warm-up before the pool (rarely needed)')
    ap.add_argument('--cos-only', action='store_true', dest='cos_only',
                    help='rebuild only COS products + build_products (skip STIS re-extraction); use when '
                         'a fix touches only the COS/coadd path and STIS output is provably unchanged')
    a = ap.parse_args()

    df = pd.read_csv(CLEAN)
    df = df.sort_values('n_spec').reset_index(drop=True)   # smallest-first so quick SNe fill early slots
    if a.only:
        want = set(a.only.split(','))
        df = df[df['name'].isin(want)]
    if a.limit:
        df = df.head(a.limit)

    rows = [dict(r) for _, r in df.iterrows()]
    if a.cos_only:                                          # only touch SNe that actually have COS data
        rows = [r for r in rows if len(cos_dirs_for(str(r['name']))) > 0]
    n_total = len(rows)

    os.makedirs(a.outroot, exist_ok=True)
    logf = open(os.path.join(a.outroot, 'run_full_catalog.log'), 'a')
    logf.write(f'\n#### run start {time.ctime()}  {n_total} SNe  workers={a.workers} ####\n')
    logf.flush()

    if a.presync:
        presync_crds(rows, logf)

    n_done = 0
    with ProcessPoolExecutor(max_workers=a.workers) as pool:
        # submit all at once; pool internally queues extras when workers are busy
        futures = {pool.submit(run_sn_worker, row, a.redo, a.outroot, a.cos_only): row['name'] for row in rows}
        for fut in as_completed(futures):
            sn_name = futures[fut]
            try:
                name, logpath, elapsed = fut.result()
                # append this SN's log to the shared log -- only the main process does this,
                # so no concurrent writes to run_full_catalog.log
                try:
                    with open(logpath) as f:
                        logf.write(f.read())
                    os.remove(logpath)
                except OSError:
                    pass
                n_done += 1
                logf.flush()
                print(f'[{n_done}/{n_total}] done {name} ({elapsed:.0f}s)', flush=True)
            except Exception as e:
                logf.write(f'\nERROR {sn_name}: {repr(e)[:300]}\n')
                logf.flush()
                n_done += 1
                print(f'[{n_done}/{n_total}] ERROR {sn_name}: {repr(e)[:200]}', flush=True)

    logf.write(f'#### run done {time.ctime()} ####\n')
    logf.close()
    print('RUNDONE')


if __name__ == '__main__':
    main()
