#!/usr/bin/env python3
# unified pipeline driver: STIS/COS reduction (WSL surf_uv) + analysis products (Windows 3.14).
# run from WORK root in Git Bash:  python3.exe scripts/rerun_all.py
# all logs go to output/logs/  -- nothing dumped to WORK root or output/ root.

import argparse, os, subprocess, sys, time, re
import paths

HERE = os.path.dirname(os.path.abspath(__file__))
WORK = os.path.dirname(HERE)
LOG_DIR = os.path.join(paths.OUT, 'logs')

# WSL mount path for WORK -- used to call run_full_catalog.py inside surf_uv.
# if the repo moves, update this one line.
WORK_WSL = '/mnt/c/Users/eluru/CALTECH_SURF_2026/WORK'


def _run(cmd, logfile, label):
    # run cmd, streaming stdout+stderr to both console and logfile in real time.
    os.makedirs(os.path.dirname(logfile), exist_ok=True)
    t0 = time.time()
    print(f'\n==> {label}', flush=True)
    with open(logfile, 'w') as lf:
        lf.write(f'$ {" ".join(str(c) for c in cmd)}\n'); lf.flush()
        proc = subprocess.Popen(
            cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
            bufsize=1, universal_newlines=True, cwd=WORK,
        )
        for line in proc.stdout:
            sys.stdout.write(line); sys.stdout.flush()
            lf.write(line); lf.flush()
        proc.wait()
    elapsed = time.time() - t0
    status = 'OK' if proc.returncode == 0 else f'FAILED (exit {proc.returncode})'
    print(f'  [{status}] {elapsed:.0f}s  log: {os.path.relpath(logfile, WORK)}', flush=True)
    return proc.returncode


def main():
    ap = argparse.ArgumentParser(
        description='run the full pipeline: STIS/COS reduction then analysis products')
    ap.add_argument('--skip-reduction', action='store_true', dest='skip_reduction',
                    help='skip the STIS/COS reduction; run analysis products only')
    ap.add_argument('--reduction-only', action='store_true', dest='reduction_only',
                    help='run only the reduction; skip analysis products')
    ap.add_argument('--workers', type=int, default=4,
                    help='parallel workers for the reduction (default 4)')
    ap.add_argument('--only', default=None, help='comma-separated SN names for a targeted rerun')
    ap.add_argument('--cos-only', action='store_true', dest='cos_only',
                    help='rebuild only COS products in the reduction step')
    a = ap.parse_args()

    os.makedirs(LOG_DIR, exist_ok=True)

    if not a.skip_reduction:
        wsl_py = f'{WORK_WSL}/scripts/run_full_catalog.py'
        wsl_cmd = (
            f'source ~/miniforge3/etc/profile.d/conda.sh && conda activate surf_uv'
            f' && PYTHONUNBUFFERED=1 python {wsl_py} --workers {a.workers}'
        )
        if a.only:
            wsl_cmd += f' --only {a.only}'
        if a.cos_only:
            wsl_cmd += ' --cos-only'

        rc = _run(
            ['wsl.exe', '-d', 'Debian', '--', 'bash', '-lc', wsl_cmd],
            os.path.join(LOG_DIR, 'reduction.log'),
            'reduction  (WSL surf_uv)',
        )
        # halt before analysis products if reduction failed (bad output would propagate)
        if rc != 0 and not a.reduction_only:
            print(f'ERROR: reduction exited {rc}. check output/logs/reduction.log')
            sys.exit(rc)

    if a.reduction_only:
        return

    py = sys.executable
    # order matters: lya_nhi writes lya_nhi_summary, which emission_products reads for F_ismcorr and
    # absorption_products joins alongside the ism CoG -- so ism (NUV) + lya_nhi run before emission and
    # absorption. ism --fuv reads the NUV ism_cog_summary for the Fe II 1608 cross-check so it runs after
    # the NUV ism pass. product_index runs LAST (scans every SN dir + sweeps stale orphans once all products
    # exist). host-sync is network-free.
    analysis = [
        ('ism.py',                       [py, os.path.join(HERE, 'ism.py')]),
        ('lya_nhi.py',                   [py, os.path.join(HERE, 'lya_nhi.py')]),
        ('emission_products.py',         [py, os.path.join(HERE, 'emission_products.py')]),
        ('absorption_products.py',       [py, os.path.join(HERE, 'absorption_products.py')]),
        ('ism.py --fuv',                 [py, os.path.join(HERE, 'ism.py'), '--fuv']),
        ('catalog_clean.py --host-sync', [py, os.path.join(HERE, 'catalog_clean.py'), '--host-sync']),
        ('product_index.py',             [py, os.path.join(HERE, 'product_index.py')]),
    ]
    any_fail = False
    for name, cmd in analysis:
        slug = re.sub(r'[^A-Za-z0-9]+', '_', name.replace('.py', '')).strip('_')
        rc = _run(cmd, os.path.join(LOG_DIR, f'{slug}.log'), name)
        if rc != 0:
            print(f'WARNING: {name} exited {rc}')
            any_fail = True

    if any_fail:
        print('\nsome analysis scripts reported errors -- check output/logs/')
    else:
        print('\nall done.  logs -> output/logs/')


if __name__ == '__main__':
    main()
