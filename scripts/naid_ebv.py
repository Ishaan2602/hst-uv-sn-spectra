import os, csv, re, time, json, datetime
import numpy as np
from scipy.integrate import trapezoid
from io import StringIO
from wiserep_api.spectra import get_target_response, get_response as _wiserep_get
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import plotstyle; plotstyle.apply()

# Na I D EW from WISeREP ground-based spectra -> Stritzinger+2018 E(B-V).
# measures the HOST component: MW Na I D is at z=0, host component is at SN redshift.
# calibration: A_V = 0.78 * EW_total(D1+D2) [A] (Stritzinger+2018, CSP-I stripped-envelope SNe)
# E(B-V) = A_V / R_V = A_V / 3.1
# WARNING: calibration scatter is large (~0.2 mag); use as supplementary, not primary reddening.

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

from paths import CATALOG as CAT, HOST_EBV, CATDIR, OUT
NAID_CACHE = os.path.join(CATDIR, ".naid_cache")
NAID_SUMMARY = os.path.join(CATDIR, "naid_ebv_summary.csv")

D2_REST = 5889.951    # Na I D2 air wavelength (A)
D1_REST = 5895.924    # Na I D1
RV = 3.1
STRITZ_SLOPE = 0.78   # A_V per angstrom of EW (Stritzinger+2018)


def get_spectra_list(iau_name):
    """return list of {url, fname, tel_inst} for available ascii spectra from wiserep."""
    try:
        resp = get_target_response(iau_name)
    except Exception:
        return []
    if resp is None:
        return []

    # URL pattern: asciifile=https%3A//www.wiserep.org/system/files/...
    # wiserep_api internally uses http:// not https:// (https returns 403)
    raw_urls = []
    for st in resp.text.split("asciifile=https%3A//"):
        url = st.split('"')[0]
        if url and "&amp" not in url and "DOCTYPE" not in url:
            raw_urls.append(url)
    raw_urls = list(dict.fromkeys(raw_urls))

    # build tel_inst lookup from the spec table
    tel_map = {}
    try:
        tbl = pd.read_html(StringIO(resp.text), match='Spec. ID')[0]
        # handle multi-index columns (wiserep table uses merged headers)
        if isinstance(tbl.columns, pd.MultiIndex):
            tbl.columns = [c[0] for c in tbl.columns]
        if 'Spectrum ascii File' in tbl.columns and 'Tel / Inst' in tbl.columns:
            for _, row in tbl.iterrows():
                fn = str(row.get('Spectrum ascii File', ''))
                ti = str(row.get('Tel / Inst', ''))
                if fn and fn != 'nan':
                    tel_map[fn.strip()] = ti.strip()
    except Exception:
        pass

    out = []
    for raw in raw_urls:
        fname = raw.split('/')[-1]
        out.append({'url': 'http://' + raw, 'fname': fname,
                    'tel_inst': tel_map.get(fname, '')})
    return out


def load_spectrum_url(url):
    """download and parse a WISeREP ascii spectrum. returns (wvl, flx, err) or None."""
    os.makedirs(NAID_CACHE, exist_ok=True)
    fname = url.split('/')[-1].split('?')[0]
    local = os.path.join(NAID_CACHE, fname)
    if os.path.exists(local):
        raw = open(local, encoding='utf-8', errors='replace').read()
    else:
        r = _wiserep_get(url)    # uses TNS bot User-Agent that WISeREP accepts
        if r is None or r.status_code != 200:
            return None
        raw = r.text
        with open(local, 'w', encoding='utf-8') as f:
            f.write(raw)
        time.sleep(0.15)

    rows = []
    for line in raw.splitlines():
        s = line.strip()
        if not s or s[0] in ('#', '\\', '|', 'W', 'w'):
            continue
        parts = s.split()
        try:
            v = [float(p) for p in parts[:3]]
            rows.append(v if len(v) == 3 else v + [np.nan])
        except ValueError:
            continue
    if len(rows) < 20:
        return None
    arr = np.array(rows)
    # sanity: wavelengths should be monotonic and in optical range
    w = arr[:, 0]
    if w[-1] < w[0]:   # descending wavelength -- flip
        arr = arr[::-1]
        w = arr[:, 0]
    if w[0] > 10000 or w[-1] < 3000:   # weird units or non-optical
        return None
    return arr[:, 0], arr[:, 1], arr[:, 2]


def estimate_resolution(wvl):
    dpx = np.median(np.diff(wvl))
    return np.median(wvl) / (2.35 * dpx) if dpx > 0 else 0


def covers_naid(wvl, z, margin=20.0):
    lam_d2 = D2_REST * (1 + z)
    lam_d1 = D1_REST * (1 + z)
    return wvl[0] < lam_d2 - margin and wvl[-1] > lam_d1 + margin


def measure_naid_ew(wvl, flx, z, cont_hw=60.0, int_hw=None):
    """
    measure HOST Na I D blended EW.
    - wvl/flx in observed frame; z = SN redshift so host D2/D1 sit at D_REST*(1+z)
    - cont_hw: half-width of each continuum window [A]
    - int_hw: half-width of integration window [A]; None = auto from resolution (2*FWHM, 2-12 A)
    returns dict: ew_a, ew_err, status, r_est, int_hw, mw_sep_a, [flags]
    """
    lam_d2 = D2_REST * (1 + z)
    lam_d1 = D1_REST * (1 + z)
    mw_sep = lam_d2 - D2_REST    # how far MW D2 is from host D2 [A]

    r_est = int(estimate_resolution(wvl))
    if int_hw is None:
        fwhm_a = lam_d2 / r_est if r_est > 0 else 10.0
        int_hw = float(np.clip(2.0 * fwhm_a, 2.0, 12.0))

    # integration window covers D2 and D1 of the host
    w_lo = lam_d2 - int_hw
    w_hi = lam_d1 + int_hw

    # continuum windows:
    # blue: blueward of MW D2 absorption (avoid both MW absorptions completely)
    # red: redward of host D1 absorption
    b_hi = D2_REST - 8.0               # stop well before MW D2
    b_lo = b_hi - cont_hw              # go 60 A further blue
    r_lo = w_hi + 5.0                  # just redward of host integration window
    r_hi = r_lo + cont_hw

    mb = (wvl >= b_lo) & (wvl <= b_hi)
    mr = (wvl >= r_lo) & (wvl <= r_hi)
    mw = (wvl >= w_lo) & (wvl <= w_hi)

    if mb.sum() < 5 or mr.sum() < 5 or mw.sum() < 5:
        return {'status': 'bad_coverage', 'ew_a': np.nan, 'ew_err': np.nan,
                'r_est': r_est, 'int_hw': int_hw, 'mw_sep_a': mw_sep}

    c = np.polyfit(np.concatenate([wvl[mb], wvl[mr]]),
                   np.concatenate([flx[mb], flx[mr]]), 1)
    cont = np.polyval(c, wvl)
    norm = flx / cont

    ew = float(trapezoid(1.0 - norm[mw], wvl[mw]))

    # error: per-pixel noise in continuum region -> propagate through integration
    cont_rms = np.std(1.0 - norm[mb | mr])
    dpx = np.median(np.diff(wvl))
    ew_err = cont_rms * np.sqrt(mw.sum()) * dpx

    if ew < -3 * ew_err:
        status = 'emission'   # He I 5876 / CSM emission raises apparent flux above continuum
    elif ew <= 0:
        status = 'upper'
    elif ew < 2 * ew_err:
        status = 'marginal'
    else:
        status = 'detected'

    flags = []
    if D1_REST > w_lo:              # MW D1 sits inside the host integration window
        flags.append('mw_close')
    if r_est < 1000:
        flags.append('low_res')   # below minimum for narrow ISM feature measurement

    return {'status': status + ('+' + ','.join(flags) if flags else ''),
            'ew_a': ew, 'ew_err': ew_err, 'r_est': r_est, 'int_hw': int_hw,
            'mw_sep_a': mw_sep, 'cont_rms': cont_rms}


def ebv_from_ew(ew_a, ew_err_a=np.nan):
    """Stritzinger+2018: E(B-V) = (A_V slope / R_V) * EW."""
    ebv = STRITZ_SLOPE / RV * ew_a
    ebv_err = STRITZ_SLOPE / RV * ew_err_a if np.isfinite(ew_err_a) else np.nan
    return ebv, ebv_err


def _naid_diag_plot(name, wvl, flx, z, meas, outdir):
    # scrutiny plot behind the Na I D E(B-V): normalized flux with host D1/D2, MW rest positions,
    # the EW integration window, and the continuum windows all marked.
    lam_d2 = D2_REST * (1 + z); lam_d1 = D1_REST * (1 + z)
    int_hw = meas.get("int_hw") or 6.0
    w_lo, w_hi = lam_d2 - int_hw, lam_d1 + int_hw
    b_hi = D2_REST - 8.0; b_lo = b_hi - 60.0; r_lo = w_hi + 5.0; r_hi = r_lo + 60.0
    mb = (wvl >= b_lo) & (wvl <= b_hi); mr = (wvl >= r_lo) & (wvl <= r_hi)
    if mb.sum() < 3 or mr.sum() < 3:
        return None
    c = np.polyfit(np.concatenate([wvl[mb], wvl[mr]]), np.concatenate([flx[mb], flx[mr]]), 1)
    norm = flx / np.polyval(c, wvl)
    sel = (wvl >= b_lo - 5) & (wvl <= r_hi + 5)
    if sel.sum() < 5:
        return None
    ewv = meas.get("ew_a", np.nan); ebvv, _ = ebv_from_ew(ewv, meas.get("ew_err", np.nan))
    fig, ax = plt.subplots(figsize=(8, 4))
    ax.plot(wvl[sel], norm[sel], 'k', lw=0.9, drawstyle='steps-mid')
    ax.axhline(1.0, color='0.6', ls=':', lw=0.8)
    ax.axvspan(w_lo, w_hi, color='gold', alpha=0.12, label='EW integration')
    ax.axvspan(b_lo, b_hi, color='tab:blue', alpha=0.07)
    ax.axvspan(r_lo, r_hi, color='tab:blue', alpha=0.07, label='continuum')
    ax.axvline(lam_d2, color='crimson', ls='--', lw=1.0); ax.axvline(lam_d1, color='crimson', ls='--', lw=1.0)
    ax.axvline(D2_REST, color='0.5', ls=':', lw=0.8); ax.axvline(D1_REST, color='0.5', ls=':', lw=0.8)
    ax.text(lam_d2, 1.06, 'host D2/D1', fontsize=7, ha='center', color='crimson')
    ax.text(D2_REST, 1.06, 'MW rest', fontsize=7, ha='center', color='0.4')
    ax.set_title(f"{name} Na I D  EW={ewv:.2f} A  E(B-V)={ebvv:.3f}  ({meas.get('status','')})", fontsize=9)
    ax.set_xlabel(r'obs $\lambda$ [$\AA$]'); ax.set_ylabel(r'$I/I_0$')
    ax.set_xlim(b_lo - 5, r_hi + 5); ax.set_ylim(-0.1, 1.3); ax.legend(fontsize=7)
    os.makedirs(outdir, exist_ok=True)
    out = os.path.join(outdir, f'{name}_naid_ew.png')
    fig.tight_layout(); fig.savefig(out, dpi=110); plt.close(fig)
    return out


def _resolve_outdir(name):
    for cand in (name, name.upper().replace(' ', '')):
        d = os.path.join(OUT, cand)
        if os.path.isdir(d):
            return d
    return None


def best_spectrum_for_naid(spectra_list, z, min_r=600):
    """pick highest-resolution spectrum that covers the Na I D region. returns dict."""
    candidates = []
    for s in spectra_list:
        result = load_spectrum_url(s['url'])
        if result is None:
            continue
        wvl, flx, err = result
        if not covers_naid(wvl, z):
            continue
        r = estimate_resolution(wvl)
        if r >= min_r:
            candidates.append({'r': r, 'wvl': wvl, 'flx': flx, 'err': err,
                                'fname': s['fname'], 'tel_inst': s['tel_inst'], 'url': s['url']})
    if not candidates:
        return None
    return max(candidates, key=lambda x: x['r'])


def survey_source(iau_name, z):
    """full pipeline for one source. returns result dict."""
    spectra_list = get_spectra_list(iau_name)
    if not spectra_list:
        return {'iau_name': iau_name, 'z': z, 'status': 'no_spectra',
                'ew_a': np.nan, 'ew_err': np.nan, 'ebv': np.nan, 'ebv_err': np.nan,
                'r_est': 0, 'mw_sep_a': np.nan, 'fname': '', 'tel_inst': '', 'n_spectra': 0}

    best = best_spectrum_for_naid(spectra_list, z)
    if best is None:
        return {'iau_name': iau_name, 'z': z, 'status': 'no_coverage',
                'ew_a': np.nan, 'ew_err': np.nan, 'ebv': np.nan, 'ebv_err': np.nan,
                'r_est': 0, 'mw_sep_a': np.nan, 'fname': '', 'tel_inst': '',
                'n_spectra': len(spectra_list)}

    meas = measure_naid_ew(best['wvl'], best['flx'], z)
    ew = meas['ew_a']
    ew_err = meas['ew_err']
    ebv, ebv_err = ebv_from_ew(ew, ew_err) if np.isfinite(ew) else (np.nan, np.nan)

    # per-SN product + EW scrutiny plot when this SN has an HST output dir (where reddening actually matters)
    outdir = _resolve_outdir(iau_name)
    diag_rel = None
    if outdir is not None:
        diag = _naid_diag_plot(iau_name, best['wvl'], best['flx'], z, meas, os.path.join(outdir, 'reddening'))
        diag_rel = os.path.relpath(diag, OUT).replace('\\', '/') if diag else None
        jrec = {"sn": iau_name, "z": z, "status": meas['status'],
                "ew_a": None if not np.isfinite(ew) else round(ew, 3),
                "ew_err": None if not np.isfinite(ew_err) else round(ew_err, 3),
                "ebv": None if not np.isfinite(ebv) else round(ebv, 4),
                "ebv_err": None if not np.isfinite(ebv_err) else round(ebv_err, 4),
                "r_est": meas['r_est'], "tel_inst": best['tel_inst'],
                "method": "Na I D EW -> Stritzinger+2018 A_V=0.78*EW, E(B-V)=A_V/3.1",
                "generated": datetime.date.today().isoformat()}
        if diag_rel:
            jrec["diag_plot"] = diag_rel
        with open(os.path.join(outdir, f"{os.path.basename(outdir)}_naid_ebv.json"), "w") as fh:
            json.dump(jrec, fh, indent=2)

    return {'iau_name': iau_name, 'z': z, 'status': meas['status'],
            'ew_a': ew, 'ew_err': ew_err, 'ebv': ebv, 'ebv_err': ebv_err,
            'r_est': meas['r_est'], 'int_hw': meas.get('int_hw', np.nan),
            'mw_sep_a': meas['mw_sep_a'],
            'fname': best['fname'], 'tel_inst': best['tel_inst'],
            'n_spectra': len(spectra_list)}


def load_catalog_targets():
    """sources with z>0.003 and no curated host E(B-V). excludes pure Ia."""
    known = set()
    with open(HOST_EBV) as f:
        for row in csv.DictReader(r for r in f if not r.startswith('#')):
            if float(row['host_ebv']) > 0:
                known.add(row['name'])

    targets = []
    with open(CAT) as f:
        for row in csv.DictReader(f):
            name = row['name']
            z = float(row.get('z') or 0)
            otype = row.get('otype', '')
            if name in known:
                continue
            if z < 0.003 or z > 0.15:
                continue
            # check both otype and classification (case-insensitive); otype is often 'SN*' even for confirmed Ia
            type_str = (otype + ' ' + row.get('classification', '')).upper()
            if ('TYPE IA' in type_str or 'SN IA' in type_str) and 'CSM' not in type_str:
                continue
            targets.append((name, z, otype))
    return targets


def run_survey(targets=None, outfile=NAID_SUMMARY, verbose=True):
    if targets is None:
        targets = load_catalog_targets()

    results = []
    for i, (name, z, otype) in enumerate(targets):
        if verbose:
            print(f'[{i+1}/{len(targets)}] {name} z={z:.4f}', end='  ', flush=True)
        res = survey_source(name, z)
        res['otype'] = otype
        results.append(res)
        if verbose:
            st = res['status']
            ew = res['ew_a']
            ebv = res['ebv']
            r = res['r_est']
            ns = res.get('n_spectra', 0)
            if np.isfinite(ew):
                print(f'{st}  EW={ew:.3f}A  E(B-V)={ebv:.3f}  R~{r}  ({ns} spectra)')
            else:
                print(f'{st}  R~{r}  ({ns} spectra)')

    if results:
        keys = ['iau_name', 'z', 'otype', 'status', 'ew_a', 'ew_err', 'ebv', 'ebv_err',
                'r_est', 'mw_sep_a', 'fname', 'tel_inst', 'n_spectra']
        with open(outfile, 'w', newline='') as f:
            w = csv.DictWriter(f, fieldnames=keys, extrasaction='ignore')
            w.writeheader()
            w.writerows(results)
        if verbose:
            print(f'\n{len(results)} sources -> {outfile}')
    return results


def _patch_host_ebv(summary_csv, write=False):
    """fill-only merge of clean Na I D detections into reference/host_ebv.csv. never overwrites a curated
    row (name already present is skipped). dry-run by default: prints the diff and writes only if write=True."""
    # read existing curated values
    existing = {}
    rows_existing = []
    with open(HOST_EBV) as f:
        rdr = csv.DictReader(r for r in f if not r.startswith('#'))
        for row in rdr:
            existing[row['name']] = row
            rows_existing.append(row)

    additions = []
    with open(summary_csv) as f:
        for row in csv.DictReader(f):
            name = row['iau_name']
            st = row['status']
            if name in existing:                         # fill-only: never touch a curated literature value
                continue
            if st == 'detected':
                ebv = float(row['ebv']); ebv_err = float(row['ebv_err'])
                src = f'NaID_Stritz18 EW={float(row["ew_a"]):.3f}A R~{row["r_est"]} {row["tel_inst"]}'
            elif st == 'upper':
                ebv, ebv_err = 0.0, 0.0                   # no Na I D detected -> assume zero host reddening
                src = f'NaID_upper R~{row["r_est"]} {row["tel_inst"]}'
            else:
                continue
            additions.append({'name': name, 'host_ebv': f'{ebv:.4f}',
                              'host_ebv_err': f'{ebv_err:.4f}', 'host_ebv_src': src})

    if not additions:
        print('no clean Na I D detections to add (all covered SNe are already curated).')
        return

    # dry-run diff
    print(f'\n{"DRY RUN -- " if not write else ""}{len(additions)} fill-only additions to {os.path.basename(HOST_EBV)}:')
    for a in additions:
        print(f'  + {a["name"]:16s} host_ebv={a["host_ebv"]} +-{a["host_ebv_err"]}  [{a["host_ebv_src"]}]')
    if not write:
        print('\nno curated row is touched. re-run with --write to apply.')
        return

    rows_existing.extend(additions)
    header_lines = []
    with open(HOST_EBV) as f:
        for line in f:
            if line.startswith('#'):
                header_lines.append(line.rstrip())
            else:
                break

    with open(HOST_EBV, 'w', newline='') as f:
        for h in header_lines:
            f.write(h + '\n')
        keys = list(rows_existing[0].keys())
        w = csv.DictWriter(f, fieldnames=keys, extrasaction='ignore')
        w.writeheader()
        w.writerows(rows_existing)

    print(f'\npatched {HOST_EBV}: added {len(additions)} sources from Na I D survey (curated rows untouched).')


if __name__ == '__main__':
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument('--test', metavar='NAME')
    ap.add_argument('--z', type=float, default=None)
    ap.add_argument('--survey', action='store_true')
    ap.add_argument('--max-n', type=int, default=None)
    ap.add_argument('--patch-ebv', action='store_true',
                    help='fill-only merge of clean detections from naid_ebv_summary.csv into host_ebv.csv (dry-run)')
    ap.add_argument('--write', action='store_true',
                    help='with --patch-ebv, actually write the additions (default is dry-run diff only)')
    args = ap.parse_args()

    if args.test:
        z = args.z or 0.01
        res = survey_source(args.test, z)
        for k, v in res.items():
            print(f'  {k}: {v}')
    elif args.survey:
        targets = load_catalog_targets()
        if args.max_n:
            targets = targets[:args.max_n]
        run_survey(targets)
    elif args.patch_ebv:
        if not os.path.exists(NAID_SUMMARY):
            print('run --survey first')
        else:
            _patch_host_ebv(NAID_SUMMARY, write=args.write)
