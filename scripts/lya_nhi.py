#!/usr/bin/env python3
# automated N(HI) measurement from damped Lya absorption for photospheric-backlight epochs.
# applies to all SNe with an early (phase < 50d) G140L spectrum bright enough to see damping wings.
# model: F = (c0 + c1*(w-1215.67) + c2*(w-1215.67)^2) * exp(-tau_Lya(w, logN, b=25, vabs=0))
# b is fixed (damping wings are b-independent), vabs fixed at 0 (host frame).
# syst_vabs = 0.5 * abs(logN(vabs=-300) - logN(vabs=+300)) captures host-ISM velocity dispersion.
#
# late-time epochs where CSM Lya emission is active need the joint emission+absorption model;
# those are NOT handled here -- only photospheric-backlight epochs (phase < 50 days).
# curated values in reference/ism_columns.csv override these automated results.
#
# usage:  python lya_nhi.py [SN ...]     (no args -> all SNe in output5 with G140L coverage)

import os, csv, glob, re, json, datetime
import numpy as np
from scipy.optimize import curve_fit
from scipy.special import wofz

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
from paths import OUT, CATALOG, LYA_NHI_SUMMARY as SUMMARY

LYA0      = 1215.67
C_KMS     = 2.99792458e5
F_LYA     = 0.4164
GAMMA_LYA = 6.265e8   # s^-1

# pipeline gates
FLUX_GATE       = 0.1e-15    # min median flux in 1225-1290 A rest (the backlight quality check)
PHASE_MAX_CLEAN = 50         # max phase (days) for a photospheric-backlight epoch
CHI2_GATE       = 0.20       # max reduced chi2
LOG_N_LO        = 18.0
LOG_N_HI        = 22.5
VABS_SYST_RANGE = 300.0      # ±km/s for the v_abs systematic (host-ISM velocity dispersion)

# --- catalog -----------------------------------------------------------------------------------
cat = {}
with open(CATALOG) as fh:
    for row in csv.DictReader(fh):
        cat[row["name"].upper()] = row

def _zof(sn):
    v = cat.get(sn.upper(), {}).get("z")
    return float(v) if v else 0.0

def _tnstype(sn):
    return cat.get(sn.upper(), {}).get("tns_type") or ""


# --- spectrum I/O -------------------------------------------------------------------------------

def load_spec(path, z):
    d = np.loadtxt(path, comments="#")
    w = d[:, 0] / (1.0 + z)
    f = d[:, 1]
    return w, f


def find_g140l(sn_dir):
    return sorted(
        [p.replace("\\", "/") for p in glob.glob(f"{OUT}/{sn_dir}/**/G140L/*_native.txt", recursive=True)
         if "epochcoadd" not in p and "/epochs/" not in p],
        key=lambda p: _phase_of(p) or 1e9
    )


def _phase_of(path):
    m = re.search(r"day([0-9.]+)", path)
    return float(m.group(1)) if m else None


def backlight_flux(path, z):
    """median flux in 1225-1290 A rest frame -- the safe Lya red wing."""
    try:
        w, f = load_spec(path, z)
        mask = (w > 1225) & (w < 1290) & np.isfinite(f) & (f > 0)
        return float(np.nanmedian(f[mask])) if mask.sum() >= 5 else None
    except Exception:
        return None


# --- Voigt tau + fit ---------------------------------------------------------------------------

def lya_tau(w, logN, b, vabs):
    lam0 = LYA0 * (1 + vabs / C_KMS)
    dlD  = lam0 * b / C_KMS
    x    = (w - lam0) / dlD
    a    = GAMMA_LYA * (lam0 * 1e-8) / (4 * np.pi * (b * 1e5))
    return 1.4973e-15 * 10**logN * F_LYA * lam0 / b * np.real(wofz(x + 1j * a))


def fit_nhi(w, f, logN0=20.5, b=25.0, vabs=0.0, wlo=1185, whi=1252):
    """
    Fit N(HI) to a photospheric-backlight Lya absorption profile.
    Returns (logN_fit, redchi2) or (None, None) on failure.
    """
    mask = (w > wlo) & (w < whi) & np.isfinite(f) & (f > 0)
    if mask.sum() < 12:
        return None, None
    wm, fm = w[mask], f[mask]
    A0 = float(np.nanmax(fm))
    if A0 <= 0:
        return None, None

    def model(w_, c0, c1, c2, logN_):
        cont = c0 + c1*(w_ - LYA0) + c2*(w_ - LYA0)**2
        return cont * np.exp(-lya_tau(w_, logN_, b, vabs))

    p0 = [np.nanmedian(fm)/A0, 0.0, 0.0, logN0]
    bounds = ([0.0, -2.0, -2.0, LOG_N_LO], [3.0, 2.0, 2.0, LOG_N_HI])
    try:
        popt, _ = curve_fit(model, wm, fm/A0, p0=p0, bounds=bounds, maxfev=80000)
    except Exception:
        return None, None

    logN_fit = popt[3]
    resid    = fm/A0 - model(wm, *popt)
    redchi2  = float(np.nansum(resid**2) / max(len(wm) - len(p0), 1))
    return logN_fit, redchi2


# --- main scan ---------------------------------------------------------------------------------

def run_catalog(names=None):
    sne_dirs = sorted(d for d in os.listdir(OUT) if os.path.isdir(os.path.join(OUT, d)) and d.upper() in cat)
    if names:
        names_up = {n.upper() for n in names}
        sne_dirs = [d for d in sne_dirs if d.upper() in names_up]

    summ = []
    for sn_dir in sne_dirs:
        z = _zof(sn_dir)
        paths = find_g140l(sn_dir)
        if not paths:
            continue
        # try epochs in phase order; take the first (earliest) that passes all gates
        result = None
        for p in paths:
            ph = _phase_of(p)
            if ph is None or ph > PHASE_MAX_CLEAN:
                continue
            bl = backlight_flux(p, z)
            if bl is None or bl < FLUX_GATE:
                continue
            w, f = load_spec(p, z)
            logN, chi2 = fit_nhi(w, f, vabs=0.0)
            if logN is None or not (LOG_N_LO < logN < LOG_N_HI):
                continue
            # reject fits stuck at either bound (sign of a failed/unconstrained fit)
            if abs(logN - LOG_N_LO) < 0.1 or abs(logN - LOG_N_HI) < 0.1:
                continue
            if chi2 is not None and chi2 > CHI2_GATE:
                continue
            # syst_vabs: half-range over ±VABS_SYST_RANGE km/s
            logN_lo, _ = fit_nhi(w, f, logN0=logN, vabs=-VABS_SYST_RANGE)
            logN_hi, _ = fit_nhi(w, f, logN0=logN, vabs=+VABS_SYST_RANGE)
            syst = 0.5 * abs((logN_lo or logN) - (logN_hi or logN))
            result = {
                "sn":          sn_dir.upper(),
                "sn_type":     _tnstype(sn_dir),
                "generated":   datetime.date.today().isoformat(),
                "grating":     "G140L",
                "phase":       ph,
                "logN_HI":     round(logN, 3),
                "logN_HI_syst_vabs": round(syst, 3),
                "backlight_flux_e15": round(bl * 1e15, 3),
                "redchi2":     round(chi2, 4) if chi2 is not None else None,
                "method":      "continuum-backlight damped Lya; vabs=0 (host frame); syst from +-300 km/s",
                "note":        "automated; curated ism_columns.csv values override this when present",
            }
            break   # take only the earliest passing epoch

        if result is None:
            continue

        sn_out = os.path.join(OUT, sn_dir, f"{sn_dir}_lya_nhi.json")
        os.makedirs(os.path.dirname(sn_out), exist_ok=True)
        with open(sn_out, "w") as fh:
            json.dump(result, fh, indent=2)
        summ.append(result)
        print(f"  {sn_dir:16s} d{result['phase']:5.0f}  logN={result['logN_HI']:.2f} ± {result['logN_HI_syst_vabs']:.2f}(syst)"
              f"  chi2={result['redchi2']}  -> {os.path.relpath(sn_out, ROOT)}")

    # write catalog summary
    if summ:
        cols = ["sn", "sn_type", "grating", "phase", "logN_HI", "logN_HI_syst_vabs",
                "backlight_flux_e15", "redchi2", "generated"]
        with open(SUMMARY, "w", newline="") as fh:
            wr = csv.DictWriter(fh, fieldnames=cols, extrasaction="ignore")
            wr.writeheader(); wr.writerows(summ)
        print(f"\nsummary -> {os.path.relpath(SUMMARY, ROOT)} ({len(summ)} sne)")
    else:
        print("no SNe passed all gates")

    return summ


if __name__ == "__main__":
    import sys
    run_catalog(sys.argv[1:] or None)
