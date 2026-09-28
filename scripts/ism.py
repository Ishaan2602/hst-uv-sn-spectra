import os, glob, csv, re, argparse, datetime, json
import numpy as np
from scipy.optimize import curve_fit, least_squares
from scipy.integrate import quad
from scipy.ndimage import median_filter
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import plotstyle; plotstyle.apply()

# ism equivalent-width -> curve-of-growth -> column density, ported from ism_ew_cog_sandbox.ipynb.
# reads the observed-frame products in output, applies z from the catalog, measures deblended EWs,
# anchors the doppler b on the Fe II series, reads off per-ion column densities.

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
from paths import OUT, CATALOG as CAT, ISM_SUMMARY
LINES = os.path.join(ROOT, "linelists", "ism_lines.csv")

C_KMS = 2.99792458e5
TAU_K = 1.4973e-15          # tau0 = TAU_K * N * f * lam / b
LIN_K = 8.853e-21           # linear CoG W/lam = LIN_K * N f lam
AOD_K = 3.7679e14           # N = AOD_K/(f lam) * int tau dv (Savage&Sembach 1991; was 1.13e17 = 300x too high)
ANCHOR_TAU = 3.0            # CoG anchor: need >=1 detected Fe II line below this tau0 (linear/transition part) or N floats

# --- single source of truth for ISM-column adoption eligibility -------------------------------------
# whether an epoch's Fe II CoG column can represent the FOREGROUND ISM (vs the SN photosphere). used by
# BOTH the per-epoch cog.csv flag and the adopted-summary (finalize_and_write) so the two never disagree.
ADOPT_MIN_FE_SNR, ADOPT_MAX_LOGN_ERR, ADOPT_MAX_LOGN = 3.0, 0.3, 16.0
ADOPT_RESOLVED_SNR = {"LMC-SN1987A-STIS-2"}

def _num(x):
    try:
        v = float(x)
        return v if np.isfinite(v) else None
    except (TypeError, ValueError):
        return None

def adopt_eligible(sn, anchored, fe_snr, b, b_err, logN, logN_err, cat, curated_sns):
    # returns (eligible, reason). the reason is written into the per-epoch product so a reader sees WHY an
    # epoch is excluded (e.g. type_Ia_photosphere) instead of silently trusting a contaminated raw column.
    if not anchored:
        return False, "not_anchored"            # every Fe II line saturated -> N floats = photosphere signature
    if sn in ADOPT_RESOLVED_SNR:
        return False, "resolved_snr"
    if "IA" in (cat.get(sn.upper(), {}).get("tns_type") or "").upper().replace(" ", ""):
        return False, "type_Ia_photosphere"     # broad SN Fe II photosphere fakes a huge fake ISM column
    fs, bb, be, ln, le = _num(fe_snr), _num(b), _num(b_err), _num(logN), _num(logN_err)
    if fs is None or fs < ADOPT_MIN_FE_SNR:
        return False, "low_fe_snr"
    if be is None or bb is None or be >= bb:
        return False, "b_unconstrained"
    if le is None or le <= 0 or le > ADOPT_MAX_LOGN_ERR:
        return False, "imprecise"
    if ln is not None and ln > ADOPT_MAX_LOGN and sn.upper() not in curated_sns:
        return False, "implausible_highN"
    return True, "ok"

NUV_GRATINGS = ("G230LB", "G230L", "G230M", "G230MB")   # where the Fe II / Mg II forest lives
FUV_GRATINGS = ("G130M", "G160M")                       # RESOLVED COS FUV -> apparent optical depth (G140L R~2000 does not resolve narrow ISM, excluded)
FUV_LINES = os.path.join(ROOT, "linelists", "ism_lines_fuv.csv")
FUV_SUMMARY = os.path.join(ROOT, "catalog", "fuv_aod_summary.csv")

ABS_TAB = os.path.join(ROOT, "reference", "ism_columns.csv")     # curated columns / N(HI) (high-column sightlines live here)
NHI_SUM = os.path.join(ROOT, "catalog", "lya_nhi_summary.csv")   # automated damped-Lya N(HI)
FEH_SOLAR = -4.50           # asplund+2009: log(Fe/H)_sun + 12 = 7.50 -> [Fe/H] = logN(Fe) - logN(HI) - (-4.50)


def _iter_csv_rows(path):
    with open(path) as fh:
        for line in fh:
            if line.strip() and not line.lstrip().startswith("#"):
                yield line


def load_curated_sns():
    # SNe with any curated row in ism_columns.csv -- these are the known high-column sightlines exempt from the
    # logN(FeII) <= 16 plausibility guard (the adopted flag should not silently drop a genuinely dense sightline).
    s = set()
    with open(ABS_TAB) as fh:
        for row in csv.reader(fh):
            if row and row[0].strip() and not row[0].startswith("#"):
                s.add(row[0].strip().upper())
    return s


def load_nhi():
    # per-SN (logN_HI, source): curated ism_columns.csv preferred, else the automated lya_nhi_summary.csv.
    # same single source of truth the emission thread's F_ismcorr uses.
    m = {}
    if os.path.exists(NHI_SUM):
        with open(NHI_SUM) as fh:
            for r in csv.DictReader(fh):
                try:
                    m[r["sn"].upper()] = (float(r["logN_HI"]), "lya_nhi_summary (automated)")
                except (ValueError, KeyError, TypeError):
                    pass
    with open(ABS_TAB) as fh:
        for row in csv.reader(fh):
            if len(row) >= 3 and not row[0].startswith("#") and row[1].strip() == "logN(HI)":
                try:
                    m[row[0].strip().upper()] = (float(row[2]), "ism_columns (curated)")
                except ValueError:
                    pass
    return m


_NHI_CACHE = None
def _nhi_map():
    global _NHI_CACHE
    if _NHI_CACHE is None:
        _NHI_CACHE = load_nhi()
    return _NHI_CACHE


def load_catalog():
    cat = {}
    with open(CAT) as fh:
        for row in csv.DictReader(fh):
            cat[row["name"].upper()] = row
    return cat


def load_lines():
    ism = []
    for row in csv.DictReader(_iter_csv_rows(LINES)):
        ism.append({"ion": row["ion"], "lam": float(row["wavelength_A"]),
                    "f": float(row["f_osc"]), "regime": row["regime"], "notes": row["notes"]})
    return ism


def load_spec(path, z):
    d = np.loadtxt(path, comments="#")
    w, f = d[:, 0], d[:, 1]
    e = d[:, 2] if d.shape[1] > 2 else np.full_like(f, np.nan)
    return w / (1.0 + z), f, e


def load_fuv_lines():
    out = []
    for row in csv.DictReader(_iter_csv_rows(FUV_LINES)):
        out.append({"ion": row["ion"], "lam": float(row["wavelength_A"]),
                    "f": float(row["f_osc"]), "regime": row["regime"], "notes": row["notes"],
                    "blend": "blend" in row["notes"].lower()})       # flag the known blends (S II 1259/Si II 1260, O I 1302/Si II 1304)
    return out


def aod_column(w, f, e, lam0, fval, vwin=150.0, cpad=250.0, sig_det=3.0):
    # apparent optical depth column for one RESOLVED line (savage&sembach 1991). w,f already in host rest frame,
    # so host ISM sits at v~0 (MW foreground is off at -cz). tau(v)=ln(Icont/Iobs); Na(v)=AOD_K*tau/(f lam);
    # N = int Na dv over +-vwin km/s. continuum = linear fit on flanks vwin..vwin+cpad. status: 'saturated' (black
    # core, lower limit >), 'detected' (integrated tau > sig_det sigma), or 'upper' (non-detection -> sig_det-sigma
    # upper limit <). the significance gate is what kills the AOD 1/f weak-line trap: a spurious low-tau blip at a
    # tiny f-value no longer masquerades as a huge detected column.
    v = (w / lam0 - 1.0) * C_KMS
    core = np.abs(v) <= vwin
    flank = (np.abs(v) > vwin) & (np.abs(v) <= vwin + cpad)
    if core.sum() < 4 or flank.sum() < 6:
        return None
    cfit = np.polyfit(v[flank], f[flank], 1)
    cont_rms = float(np.nanstd(f[flank] - np.polyval(cfit, v[flank])))
    vc = v[core]; fc = f[core]; ec = e[core]; contc = np.polyval(cfit, vc)
    if np.nanmedian(contc) <= 0:
        return None
    contc = np.clip(contc, 1e-6 * float(np.nanmax(contc)), None)     # flank fit can extrapolate negative; floor it
    floor = np.maximum(1e-3 * contc, np.where(ec > 0, 0.5 * ec, 1e-3 * contc))
    fsafe = np.clip(fc, floor, None)                            # keep ln finite; the sat flag records the clip
    dv = np.gradient(vc)
    tau = np.log(contc / fsafe)
    tint = float(np.sum(tau * dv))
    stau = np.where(fc > 0, ec / np.clip(fc, ec, None), 0.0)    # per-pixel photon sigma on tau = df/f
    st_tau = float(np.sqrt(np.sum((stau * dv) ** 2)))          # photon sigma on the tau integral
    thi = np.log(np.clip((contc + cont_rms) / fsafe, 1e-6, None))    # continuum-placement +-1 RMS
    tlo = np.log(np.clip((contc - cont_rms) / fsafe, 1e-6, None))
    sc_tau = 0.5 * abs(float(np.sum((thi - tlo) * dv)))
    k = AOD_K / (fval * lam0)
    sig = tint / st_tau if st_tau > 0 else 0.0                 # detection significance (photon only; continuum is systematic)
    # saturation is a SUB-CASE of detection: a genuinely black narrow core (|v|<60, 5px-smoothed min below 2 sigma).
    # a non-detection (low/negative sig = net emission or continuum slope) can never be "saturated".
    narrow = np.abs(vc) <= min(60.0, vwin)
    fsm = median_filter(fc, size=5)
    eloc = float(np.nanmedian(ec[ec > 0])) if np.any(ec > 0) else 0.0
    black_core = bool(narrow.any() and np.min(fsm[narrow]) <= max(0.0, 2.0 * eloc))
    if sig < sig_det:
        status, flag, sat = "upper", "<", False        # non-detection -> report the sig_det-sigma upper bound
        N = k * sig_det * st_tau; sN = np.nan
    elif black_core:
        status, flag, sat = "saturated", ">", True     # detected AND black core -> AOD undercounts, lower limit
        N = k * tint; sN = k * float(np.hypot(st_tau, sc_tau))
    else:
        status, flag, sat = "detected", "", False
        N = k * tint; sN = k * float(np.hypot(st_tau, sc_tau))
    if N <= 0:
        return None
    return {"logN": float(np.log10(N)), "logN_err": float(sN / (N * np.log(10))) if np.isfinite(sN) else np.nan,
            "status": status, "flag": flag, "sig": float(sig), "sat": sat, "tau_int": tint,
            "npix": int(core.sum()), "cont_rms": cont_rms, "v": vc, "f": fc, "cont": contc}


def analyze_fuv(w, f, e, fuv_lines, vwin=150.0):
    # AOD column for every covered FUV ISM line in one resolved COS FUV product (host rest frame already applied).
    ok = np.isfinite(f)
    w, f, e = w[ok], f[ok], e[ok]
    med_e = np.nanmedian(e[np.isfinite(e) & (e > 0)]) if np.any(np.isfinite(e) & (e > 0)) else 0.0
    e = np.where(np.isfinite(e) & (e > 0), e, med_e)
    out = {}
    for L in fuv_lines:
        if not (w.min() + 1 < L["lam"] < w.max() - 1):
            continue
        res = aod_column(w, f, e, L["lam"], L["f"], vwin=vwin)
        if res is not None:
            res["blend"] = L["blend"]
            out[(L["ion"], round(L["lam"], 3))] = res
    return out



def _despike_up(f, size=5, nsig=6.0):
    # clip isolated UPWARD spikes only (hot px / CR); absorption-safe. local to ism (distinct from coadd.despike,
    # which is the positive+DQ-aware product net).
    f = np.asarray(f, float).copy()
    med = median_filter(f, size=size)
    resid = f - med
    mad = np.nanmedian(np.abs(resid[np.isfinite(resid)])) or 0.0
    hot = resid > nsig * 1.4826 * mad
    f[hot] = med[hot]
    return f


def lines_in(ism, w, pad=5.0):
    # return line list entries whose rest wavelength falls at least pad A inside the spectrum coverage
    return [L for L in ism if w.min() + pad < L["lam"] < w.max() - pad]


def group_lines(inr, link=12.0):
    # cluster lines into blended groups: any two lines within `link` A of each other share a group
    ls = sorted(inr, key=lambda L: L["lam"])
    grps = [[ls[0]]]
    for L in ls[1:]:
        if L["lam"] - grps[-1][-1]["lam"] <= link:
            grps[-1].append(L)
        else:
            grps.append([L])
    return grps


def _mg(x, sig, amps, lams):
    # superposition of Gaussians all sharing the same sigma (the instrumental + thermal broadening)
    y = np.zeros_like(x, float)
    for a, mu in zip(amps, lams):
        y += a * np.exp(-0.5 * ((x - mu) / sig) ** 2)
    return y


def fit_group(w, f, group, all_lams, pad=16.0, cont_shift=0.0):
    lams = [L["lam"] for L in group]
    win = (w > min(lams) - pad) & (w < max(lams) + pad)
    if win.sum() < 10:
        return None
    wl, fl = w[win], f[win]
    fk = np.ones_like(wl, bool)
    for lm in all_lams:
        fk &= np.abs(wl - lm) > 4.0
    if fk.sum() < 4:
        return None
    cfit = np.polyval(np.polyfit(wl[fk], fl[fk], 1), wl)
    cont_rms = float(np.nanstd(fl[fk] - cfit[fk]))          # flank scatter = continuum-placement uncertainty
    cont = cfit + cont_shift                                # additive +-1/3-RMS shift (sembach&savage 1992)
    if np.nanmedian(cont) < 1e-20:
        return None
    absn = 1.0 - fl / cont
    n = len(lams)
    p0 = [1.4] + [float(np.clip(np.interp(mu, wl, absn), 0.02, 1.0)) for mu in lams]
    try:
        popt, _ = curve_fit(lambda x, sig, *a: _mg(x, sig, a, lams), wl, absn,
                            p0=p0, bounds=([0.6] + [0.0] * n, [5.0] + [1.6] * n), maxfev=30000)
    except Exception:
        return None
    sig, amps = popt[0], np.array(popt[1:])
    return {"lams": lams, "sig": sig, "amps": amps, "cont": cont, "wl": wl, "win": win, "cont_rms": cont_rms,
            "ews": {round(mu, 3): a * sig * np.sqrt(2 * np.pi) for mu, a in zip(lams, amps)}}


def deblend_ews(w, f, e, inr, n_mc=200, seed=42, link=12.0, return_fits=False):
    all_lams = [L["lam"] for L in inr]
    rng = np.random.default_rng(seed)
    have_e = np.any(np.isfinite(e)) and np.any(np.asarray(e) > 0)
    out = {}
    cerr = {}                                    # per-line continuum-placement EW error (sembach&savage)
    fits = {}
    for gi, g in enumerate(group_lines(inr, link)):
        fit = fit_group(w, f, g, all_lams)
        fits[gi] = fit
        if fit is None:
            for L in g:
                out[round(L["lam"], 3)] = (np.nan, np.nan, np.nan); cerr[round(L["lam"], 3)] = np.nan
            continue
        if have_e:
            en = np.where(np.isfinite(e) & (e > 0), e, 0.0)
        else:
            wl, cont, win = fit["wl"], fit["cont"], fit["win"]
            fk = np.ones_like(wl, bool)
            for lm in all_lams:
                fk &= np.abs(wl - lm) > 4.0
            en = np.full_like(f, np.nanstd((f[win] - cont)[fk]) or np.nanmedian(np.abs(f[win])))
        rms = fit.get("cont_rms", 0.0)           # re-fit with the flank continuum shifted +-1/3 RMS
        fhi = fit_group(w, f, g, all_lams, cont_shift=+rms / 3.0) if rms > 0 else None
        flo = fit_group(w, f, g, all_lams, cont_shift=-rms / 3.0) if rms > 0 else None
        draws = {round(L["lam"], 3): [] for L in g}
        for _ in range(n_mc):
            fm = fit_group(w, f + rng.normal(0.0, en), g, all_lams)
            if fm is None:
                continue
            for lm, val in fm["ews"].items():
                if lm in draws:
                    draws[lm].append(val)
        for L in g:
            k = round(L["lam"], 3)
            arr = np.array(draws[k])
            if len(arr) >= 20:
                med = np.median(arr)
                out[k] = (fit["ews"][k], med - np.percentile(arr, 16), np.percentile(arr, 84) - med)
            else:
                out[k] = (fit["ews"].get(k, np.nan), np.nan, np.nan)
            vals = [fit["ews"].get(k)] + ([fhi["ews"].get(k)] if fhi else []) + ([flo["ews"].get(k)] if flo else [])
            vals = [v for v in vals if v is not None and np.isfinite(v)]
            cerr[k] = 0.5 * (max(vals) - min(vals)) if len(vals) > 1 else 0.0
    if return_fits:
        return out, cerr, fits
    return out, cerr


_ltau = np.linspace(-4.0, 6.0, 600)
_Fint = np.array([quad(lambda x, t=10.0 ** lt: 1.0 - np.exp(-t * np.exp(-x * x)), 0, 20)[0] for lt in _ltau])


def cog_F(tau0):
    # interpolate the pre-tabulated CoG integral F(tau0) = int_0^inf [1 - exp(-tau0 * exp(-u^2))] du
    return np.interp(np.log10(np.clip(tau0, 1e-4, 1e6)), _ltau, _Fint)


def red_ew(N, f, lam, b):
    # reduced equivalent width W/lam = (2b/c) * F(tau0);  tau0 = TAU_K * N * f * lam / b
    return (2.0 * b / C_KMS) * cog_F(TAU_K * N * f * lam / b)


def collect(inr, deb, ion=None, snr_min=2.0, cerr=None):
    # detection gate + CoG points. default gates on the MC photon EW error; if cerr (per-line continuum-placement EW
    # error) is given, gate AND weight on the TOTAL error hypot(photon, continuum) at snr_min sigma. that way a line
    # whose EW is mostly continuum-placement artifact no longer counts
    # as a detection; the Fe II anchor is essentially unchanged, only the marginal weak-line detections are cleaned).
    lam, fo, y, ye = [], [], [], []
    for L in inr:
        if ion and L["ion"] != ion:
            continue
        ew, lo, hi = deb[round(L["lam"], 3)]
        if not np.isfinite(ew) or ew <= 0.02 or not (hi > 0):
            continue
        err = 0.5 * (lo + hi)
        if cerr is not None:
            err = float(np.hypot(err, cerr.get(round(L["lam"], 3), 0.0) or 0.0))    # total EW error
        if err > 0 and ew / err < snr_min:      # drop non-detections (2-sigma EW gate)
            continue
        lam.append(L["lam"]); fo.append(L["f"]); y.append(ew / L["lam"]); ye.append(err / L["lam"])
    return map(np.array, (lam, fo, y, ye))


def fit_b_N(lam, fo, y, ye, b0=60.0, logN0=14.5):
    # joint Fe II fit: minimize log-EW residuals in log space (compressed dynamic range across the CoG branches)
    def resid(p):
        pred = red_ew(10.0 ** p[1], fo, lam, p[0])
        return (np.log10(y) - np.log10(np.clip(pred, 1e-30, None))) / (ye / (y * np.log(10)) + 1e-3)
    return least_squares(resid, [b0, logN0], bounds=([5, 10], [300, 20]), max_nfev=5000).x


def fit_N(lam, fo, y, ye, b, logN0=14.0):
    # logN with b fixed; for all non-Fe ions (Fe II sets b, every other ion inherits it)
    def resid(p):
        return np.log10(y) - np.log10(np.clip(red_ew(10.0 ** p[0], fo, lam, b), 1e-30, None))
    return least_squares(resid, [logN0], bounds=([9], [21]), max_nfev=3000).x[0]


def _total_ew_err(inr, deb, cerr, ion):
    # per-detected-line TOTAL EW error = quadrature(photon MC, continuum placement), returned in W/lam units.
    # collect(cerr=...) already gates + returns the total error, so this is just the option-A collect for `ion`.
    return collect(inr, deb, ion=ion, cerr=cerr)


def _bN_mc(inr, deb, cerr, n=15, seed=1):
    # b + logN(FeII) uncertainty: MC over the Fe II points perturbed by their TOTAL EW errors, refit each draw
    lam, fo, y, ye = _total_ew_err(inr, deb, cerr, "Fe II")
    if len(lam) < 3:
        return np.nan, np.nan
    rng = np.random.default_rng(seed); bs, Ns = [], []
    for _ in range(n):
        yp = np.clip(y + rng.normal(0, ye), 1e-12, None)
        bb, NN = fit_b_N(lam, fo, yp, ye)
        bs.append(bb); Ns.append(NN)
    return float(np.std(bs)), float(np.std(Ns))


def _N_mc(inr, deb, cerr, ion, b, n=15, seed=2):
    lam, fo, y, ye = _total_ew_err(inr, deb, cerr, ion)
    if len(lam) == 0:
        return np.nan
    rng = np.random.default_rng(seed); Ns = []
    for _ in range(n):
        yp = np.clip(y + rng.normal(0, ye), 1e-12, None)
        Ns.append(fit_N(lam, fo, yp, ye, b))
    return float(np.std(Ns))


def analyze(w, f, e, ism, n_mc=200, clean_spikes=True):
    ok = np.isfinite(f) & (f > 0)
    w, f, e = w[ok], f[ok], e[ok]
    if clean_spikes:
        f = _despike_up(f)
    inr = lines_in(ism, w)
    if len(inr) < 4:
        return None
    deb, cerr, fits = deblend_ews(w, f, e, inr, n_mc=n_mc, return_fits=True)
    fl, ff, fy, fye = collect(inr, deb, ion="Fe II", cerr=cerr)    # option A: total-error 2-sigma gate
    if len(fl) < 3:
        return None
    b, logN_fe = fit_b_N(fl, ff, fy, fye)
    # physical plausibility gate: ISM b ~ 10-150 km/s, Fe II logN ~ 13-17
    # values outside this range mean the fit is chasing noise/CSM emission
    if not (8.0 < b < 250.0) or not (13.0 < logN_fe < 17.5):
        return None
    fe_snr = float(np.median(fy / fye)) if len(fy) else np.nan   # per-epoch quality: median Fe II detection SNR
    # CoG anchor: at least one detected Fe II line on the linear/transition part (tau0 < ANCHOR_TAU) so N is
    # pinned. if every Fe II line is saturated (all tau0 high) b sets the flat height but N floats up -> that is
    # the SN iron-photosphere contamination signature (deep photospheric troughs read as a huge fake ISM column).
    anchored = bool(np.min(TAU_K * 10.0 ** logN_fe * ff * fl / b) < ANCHOR_TAU)
    b_err, logN_fe_err = _bN_mc(inr, deb, cerr)      # b + logN(FeII) uncertainty (photon MC + continuum placement)
    Ncol = {}
    for ion in sorted(set(L["ion"] for L in inr)):
        lam, fo, y, ye = collect(inr, deb, ion=ion, cerr=cerr)     # option A: total-error gate for every ion
        if len(lam) == 0:
            continue
        ln = logN_fe if ion == "Fe II" else fit_N(lam, fo, y, ye, b)
        ln_err = logN_fe_err if ion == "Fe II" else _N_mc(inr, deb, cerr, ion, b)
        # lower-limit flag: N only a bound if even the weakest detected line is saturated (tau0 > 5)
        is_limit = bool(np.min(TAU_K * 10 ** ln * fo * lam / b) > 5.0)
        Ncol[ion] = (ln, len(lam), is_limit, ln_err)
    return {"inr": inr, "deb": deb, "cerr": cerr, "fits": fits, "w_used": w, "f_used": f,
            "b": b, "b_err": b_err, "logN_fe": logN_fe, "logN_fe_err": logN_fe_err,
            "Ncol": Ncol, "n_fe": len(fl), "fe_snr": fe_snr, "anchored": anchored}


def _write_diag(sn, grating, ph, w, f, inr, deb, b, logN_fe, Ncol, outdir, fits=None):
    # per-epoch absorption scrutiny plots: continuum panels + CoG.
    # w, f are the SAME despiked arrays that were fit; fits are the actual per-group fits (reused, not re-fit).
    all_lams = [L["lam"] for L in inr]
    grps = group_lines(inr)

    # 1. per-group continuum-subtraction panels
    n_grp = len(grps)
    ncols = min(n_grp, 4); nrows = (n_grp + ncols - 1) // ncols
    fig, axes = plt.subplots(nrows, ncols, figsize=(4 * ncols, 3 * nrows), squeeze=False)
    for gi, g in enumerate(grps):
        ax = axes[gi // ncols][gi % ncols]
        fit = fits.get(gi) if fits is not None else fit_group(w, f, g, all_lams)
        lams = [L["lam"] for L in g]
        lo, hi = min(lams) - 16, max(lams) + 16
        sel = (w > lo) & (w < hi)
        ax.plot(w[sel], f[sel] * 1e15, color="0.35", lw=0.9)
        if fit is not None:
            ax.plot(fit["wl"], fit["cont"] * 1e15, color="crimson", lw=1.2, ls="--")
        for L in g:
            k = round(L["lam"], 3)
            ew, elo, ehi = deb.get(k, (np.nan, np.nan, np.nan))
            err = 0.5 * (elo + ehi) if np.isfinite(elo) else np.nan
            det = np.isfinite(ew) and err > 0 and ew / err >= 2.0     # 2-sigma detection
            col = "steelblue" if det else "0.6"
            ax.axvline(L["lam"], color=col, lw=0.8, ls=":")
            lab = f"{L['ion']} {L['lam']:.0f}"
            if np.isfinite(ew):
                lab += f"\nEW={ew*1e3:.0f}+-{err*1e3:.0f} mA" + ("" if det else " n.d.")
            ax.text(L["lam"], ax.get_ylim()[1] if ax.get_ylim()[1] != 1.0 else 1.0,
                    lab, fontsize=6, ha="center", va="bottom", rotation=90, color=col)
        ax.set_xlabel("rest wvl (A)", fontsize=7); ax.set_ylabel("flux (1e-15)", fontsize=7)
        ax.tick_params(labelsize=6)
    for gi in range(n_grp, nrows * ncols):
        axes[gi // ncols][gi % ncols].set_visible(False)
    fig.suptitle(f"{sn}  {grating}  day{ph:.0f}  - continuum fits per line group", fontsize=8)
    fig.tight_layout()
    fig.savefig(os.path.join(outdir, f"{sn}_{grating}_day{ph:.0f}_cont.png"), dpi=110)
    plt.close(fig)

    # 2. curve-of-growth: observed W/lam vs the fitted curve, colored by ion
    lam_fe, fo_fe, y_fe, ye_fe = collect(inr, deb, ion="Fe II")
    if len(lam_fe) < 2:
        return
    tau0_fe = TAU_K * 10 ** logN_fe * fo_fe * lam_fe / b
    cog_x = np.logspace(-3, 4, 400)
    cog_y = np.array([cog_F(t) for t in cog_x])
    fig, ax = plt.subplots(figsize=(5.5, 4.2))
    ax.loglog(cog_x, 2 * b / C_KMS * np.array([cog_F(t) for t in cog_x]), color="0.6", lw=1.3, label="theory")
    ion_colors = {"Fe II": "steelblue", "Mg II": "darkorange", "Cr II": "green",
                  "Zn II": "purple", "Mn II": "crimson"}
    for ion in sorted(set(L["ion"] for L in inr)):
        lm, fo, y, ye = collect(inr, deb, ion=ion)
        if len(lm) == 0:
            continue
        is_lim = bool(Ncol[ion][2]) if (ion in Ncol and len(Ncol[ion]) > 2) else False
        ax.errorbar(TAU_K * 10 ** (Ncol[ion][0] if ion in Ncol else logN_fe) * fo * lm / b,
                    y, yerr=[ye, ye], fmt=">" if is_lim else "o", ms=5, capsize=3,
                    color=ion_colors.get(ion, "gray"), label=ion + (" (lim)" if is_lim else ""))
    ax.set_xlabel(r"optical depth $\tau_0$"); ax.set_ylabel(r"$W/\lambda$")
    ax.grid(True, which="major", alpha=0.45, lw=0.7)
    ax.grid(True, which="minor", alpha=0.2, lw=0.4)
    ax.set_title(f"{sn}  {grating}  day{ph:.0f}  b={b:.0f}  logN(FeII)={logN_fe:.2f}", fontsize=8)
    ax.legend(fontsize=7, loc="upper left")
    fig.tight_layout()
    fig.savefig(os.path.join(outdir, f"{sn}_{grating}_day{ph:.0f}_cog.png"), dpi=110)
    plt.close(fig)


def phase_of(path, cat=None):
    m = re.search(r"day([0-9.]+)", path)
    if m:
        return float(m.group(1))
    # date-only filenames: compute from catalog discovery MJD
    m_date = re.search(r"(\d{4}-\d{2}-\d{2})", path)
    if m_date and cat:
        sn = os.path.basename(path).split("_")[0].upper()
        row = cat.get(sn, {})
        mjd_str = str(row.get("tns_disc_mjd", "")).strip()
        if mjd_str and mjd_str not in ("", "nan", "None"):
            obs = datetime.date.fromisoformat(m_date.group(1))
            disc = datetime.date(1858, 11, 17) + datetime.timedelta(days=float(mjd_str))
            return float((obs - disc).days)
    return np.nan


def _write_ism_csv(sn, g, ph, r, ismdir, cat=None, curated_sns=None):
    # per-ion columns (+ logN error) and a rich per-line table: EW, the photon/continuum/total error budget,
    # detection flag, optical depth, saturation. this is the B4 rich output (was queued, now wired).
    # the per-ion table also carries adopted_eligible/exclude_reason so a reader never mistakes a raw
    # photospheric Fe II column for a foreground ISM measurement (only the eligible rows feed the adopted value).
    cat = cat if cat is not None else load_catalog()
    curated_sns = curated_sns if curated_sns is not None else load_curated_sns()
    elig, reason = adopt_eligible(sn, bool(r["anchored"]), r["fe_snr"], r["b"], r["b_err"],
                                  r["logN_fe"], r["logN_fe_err"], cat, curated_sns)
    with open(os.path.join(ismdir, f"{sn}_{g}_day{ph:.0f}_cog.csv"), "w", newline="") as fh:
        wri = csv.writer(fh); wri.writerow(["ion", "logN", "logN_err", "n_lines", "limit", "adopted_eligible", "exclude_reason"])
        for ion, (ln, nl, lim, lnerr) in r["Ncol"].items():
            # eligibility is a Fe II-anchor property; carry it on the Fe II row (blank for the other ions).
            ec = ("yes" if elig else "no") if ion == "Fe II" else ""
            rc = ("" if elig else reason) if ion == "Fe II" else ""
            wri.writerow([ion, f"{ln:.3f}", f"{lnerr:.3f}" if np.isfinite(lnerr) else "", nl, ">" if lim else "", ec, rc])
        hi = _nhi_map().get(sn.upper())          # sightline N(HI) from damped Lya (same every epoch); the H the cog output was always supposed to carry
        if hi:
            wri.writerow(["H I", f"{hi[0]:.3f}", "", "", "", "", ""])
    deb, cerr, b = r["deb"], r["cerr"], r["b"]
    ioncol = {ion: v[0] for ion, v in r["Ncol"].items()}
    with open(os.path.join(ismdir, f"{sn}_{g}_day{ph:.0f}_lines.csv"), "w", newline="") as fh:
        wri = csv.writer(fh)
        wri.writerow(["ion", "lam", "f", "EW_A", "EW_err_mc", "EW_err_cont", "EW_err_tot", "detected", "tau0", "saturated"])
        for L in sorted(r["inr"], key=lambda L: L["lam"]):
            k = round(L["lam"], 3)
            ew, lo, hi = deb.get(k, (np.nan, np.nan, np.nan))
            mc = 0.5 * (lo + hi) if np.isfinite(lo) else np.nan
            ce = cerr.get(k, np.nan)
            tot = np.hypot(mc, ce) if (np.isfinite(mc) and np.isfinite(ce)) else (mc if np.isfinite(mc) else ce)
            det = bool(np.isfinite(ew) and np.isfinite(tot) and tot > 0 and ew / tot >= 2.0)   # option A: total-error 2-sigma gate
            ln = ioncol.get(L["ion"], r["logN_fe"])
            tau0 = TAU_K * 10 ** ln * L["f"] * L["lam"] / b
            wri.writerow([L["ion"], f"{L['lam']:.3f}", f"{L['f']:.4f}",
                          f"{ew:.4f}" if np.isfinite(ew) else "", f"{mc:.4f}" if np.isfinite(mc) else "",
                          f"{ce:.4f}" if np.isfinite(ce) else "", f"{tot:.4f}" if np.isfinite(tot) else "",
                          "yes" if det else "no", f"{tau0:.2f}", "yes" if tau0 > 5 else "no"])


def finalize_and_write(summ, cat):
    # decide the adopted foreground column per SN and write the master summary. the adopted flag
    # now means what a user thinks it means, with four guards on top of the anchor: (1) precision -- logN(FeII)_err
    # <= 0.3 dex; (2) plausibility -- logN(FeII) <= 16 unless the SN is a curated high-column sightline
    # (ism_columns.csv); (3) never a Type Ia (the broad SN Fe II photosphere fakes an enormous ISM column);
    # (4) fe_snr >= 3, b_err < b, logN_err > 0. when >= 2 anchored gated epochs survive, the column CANNOT change
    # between epochs (the gas is light-years away), so we report the MEAN over those epochs with the epoch scatter as
    # the honest error (the per-epoch photon error understates the truth ~10x); a single surviving epoch keeps its own
    # value/error. N(HI) + [Fe/H]_gas are joined from the single N(HI) source of truth (curated ism_columns.csv, else
    # automated lya_nhi_summary.csv). this runs on the per-epoch rows only, so it can be re-applied without the
    # (slow) CoG rerun via `--finalize`.
    curated_sns = load_curated_sns()
    nhi_map = load_nhi()
    new_cols = ["logN_FeII_adopted", "logN_FeII_adopted_err", "b_adopted", "n_epochs_adopted",
                "logN_HI", "logN_HI_src", "FeH_gas"]
    gated = {}
    for row in summ:
        sn = row["sn"]
        elig, _ = adopt_eligible(sn, row["anchored"] == "yes", row["fe_snr"], row["b"],
                                 row["b_err"], row["logN_FeII"], row["logN_FeII_err"], cat, curated_sns)
        if elig:
            gated.setdefault(sn, []).append(row)
    for row in summ:                                # default the new columns blank
        row["adopted"] = ""
        for c in new_cols:
            row[c] = ""
    for sn, rows in gated.items():
        rep = max(rows, key=lambda r: r["fe_snr"])
        logNs = [r["logN_FeII"] for r in rows]; bs = [r["b"] for r in rows]; n = len(rows)
        if n >= 2:
            mean = float(np.mean(logNs)); bmean = float(np.mean(bs))
            scat = float(np.std(logNs, ddof=1)) if n >= 3 else abs(logNs[0] - logNs[1]) / 2.0
            rep["adopted"] = f"mean({n})"
        else:
            mean, bmean, scat = logNs[0], bs[0], float(rep["logN_FeII_err"]); rep["adopted"] = "yes"
        rep["logN_FeII_adopted"] = round(mean, 3)
        rep["logN_FeII_adopted_err"] = round(scat, 3)
        rep["b_adopted"] = round(bmean, 1)
        rep["n_epochs_adopted"] = n
        hi = nhi_map.get(sn.upper())
        if hi:
            rep["logN_HI"] = round(hi[0], 3); rep["logN_HI_src"] = hi[1]
            rep["FeH_gas"] = round(mean - hi[0] - FEH_SOLAR, 2)   # gas-phase [Fe/H]; iron is depleted onto dust so this is a LOWER limit on Z
    with open(ISM_SUMMARY, "w", newline="") as fh:
        wri = csv.DictWriter(fh, restval="", fieldnames=["sn", "grating", "phase", "z", "b", "b_err",
            "logN_FeII", "logN_FeII_err", "n_fe", "fe_snr", "anchored", "adopted"] + new_cols)
        wri.writeheader()
        for row in summ:
            wri.writerow(row)


def _load_summ_csv():
    # re-load the per-epoch rows from an existing ism_cog_summary.csv (coerce the numeric fields back) so the adopted
    # finalization can be re-applied without the slow CoG rerun (the per-epoch numbers do not change).
    rows = []
    with open(ISM_SUMMARY) as fh:
        for r in csv.DictReader(fh):
            for k in ("z", "b", "logN_FeII", "fe_snr", "phase"):
                r[k] = float(r[k])
            r["n_fe"] = int(r["n_fe"])
            for k in ("b_err", "logN_FeII_err"):
                r[k] = float(r[k]) if r.get(k) not in ("", None) else ""
            rows.append({k: r[k] for k in ("sn", "grating", "phase", "z", "b", "b_err", "logN_FeII",
                                            "logN_FeII_err", "n_fe", "fe_snr", "anchored")})
    return rows


def run_catalog(n_mc=150, min_fe=3):
    # loop over every NUV per-grating product, write per-epoch cog csv + a master summary
    cat = load_catalog()
    ism = load_lines()
    curated_sns = load_curated_sns()
    summ = []
    sne = sorted(d for d in os.listdir(OUT) if os.path.isdir(os.path.join(OUT, d)) and d.upper() in cat)
    for sn in sne:
        z = float(cat[sn.upper()]["z"])
        ismdir = os.path.join(OUT, sn, "absorption")
        for old in glob.glob(os.path.join(ismdir, f"{sn}_*_day*_cog.csv")) + \
                   glob.glob(os.path.join(ismdir, f"{sn}_*_day*_lines.csv")) + \
                   glob.glob(os.path.join(ismdir, f"{sn}_*_day*_cog.png")) + \
                   glob.glob(os.path.join(ismdir, f"{sn}_*_day*_cont.png")):
            os.remove(old)                       # drop stale-phase orphans so only current epochs ship
        prods = []
        for g in NUV_GRATINGS:
            prods += glob.glob(f"{OUT}/{sn}/**/{g}/{sn}_*_{g}_native.txt", recursive=True)
        for p in sorted(set(prods)):
            g = next((g for g in NUV_GRATINGS if f"_{g}_native" in p), "?")
            try:
                r = analyze(*load_spec(p, z), ism, n_mc=n_mc)
            except Exception as ex:
                print(f"  {sn} {os.path.basename(p)}: {ex}")
                continue
            if r is None:
                continue
            ph = phase_of(p, cat)
            ismdir = os.path.join(OUT, sn, "absorption")
            os.makedirs(ismdir, exist_ok=True)
            _write_ism_csv(sn, g, ph, r, ismdir, cat, curated_sns)
            try:
                _write_diag(sn, g, ph, r["w_used"], r["f_used"], r["inr"], r["deb"], r["b"], r["logN_fe"], r["Ncol"], ismdir, r["fits"])
            except Exception as ex:
                print(f"  WARNING diag {sn} {g} day{ph:.0f}: {ex}")
            summ.append({"sn": sn, "grating": g, "phase": ph, "z": z, "b": round(r["b"], 1),
                         "b_err": round(r["b_err"], 1) if np.isfinite(r["b_err"]) else "",
                         "logN_FeII": round(r["logN_fe"], 3),
                         "logN_FeII_err": round(r["logN_fe_err"], 3) if np.isfinite(r["logN_fe_err"]) else "",
                         "n_fe": r["n_fe"], "fe_snr": round(r["fe_snr"], 1), "anchored": "yes" if r["anchored"] else "no"})
            print(f"  {sn:14} {g:7} day{ph:6.1f}  b={r['b']:5.1f}  logN(FeII)={r['logN_fe']:.2f}  feSNR={r['fe_snr']:.1f}  anchored={r['anchored']}")
    finalize_and_write(summ, cat)
    print(f"\n{len(summ)} epochs measured -> {ISM_SUMMARY}")
    return summ


def _nuv_feii_lookup():
    # per-SN NUV Fe II column for the same-ion cross-check: prefer the adopted (mean) value, else the best single
    # epoch (highest fe_snr). this is what the FUV AOD Fe II 1608 gets compared against.
    m = {}; best = {}
    if not os.path.exists(ISM_SUMMARY):
        return m
    with open(ISM_SUMMARY) as fh:
        for r in csv.DictReader(fh):
            sn = r["sn"].upper()
            if r.get("logN_FeII_adopted") not in ("", None):
                le = r.get("logN_FeII_adopted_err")
                m[sn] = (float(r["logN_FeII_adopted"]), float(le) if le not in ("", None) else np.nan, r.get("adopted") or "adopted")
            else:
                try:
                    fs = float(r["fe_snr"])
                except (ValueError, KeyError, TypeError):
                    continue
                if sn not in best or fs > best[sn][0]:
                    le = r.get("logN_FeII_err")
                    best[sn] = (fs, float(r["logN_FeII"]), float(le) if le not in ("", None) else np.nan)
    for sn, (fs, ln, le) in best.items():
        m.setdefault(sn, (ln, le, "single-epoch NUV"))
    return m


def _fuv_diag_plot(sn, g, ph, res, outdir):
    # scrutiny plot behind the FUV AOD columns: normalized flux vs velocity for every measured line,
    # with the +-vwin integration window shaded so the continuum + core choice is visible.
    items = sorted(res.items(), key=lambda kv: kv[0][1])
    n = len(items)
    if n == 0:
        return None
    ncol = min(3, n); nrow = int(np.ceil(n / ncol))
    fig, axes = plt.subplots(nrow, ncol, figsize=(3.6*ncol, 2.4*nrow), squeeze=False)
    for ax, ((ion, lam), d) in zip(axes.ravel(), items):
        v, fc, cont = d["v"], d["f"], d["cont"]
        norm = fc / np.where(cont > 0, cont, np.nan)
        ax.axhline(1.0, color='0.6', ls=':', lw=0.8); ax.axvline(0.0, color='0.6', ls=':', lw=0.8)
        ax.plot(v, norm, 'k', lw=0.9, drawstyle='steps-mid')
        ax.axvspan(-150, 150, color='gold', alpha=0.08)
        tag = d["flag"] or '='
        ax.set_title(f"{ion} {lam:.0f}  {tag}logN={d['logN']:.2f} ({d['status']}, {d['sig']:.0f}$\\sigma$)", fontsize=7)
        ax.set_xlim(-400, 400); ax.set_ylim(-0.1, 1.6)
        ax.set_xlabel(r'$v$ [km s$^{-1}$]', fontsize=7); ax.set_ylabel(r'$I/I_0$', fontsize=7)
        ax.tick_params(labelsize=6)
    for ax in axes.ravel()[n:]:
        ax.set_visible(False)
    fig.suptitle(f'{sn} {g} d{ph:.0f} FUV AOD', fontsize=9)
    os.makedirs(outdir, exist_ok=True)
    out = os.path.join(outdir, f'{sn}_{g}_day{ph:.0f}_aod.png')
    fig.tight_layout(); fig.savefig(out, dpi=110); plt.close(fig)
    return out


def run_fuv(vwin=150.0):
    # loop every resolved COS FUV product, measure the FUV ISM lines by AOD, and cross-check Fe II 1608 (AOD)
    # against the NUV Fe II CoG for the same sightline (same ion, independent method + spectrograph).
    cat = load_catalog()
    fuv_lines = load_fuv_lines()
    nuv_fe = _nuv_feii_lookup()
    sne = sorted(d for d in os.listdir(OUT) if os.path.isdir(os.path.join(OUT, d)) and d.upper() in cat)
    rows = []
    for sn in sne:
        z = float(cat[sn.upper()]["z"])
        prods = []
        for g in FUV_GRATINGS:
            prods += glob.glob(f"{OUT}/{sn}/**/{g}/{sn}_*_{g}_native.txt", recursive=True)
        for p in sorted(set(prods)):
            g = next((g for g in FUV_GRATINGS if f"_{g}_native" in p), "?")
            try:
                res = analyze_fuv(*load_spec(p, z), fuv_lines, vwin=vwin)
            except Exception as ex:
                print(f"  {sn} {os.path.basename(p)}: {ex}")
                continue
            if not res:
                continue
            ph = phase_of(p, cat)
            ismdir = os.path.join(OUT, sn, "absorption")
            os.makedirs(ismdir, exist_ok=True)
            with open(os.path.join(ismdir, f"{sn}_{g}_day{ph:.0f}_aod.csv"), "w", newline="") as fh:
                wri = csv.writer(fh); wri.writerow(["ion", "lam", "logN", "logN_err", "status", "flag", "sig", "blend", "tau_int", "npix"])
                for (ion, lam), d in sorted(res.items(), key=lambda kv: kv[0][1]):
                    wri.writerow([ion, f"{lam:.3f}", f"{d['logN']:.3f}",
                                  f"{d['logN_err']:.3f}" if np.isfinite(d["logN_err"]) else "", d["status"], d["flag"],
                                  f"{d['sig']:.1f}", "yes" if d["blend"] else "", f"{d['tau_int']:.2f}", d["npix"]])
            png = _fuv_diag_plot(sn, g, ph, res, ismdir)
            # per-SN JSON so the FUV AOD columns are surfaced next to the emission/absorption products
            jlines = [{"ion": ion, "lam": round(lam, 3), "logN_aod": round(d["logN"], 3),
                       "logN_aod_err": round(d["logN_err"], 3) if np.isfinite(d["logN_err"]) else None,
                       "status": d["status"], "flag": d["flag"], "sig": round(d["sig"], 1),
                       "blend": bool(d["blend"])} for (ion, lam), d in sorted(res.items(), key=lambda kv: kv[0][1])]
            jrec = {"sn": sn.upper(), "grating": g, "phase": round(ph, 1),
                    "generated": datetime.date.today().isoformat(), "method": "apparent optical depth (Savage&Sembach 1991), host rest frame",
                    "lines": jlines}
            if png:
                jrec["diag_plot"] = os.path.relpath(png, OUT).replace("\\", "/")
            with open(os.path.join(OUT, sn, f"{sn}_fuv_aod.json"), "w") as fh:
                json.dump(jrec, fh, indent=2)
            for (ion, lam), d in sorted(res.items(), key=lambda kv: kv[0][1]):
                rows.append({"sn": sn, "grating": g, "phase": round(ph, 1), "ion": ion, "lam": round(lam, 3),
                             "logN_aod": round(d["logN"], 3),
                             "logN_aod_err": round(d["logN_err"], 3) if np.isfinite(d["logN_err"]) else "",
                             "status": d["status"], "flag": d["flag"], "sig": round(d["sig"], 1),
                             "blend": "yes" if d["blend"] else "", "npix": d["npix"]})
            fe = res.get(("Fe II", 1608.451))
            if fe is not None:
                nu = nuv_fe.get(sn.upper())
                tag = (f" | NUV(FeII)={nu[0]:.2f}+-{nu[1]:.2f} [{nu[2]}] d={fe['logN']-nu[0]:+.2f}") if nu else " | NUV=n/a"
                print(f"  {sn:14} {g:6} day{ph:6.1f}  FeII1608 {fe['flag'] or '='}{fe['logN']:5.2f} ({fe['status']}, sig={fe['sig']:.1f}){tag}")
            else:
                det = sum(1 for d in res.values() if d["status"] == "detected")
                print(f"  {sn:14} {g:6} day{ph:6.1f}  {len(res)} FUV lines ({det} detected), no Fe II 1608")
    with open(FUV_SUMMARY, "w", newline="") as fh:
        wri = csv.DictWriter(fh, fieldnames=["sn", "grating", "phase", "ion", "lam", "logN_aod", "logN_aod_err",
                                             "status", "flag", "sig", "blend", "npix"])
        wri.writeheader()
        for r in rows:
            wri.writerow(r)
    ndet = sum(1 for r in rows if r["status"] == "detected")
    print(f"\n{len(rows)} FUV AOD line-columns ({ndet} detected, rest upper/lower limits) -> {FUV_SUMMARY}")
    return rows


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--nmc", type=int, default=150)
    ap.add_argument("--min-fe", type=int, default=3)
    ap.add_argument("--sn", default=None, help="run a single SN (uppercase dir name)")
    ap.add_argument("--finalize", action="store_true", help="re-derive the adopted flag + N(HI) join from the existing summary CSV (no CoG rerun)")
    ap.add_argument("--fuv", action="store_true", help="run the resolved COS FUV AOD pass + Fe II 1608 cross-check vs the NUV CoG")
    a = ap.parse_args()
    if a.finalize:
        finalize_and_write(_load_summ_csv(), load_catalog())
        print(f"re-finalized -> {ISM_SUMMARY}")
    elif a.fuv:
        run_fuv()
    elif a.sn:
        cat = load_catalog(); ism = load_lines(); z = float(cat[a.sn.upper()]["z"]); curated_sns = load_curated_sns()
        ismdir = os.path.join(OUT, a.sn, "absorption")
        for old in glob.glob(os.path.join(ismdir, f"{a.sn}_*_day*_cog.csv")) + \
                   glob.glob(os.path.join(ismdir, f"{a.sn}_*_day*_lines.csv")) + \
                   glob.glob(os.path.join(ismdir, f"{a.sn}_*_day*_cog.png")) + \
                   glob.glob(os.path.join(ismdir, f"{a.sn}_*_day*_cont.png")):
            os.remove(old)                       # drop stale-phase orphans so only current epochs ship
        for g in NUV_GRATINGS:
            for p in sorted(glob.glob(f"{OUT}/{a.sn}/**/{g}/{a.sn}_*_{g}_native.txt", recursive=True)):
                r = analyze(*load_spec(p, z), ism, n_mc=a.nmc)
                if r:
                    ph = phase_of(p, cat)
                    print(f"{a.sn} {g} day{ph:.0f}: b={r['b']:.1f} logN(FeII)={r['logN_fe']:.2f} feSNR={r['fe_snr']:.1f}")
                    ismdir = os.path.join(OUT, a.sn, "absorption")
                    os.makedirs(ismdir, exist_ok=True)
                    _write_ism_csv(a.sn, g, ph, r, ismdir, cat, curated_sns)
                    # write diagnostic plots (reuse the actual despiked flux + fits)
                    try:
                        _write_diag(a.sn, g, ph, r["w_used"], r["f_used"], r["inr"], r["deb"], r["b"], r["logN_fe"], r["Ncol"], ismdir, r["fits"])
                    except Exception as ex:
                        print(f"  WARNING diag {g} day{ph:.0f}: {ex}")
    else:
        run_catalog(n_mc=a.nmc, min_fe=a.min_fe)
