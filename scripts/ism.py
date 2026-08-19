import os, glob, csv, re, argparse, datetime
import numpy as np
from scipy.optimize import curve_fit, least_squares
from scipy.integrate import quad
from scipy.ndimage import median_filter
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# ism equivalent-width -> curve-of-growth -> column density, ported from ism_ew_cog_sandbox.ipynb.
# reads the observed-frame products in output5, applies z from the catalog, measures deblended EWs,
# anchors the doppler b on the Fe II series, reads off per-ion column densities. see docs/ISM_work1.md.

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
from paths import OUT, CATALOG as CAT, ISM_SUMMARY
LINES = os.path.join(ROOT, "linelists", "ism_lines.csv")

C_KMS = 2.99792458e5
TAU_K = 1.4973e-15          # tau0 = TAU_K * N * f * lam / b
LIN_K = 8.853e-21           # linear CoG W/lam = LIN_K * N f lam
AOD_K = 3.7679e14           # N = AOD_K/(f lam) * int tau dv (Savage&Sembach 1991; was 1.13e17 = 300x too high)
ANCHOR_TAU = 3.0            # CoG anchor: need >=1 detected Fe II line below this tau0 (linear/transition part) or N floats
NUV_GRATINGS = ("G230LB", "G230L", "G230M", "G230MB")   # where the Fe II / Mg II forest lives


def load_catalog():
    cat = {}
    with open(CAT) as fh:
        for row in csv.DictReader(fh):
            cat[row["name"].upper()] = row
    return cat


def load_lines():
    ism = []
    with open(LINES) as fh:
        for row in csv.DictReader(r for r in fh if not r.startswith("#")):
            ism.append({"ion": row["ion"], "lam": float(row["wavelength_A"]),
                        "f": float(row["f_osc"]), "regime": row["regime"], "notes": row["notes"]})
    return ism


def load_spec(path, z):
    d = np.loadtxt(path, comments="#")
    w, f = d[:, 0], d[:, 1]
    e = d[:, 2] if d.shape[1] > 2 else np.full_like(f, np.nan)
    return w / (1.0 + z), f, e


def despike(f, size=5, nsig=6.0):
    # clip isolated UPWARD spikes only (hot px / CR); absorption-safe
    f = np.asarray(f, float).copy()
    med = median_filter(f, size=size)
    resid = f - med
    mad = np.nanmedian(np.abs(resid[np.isfinite(resid)])) or 0.0
    hot = resid > nsig * 1.4826 * mad
    f[hot] = med[hot]
    return f


def lines_in(ism, w, pad=5.0):
    return [L for L in ism if w.min() + pad < L["lam"] < w.max() - pad]


def group_lines(inr, link=12.0):
    ls = sorted(inr, key=lambda L: L["lam"])
    grps = [[ls[0]]]
    for L in ls[1:]:
        if L["lam"] - grps[-1][-1]["lam"] <= link:
            grps[-1].append(L)
        else:
            grps.append([L])
    return grps


def _mg(x, sig, amps, lams):
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
    return np.interp(np.log10(np.clip(tau0, 1e-4, 1e6)), _ltau, _Fint)


def red_ew(N, f, lam, b):
    return (2.0 * b / C_KMS) * cog_F(TAU_K * N * f * lam / b)


def collect(inr, deb, ion=None, snr_min=2.0):
    lam, fo, y, ye = [], [], [], []
    for L in inr:
        if ion and L["ion"] != ion:
            continue
        ew, lo, hi = deb[round(L["lam"], 3)]
        if not np.isfinite(ew) or ew <= 0.02 or not (hi > 0):
            continue
        err = 0.5 * (lo + hi)
        if err > 0 and ew / err < snr_min:      # drop non-detections (2-sigma EW gate)
            continue
        lam.append(L["lam"]); fo.append(L["f"]); y.append(ew / L["lam"]); ye.append(err / L["lam"])
    return map(np.array, (lam, fo, y, ye))


def fit_b_N(lam, fo, y, ye, b0=60.0, logN0=14.5):
    def resid(p):
        pred = red_ew(10.0 ** p[1], fo, lam, p[0])
        return (np.log10(y) - np.log10(np.clip(pred, 1e-30, None))) / (ye / (y * np.log(10)) + 1e-3)
    return least_squares(resid, [b0, logN0], bounds=([5, 10], [300, 20]), max_nfev=5000).x


def fit_N(lam, fo, y, ye, b, logN0=14.0):
    def resid(p):
        return np.log10(y) - np.log10(np.clip(red_ew(10.0 ** p[0], fo, lam, b), 1e-30, None))
    return least_squares(resid, [logN0], bounds=([9], [21]), max_nfev=3000).x[0]


def _total_ew_err(inr, deb, cerr, ion):
    # per-detected-line TOTAL EW error = quadrature(photon MC, continuum placement), returned in W/lam units
    lam, fo, y, ye = collect(inr, deb, ion=ion)               # ye = MC err in W/lam
    tot = np.array([np.hypot(ye[i], (cerr.get(round(lam[i], 3), 0.0) or 0.0) / lam[i]) for i in range(len(lam))])
    return lam, fo, y, tot


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
        f = despike(f)
    inr = lines_in(ism, w)
    if len(inr) < 4:
        return None
    deb, cerr, fits = deblend_ews(w, f, e, inr, n_mc=n_mc, return_fits=True)
    fl, ff, fy, fye = collect(inr, deb, ion="Fe II")
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
        lam, fo, y, ye = collect(inr, deb, ion=ion)
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
    ax.set_xlabel("optical depth tau0"); ax.set_ylabel("W/lambda")
    ax.grid(True, which="both", alpha=0.2)
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


def _write_ism_csv(sn, g, ph, r, ismdir):
    # per-ion columns (+ logN error) and a rich per-line table: EW, the photon/continuum/total error budget,
    # detection flag, optical depth, saturation. this is the B4 rich output (was queued, now wired).
    with open(os.path.join(ismdir, f"{sn}_{g}_day{ph:.0f}_cog.csv"), "w", newline="") as fh:
        wri = csv.writer(fh); wri.writerow(["ion", "logN", "logN_err", "n_lines", "limit"])
        for ion, (ln, nl, lim, lnerr) in r["Ncol"].items():
            wri.writerow([ion, f"{ln:.3f}", f"{lnerr:.3f}" if np.isfinite(lnerr) else "", nl, ">" if lim else ""])
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
            det = bool(np.isfinite(ew) and np.isfinite(mc) and mc > 0 and ew / mc >= 2.0)
            ln = ioncol.get(L["ion"], r["logN_fe"])
            tau0 = TAU_K * 10 ** ln * L["f"] * L["lam"] / b
            wri.writerow([L["ion"], f"{L['lam']:.3f}", f"{L['f']:.4f}",
                          f"{ew:.4f}" if np.isfinite(ew) else "", f"{mc:.4f}" if np.isfinite(mc) else "",
                          f"{ce:.4f}" if np.isfinite(ce) else "", f"{tot:.4f}" if np.isfinite(tot) else "",
                          "yes" if det else "no", f"{tau0:.2f}", "yes" if tau0 > 5 else "no"])


def run_catalog(n_mc=150, min_fe=3):
    # loop over every NUV per-grating product, write per-epoch cog csv + a master summary
    cat = load_catalog()
    ism = load_lines()
    summ = []
    sne = sorted(d for d in os.listdir(OUT) if os.path.isdir(os.path.join(OUT, d)) and d.upper() in cat)
    for sn in sne:
        z = float(cat[sn.upper()]["z"])
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
            _write_ism_csv(sn, g, ph, r, ismdir)
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
    sfile = ISM_SUMMARY
    # adopted = highest-SNR epoch per SN AMONG ANCHORED epochs (screen is constant -> best anchored epoch wins).
    # unanchored epochs (fake Ia iron-photosphere columns) are reported but never adopted.
    best = {}
    for row in summ:
        if row["anchored"] != "yes":
            continue
        if row["sn"] not in best or row["fe_snr"] > best[row["sn"]]["fe_snr"]:
            best[row["sn"]] = row
    for row in summ:
        row["adopted"] = "yes" if best.get(row["sn"]) is row else ""
    with open(sfile, "w", newline="") as fh:
        wri = csv.DictWriter(fh, fieldnames=["sn", "grating", "phase", "z", "b", "b_err", "logN_FeII", "logN_FeII_err", "n_fe", "fe_snr", "anchored", "adopted"])
        wri.writeheader()
        for row in summ:
            wri.writerow(row)
    print(f"\n{len(summ)} epochs measured -> {sfile}")
    return summ


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--nmc", type=int, default=150)
    ap.add_argument("--min-fe", type=int, default=3)
    ap.add_argument("--sn", default=None, help="run a single SN (uppercase dir name)")
    a = ap.parse_args()
    if a.sn:
        cat = load_catalog(); ism = load_lines(); z = float(cat[a.sn.upper()]["z"])
        for g in NUV_GRATINGS:
            for p in sorted(glob.glob(f"{OUT}/{a.sn}/**/{g}/{a.sn}_*_{g}_native.txt", recursive=True)):
                r = analyze(*load_spec(p, z), ism, n_mc=a.nmc)
                if r:
                    ph = phase_of(p, cat)
                    print(f"{a.sn} {g} day{ph:.0f}: b={r['b']:.1f} logN(FeII)={r['logN_fe']:.2f} feSNR={r['fe_snr']:.1f}")
                    ismdir = os.path.join(OUT, a.sn, "absorption")
                    os.makedirs(ismdir, exist_ok=True)
                    _write_ism_csv(a.sn, g, ph, r, ismdir)
                    # write diagnostic plots (reuse the actual despiked flux + fits)
                    try:
                        _write_diag(a.sn, g, ph, r["w_used"], r["f_used"], r["inr"], r["deb"], r["b"], r["logN_fe"], r["Ncol"], ismdir, r["fits"])
                    except Exception as ex:
                        print(f"  WARNING diag {g} day{ph:.0f}: {ex}")
    else:
        run_catalog(n_mc=a.nmc, min_fe=a.min_fe)
