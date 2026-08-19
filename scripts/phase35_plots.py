#!/usr/bin/env python3
# phase-3.5 mentor/paper plots. short titles, no giant headers. saves to plot_outputs/phase_3/.
# 1) joint absorption+emission fit on one epoch (2023ixf d66)
# 2) flux uncertainty budget (photon MC vs the continuum/window/smoothing systematics)
# 3) emission model selection across the catalog (which profile wins, by type)
import os, sys, glob, json, collections
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import paths, ism
import emission_products as ep

DEST = os.path.join(paths.ROOT, "plot_outputs", "phase_3")
os.makedirs(DEST, exist_ok=True)
C = 2.99792458e5
cat = ism.load_catalog(); lines = ism.load_lines()


def _ixf_native(day_tag, grating="G230LB"):
    for case in ("SN2023ixf", "SN2023IXF"):
        hits = glob.glob(f"{paths.OUT}/{case}/**/{grating}/*{day_tag}*native.txt", recursive=True)
        if hits:
            return hits[0]
    return None


def plot_joint():
    # left panel: ISM absorption at d14 (emission hasn't wrecked the continuum yet)
    # right panel: Mg II CSM emission at d66
    p14 = _ixf_native("day14"); p66 = _ixf_native("day66")
    if not p14 or not p66:
        print("missing 2023ixf files for joint plot"); return

    z = float(cat["SN2023IXF"]["z"]); mw = float(cat["SN2023IXF"]["ebv"]); host = float(cat["SN2023IXF"].get("host_ebv") or 0)
    fig, (a0, a1) = plt.subplots(1, 2, figsize=(13, 4.4))

    # left: Fe II ISM forest at d14 (continuum-normalized), line positions labeled
    w, f, e = ep.load_spec(p14, z)
    f = f * ep.deredden(w, z, mw, host)
    ok = np.isfinite(f); w, f = w[ok], f[ok]
    m = (w > 2310) & (w < 2650)
    from scipy.ndimage import median_filter as _mf
    cont_sm = _mf(f, size=151)           # wide median filter tracks the broad continuum shape
    norm = np.where(cont_sm > 0, f / cont_sm, np.nan)
    a0.plot(w[m], norm[m], color="0.3", lw=0.7)
    a0.axhline(1, color="0.8", lw=0.6)
    for L in ism.lines_in(ism.load_lines(), w):
        if 2310 < L["lam"] < 2650:
            a0.axvline(L["lam"], color="steelblue", lw=0.5, ls=":", alpha=0.8)
            a0.text(L["lam"], 1.10, L["ion"].split()[0], rotation=90, fontsize=5.5,
                    color="steelblue", va="bottom", ha="center")
    a0.set_ylim(0.3, 1.25); a0.set_xlabel("rest wavelength (A)"); a0.set_ylabel("flux / continuum")
    # b=61.7 km/s, logN(FeII)=15.28 from the ISM summary (adopted epoch)
    a0.set_title("2023ixf d14  ISM absorption  b=62 km/s  logN(FeII)=15.3", fontsize=9)

    # right: Mg II emission at d66 + best-fit model
    w, f, e = ep.load_spec(p66, z)
    f = f * ep.deredden(w, z, mw, host)
    ok = np.isfinite(f); w, f = w[ok], f[ok]
    bf = ep.bostroem_flux(w, f, e=e, lam0=ep.MG2, vline=ep._MG2_DEFAULT, deg=1)
    prof = ep._fit_profiles(bf)
    keep = bf["emis"] & ~bf["notch"]; vv, ff = bf["v"][keep], bf["fc"][keep]
    o = np.argsort(vv); a1.plot(vv[o], ff[o] * 1e15, color="k", lw=1.1, label="data")
    a1.axhline(0, color="0.8", lw=0.6)
    mcol = {"gaussian": "tab:green", "lorentzian": "tab:purple", "skew": "tab:orange", "kwok": "crimson"}
    best = prof["best_model"] if prof else None; vg = np.linspace(vv.min(), vv.max(), 500)
    for name, (fn, pp) in ((prof or {}).get("_fits", {}) or {}).items():
        a1.plot(vg, fn(vg, *pp) * 1e15, color=mcol.get(name, "0.5"),
                lw=2.4 if name == best else 0.9, ls="-" if name == best else "--",
                label=name + (" *" if name == best else ""))
    a1.set_xlabel("velocity (km/s)"); a1.set_ylabel("flux (1e-15)")
    a1.set_title(f"2023ixf d66  Mg II emission  (best: {best})", fontsize=9); a1.legend(fontsize=7)
    fig.tight_layout(); fig.savefig(f"{DEST}/joint_absorption_emission.png", dpi=130); plt.close(fig)
    print("saved joint_absorption_emission.png")


def _syst_components(w, f, lam0, vline, vcont, deg, baseF, sigphot):
    # per-source error contributions (1e-15): photon MC + the three systematics separately
    d_int = 0.0
    for dv in (-1500, 1500):
        vl = (vline[0] + dv, vline[1] - dv)
        Fv = ep.bostroem_flux(w, f, lam0=lam0, vline=vl, vcont=vcont, deg=deg)["F"]
        d_int = max(d_int, abs(Fv - baseF))
    d_cont = 0.0
    for dc in (-1500, 1500):
        Fv = ep.bostroem_flux(w, f, lam0=lam0, vline=vline, vcont=(vcont[0] - dc, vcont[1] + dc), deg=deg)["F"]
        d_cont = max(d_cont, abs(Fv - baseF))
    d_sm = abs(ep.bostroem_flux(w, f, lam0=lam0, vline=vline, vcont=vcont, deg=deg, smooth=5)["F"] - baseF)
    return np.array([sigphot, d_int, d_cont, d_sm]) * 1e15


def plot_uncertainty():
    # read the already-computed systematic errors from the emission.json (avoids re-run + duplicate paths)
    import json as _json
    ixf_json = glob.glob(os.path.join(paths.OUT, "SN2023ixf", "SN2023IXF_emission.json")) or \
               glob.glob(os.path.join(paths.OUT, "SN2023IXF", "SN2023IXF_emission.json"))
    if not ixf_json:
        print("no 2023ixf emission json found"); return
    d = _json.load(open(ixf_json[0]))
    rows = [(r["phase"], r.get("flux_err_photon", 0) or 0, r.get("flux_err", 0) or 0, r["flux"])
            for r in d["emission"]["mg2"]
            if r.get("flux_reliable") and r.get("flux") and r["flux"] > 0]
    if not rows:
        print("no reliable rows"); return
    rows.sort(key=lambda x: x[0])
    labs = [f"d{int(ph)}" for ph, *_ in rows]
    # flux/flux_err already stored in 1e-15 units in the json
    phot = np.array([r[1] for r in rows])
    syst = np.array([r[2] for r in rows])
    flux = np.array([r[3] for r in rows])
    x = np.arange(len(labs)); wbar = 0.35
    fig, (a0, a1) = plt.subplots(1, 2, figsize=(13, 4.4))
    a0.bar(x - wbar / 2, phot, wbar, label="photon MC", color="0.6")
    a0.bar(x + wbar / 2, syst, wbar, label="systematic (continuum/window/smooth)", color="tab:blue")
    a0.set_xticks(x); a0.set_xticklabels(labs, fontsize=8)
    a0.set_ylabel("flux error (1e-15)"); a0.legend(fontsize=8)
    a0.set_title("2023ixf Mg II: systematic vs photon error per epoch", fontsize=9)
    a1.plot(flux, syst / flux * 100, "o", color="tab:blue", label="systematic")
    a1.plot(flux, phot / flux * 100, "s", color="0.5", label="photon MC")
    a1.set_xlabel("Mg II flux (1e-15)"); a1.set_ylabel("fractional error (%)")
    a1.legend(fontsize=8)
    a1.set_title("fractional error vs flux: syst ~10-30%, photon <5%", fontsize=9)
    fig.tight_layout(); fig.savefig(f"{DEST}/uncertainty_budget.png", dpi=130); plt.close(fig)
    print("saved uncertainty_budget.png")


def plot_model_selection():
    def _t(typ):
        t = (typ or "?").upper()
        return "IIn" if "IIN" in t else "II" if "II" in t else "Ia" if "IA" in t else \
               "Ibc" if ("IB" in t or "IC" in t) else "SLSN" if "SLSN" in t else "other"
    by_t = collections.defaultdict(collections.Counter)
    for p in glob.glob(os.path.join(paths.OUT, "*", "*_emission.json")):
        d = json.load(open(p)); typ = _t(d["sn_type"] or "?")
        for r in d["emission"]["mg2"]:
            if r.get("best_model"):
                by_t[typ][r["best_model"]] += 1
    order = [t for t in ["II", "IIn", "Ia", "Ibc", "SLSN", "other"] if t in by_t]
    models = ["gaussian", "lorentzian", "skew", "kwok"]
    mcol = {"gaussian": "tab:green", "lorentzian": "tab:purple", "skew": "tab:orange", "kwok": "crimson"}
    fig, ax = plt.subplots(figsize=(8.5, 4.4))
    bottom = np.zeros(len(order))
    for m in models:
        vals = np.array([by_t[t].get(m, 0) for t in order], float)
        ax.bar(order, vals, bottom=bottom, label=m, color=mcol[m])
        bottom += vals
    ax.set_ylabel("clean-emission epochs"); ax.set_xlabel("SN type")
    ax.set_title("Mg II best-fit profile by SN type (catalog-wide)", fontsize=9)
    ax.legend(fontsize=8)
    fig.tight_layout(); fig.savefig(f"{DEST}/model_selection_by_type.png", dpi=130); plt.close(fig)
    print("saved model_selection_by_type.png")


if __name__ == "__main__":
    plot_joint()
    plot_uncertainty()
    plot_model_selection()
    print("DONE")
