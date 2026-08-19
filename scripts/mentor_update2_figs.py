#!/usr/bin/env python3
# generate the 4 figures for docs/mentor_update2/ from the validated machinery.
import os, sys, json
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__))))
import emission_products as ep

ROOT = ep.ROOT
OUT = ep.OUT
DST = os.path.join(ROOT, "docs", "mentor_update2")
os.makedirs(DST, exist_ok=True)
MG2 = ep.MG2


def prof(sn, rel):
    z, mw, host = ep.zof(sn), ep.ebvof(sn), ep.hostof(sn)
    w, f, e = ep.load_spec(os.path.join(OUT, sn, rel), z)
    fac = ep.deredden(w, z, mw, host); f = f * fac
    ok = np.isfinite(f); w, f = w[ok], f[ok]
    bf = ep.bostroem_flux(w, f, lam0=MG2, vline=ep.MG2_WINDOW.get(sn.upper(), ep._MG2_DEFAULT))
    keep = bf["emis"] & ~bf["notch"]
    return w, f, bf, bf["v"][keep], bf["fc"][keep]


# --- plot 1: emission model selection across regimes ------------------------------------------------
def plot1():
    tests = [("SN2023IXF", "STIS/CCD/2023-07-24_day66/G230LB/SN2023IXF_2023-07-24_day66_G230LB_native.txt", "2023ixf d66 (broad shell)"),
             ("SN-2005IP", "STIS/MAMA/2014-03-28_day3065/G230L/SN-2005IP_2014-03-28_day3065_G230L_native.txt", "2005ip d3065 (narrow IIn)"),
             ("SN2010JL", "STIS/CCD/2011-01-23_day82/G230LB/SN2010JL_2011-01-23_day82_G230LB_native.txt", "2010JL d82 (narrow IIn)")]
    col = {"gaussian": "tab:green", "lorentzian": "tab:purple", "skew": "tab:orange", "kwok": "crimson"}
    fig, axes = plt.subplots(1, 3, figsize=(15, 4.3))
    for ax, (sn, rel, lab) in zip(axes, tests):
        w, f, bf, vv, ff = prof(sn, rel)
        sig = float(np.nanstd(bf["fc"][bf["contfit"]])) or 0.05 * np.nanmax(ff)
        fits = ep._fit_models(vv, ff)
        ics = {n: ep._model_ic(vv, ff, fn, p, sig) for n, (fn, p) in fits.items()}
        best = min(ics, key=lambda n: ics[n]["bic"])
        o = np.argsort(vv); ax.axhline(0, color="0.8", lw=0.6)
        ax.plot(vv[o], ff[o] * 1e15, color="k", lw=1.1, label="data")
        vg = np.linspace(vv.min(), vv.max(), 500)
        for n, (fn, p) in fits.items():
            ax.plot(vg, fn(vg, *p) * 1e15, color=col[n], lw=2.4 if n == best else 1.0,
                    ls="-" if n == best else "--", label=n + (" *" if n == best else ""))
        ax.set_title(f"{lab}\nBIC best: {best}", fontsize=10)
        ax.set_xlabel("velocity [km/s]"); ax.set_ylabel("flux [1e-15]"); ax.legend(fontsize=7)
    fig.suptitle("Mg II emission profile: broad CSM shells -> skew/kwok, narrow-line IIn -> Lorentzian", fontsize=11)
    fig.tight_layout(); fig.savefig(os.path.join(DST, "1_emission_models.png"), dpi=120); plt.close(fig)


# --- plot 2: P-Cygni (degenerate cont-sub vs normalized v_phot) -------------------------------------
def pcyg2(v, Ae, mue, sige, Aa, mua, siga):
    return Ae * np.exp(-0.5 * ((v - mue) / sige) ** 2) - Aa * np.exp(-0.5 * ((v - mua) / siga) ** 2)


def plot2():
    sn, rel = "PTF11KLY", "STIS/MAMA/2011-10-07_day44/G230L/PTF11KLY_2011-10-07_day44_G230L_native.txt"
    w, f, bf, vv, ff = prof(sn, rel)
    fig, (a0, a1) = plt.subplots(1, 2, figsize=(13, 4.5))
    A = np.nanmax(ff)
    p = curve_fit(pcyg2, vv, ff / A, p0=[1, -1000, 4000, 1, -5000, 3000],
                  bounds=([0, -6000, 500, 0, -12000, 500], [5, 6000, 12000, 5, 2000, 12000]), maxfev=40000)[0]
    Ae, mue, sige, Aa, mua, siga = p; Ae *= A; Aa *= A
    o = np.argsort(vv); vg = np.linspace(vv.min(), vv.max(), 500)
    a0.axhline(0, color="0.8", lw=0.6); a0.plot(vv[o], ff[o] * 1e15, color="k", lw=1.1, label="data (cont-subtracted)")
    a0.plot(vg, Ae * np.exp(-0.5 * ((vg - mue) / sige) ** 2) * 1e15, color="tab:green", ls="--", label="emission comp")
    a0.plot(vg, -Aa * np.exp(-0.5 * ((vg - mua) / siga) ** 2) * 1e15, color="tab:red", ls="--", label="absorption comp")
    a0.plot(vg, pcyg2(vg, Ae, mue, sige, Aa, mua, siga) * 1e15, color="tab:blue", lw=2, label="sum")
    a0.set_title("continuum-subtracted P-Cygni fit = DEGENERATE\n(two big components cancel; velocity meaningless)", fontsize=9)
    a0.set_xlabel("velocity [km/s]"); a0.set_ylabel("flux [1e-15]"); a0.legend(fontsize=7)
    # right: normalized spectrum + absorption minimum
    z, mw, host = ep.zof(sn), ep.ebvof(sn), ep.hostof(sn)
    v = (w - MG2) / MG2 * ep.C_KMS
    reg = (v > -22000) & (v < 12000); feat = (v > -16000) & (v < 7000); cf = reg & ~feat
    cont = np.polyval(np.polyfit(w[cf], f[cf], 2), w); fn = f / cont
    vphot = ep._vphot_normalized(w, f, MG2)
    a1.axhline(1, color="0.8", lw=0.6); a1.plot(v[reg], fn[reg], color="k", lw=1.0)
    a1.axvline(vphot, color="crimson", ls="--", lw=1.3, label=f"v_phot = {vphot} km/s")
    a1.set_ylim(0, 2.2); a1.set_xlim(-22000, 12000)
    a1.set_title("continuum-NORMALIZED spectrum: photospheric velocity\nfrom the blueshifted absorption minimum (rough)", fontsize=9)
    a1.set_xlabel("velocity [km/s]"); a1.set_ylabel("f / pseudo-continuum"); a1.legend(fontsize=8)
    fig.suptitle("PTF11KLY (SN2011fe) day 44 -- how to characterize a photospheric P-Cygni epoch", fontsize=11)
    fig.tight_layout(); fig.savefig(os.path.join(DST, "2_pcygni.png"), dpi=120); plt.close(fig)


# --- plot 3: ISM CoG anchor -- clean vs Ia contamination (compose the generated diagnostics) --------
def plot3():
    clean = os.path.join(OUT, "SN2023IXF", "absorption", "SN2023IXF_G230LB_day14_cog.png")
    dirty = os.path.join(OUT, "SN2017CBV", "absorption", "SN2017CBV_G230L_day20_cog.png")
    fig, (a0, a1) = plt.subplots(1, 2, figsize=(13, 5))
    for ax, img, lab in ((a0, clean, "clean backlight (2023ixf): anchored, logN=15.3"),
                         (a1, dirty, "Ia contamination (2017CBV): all saturated, logN floats to 16.5")):
        if os.path.exists(img):
            ax.imshow(plt.imread(img)); ax.set_title(lab, fontsize=10)
        ax.axis("off")
    fig.suptitle("Curve-of-growth anchor test: a real absorber has a low-tau0 line pinning N; the Ia does not", fontsize=11)
    fig.tight_layout(); fig.savefig(os.path.join(DST, "3_ism_anchor.png"), dpi=120); plt.close(fig)


# --- plot 4: 2023ixf flux vs phase, reliable vs pcygni, vs Bostroem+2026 -------------------------
def plot4():
    d = json.load(open(os.path.join(OUT, "SN2023IXF", "SN2023IXF_emission.json")))
    # Bostroem+2026 Table 4 values (1e-15 erg/s/cm2, MW+host dered)
    bost_mg2 = {66: 367, 199: 173, 311: 143, 619: 56.5, 722: 47.6}
    bost_lya  = {199: 214, 311: 249, 619: 132, 722: 128}
    fig, ax = plt.subplots(figsize=(7, 4.6))
    for line, col, lab in (("mg2", "darkorange", "Mg II 2800"), ("lya", "steelblue", "Lya 1216")):
        recs = d["emission"][line]
        rel = [r for r in recs if r.get("flux_reliable", True)]
        unr = [r for r in recs if not r.get("flux_reliable", True)]
        if rel:
            ax.errorbar([r["phase"] for r in rel], [r["flux"] for r in rel], yerr=[r["flux_err"] or 0 for r in rel],
                        fmt="o-", color=col, ms=5, capsize=2, label=lab)
        if unr:
            ax.plot([r["phase"] for r in unr], [r["flux"] for r in unr], "x", color=col, ms=7, alpha=0.6,
                    label=f"{lab} pcygni (unreliable)")
    bost_lines = [(bost_mg2, "darkorange"), (bost_lya, "steelblue")]
    first = True
    for bd, col in bost_lines:
        ax.plot(list(bd.keys()), list(bd.values()), "*", color=col, ms=10, ls="none",
                label="Bostroem+2026" if first else "_")
        first = False
    ax.legend(fontsize=8)
    ax.set_yscale("log"); ax.set_xlabel("phase [days]"); ax.set_ylabel("line flux [1e-15 erg/s/cm2]")
    ax.set_title("SN2023ixf UV emission-line flux vs phase\n(reliable epochs solid; P-Cygni flagged; stars = Bostroem+2026)", fontsize=9)
    fig.tight_layout(); fig.savefig(os.path.join(DST, "4_flux_phase.png"), dpi=120); plt.close(fig)


if __name__ == "__main__":
    plot1(); print("plot1 done")
    plot2(); print("plot2 done")
    plot3(); print("plot3 done")
    plot4(); print("plot4 done")
