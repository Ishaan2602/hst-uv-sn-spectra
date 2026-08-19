#!/usr/bin/env python3
# generates 4 clean summary plots for the mentor progress update email.
# run from the repo root: python scripts/make_mentor_plots.py

import os, glob, sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy.special import wofz
from scipy.signal import savgol_filter

sys.path.insert(0, os.path.join(os.path.dirname(__file__)))
import emission_products as ep

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
from paths import OUT as OUT5
DEST = os.path.join(ROOT, "docs", "mentor_update")
os.makedirs(DEST, exist_ok=True)

plt.rcParams.update({"font.size": 11, "axes.spines.top": False, "axes.spines.right": False,
                     "figure.dpi": 150, "savefig.dpi": 150, "savefig.bbox": "tight"})

# ---- helpers (same as investigation notebooks) --------------------------------------------------
C_KMS = 2.99792458e5
MG2, LYA = 2799.94, 1215.67
F_LYA, GAMMA_LYA = 0.4164, 6.265e8
TAU_K = 1.4973e-15

def load_spec(path, z):
    d = np.loadtxt(path, comments="#")
    w, f = d[:, 0], d[:, 1]
    return w / (1 + z), f

def find1(pattern):
    hits = glob.glob(os.path.join(OUT5, pattern), recursive=True)
    return hits[0] if hits else None

def kwok_shell(v, A, mu, fwhm, vc, vin):
    sig = fwhm / (2 * np.sqrt(2 * np.log(2)))
    g = np.exp(-0.5 * ((v - mu) / sig) ** 2)
    vh = mu + vc
    supp = np.ones_like(v, float)
    inside = np.abs(v - vh) < vin
    supp[inside] = np.exp(-0.5 * (vin**2 - (v[inside] - vh)**2) / sig**2)
    return A * g * supp

def voigt_tau(w_A, logN, b_kms):
    # Voigt optical depth at rest LYA (Milky-Way + host frame at v=0)
    b = b_kms * 1e5
    lam0 = LYA * 1e-8
    a = GAMMA_LYA * lam0 / (4 * np.pi * b)
    x = (w_A * 1e-8 - lam0) * C_KMS * 1e5 / (b * lam0)
    return TAU_K * 10**logN * F_LYA * lam0 * 1e8 / (b_kms) * np.real(wofz(x + 1j * a))

# =================================================================================================
# PLOT 1 — Mg II day-66 profile: flat-top shell vs. Kwok asymmetric shell
# =================================================================================================
sn = "SN2023IXF"; z = ep.zof(sn); mw = ep.ebvof(sn); host = ep.hostof(sn)
p66 = find1("SN2023IXF/**/G230LB/*day66*_native.txt")
if p66:
    w, f = load_spec(p66, z)
    f = f * ep.deredden(w, z, mw, host)
    ok = np.isfinite(f); w, f = w[ok], f[ok]
    v = (w - MG2) / MG2 * C_KMS
    m = (v > -14000) & (v < 10000)
    vv, ff = v[m], f[m]
    # continuum + cont-sub
    cont_m = (v > 10500) & (v < 14000)
    cont = np.polyval(np.polyfit(v[cont_m], f[cont_m], 1), v)
    fc = f - cont; fcc = fc[m]
    # fit Kwok
    fitm = m & (np.abs(v) > 350)
    vfit, ffit = v[fitm], fc[fitm]
    A0 = np.nanmax(ffit)
    lo, hi = [0,-8000,3000,-3000,1000],[3,2000,15000,6000,9000]
    pk = curve_fit(kwok_shell, vfit, ffit/A0, p0=[1,-2995,7580,1750,5000],
                   bounds=(lo,hi), maxfev=30000)[0]
    pk[0] *= A0
    # flat-top shell from phase-2 (known over-shooter)
    from scipy.special import erf as _erf
    def shell_flat(v, A, v0, vout, tau, edge):
        x = np.abs(v - v0)
        prof = 0.5 * (1 - _erf((x - vout) / edge))
        atten = np.where(v > v0, np.exp(-tau * (v - v0) / vout), 1.0)
        return A * prof * atten
    fitm2 = m & (np.abs(v) > 350)
    vfit2, ffit2 = v[fitm2], fc[fitm2]; A2 = np.nanmax(ffit2)
    try:
        ps = curve_fit(shell_flat, vfit2, ffit2, p0=[A2,-2000,6000,1,1500],
                       bounds=([0,-6000,3000,0,300],[2*A2,2000,11000,6,4000]), maxfev=20000)[0]
    except Exception:
        ps = None
    vg = np.linspace(-14000, 10000, 1500)
    fig, ax = plt.subplots(figsize=(8, 4.5))
    ax.plot(vv, fcc*1e15, color="0.2", lw=1.2, label="data (cont-sub)")
    if ps is not None:
        ax.plot(vg, shell_flat(vg, *ps)*1e15, color="crimson", lw=1.5, ls="--",
                label=f"flat-top shell (old; overshoots)")
    ax.plot(vg, kwok_shell(vg, *pk)*1e15, color="navy", lw=2, label="Kwok off-center shell (fits)")
    for lam, c in [(2796.35, "steelblue"), (2803.53, "steelblue")]:
        ax.axvline((lam-MG2)/MG2*C_KMS, color=c, ls=":", lw=1, alpha=0.6)
    ax.axvline(pk[1]+pk[3]-pk[4], color="gray", ls=":", lw=1, alpha=0.6, label=f"shell blue edge {pk[1]+pk[3]-pk[4]:.0f} km/s")
    ax.axhline(0, color="0.7", lw=0.7)
    ax.set_xlabel("velocity from Mg II 2800 (km/s)")
    ax.set_ylabel("flux (1e-15 erg/s/cm²/Å)")
    ax.set_title("SN2023IXF day 66 — Mg II CSM shell")
    ax.legend(fontsize=9)
    ax.set_xlim(-14000, 9000)
    fig.tight_layout()
    fig.savefig(os.path.join(DEST, "1_mg2_shell_fit.png"))
    plt.close()
    print("plot 1 done")
else:
    print("plot 1 SKIP: day66 G230LB not found")

# =================================================================================================
# PLOT 2 — kinematic evolution: inner velocity + flux vs phase
# =================================================================================================
# read from the already-generated emission.json so we don't re-fit everything
import json
jp = os.path.join(OUT5, "SN2023IXF", "SN2023IXF_emission.json")
d = json.load(open(jp))
mg2_epochs = d["emission"]["mg2"]
lya_epochs  = d["emission"]["lya"]

ph_mg  = np.array([r["phase"] for r in mg2_epochs])
fl_mg  = np.array([r["flux"]  for r in mg2_epochs])
iv_mg  = np.array([r["shell"]["inner_v"] for r in mg2_epochs if "shell" in r])
ph_iv  = np.array([r["phase"] for r in mg2_epochs if "shell" in r])
ph_ly  = np.array([r["phase"] for r in lya_epochs])
fl_ly  = np.array([r["flux"]  for r in lya_epochs])

BOST_MG2 = {66:367, 199:173, 311:143, 619:56.5, 722:47.6}
BOST_LYA = {199:214, 311:249, 619:132, 722:128}

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10.5, 4.5))

ax1.plot(ph_mg, fl_mg, "o-", color="darkorange", ms=6, lw=1.3, label="Mg II (ours, MW+host dered)")
ax1.plot(ph_ly, fl_ly, "s-", color="steelblue",  ms=6, lw=1.3, label="Lya (ours)")
ax1.plot(list(BOST_MG2), list(BOST_MG2.values()), "*", color="saddlebrown", ms=12, ls="none",
         label="Mg II Bostroem+2026")
ax1.plot(list(BOST_LYA), list(BOST_LYA.values()), "P", color="navy", ms=9, ls="none",
         label="Lya Bostroem+2026")
ax1.set_xscale("log"); ax1.set_yscale("log")
ax1.set_xlabel("phase (days)"); ax1.set_ylabel("line flux (1e-15 erg/s/cm²)")
ax1.set_title("CSM line flux vs time")
ax1.legend(fontsize=8.5)

ax2.plot(ph_iv, iv_mg, "o-", color="purple", ms=6, lw=1.3)
ax2.axhline(-3500, color="gray", ls="--", lw=1, alpha=0.6)
ax2.set_xlabel("phase (days)"); ax2.set_ylabel("shell blue-edge velocity (km/s)")
ax2.set_title("Mg II shell inner velocity (blue-side asymmetry fading)")
ax2.set_xscale("log")

fig.tight_layout()
fig.savefig(os.path.join(DEST, "2_flux_and_kinematics.png"))
plt.close()
print("plot 2 done")

# =================================================================================================
# PLOT 3 — N(HI) Lya fit at day 311
# =================================================================================================
p_lya = find1("SN2023IXF/**/G140L/*day311*_native.txt")
if not p_lya:
    p_lya = find1("SN2023IXF/**/G140L/*day308*_native.txt")
if p_lya:
    w, f = load_spec(p_lya, z)
    f = f * ep.deredden(w, z, mw, 0.0)   # MW-only at Lya; host dered makes it worse (host Lya fully absorbed)
    ok = np.isfinite(f); w, f = w[ok], f[ok]
    m = (w > 1180) & (w < 1270)
    wm, fm = w[m], f[m]
    # adopted best-fit params from ism_investigation2
    logN, b_kms, vabs = 20.83, 30.0, 0.0
    # simple continuum + kwok emission background (approximate)
    contm = ((wm > 1180) & (wm < 1205)) | ((wm > 1232) & (wm < 1270))
    if contm.sum() > 3:
        cp = np.polyfit(wm[contm], fm[contm], 1)
    else:
        cp = np.polyfit([wm[0], wm[-1]], [np.nanmedian(fm[:10]), np.nanmedian(fm[-10:])], 1)
    cont = np.polyval(cp, wm)
    # model: cont * exp(-tau)
    tau = voigt_tau(wm, logN, b_kms)
    model = cont * np.exp(-tau)

    fig, ax = plt.subplots(figsize=(8, 4.5))
    ax.plot(wm, fm*1e15, color="0.2", lw=1.2, label="data (MW-dered)")
    ax.plot(wm, model*1e15, color="firebrick", lw=2, label=f"Voigt fit  logN(HI) = {logN}")
    ax.plot(wm, cont*1e15, color="steelblue", lw=1.2, ls="--", label="continuum + emission")
    # mark the absorber
    ax.axvline(LYA, color="gray", ls=":", lw=1, alpha=0.7, label="Lya rest 1216 Å")
    # geocoronal mask region
    ax.axvspan(1213, 1218, color="gold", alpha=0.15, label="geocoronal mask")
    ax.set_xlabel("rest wavelength (Å)")
    ax.set_ylabel("flux (1e-15 erg/s/cm²/Å)")
    ax.set_title("SN2023IXF damped Lya — N(HI) measurement")
    ax.legend(fontsize=9)
    ax.set_xlim(1182, 1265)
    fig.tight_layout()
    fig.savefig(os.path.join(DEST, "3_lya_NHI_fit.png"))
    plt.close()
    print("plot 3 done")
else:
    print("plot 3 SKIP: G140L day311/308 not found")

# =================================================================================================
# PLOT 4 — metallicity bar chart
# =================================================================================================
solar = {"Fe": -4.50, "Cr": -6.36, "Zn": -7.44, "Mn": -6.57}
logN  = {"Fe": 15.16, "Cr": 13.9, "Zn": 13.3, "Mn": 13.1}
logNH = 20.83
ions  = ["Zn", "Cr", "Mn", "Fe"]
xH    = {ion: logN[ion] - logNH - solar[ion] for ion in ions}
colors = {"Zn": "steelblue", "Cr": "darkorange", "Mn": "purple", "Fe": "crimson"}
labels = {"Zn": "Zn II\n(undepleted)", "Cr": "Cr II", "Mn": "Mn II", "Fe": "Fe II\n(dust-depleted)"}

fig, ax = plt.subplots(figsize=(7, 4.5))
xs = np.arange(len(ions))
bars = ax.bar(xs, [xH[ion] for ion in ions], color=[colors[ion] for ion in ions],
              width=0.55, edgecolor="white", linewidth=0.5)
ax.axhline(0, color="0.3", lw=1.2, ls="--", label="solar ([X/H] = 0)")
ax.set_xticks(xs); ax.set_xticklabels([labels[ion] for ion in ions])
ax.set_ylabel("[X/H]  (log, solar-scaled)")
ax.set_title("Foreground gas metallicity — SN2023IXF sightline")
# annotate
for bar, ion in zip(bars, ions):
    val = xH[ion]
    if val < -0.25:
        ax.text(bar.get_x() + bar.get_width()/2, val/2,
                f"{val:+.2f}", ha="center", va="center", fontsize=9.5, color="white")
    else:
        ax.text(bar.get_x() + bar.get_width()/2, 0.04,
                f"{val:+.2f}", ha="center", va="bottom", fontsize=9.5, color="0.3")
ax.set_ylim(-1.5, 0.5)
ax.legend(fontsize=9)
ax.text(0.98, 0.05, f"logN(HI) = {logNH}  (damped Lya)", transform=ax.transAxes,
        ha="right", fontsize=9, color="0.4")
fig.tight_layout()
fig.savefig(os.path.join(DEST, "4_metallicity.png"))
plt.close()
print("plot 4 done")

print(f"\nplots written to {DEST}")
