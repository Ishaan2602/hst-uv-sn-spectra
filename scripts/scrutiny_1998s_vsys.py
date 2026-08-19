import os, glob, sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from paths import OUT

C = 2.99792458e5
ROOT = os.path.join(os.path.dirname(__file__), "..")
OUTDIR = os.path.join(ROOT, "plot_outputs", "phase_2")
EPO = os.path.join(OUT, "SN1998S/STIS/ECHELLE/1998-04-04_day34/E230M")
native = os.path.join(EPO, "SN1998S_1998-04-04_day34_E230M_native.txt")

d = np.loadtxt(native, comments="#")
w_obs, f = d[:, 0], d[:, 1]
ok = np.isfinite(f) & (f > 0)
w_obs, f = w_obs[ok], f[ok]

# candidate systemics (heliocentric km/s): catalog z=0.002729 -> 818 ; NED galaxy ~895 ; Bowen 910
cands = {"catalog 818": 818.0, "NED gal 895": 895.0, "Bowen 910": 910.0}
# Bowen host ISM components span dv=-102..+9 rel to systemic 910 -> helio 808..919, strongest mid (~-60)

def helio_v(w, lam0):
    return (w - lam0) / lam0 * C

def abs_centroid(lam0, vwin=(600, 1120), contpad=260):
    v = helio_v(w_obs, lam0)
    m = (v > vwin[0] - contpad) & (v < vwin[1] + contpad)
    vv, ff = v[m], f[m]
    fk = (vv < vwin[0]) | (vv > vwin[1])          # continuum flanks outside the trough window
    if fk.sum() < 6:
        return None
    cont = np.polyval(np.polyfit(vv[fk], ff[fk], 1), vv)
    dec = 1 - ff / cont
    tr = (vv > vwin[0]) & (vv < vwin[1]) & (dec > 0)
    if dec[tr].sum() <= 0:
        return None
    cen = np.sum(vv[tr] * dec[tr]) / np.sum(dec[tr])
    return dict(v=vv, ff=ff, cont=cont, cen=cen)

lines = [(2600.173, "Fe II 2600"), (2586.650, "Fe II 2586"), (2374.461, "Fe II 2374"), (2796.352, "Mg II 2796")]
fig, axes = plt.subplots(1, 4, figsize=(18, 4.6))
print(f"{'line':13}{'abs centroid v_helio (km/s)':>28}")
for ax, (lam0, lab) in zip(axes, lines):
    v = helio_v(w_obs, lam0)
    m = (v > 400) & (v < 1300)
    ax.plot(v[m], f[m] / np.nanpercentile(f[m], 90), "k", lw=0.8)
    ax.axvspan(808, 919, color="crimson", alpha=0.12, label="Bowen comps 808-919")
    for lb, vs in cands.items():
        ax.axvline(vs, ls="--", lw=1.1, label=lb)
    r = abs_centroid(lam0) if "Mg II" not in lab else None
    if r is not None:
        ax.axvline(r["cen"], color="magenta", lw=1.6, label=f"centroid {r['cen']:.0f}")
        print(f"{lab:13}{r['cen']:>28.0f}")
    else:
        print(f"{lab:13}{'(saturated/skip centroid)':>28}")
    ax.set_xlim(400, 1300); ax.set_xlabel("v_helio (km/s)"); ax.set_title(lab, fontsize=9)
axes[0].set_ylabel("norm flux"); axes[0].legend(fontsize=6.5, loc="lower left")
plt.suptitle("SN1998S E230M day34: host ISM absorption in HELIOCENTRIC velocity. "
             "Which systemic (dashed) sits where Bowen's comps (red band) + measured centroid (magenta) are?", fontsize=9.5)
plt.tight_layout(); plt.savefig(os.path.join(OUTDIR, "1998s_echelle_vsys_helio.png"), dpi=110); plt.close()
print("saved 1998s_echelle_vsys_helio.png")
