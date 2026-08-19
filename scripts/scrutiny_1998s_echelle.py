import os, glob, sys
import numpy as np
from astropy.io import fits
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from paths import OUT

C = 2.99792458e5
ROOT = os.path.join(os.path.dirname(__file__), "..")
OUTDIR = os.path.join(ROOT, "plot_outputs", "phase_2")
os.makedirs(OUTDIR, exist_ok=True)

EPO = os.path.join(OUT, "SN1998S/STIS/ECHELLE/1998-04-04_day34/E230M")
x1ds = sorted(glob.glob(os.path.join(EPO, "o4sz*_x1d.fits")))
native = os.path.join(EPO, "SN1998S_1998-04-04_day34_E230M_native.txt")

# ---- A3: per-order background fraction (Bowen: bg count rate a significant fraction of total) ----
h = fits.open(x1ds[0])[1].data
frac, wctr = [], []
for i in range(len(h)):
    g, bkg = h["GROSS"][i], h["BACKGROUND"][i]
    w = h["WAVELENGTH"][i]
    good = np.isfinite(g) & (g > 0)
    if good.sum() < 10:
        frac.append(np.nan); wctr.append(np.nanmedian(w)); continue
    frac.append(np.nanmedian(bkg[good] / g[good]))
    wctr.append(np.nanmedian(w))
frac = np.array(frac); wctr = np.array(wctr)

fig, (a0, a1) = plt.subplots(1, 2, figsize=(15, 4.6))
a0.plot(wctr, 100 * frac, "o-", color="firebrick")
a0.axhline(100 * np.nanmedian(frac), color="gray", ls="--", label=f"median {100*np.nanmedian(frac):.0f}%")
a0.set_xlabel("order central wvl (A)"); a0.set_ylabel("median BACKGROUND / GROSS  (%)")
a0.set_title("A3: echelle background as fraction of gross (per order)"); a0.legend()

# zoom on the Mg II 2796 order (well-centered, sporder 73)
oi = next(i for i in range(len(h)) if h["WAVELENGTH"][i].min() < 2796 < h["WAVELENGTH"][i].max()
          and abs((2796 - h["WAVELENGTH"][i].min())/(h["WAVELENGTH"][i].max()-h["WAVELENGTH"][i].min()) - 0.5) < 0.35)
w = h["WAVELENGTH"][oi]
a1.plot(w, h["GROSS"][oi], "k", lw=0.7, label="GROSS")
a1.plot(w, h["BACKGROUND"][oi], "firebrick", lw=0.9, label="BACKGROUND")
a1.plot(w, h["NET"][oi], "steelblue", lw=0.7, label="NET (gross-bg)")
a1.axvline(2796.35, color="orange", ls=":"); a1.axvline(2803.53, color="orange", ls=":")
a1.set_xlim(2785, 2812); a1.set_xlabel("wvl (A)"); a1.set_ylabel("counts")
a1.set_title(f"Mg II order (sporder {h['SPORDER'][oi]}): gross vs background"); a1.legend(fontsize=8)
plt.tight_layout(); plt.savefig(os.path.join(OUTDIR, "1998s_echelle_A3_background.png"), dpi=110); plt.close()
print("A3 median bg/gross fraction:", f"{100*np.nanmedian(frac):.1f}%", " range",
      f"{100*np.nanmin(frac):.0f}-{100*np.nanmax(frac):.0f}%")
print("Mg II order extrsize/bk:", h["EXTRSIZE"][oi], h["BK1SIZE"][oi], h["BK2SIZE"][oi],
      "offsets", h["BK1OFFST"][oi], h["BK2OFFST"][oi])

# ---- manual velocity-shift exploration (address z concern empirically) ----
d = np.loadtxt(native, comments="#")
w_obs, f = d[:, 0], d[:, 1]
ok = np.isfinite(f) & (f > 0)
w_obs, f = w_obs[ok], f[ok]
bowen_dv = np.array([-102.1, -90.1, -78.7, -71.3, -60.1, -50.5, -40.9, -24.9, 9.4])
vsys_try = [818.0, 850.0, 895.0, 910.0]     # catalog(818) .. Bowen(910)

lines = [(2600.173, "Fe II 2600"), (2586.650, "Fe II 2586"), (2796.352, "Mg II 2796")]
fig, axes = plt.subplots(1, 3, figsize=(16, 4.6))
for ax, (lam0, lab) in zip(axes, lines):
    for vs in vsys_try:
        zc = vs / C
        wr = w_obs / (1 + zc)
        v = (wr - lam0) / lam0 * C
        m = np.abs(v) < 260
        ax.plot(v[m], f[m] / np.nanpercentile(f[m], 90), lw=0.8, alpha=0.8, label=f"vsys={vs:.0f}")
    ax.axvline(0, color="gray", ls=":")
    # bowen components relative to systemic -> at v=dv when the frame IS the systemic
    for dv in bowen_dv:
        ax.axvline(dv, color="crimson", lw=0.5, alpha=0.4)
    ax.set_xlim(-260, 260); ax.set_xlabel("v rel line (km/s)"); ax.set_title(lab, fontsize=9)
axes[0].set_ylabel("norm flux"); axes[0].legend(fontsize=7, title="rest-frame anchor")
plt.suptitle("SN1998S E230M day34: which systemic centers the host ISM trough? "
             "(red = Bowen host comps at their systemic-relative v)", fontsize=10)
plt.tight_layout(); plt.savefig(os.path.join(OUTDIR, "1998s_echelle_vshift.png"), dpi=110); plt.close()
print("saved vshift + A3 plots")
