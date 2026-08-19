import os, glob, sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from paths import OUT

C = 2.99792458e5
OUTDIR = os.path.join(os.path.dirname(__file__), "..", "plot_outputs", "phase_2")
os.makedirs(OUTDIR, exist_ok=True)

z_cat = 0.002729          # catalog (NED, from the SN) -> 818 km/s
v_sys = 910.0             # Bowen host systemic (heliocentric), km/s
z_host = v_sys / C        # 0.003035

p = sorted(glob.glob(f"{OUT}/SN1998S/**/E230M/*day34*native.txt", recursive=True))[0]
d = np.loadtxt(p, comments="#")
w_obs, f = d[:, 0], d[:, 1]
ok = np.isfinite(f) & (f > 0)
w_obs, f = w_obs[ok], f[ok]

# bowen table 3 host components (delta v relative to systemic)
bowen_dv = np.array([-102.1, -90.1, -78.7, -71.3, -60.1, -50.5, -40.9, -24.9, 9.4])

lines = [(2600.173, "Fe II 2600"), (2382.765, "Fe II 2382"), (2796.352, "Mg II 2796")]

fig, axes = plt.subplots(1, 3, figsize=(16, 4.6))
for ax, (lam0, lab) in zip(axes, lines):
    # velocity in the NOTEBOOK rest frame (z_cat) relative to lam0
    wr_cat = w_obs / (1.0 + z_cat)
    v_cat = (wr_cat - lam0) / lam0 * C
    m = np.abs(v_cat) < 320
    ax.plot(v_cat[m], f[m], "k", lw=0.8, label="flux (z_cat frame)")

    # AOD windows used in the notebook (cell 24): core |v|<80, flank 100-260
    ax.axvspan(-80, 80, color="steelblue", alpha=0.12, label="AOD core |v|<80")
    for s in (-1, 1):
        ax.axvspan(s * 100, s * 260, color="orange", alpha=0.10)

    # where the host absorption ACTUALLY is in the z_cat frame:
    # v_notebook = v_abs - v_zcat, with v_abs = v_sys + dv, v_zcat = z_cat*C
    v_host_in_cat = (v_sys + bowen_dv) - z_cat * C
    for vv in v_host_in_cat:
        ax.axvline(vv, color="crimson", ls="-", lw=0.7, alpha=0.6)
    ax.axvline(v_host_in_cat.mean(), color="crimson", ls="--", lw=0.0)

    # MW version of this line (heliocentric ~0) in the z_cat frame, if it lands nearby
    for mw_lam, mw_lab in [(lam0, "MW self")] + ([(2803.531, "MW MgII2803")] if "Mg II" in lab else []):
        v_mw = (mw_lam / (1.0 + z_cat) - lam0) / lam0 * C
        if abs(v_mw) < 320:
            ax.axvline(v_mw, color="green", ls=":", lw=1.2)
            ax.text(v_mw, ax.get_ylim()[1] * 0.9, mw_lab, color="green", fontsize=6, rotation=90, va="top")

    ax.axvline(0, color="gray", ls=":", lw=0.8)
    ax.set_title(f"{lab}\nred=Bowen host comps (in z_cat frame), blue=AOD core", fontsize=8)
    ax.set_xlabel("v in z_cat rest frame (km/s)")
axes[0].set_ylabel("flux")
off = (v_sys - z_cat * C)
plt.suptitle(f"SN1998S E230M day34: z_cat=0.002729 (818 km/s) vs host systemic 910 km/s -> "
             f"host ISM sits at ~+{off:.0f} km/s, NOT v=0. AOD core |v|<80 mis-centered.", fontsize=10)
plt.tight_layout()
outp = os.path.join(OUTDIR, "1998s_zcheck.png")
plt.savefig(outp, dpi=110)
print("saved", outp)
print(f"z_cat={z_cat} -> {z_cat*C:.0f} km/s ; host systemic={v_sys:.0f} km/s ; offset={off:.0f} km/s")
print(f"host comps in z_cat frame span {(v_sys+bowen_dv-z_cat*C).min():.0f} .. {(v_sys+bowen_dv-z_cat*C).max():.0f} km/s")
print(f"MW MgII2803 in host-MgII2796 z_cat frame: {(2803.531/(1+z_cat)-2796.352)/2796.352*C:.0f} km/s")
