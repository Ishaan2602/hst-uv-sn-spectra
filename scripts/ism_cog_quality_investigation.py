#!/usr/bin/env python3
# scratch: does the curve-of-growth coherence (scatter of the Fe II points about the fitted curve)
# separate a real ISM absorber from SN iron-photosphere contamination? compare known-clean backlights
# vs the Ia false positives that are marked adopted in ism_cog_summary.csv.
import os, sys, glob
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import ism

OUT = ism.OUT
cat = ism.load_catalog()
lines = ism.load_lines()


def cog_resid(inr, deb, b, logN_fe):
    # rms scatter (dex) of the Fe II W/lam points about the single fitted curve of growth,
    # plus the minimum tau0 over the detected Fe II lines (the linear-part ANCHOR: <~1 means a
    # weak line pins N; all lines saturated -> N floats = the contamination signature).
    lam, fo, y, ye = ism.collect(inr, deb, ion="Fe II")
    if len(lam) < 3:
        return np.nan, 0, np.nan
    pred = ism.red_ew(10.0 ** logN_fe, fo, lam, b)
    r = np.log10(y) - np.log10(np.clip(pred, 1e-30, None))
    tau0 = ism.TAU_K * 10.0 ** logN_fe * fo * lam / b
    return float(np.sqrt(np.nanmean(r ** 2))), len(lam), float(np.min(tau0))


CLEAN = [("SN2023IXF", 14), ("PTF11KLY", 17), ("SN1998S", 29), ("SN2017EGM", 30), ("SN2022WSP", 9)]
DIRTY = [("SN2017CBV", 20), ("SN2012CG", 32), ("AT2021J", 19), ("ASASSN14LP", 20), ("SN2025BCO", 5), ("SN2013DY", 17)]


def run(sn, near_day):
    z = float(cat[sn.upper()]["z"])
    best = None
    for g in ism.NUV_GRATINGS:
        for p in glob.glob(f"{OUT}/{sn}/**/{g}/{sn}_*_{g}_native.txt", recursive=True):
            ph = ism.phase_of(p, cat)
            if not np.isfinite(ph) or abs(ph - near_day) > 12:
                continue
            r = ism.analyze(*ism.load_spec(p, z), lines, n_mc=60)
            if r is None:
                continue
            rms, nfe, tau0min = cog_resid(r["inr"], r["deb"], r["b"], r["logN_fe"])
            best = (ph, r["b"], r["logN_fe"], nfe, r["fe_snr"], rms, tau0min)
    return best


print(f"{'SN':14} {'type':16} {'day':>5} {'b':>6} {'logN':>6} {'nfe':>4} {'feSNR':>6} {'cogRMS':>7} {'tau0min':>8}")
for grp, tag in ((CLEAN, "CLEAN"), (DIRTY, "DIRTY")):
    print(f"--- {tag} ---")
    for sn, d in grp:
        try:
            r = run(sn, d)
        except Exception as e:
            print(f"{sn:14} FAIL {e}"); continue
        if r is None:
            print(f"{sn:14} no passing epoch near day{d}"); continue
        ph, b, ln, nfe, snr, rms, tau0min = r
        typ = (cat[sn.upper()].get("tns_type") or "?")[:15]
        print(f"{sn:14} {typ:16} {ph:5.0f} {b:6.1f} {ln:6.2f} {nfe:4d} {snr:6.1f} {rms:7.3f} {tau0min:8.2f}")
