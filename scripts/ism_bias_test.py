#!/usr/bin/env python3
# phase-3.5 bias test: does anchoring b on Fe II ALONE (current) bias the columns vs a joint fit of
# ONE shared b + one logN per ion (Zimmerman's actual method)? run both on the anchored multi-ion SNe.
import os, sys, glob
import numpy as np
from scipy.optimize import least_squares
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import ism

cat = ism.load_catalog()
lines = ism.load_lines()
OUT = ism.OUT


def fit_joint(inr, deb, ions):
    # one shared b + one logN per ion, all detected lines on the universal F(tau0)
    L, F, Y, YE, IX = [], [], [], [], []
    for ii, ion in enumerate(ions):
        lam, fo, y, ye = ism.collect(inr, deb, ion=ion)
        for j in range(len(lam)):
            L.append(lam[j]); F.append(fo[j]); Y.append(y[j]); YE.append(ye[j]); IX.append(ii)
    if len(L) < len(ions) + 2:
        return None
    L, F, Y, YE, IX = map(np.array, (L, F, Y, YE, IX))
    def resid(p):
        b = p[0]; logN = p[1:]
        pred = np.array([ism.red_ew(10.0 ** logN[IX[k]], F[k], L[k], b) for k in range(len(L))])
        return (np.log10(Y) - np.log10(np.clip(pred, 1e-30, None))) / (YE / (Y * np.log(10)) + 1e-3)
    p0 = [60.0] + [14.5] * len(ions)
    res = least_squares(resid, p0, bounds=([5] + [10] * len(ions), [300] + [20] * len(ions)), max_nfev=8000)
    return res.x[0], dict(zip(ions, res.x[1:]))


def run(sn, near_day):
    z = float(cat[sn.upper()]["z"])
    for g in ism.NUV_GRATINGS:
        for p in glob.glob(f"{OUT}/{sn}/**/{g}/{sn}_*_{g}_native.txt", recursive=True):
            ph = ism.phase_of(p, cat)
            if not np.isfinite(ph) or abs(ph - near_day) > 8:
                continue
            r = ism.analyze(*ism.load_spec(p, z), lines, n_mc=60)
            if r is None:
                continue
            inr, deb = r["inr"], r["deb"]
            def ndet(ion):
                lam, _, _, _ = ism.collect(inr, deb, ion=ion); return len(lam)
            ions = [ion for ion in sorted(set(L["ion"] for L in inr))
                    if ndet(ion) >= 1 and ion != "Mg II"]  # Mg II saturated, skip
            jb = fit_joint(inr, deb, ions)
            print(f"\n== {sn} {g} d{ph:.0f}: FeII-only b={r['b']:.1f}  |  joint b={jb[0]:.1f}" if jb else f"\n== {sn} {g} d{ph:.0f}: joint failed")
            if jb:
                jbb, jN = jb
                print(f"   {'ion':6} {'FeII-anchor':>11} {'joint':>7} {'d(dex)':>7}")
                for ion in ions:
                    a = r["Ncol"][ion][0]; b2 = jN[ion]
                    print(f"   {ion:6} {a:11.3f} {b2:7.3f} {b2 - a:+7.3f}")
            return


for sn, d in [("SN2023IXF", 14), ("SN2017EGM", 30), ("SN1998S", 29), ("PTF11KLY", 17), ("SN2010JL", 9)]:
    try:
        run(sn, d)
    except Exception as ex:
        print(f"{sn}: FAIL {ex}")
