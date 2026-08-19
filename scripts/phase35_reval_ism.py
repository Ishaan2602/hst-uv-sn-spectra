#!/usr/bin/env python3
# phase-3.5 held-out re-validation (ISM anchor): ANCHOR_TAU=3.0 was set on ~10 SNe. does the clean/contaminated
# split hold on a HELD-OUT set (different SNe), with no boundary (tau0min ~ 2.5-3.5) misclassifications?
import os, sys, glob
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import ism

sys.stdout = open(os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "output", "_reval_ism.txt"), "w")

cat = ism.load_catalog(); lines = ism.load_lines(); OUT = ism.OUT
DEVISED = {"SN2023IXF", "PTF11KLY", "SN1998S", "SN2017EGM", "SN2022WSP",
           "SN2017CBV", "SN2012CG", "AT2021J", "ASASSN14LP", "SN2025BCO", "SN2013DY"}
# held-out SNe with NUV ISM coverage, NOT in the devising set (mix of expected-clean + possible-contaminated)
HELDOUT = ["SN2010JL", "SN2010AL", "SN1999EM", "SN2018BSZ", "SN2021YJA", "AT2022ACKO",
           "AT2023IUC", "SNNGC4414", "SN2024GGI", "SN2022HRS-HV", "SN2012CG", "SN2017CBV"]

print(f"{'SN':16} {'type':12} {'day':>5} {'b':>6} {'logN':>6} {'tau0min':>8} {'anchored':>9} note")
rows = []
for sn in HELDOUT:
    if sn.upper() not in cat:
        continue
    z = float(cat[sn.upper()]["z"])
    for g in ism.NUV_GRATINGS:
        for p in sorted(glob.glob(f"{OUT}/{sn}/**/{g}/{sn}_*_{g}_native.txt", recursive=True)):
            try:
                r = ism.analyze(*ism.load_spec(p, z), lines, n_mc=50)
            except Exception:
                continue
            if r is None:
                continue
            fl, ff, fy, fye = ism.collect(r["inr"], r["deb"], ion="Fe II")
            tau0 = ism.TAU_K * 10.0 ** r["logN_fe"] * ff * fl / r["b"]
            t0m = float(np.min(tau0))
            ph = ism.phase_of(p, cat)
            typ = (cat[sn.upper()].get("tns_type") or "?")[:11]
            note = "<< BOUNDARY" if 2.0 < t0m < 4.0 else ("<< high logN + anchored?" if (r["anchored"] and r["logN_fe"] > 16.0) else "")
            print(f"{sn:16} {typ:12} {ph:5.0f} {r['b']:6.1f} {r['logN_fe']:6.2f} {t0m:8.2f} {str(r['anchored']):>9} {note}")
            rows.append((sn, r["logN_fe"], t0m, r["anchored"]))
# summary: does anchored track low tau0min / sensible logN?
anc = [x for x in rows if x[3]]; unanc = [x for x in rows if not x[3]]
print(f"\nanchored:   n={len(anc)}  logN median={np.median([x[1] for x in anc]):.2f}  tau0min max={max(x[2] for x in anc):.2f}")
print(f"unanchored: n={len(unanc)} logN median={np.median([x[1] for x in unanc]):.2f} tau0min min={min((x[2] for x in unanc), default=np.nan):.2f}")
