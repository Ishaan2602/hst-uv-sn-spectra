#!/usr/bin/env python3
# phase-3.5 held-out re-validation (EMISSION): confirm the fixes generalize to sources they were NOT
# devised on. the model selection was devised on a 6-source montage; here we look at the WHOLE catalog.
import json, glob, collections, os, sys
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from paths import OUT

DEVISED = {"SN2023IXF", "SN2024GGI", "SN1998S", "SN2010JL", "SN-2005IP", "PTF11KLY", "SN2017EGM"}
bm = collections.Counter()
by_type = collections.defaultdict(collections.Counter)
ferr_ratio = []           # systematic flux_err / flux for reliable epochs
heldout_shape = []        # (sn, type, phase, best_model) for reliable clean-emission epochs NOT in DEVISED
pcyg_vphot = []

for p in glob.glob(os.path.join(OUT, "*", "*_emission.json")):
    d = json.load(open(p)); sn = d["sn"]; typ = (d["sn_type"] or "?")
    for r in d["emission"]["mg2"]:
        b = r.get("best_model"); bm[b] += 1
        if r.get("flux_reliable") and r.get("flux_err") and r["flux"] > 0:
            ferr_ratio.append(r["flux_err"] / r["flux"])
        if r.get("pcygni") and r.get("v_phot_kms") is not None:
            pcyg_vphot.append(r["v_phot_kms"])
        if b is not None and r.get("flux_reliable") and sn not in DEVISED:
            heldout_shape.append((sn, typ, r["phase"], b))

def _t(typ):
    t = typ.upper()
    if "IIN" in t: return "IIn"
    if "II" in t: return "II"
    if "IA" in t: return "Ia"
    if "IB" in t or "IC" in t: return "Ibc"
    if "SLSN" in t: return "SLSN"
    return "other"

by_t = collections.defaultdict(collections.Counter)
for p in glob.glob(os.path.join(OUT, "*", "*_emission.json")):
    d = json.load(open(p)); typ = _t(d["sn_type"] or "?")
    for r in d["emission"]["mg2"]:
        if r.get("best_model"): by_t[typ][r["best_model"]] += 1

print("=== catalog-wide best_model distribution ===")
print(dict(bm))
print("\n=== best_model by SN type (clean-emission epochs) ===")
for t in sorted(by_t):
    print(f"  {t:6}: {dict(by_t[t])}")
print("\n=== systematic flux_err / flux (reliable epochs) ===")
fr = np.array(ferr_ratio)
print(f"  n={len(fr)}  median={np.median(fr):.2f}  16-84%={np.percentile(fr,16):.2f}-{np.percentile(fr,84):.2f}  (photon-only was ~0.01-0.03)")
print("\n=== HELD-OUT clean-emission shapes (SNe NOT in the 6-regime devising set) ===")
for sn, typ, ph, b in sorted(heldout_shape):
    print(f"  {sn:16} {typ:14} d{ph:.0f}  -> {b}")
print(f"\n=== pcygni v_phot (km/s) n={len(pcyg_vphot)}  median={np.median(pcyg_vphot):.0f}  range={min(pcyg_vphot)}..{max(pcyg_vphot)}")
