#!/usr/bin/env python3
# per-source ABSORPTION / ISM products for the public repo (foreground metal columns + metallicity).
# assembles the curated absorption table (reference/ism_columns.csv) into a per-SN JSON sidecar:
# the b-value, the Fe/Cr/Zn/Mn columns, N(H I) from the damped Lya fit, the [X/H] metallicities and
# the [Fe/Zn] depletion. one summary CSV across the catalog.
#
# this is the ABSORPTION / ISM thread only (foreground gas absorbing the SN light). the late-time
# CSM emission (Mg II / Lya line flux + shell) lives in emission_products.py - same 2800A wavelength,
# opposite physics, kept separate on purpose. the per-epoch EW -> curve-of-growth engine is ism.py.
#
# usage:  python absorption_products.py [SN ...]     (no args -> every SN curated in ism_columns.csv)

import os, csv, re, json, datetime

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
from paths import OUT, CATALOG
ABS_TAB = os.path.join(ROOT, "reference", "ism_columns.csv")
SUMMARY = os.path.join(OUT, "absorption_summary.csv")

# --- catalog ---------------------------------------------------------------------------------------
cat = {}
with open(CATALOG) as fh:
    for row in csv.DictReader(fh):
        cat[row["name"].upper()] = row
def _catf(sn, key):
    v = cat[sn.upper()].get(key)
    return float(v) if v not in (None, "") else 0.0
zof = lambda sn: _catf(sn, "z")
tnstype = lambda sn: cat[sn.upper()].get("tns_type") or ""


def load_absorption_table():
    # long format: quantity is b_kms | logN(<ion>) | logN(HI) | [X/H] | [Fe/Zn]
    d = {}
    with open(ABS_TAB) as fh:
        lines = [ln for ln in fh if not ln.lstrip().startswith("#")]   # drop leading comment block
    for row in csv.DictReader(lines):
        d.setdefault(row["sn"].upper(), []).append(row)
    return d


def absorption_block(rows):
    if not rows:
        return None
    cols, met = {}, {}
    b = nhi = dep = None
    for r in rows:
        q = r["quantity"].strip()
        val = float(r["value"]) if r["value"] else None
        err = float(r["err"]) if r["err"] else None
        note = r.get("note", "").strip()
        meth = r.get("method", "").strip()
        m = re.fullmatch(r"logN\((.+)\)", q)
        if q == "b_kms":
            b = val
        elif q == "logN(HI)":
            nhi = {"logN": val, "err": err, "method": meth, "note": note}
        elif q == "[Fe/Zn]":
            dep = {"value": val, "method": meth, "note": note}
        elif m and m.group(1) == "HI":
            nhi = {"logN": val, "err": err, "method": meth, "note": note}
        elif m:
            cols[m.group(1)] = {"logN": val, "err": err, "method": meth, "note": note}
        elif q.startswith("[") and q.endswith("]"):
            met[q] = {"value": val, "err": err, "method": meth, "note": note}
    return {
        "note": "foreground MW+host ISM screen, blended at low-res (see docs/analysis_phase3.md sec 3)",
        "b_kms": b, "columns": cols, "N_HI": nhi, "metallicity": met, "depletion": dep,
    }


# --- assemble + write ------------------------------------------------------------------------------
def build_absorption(sn, tab):
    block = absorption_block(tab.get(sn.upper(), []))
    if block is None:
        return None, None
    prod = {
        "sn": sn.upper(),
        "sn_type": tnstype(sn),
        "generated": datetime.date.today().isoformat(),
        "provenance": {
            "z": zof(sn),
            "method": "EW -> curve-of-growth metal columns (early-epoch anchor); N(HI) from the damped "
                      "Lya joint emission+absorption fit; [X/H] vs Asplund09 solar. see docs/analysis_phase3.md sec 3.",
            "screen": "foreground MW+host ISM, blended at low-res (host-dominated where noted)",
            "cog_engine": "scripts/ism.py (per-epoch EW/CoG batch)",
        },
        "absorption": block,
    }
    dst = os.path.join(OUT, sn, f"{sn}_absorption.json")
    os.makedirs(os.path.dirname(dst), exist_ok=True)
    with open(dst, "w") as fh:
        json.dump(prod, fh, indent=2)
    return prod, dst


def summary_row(prod):
    ab = prod["absorption"] or {}
    nhi = (ab.get("N_HI") or {}).get("logN")
    zh = ((ab.get("metallicity") or {}).get("[Zn/H]") or {}).get("value")
    dep = (ab.get("depletion") or {}).get("value")
    return {
        "sn": prod["sn"], "sn_type": prod["sn_type"], "z": prod["provenance"]["z"],
        "b_kms": ab.get("b_kms"), "n_columns": len(ab.get("columns") or {}),
        "logN_HI": nhi, "Zn_H": zh, "Fe_Zn": dep,
    }


def main(names):
    tab = load_absorption_table()
    if not names:
        names = sorted(tab.keys())
    summ = []
    for sn in names:
        if sn.upper() not in cat:
            print(f"  skip {sn}: not in catalog"); continue
        prod, dst = build_absorption(sn, tab)
        if prod is None:
            print(f"  skip {sn}: no curated absorption rows"); continue
        summ.append(summary_row(prod))
        ab = prod["absorption"]
        print(f"  {sn:14s} b={ab.get('b_kms')} cols={len(ab.get('columns') or {})}  -> {os.path.relpath(dst, ROOT)}")
    if summ:
        with open(SUMMARY, "w", newline="") as fh:
            wr = csv.DictWriter(fh, fieldnames=list(summ[0].keys()))
            wr.writeheader(); wr.writerows(summ)
        print(f"\nsummary -> {os.path.relpath(SUMMARY, ROOT)} ({len(summ)} sne)")


if __name__ == "__main__":
    import sys
    main(sys.argv[1:])
