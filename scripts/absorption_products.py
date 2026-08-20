#!/usr/bin/env python3
# per-source ABSORPTION / ISM products for the public repo (foreground metal columns + metallicity).
# combines three data sources at different depth levels:
#
#   1. AUTOMATED CoG  -- catalog/ism_cog_summary.csv + output5/<SN>/absorption/*_cog.csv
#      produced by ism.py for all SNe with NUV spectra where >=3 Fe II lines clear the 2-sigma
#      detection gate. provides b, per-ion logN, quality metrics.
#
#   2. AUTOMATED N(HI) -- catalog/lya_nhi_summary.csv
#      produced by lya_nhi.py for all SNe with an early photospheric G140L backlight.
#      enables automated [X/H] computation when combined with source 1 above.
#
#   3. CURATED DEEP   -- reference/ism_columns.csv (hand-maintained per-notebook)
#      N(HI) from careful damped Lya fits, absolute metallicities [X/H], [Fe/Zn] depletion.
#      only available for SNe where the full analysis notebook has been run.
#      curated values override automated where both exist.
#
# usage:  python absorption_products.py [SN ...]     (no args -> all SNe with any absorption data)

import os, csv, glob, re, json, datetime

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
from paths import OUT, CATALOG, ISM_SUMMARY
from paths import ABSORPTION_SUMMARY as SUMMARY
ABS_TAB  = os.path.join(ROOT, "reference", "ism_columns.csv")
NHI_SUM  = os.path.join(ROOT, "catalog",   "lya_nhi_summary.csv")

# only dominant-ionization ISM species valid for [X/H]; Mg I is a minority fraction, Al III is not cold-ISM
_SOLAR_REF = {"Fe II": -4.50, "Zn II": -7.44, "Mn II": -6.57, "Cr II": -6.36}
_ION_ERR_GATE = 0.5   # skip [X/H] for ions with logN_err >= this

# --- catalog ---------------------------------------------------------------------------------------
cat = {}
with open(CATALOG) as fh:
    for row in csv.DictReader(fh):
        cat[row["name"].upper()] = row

def _catf(sn, key):
    v = cat[sn.upper()].get(key)
    return float(v) if v not in (None, "") else 0.0

zof     = lambda sn: _catf(sn, "z")
tnstype = lambda sn: cat[sn.upper()].get("tns_type") or ""


# --- data loaders ----------------------------------------------------------------------------------

def load_nhi_summary():
    # return {SN_UPPER: nhi_row} for every SN with an automated N(HI) measurement.
    d = {}
    try:
        with open(NHI_SUM) as fh:
            for row in csv.DictReader(fh):
                d[row["sn"].upper()] = row
    except FileNotFoundError:
        pass
    return d


def load_cog_summary():
    # return {sn_dir: adopted_row} for every SN with an adopted CoG epoch.
    # keys are the SN directory names as ism.py wrote them (from os.listdir output5).
    adopted = {}
    try:
        with open(ISM_SUMMARY) as fh:
            for row in csv.DictReader(fh):
                if row.get("adopted") == "yes":
                    adopted[row["sn"]] = row
    except FileNotFoundError:
        pass
    return adopted


def load_cog_ions(sn_dir, grating, phase):
    # read per-epoch *_cog.csv for the adopted epoch; return {ion: {logN, logN_err, n_lines, limit}}.
    ismdir = os.path.join(OUT, sn_dir, "absorption")
    ph_int = int(round(float(phase)))
    for pat in [
        os.path.join(ismdir, f"{sn_dir}_{grating}_day{ph_int}_cog.csv"),
        os.path.join(ismdir, f"{sn_dir}_{grating}_day{float(phase):.1f}_cog.csv"),
    ]:
        if os.path.exists(pat):
            break
    else:
        candidates = glob.glob(os.path.join(ismdir, f"{sn_dir}_{grating}_day*_cog.csv"))
        if not candidates:
            return {}
        def _ph(p):
            m = re.search(r"_day([0-9.]+)_cog", p)
            return abs(float(m.group(1)) - float(phase)) if m else 1e9
        pat = min(candidates, key=_ph)
    ions = {}
    with open(pat) as fh:
        for row in csv.DictReader(fh):
            ions[row["ion"]] = {
                "logN":     float(row["logN"])     if row.get("logN")     else None,
                "logN_err": float(row["logN_err"]) if row.get("logN_err") else None,
                "n_lines":  int(row["n_lines"])    if row.get("n_lines")  else 0,
                "limit":    row.get("limit", "").strip() == ">",
            }
    return ions


def load_absorption_table():
    # long format: quantity is b_kms | logN(<ion>) | logN(HI) | [X/H] | [Fe/Zn]
    d = {}
    with open(ABS_TAB) as fh:
        lines = [ln for ln in fh if not ln.lstrip().startswith("#")]   # drop leading comment block
    for row in csv.DictReader(lines):
        d.setdefault(row["sn"].upper(), []).append(row)
    return d


# --- block builders --------------------------------------------------------------------------------

def _curated_block(rows):
    if not rows:
        return None
    cols, met = {}, {}
    b = nhi = dep = None
    for r in rows:
        q    = r["quantity"].strip()
        val  = float(r["value"]) if r.get("value") else None
        err  = float(r["err"])   if r.get("err")   else None
        note = r.get("note", "").strip()
        meth = r.get("method", "").strip()
        m = re.fullmatch(r"logN\((.+)\)", q)
        if q == "b_kms":
            b = val
        elif q == "logN(HI)" or (m and m.group(1) == "HI"):
            nhi = {"logN": val, "err": err, "method": meth, "note": note}
        elif q == "[Fe/Zn]":
            dep = {"value": val, "method": meth, "note": note}
        elif m:
            cols[m.group(1)] = {"logN": val, "err": err, "method": meth, "note": note}
        elif q.startswith("[") and q.endswith("]"):
            met[q] = {"value": val, "err": err, "method": meth, "note": note}
    return {"b_kms": b, "columns": cols or None, "N_HI": nhi, "metallicity": met or None, "depletion": dep}


# --- main assembler --------------------------------------------------------------------------------

def build_absorption(sn_dir, cog_row, cog_ions, curated_rows, nhi_row=None):
    """
    sn_dir       -- exact SN directory name in output5/ (used for paths + provenance)
    cog_row      -- adopted row from ism_cog_summary.csv, or None
    cog_ions     -- {ion: ...} from per-epoch *_cog.csv, or {}
    curated_rows -- list of rows from ism_columns.csv for this SN, or []
    nhi_row      -- row from lya_nhi_summary.csv, or None
    """
    has_cog     = cog_row is not None
    curated     = _curated_block(curated_rows) if curated_rows else None
    has_curated = curated is not None
    has_nhi_auto = nhi_row is not None

    if not has_cog and not has_curated and not has_nhi_auto:
        return None, None

    block = {"note": "foreground MW+host ISM screen, blended at low-res"}

    if has_cog:
        block["b_kms"]    = round(float(cog_row["b"]), 1)
        block["b_err"]    = round(float(cog_row["b_err"]), 1) if cog_row.get("b_err") else None
        block["anchored"] = True   # adopted rows are always anchored
        block["columns"]  = {
            ion: {
                "logN":     round(v["logN"],     3) if v["logN"]     is not None else None,
                "logN_err": round(v["logN_err"], 3) if v["logN_err"] is not None else None,
                "n_lines":  v["n_lines"],
                "limit":    v["limit"],
            }
            for ion, v in cog_ions.items() if v["logN"] is not None
        }

    if has_curated:
        # curated b/columns override automated (hand-verified in notebooks)
        if curated["b_kms"] is not None:
            block["b_kms"] = curated["b_kms"]
        if curated["columns"]:
            block["columns"] = curated["columns"]
        if curated["N_HI"]:
            block["N_HI"] = curated["N_HI"]
        if curated["metallicity"]:
            block["metallicity"] = curated["metallicity"]
        if curated["depletion"]:
            block["depletion"] = curated["depletion"]

    # automated N(HI) from lya_nhi.py -- add if not already set by curated
    if has_nhi_auto and "N_HI" not in block:
        block["N_HI"] = {
            "logN":             float(nhi_row["logN_HI"]),
            "logN_HI_syst_vabs": float(nhi_row["logN_HI_syst_vabs"]),
            "method":  nhi_row.get("method", "automated damped Lya"),
            "source":  "automated",
        }
        # compute [X/H] from automated N(HI) + CoG columns when ions are reliable
        cols = block.get("columns") or {}
        met = {}
        logN_HI = float(nhi_row["logN_HI"])
        for ion, vals in cols.items():
            if ion not in _SOLAR_REF:
                continue
            logN_ion = vals.get("logN")
            err_ion  = vals.get("logN_err") or vals.get("err")
            if logN_ion is None:
                continue
            if vals.get("limit"):
                continue   # lower limit -- [X/H] would be a lower bound, not meaningful
            if err_ion is not None and float(err_ion) >= _ION_ERR_GATE:
                continue
            xh = round(logN_ion - logN_HI - _SOLAR_REF[ion], 2)
            err_xh = round(float(err_ion or 0.0), 2) if err_ion else None  # does not include N_HI syst
            met[f"[{ion.split()[0]}/H]"] = {"value": xh, "err": err_xh,
                                              "method": "logN(X)-logN(HI)-solar; N(HI) automated",
                                              "note": "see logN_HI_syst_vabs for dominant N(HI) uncertainty"}
        if met:
            block["metallicity"] = met

    has_deep = has_curated and (curated.get("N_HI") or curated.get("metallicity"))
    has_nhi  = "N_HI" in block
    if has_cog and has_deep:
        tier = "cog_plus_metallicity"
    elif has_cog and has_nhi:
        tier = "cog_plus_nhi_auto"
    elif has_cog:
        tier = "cog_automated"
    elif has_curated:
        tier = "curated_only"    # curated ism_columns.csv data but no adopted CoG
    elif has_nhi_auto:
        tier = "nhi_auto_only"   # only automated N(HI) from lya_nhi.py
    else:
        tier = "curated_only"    # fallback (shouldn't reach here given the guard above)

    prov = {
        "z": zof(sn_dir),
        "method": "EW -> curve-of-growth metal columns; automated pipeline (ism.py), adopted epoch",
        "screen": "foreground MW+host ISM, blended at low-res (host-dominated where noted)",
        "cog_engine": "scripts/ism.py (per-epoch EW/CoG batch)",
        "data_tier": tier,
    }
    if has_cog:
        prov["adopted_epoch"] = {
            "grating": cog_row["grating"],
            "phase":   float(cog_row["phase"]),
            "fe_snr":  float(cog_row["fe_snr"]) if cog_row.get("fe_snr") else None,
            "n_fe":    int(cog_row["n_fe"])      if cog_row.get("n_fe")   else None,
        }
    if has_deep:
        prov["metallicity_method"] = "N(HI) from damped Lya joint fit; [X/H] vs Asplund09 solar"
    if has_nhi_auto and not has_deep:
        prov["nhi_auto_epoch"] = {
            "grating": nhi_row["grating"],
            "phase":   float(nhi_row["phase"]),
        }

    prod = {
        "sn":        sn_dir.upper(),
        "sn_type":   tnstype(sn_dir),
        "generated": datetime.date.today().isoformat(),
        "provenance": prov,
        "absorption": block,
    }

    dst = os.path.join(OUT, sn_dir, f"{sn_dir}_absorption.json")
    os.makedirs(os.path.dirname(dst), exist_ok=True)
    with open(dst, "w") as fh:
        json.dump(prod, fh, indent=2)
    return prod, dst


# --- summary row -----------------------------------------------------------------------------------

def summary_row(prod):
    ab   = prod["absorption"] or {}
    prov = prod["provenance"] or {}
    ep   = prov.get("adopted_epoch") or {}
    cols = ab.get("columns") or {}
    fe   = cols.get("Fe II") or {}
    nhi  = (ab.get("N_HI")         or {}).get("logN")
    zh   = ((ab.get("metallicity") or {}).get("[Zn/H]") or {}).get("value")
    dep  = (ab.get("depletion")    or {}).get("value")
    return {
        "sn":              prod["sn"],
        "sn_type":         prod["sn_type"],
        "z":               prov.get("z"),
        "data_tier":       prov.get("data_tier"),
        "b_kms":           ab.get("b_kms"),
        "b_err":           ab.get("b_err"),
        "logN_FeII":       fe.get("logN"),
        "logN_FeII_err":   fe.get("logN_err") or fe.get("err"),
        "n_columns":       len(cols),
        "fe_snr":          ep.get("fe_snr"),
        "adopted_grating": ep.get("grating"),
        "adopted_phase":   ep.get("phase"),
        "logN_HI":         nhi,
        "Zn_H":            zh,
        "Fe_Zn":           dep,
    }


# --- main ------------------------------------------------------------------------------------------

def main(names):
    cog_summary = load_cog_summary()       # {sn_dir: adopted_row}
    nhi_summary = load_nhi_summary()       # {SN_UPPER: nhi_row}
    tab         = load_absorption_table()  # {SN_UPPER: [curated rows]}

    # resolve all candidate SN directory names from all three sources
    dir_map = {d.upper(): d for d in os.listdir(OUT) if os.path.isdir(os.path.join(OUT, d))}

    if not names:
        cog_keys = set(cog_summary.keys())
        cur_keys = {dir_map[k] for k in tab      if k in dir_map}
        nhi_keys = {dir_map[k] for k in nhi_summary if k in dir_map}
        candidates_dirs = sorted(cog_keys | cur_keys | nhi_keys)
    else:
        candidates_dirs = []
        for n in names:
            d = dir_map.get(n.upper())
            if d:
                candidates_dirs.append(d)
            else:
                print(f"  skip {n}: directory not found in {OUT}")

    summ = []
    for sn_dir in candidates_dirs:
        if sn_dir.upper() not in cat:
            print(f"  skip {sn_dir}: not in catalog"); continue

        cog_row  = cog_summary.get(sn_dir)
        cog_ions = load_cog_ions(sn_dir, cog_row["grating"], cog_row["phase"]) if cog_row else {}
        cur_rows = tab.get(sn_dir.upper(), [])
        nhi_row  = nhi_summary.get(sn_dir.upper())

        prod, dst = build_absorption(sn_dir, cog_row, cog_ions, cur_rows, nhi_row)
        if prod is None:
            print(f"  skip {sn_dir}: no absorption data"); continue

        summ.append(summary_row(prod))
        ab   = prod["absorption"]
        tier = prod["provenance"]["data_tier"]
        print(f"  {sn_dir:16s} b={ab.get('b_kms')} cols={len(ab.get('columns') or {})} [{tier}]  -> {os.path.relpath(dst, ROOT)}")

    if summ:
        with open(SUMMARY, "w", newline="") as fh:
            wr = csv.DictWriter(fh, fieldnames=list(summ[0].keys()))
            wr.writeheader(); wr.writerows(summ)
        print(f"\nsummary -> {os.path.relpath(SUMMARY, ROOT)} ({len(summ)} sne)")


if __name__ == "__main__":
    import sys
    main(sys.argv[1:])
