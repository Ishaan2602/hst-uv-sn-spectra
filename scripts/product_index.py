# per-SN product index: scans each output/<SN>/ dir and writes <SN>_products.json listing every product
# that exists (data-analysis threads + reduced spectra) and its paths, so a repo user reads ONE file to
# know what is available for a SN instead of learning the naming grammar and globbing. also sweeps the
# casing orphans that older reruns left behind (writers overwrite by name but did not clear old-cased names).
import os, glob, json, csv, datetime, argparse
from paths import OUT


def _rel(p):
    return os.path.relpath(p, OUT).replace("\\", "/")


def _exists(sn_dir, name):
    p = os.path.join(sn_dir, name)
    return _rel(p) if os.path.exists(p) else None


def index_one(sn_dir):
    sn = os.path.basename(sn_dir)
    ab = os.path.join(sn_dir, "absorption")
    em = os.path.join(sn_dir, "emission")
    rd = os.path.join(sn_dir, "reddening")
    prods = {
        "reduced_spectra": {
            "epoch_native": sorted(_rel(p) for p in glob.glob(f"{sn_dir}/epochs/*_native.txt")),
            "epoch_resel":  sorted(_rel(p) for p in glob.glob(f"{sn_dir}/epochs/*_resel.txt")),
        },
        "emission": {
            "json":  _exists(sn_dir, f"{sn}_emission.json"),
            "plots": sorted(_rel(p) for p in glob.glob(f"{em}/*.png")),
        },
        "absorption_ism": {
            "json":     _exists(sn_dir, f"{sn}_absorption.json"),
            "cog_csv":  sorted(_rel(p) for p in glob.glob(f"{ab}/*_cog.csv")),
            "line_csv": sorted(_rel(p) for p in glob.glob(f"{ab}/*_lines.csv")),
            "plots":    sorted(_rel(p) for p in glob.glob(f"{ab}/*_cog.png") + glob.glob(f"{ab}/*_cont.png")),
        },
        "fuv_aod": {
            "json":    _exists(sn_dir, f"{sn}_fuv_aod.json"),
            "aod_csv": sorted(_rel(p) for p in glob.glob(f"{ab}/*_aod.csv")),
            "plots":   sorted(_rel(p) for p in glob.glob(f"{ab}/*_aod.png")),
        },
        "lya_nhi": {
            "json":  _exists(sn_dir, f"{sn}_lya_nhi.json"),
            "plots": sorted(_rel(p) for p in glob.glob(f"{ab}/*_lya_nhi_fit.png")),
        },
        "naid_reddening": {
            "json":  _exists(sn_dir, f"{sn}_naid_ebv.json"),
            "plots": sorted(_rel(p) for p in glob.glob(f"{rd}/*_naid_ew.png")),
        },
        "manifest": {
            "combined": _exists(sn_dir, f"{sn}_manifest.json"),
            "stis":     _exists(sn_dir, f"{sn}_stis_manifest.json"),
            "cos":      _exists(sn_dir, f"{sn}_cos_manifest.json"),
        },
        "timeseries_plot": _exists(sn_dir, f"{sn}_timeseries.png"),
    }
    # a compact "available" list: which analysis threads actually produced a product for this SN
    available = [k for k in ("emission", "absorption_ism", "fuv_aod", "lya_nhi", "naid_reddening")
                 if prods[k].get("json")]
    return {"sn": sn, "generated": datetime.date.today().isoformat(),
            "available_threads": available, "products": prods}


def sweep_orphans(sn_dir):
    # fix product JSONs whose SN-prefix casing does not match the directory name (old lowercase names left by
    # pre-uppercase reruns). the filesystem here is case-INSENSITIVE, so os.path.exists cannot tell the two
    # casings apart; use the EXACT-case names from os.listdir instead. rule: if the correctly-cased sibling
    # ALSO exists (exact case) the miscased one is a true orphan -> remove; otherwise it is the LIVE product
    # -> rename via a temp name (a direct rename to a case-only-different target is a no-op on this FS).
    sn = os.path.basename(sn_dir)
    listing = set(os.listdir(sn_dir))               # exact casing as stored on disk
    removed = []; renamed = []
    for b in list(listing):
        if not b.endswith(".json"):
            continue
        if b.startswith(sn + "_") or not b.lower().startswith(sn.lower() + "_"):
            continue                                # correctly cased, or unrelated
        correct = sn + "_" + b[len(sn) + 1:]
        if correct in listing:
            os.remove(os.path.join(sn_dir, b)); removed.append(b)     # true orphan
        else:
            tmp = os.path.join(sn_dir, b + ".casefix")
            os.rename(os.path.join(sn_dir, b), tmp)
            os.rename(tmp, os.path.join(sn_dir, correct)); renamed.append((b, correct))
    return removed, renamed


def run(names=None, sweep=True):
    dirs = sorted(d for d in os.listdir(OUT) if os.path.isdir(os.path.join(OUT, d)))
    if names:
        want = {n.upper() for n in names}
        dirs = [d for d in dirs if d.upper() in want]
    n_orphans = 0
    for d in dirs:
        sn_dir = os.path.join(OUT, d)
        if sweep:
            orph, ren = sweep_orphans(sn_dir)
            if orph:
                n_orphans += len(orph); print(f"  swept {d}: removed {orph}")
            if ren:
                print(f"  fixed {d}: renamed {ren}")
        idx = index_one(sn_dir)
        with open(os.path.join(sn_dir, f"{d}_products.json"), "w") as fh:
            json.dump(idx, fh, indent=2)
    print(f"wrote {len(dirs)} product indexes; removed {n_orphans} true orphans")


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("sn", nargs="*")
    ap.add_argument("--no-sweep", action="store_true")
    a = ap.parse_args()
    run(a.sn or None, sweep=not a.no_sweep)
