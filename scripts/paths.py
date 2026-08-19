import os

# single source of truth for where things live. bump CANONICAL to promote a newly-validated reduction
# (output5 -> output6 -> ...); every analysis script reads OUT so they all follow automatically.
# reduction scripts still take --outroot (target a fresh candidate tree each iteration) but default to
# OUT so a bare run can never land in a stale dir.

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

CANONICAL = "output5"                         # <-- the ONE line to bump when promoting a new reduction
OUT = os.path.join(ROOT, CANONICAL)           # current reduced-spectra tree; analysis READS this

# catalog + catalog-level products live outside any versioned tree so they survive tree bumps and can't
# be nuked with a stale output dir.
CATDIR = os.path.join(ROOT, "catalog")
CATALOG = os.path.join(CATDIR, "uv_sn_catalog_clean.csv")
CATALOG_RAW = os.path.join(CATDIR, "uv_sn_catalog.csv")
ISM_SUMMARY = os.path.join(CATDIR, "ism_cog_summary.csv")
XMATCH_CACHE = os.path.join(CATDIR, ".xmatch_cache")
TNS_CACHE = os.path.join(CATDIR, ".tns_cache")

# curated inputs (hand-maintained, read by scripts) already live in dedicated top-level dirs
REFERENCE = os.path.join(ROOT, "reference")
LINELISTS = os.path.join(ROOT, "linelists")
HOST_EBV = os.path.join(REFERENCE, "host_ebv.csv")     # AUTHORITATIVE host reddening (the catalog only mirrors it)


def host_ebv_map():
    # read the AUTHORITATIVE curated host reddening straight from reference/host_ebv.csv. products call this
    # instead of trusting the catalog's host_ebv column, which is a derived mirror that goes stale whenever the
    # (network-heavy) catalog_clean.py isn't re-run after an edit here. reading the source removes that failure mode.
    import csv
    out = {}
    if not os.path.exists(HOST_EBV):
        return out
    with open(HOST_EBV) as fh:
        for row in csv.DictReader(ln for ln in fh if not ln.lstrip().startswith("#")):
            try:
                out[row["name"].upper()] = (float(row["host_ebv"]),
                                            float(row["host_ebv_err"]) if row.get("host_ebv_err") else None,
                                            row.get("host_ebv_src") or "")
            except (KeyError, ValueError):
                pass
    return out
