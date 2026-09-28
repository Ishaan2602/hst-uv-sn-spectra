# MAST discovery diff: query HST STIS+COS supernova spectroscopy, aggregate by target, and diff against
# the current raw catalog by coordinate. read-only - reports new targets, writes nothing to the catalogs.
import os, sys
import numpy as np
import pandas as pd
from astroquery.mast import Observations
from astropy.coordinates import SkyCoord
import astropy.units as u
from paths import CATALOG_RAW

print("querying MAST for HST STIS+COS supernova spectroscopy ...", flush=True)
obs = Observations.query_criteria(
    obs_collection="HST",
    instrument_name=["STIS", "COS", "STIS/CCD", "STIS/NUV-MAMA", "STIS/FUV-MAMA", "COS/FUV", "COS/NUV"],
    dataproduct_type="spectrum",
    target_classification="*upernova*",
)
df = obs.to_pandas()
print(f"raw MAST rows: {len(df)}", flush=True)
if len(df) == 0:
    sys.exit("no rows returned")

# aggregate by target position (0.02 deg tolerance clusters the same SN across programs/epochs)
df = df.dropna(subset=["s_ra", "s_dec"])
coords = SkyCoord(df["s_ra"].values * u.deg, df["s_dec"].values * u.deg)
used = np.zeros(len(df), bool)
targets = []
for i in range(len(df)):
    if used[i]:
        continue
    sep = coords[i].separation(coords).deg
    grp = sep < 0.02
    used |= grp
    sub = df[grp]
    targets.append({
        "target": str(sub["target_name"].iloc[0]),
        "ra": float(np.median(sub["s_ra"])), "dec": float(np.median(sub["s_dec"])),
        "n_spec": int(grp.sum()),
        "instruments": ",".join(sorted(set(str(x).split("/")[0] for x in sub["instrument_name"]))),
        "t_min": float(np.nanmin(sub["t_min"])) if "t_min" in sub else np.nan,
        "t_max": float(np.nanmax(sub["t_max"])) if "t_max" in sub else np.nan,
    })
T = pd.DataFrame(targets)
print(f"unique MAST supernova targets (0.02deg cluster): {len(T)}", flush=True)

cur = pd.read_csv(CATALOG_RAW)
ccoord = SkyCoord(cur["ra"].values * u.deg, cur["dec"].values * u.deg)
tcoord = SkyCoord(T["ra"].values * u.deg, T["dec"].values * u.deg)
idx, sep2d, _ = tcoord.match_to_catalog_sky(ccoord)
T["match_sep_arcsec"] = sep2d.arcsec
T["in_catalog"] = T["match_sep_arcsec"] < 5.0
new = T[~T["in_catalog"]].sort_values("t_max", ascending=False)
print(f"\ncatalog has {len(cur)} targets; MAST now has {len(T)}; NEW (not within 5 arcsec of catalog): {len(new)}")
from astropy.time import Time
for _, r in new.iterrows():
    tmax = Time(r["t_max"], format="mjd").iso[:10] if np.isfinite(r["t_max"]) else "?"
    print(f"  {r['target']:22s} ra={r['ra']:.4f} dec={r['dec']:+.4f}  n_spec={r['n_spec']:3d}  {r['instruments']:10s} last={tmax}")
new.to_csv(os.path.join(os.path.dirname(CATALOG_RAW), "_discovery_new.csv"), index=False)
print(f"\nnew-target list -> catalog/_discovery_new.csv (for review before folding in)")
