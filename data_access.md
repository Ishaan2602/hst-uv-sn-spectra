# Data Access & Software Tools

How HST UV spectra get from the archive to our notebook, and how all the
acronyms relate. Written so future-me stops being confused about MAST vs HSLA
vs HASP vs astroquery.

## The one thing to remember
Everything lives in **MAST** (Mikulski Archive for Space Telescopes). MAST Portal,
the HST Search Form, astroquery, HASP, HSLA are NOT separate archives. They are
either different doors into MAST, or different processing levels of the same data.

So there are really two separate questions:
1. HOW do I fetch the files? (Portal / HST Search Form / astroquery)
2. WHAT processing level do I want? (raw 2D / pipeline x1d / HASP coadd / HSLA abutted)

---

## 1. Access tools (the "how do I fetch it" doors)

| Tool | What it is | Best for |
|------|-----------|----------|
| MAST Portal (`mast.stsci.edu`) | Interactive web GUI. Advanced search, set Project=HSLA to filter to legacy products. | Eyeballing a handful of targets, or one class of object. |
| HST Search Form in MAST (`mast.stsci.edu/search/ui/#hst`) | HST-specific web search form. HASP `cspec.fits` and HSLA `cspec`/`aspec` show up here as default downloadable products. | A handful of targets when you want the higher-level coadds too. This is the "archive search" form. |
| `astroquery.mast` (`Observations`) | Python API into MAST. HASP/HSLA products are in MAST so they come through here too. | Scripted, reproducible, bulk downloads. This is our main path. |

All three hit the same MAST back end. Manual download from the Portal/Search Form
and a scripted `astroquery` pull return the same files. We prefer `astroquery` for
reproducibility, but manual download is fine for one-off grabbing.

---

## 2. Processing levels (the "what do I actually want" ladder)

From least to most processed:

1. **Raw / 2D frames**
   - STIS: `RAW` -> `FLT` (flat-fielded 2D) -> `CRJ` (cosmic-ray-rejected 2D).
   - COS: `rawtag` / `corrtag` TIME-TAG photon event lists (COS is photon-counting,
     not a simple 2D image).
2. **Pipeline 1D extraction (`x1d`)**
   - STIS pipeline = `calstis`, COS pipeline = `calcos`. Both spit out an `x1d`
     (1D wavelength/flux/error) using DEFAULT extraction params and DEFAULT
     background regions.
   - This is what we grabbed on Thursday for SN2020fqv. It is MAST's default
     reduction, not ours.
3. **HASP coadds (`cspec.fits`)** — see below.
4. **HSLA abutted spectra (`aspec.fits`) + classification** — see below.

---

## 3. HASP (Hubble Advanced Spectral Products)

- Automated coaddition pipeline that combines the pipeline `x1d` files into
  higher quality 1D spectra, for nearly every COS and STIS spectrum in MAST
  (3200+ programs, 64000+ datasets).
- Two product types, each at visit level and program level:
  - **coadd**: combine exposures from the SAME grating.
  - **abutment**: stitch together DIFFERENT gratings / instruments.
- Output naming: `hst_<PID>_<instrument>_<target>_<opt_elem>_<ippp>_cspec.fits`
  (instrument is one of COS, STIS, or COS-STIS).
- Public Python coadd script + Jupyter notebooks let you run CUSTOM coadds
  (Setup, CoaddTutorial, FluxScaleTutorial, DataDiagnostic, WavelengthAdjustment).
- **The caveat that matters for us:** HASP coadds are built from the STANDARD
  pipeline `x1d`. It does NOT redo the trace, background, defringing, or cosmic-ray
  handling. "Users must generate custom coadds to address these issues."

## 4. HSLA (Hubble Spectroscopic Legacy Archive)

- Sits ON TOP of HASP. Two extra steps:
  1. Combines data across MULTIPLE programs/instruments/gratings per target into
     the highest-SNR product, splicing COS+STIS into a single "abutted spectrum"
     (`aspec.fits`) spanning the full wavelength range.
  2. Automatic target classification (star/galaxy/etc.) via SIMBAD/NED, so you can
     search by object type, not just name.
- History: original HSLA (2016-2018, COS only, Peeples et al. 2017) needed a lot of
  manual work. The new HSLA is fully automated and re-updates when a target is
  re-observed or recalibrated.
- Same caveat as HASP: derived from standard pipeline `x1d`, no custom reduction.
  Abutted spectra also have discontinuous resolution / S/N across the splice points,
  and COS LSF varies with cenwave + lifetime position.

---

## 5. What this means for OUR project

The whole point (per PI) is to do our OWN uniform reduction: custom trace, custom
background regions, our own cosmic-ray rejection. HASP and HSLA both start from the
default pipeline `x1d`, so their coadds inherit whatever the default trace and
background got wrong. For faint SNe sitting in extended host galaxies, the default
trace/background can be plain wrong, and SNe are variable so HASP/HSLA auto
associations and classifications are not reliable for us.

So we use them like this:
- **Discovery / quick look:** HSLA + Portal to see what exists for a target and
  sanity-check shapes.
- **Reuse the coadd code:** once WE produce custom `x1d` files (our trace/bg/CR),
  we can feed them to the HASP coadd script to do the combine step, instead of
  writing coaddition from scratch. The `combine_spectra_reference.ipynb` from the
  PI does the equivalent with `specutils` resampling.
- **Validation:** compare our reduced spectrum against the HASP/HSLA product to
  confirm we did not break flux calibration.
- We do NOT ship HASP/HSLA default coadds as our final science products.

---

## 6. STIS vs COS extraction (they are NOT the same)

This trips people up. "Do our own extraction" means different things per instrument.

**STIS (CCD gratings G230LB/G430L/G750L, MAMA G140L/G230L):**
- 2D image data. Workflow we already have in `1dspectrumreference.ipynb`:
  1. `stistools.ocrreject` on the `FLT` -> a cleaner `CRJ` (cosmic-ray reject).
  2. `stistools.x1d` on the `CRJ` with custom `extrsize`, `bk1offst`, `bk2offst`,
     and optionally a manual center. This puts the extraction aperture on the
     trace and the background apertures off to the sides.
- This is the red/orange-lines picture the PI shared: solid lines = extraction apertures OR bg apertures,
  dashed line = center of trace.

**COS (FUV G130M/G160M/G140L, NUV G185M/G225M/G285M/G230L):**
- Photon-counting TIME-TAG data, no simple 2D trace like STIS CCD.
- Pipeline = `calcos`. Custom extraction = edit the BOXCAR extraction boxes in the
  `XTRACTAB` reference table (COS "Extract" notebook), then re-run `calcos`.
- Other COS-specific knobs: `SplitTag` (split TIME-TAG into sub-exposures),
  DayNight filtering, association (`asn`) file editing.
- So COS "full extraction" is a `calcos` + `XTRACTAB` job, not a `stistools.x1d`
  job. Plan COS work separately from STIS.

---

## 7. Running the reduction needs Linux/WSL (what happened 6/21)

When we tried to actually run our own STIS extraction on native Windows, it hit a
hard wall. Logging it here so we do not repeat the dead end.

**What we did (sandbox, SN2024iss of8b02010, G230LB):**
1. `pip install stistools crds` into the Windows Python 3.14 kernel.
2. astroquery: resolved SN2024iss by name, downloaded `flt`, `crj`, `sx1` into
   `./Data/2024iss/...`.
3. Plotted the pipeline `sx1` (default 1D) fine.
4. Set `CRDS_PATH`/`CRDS_SERVER_URL`, pre-created `crds_cache/{mappings,references,config}/hst`
   (CRDS errors on a fresh Windows cache if those dirs do not exist), then
   `crds.assign_bestrefs([...], sync_references=True)` pulled 21 STIS reference
   files into `crds_cache/references/hst/stis/` and we pointed `oref` at that dir.
5. Called `stistools.ocrreject.ocrreject(flt, ...)`.

**What happened:** `FileNotFoundError: [WinError 2] The system cannot find the file
specified`. `stistools` does not do the math in Python, it shells out to the
compiled HSTCAL C binaries (`cs2.e` = ocrreject, `cs6.e` = x1d, `cs0.e` = calstis).
Those binaries are not on the machine.

**Why:** HSTCAL (and COS `calcos`) are compiled C programs distributed through
conda-forge, and historically built for Linux/macOS only. There is no pip wheel
and effectively no native Windows build. So `ocrreject`, `x1d`, `calstis`, and
`calcos` cannot run on native Windows, no matter what we pip install. The only
Windows-side binary present is `ocrreject_exam.exe`, which is a pure-Python
diagnostic, not the calibration step.

**What this means, split by where it can run:**
- Works on native Windows now: astroquery search/download, reading `x1d`/`sx1`,
  DQ bitmask filtering, plotting, the 2D trace + extraction-region visualization
  (`show_extraction_regions`), and specutils coaddition.
- Needs Linux: `ocrreject`, `x1d` re-extraction, `calstis`, and all COS `calcos`
  custom extraction. This is the part the PI actually cares about (our own trace,
  background, CR rejection).

**Paths forward (needs a decision):**
- WSL: this machine already has WSL with Debian and Ubuntu. Install `stenv` or
  `conda install -c conda-forge hstcal` inside the WSL distro, mount the project
  at `/mnt/c/Users/eluru/CALTECH_SURF_2026/WORK`, run the compiled steps there.
- Remote Linux: run the reduction on an astro server over SSH, pull products back.
- Native Windows conda: not viable, conda-forge `hstcal` has no win-64 build.

## 8. Target naming in MAST (why the cone search happened)

- MAST `target_name` is whatever the PROPOSER typed. SN2020fqv was logged as
  `TESS-SN`, which is why `query_criteria(target_name="SN2020fqv")` returned zero
  and we fell back to a coordinate cone search.
- Cleaner fix: use the name resolver instead of a hand-typed cone search. Either
  `Observations.query_criteria(objectname="SN 2020fqv", radius="...", ...)` or
  `Mast.resolve_object("SN 2020fqv")` to get coords. The resolver runs the IAU /
  SIMBAD / NED lookup for us, then does the cone search under the hood. (To be
  tested in the sandbox.)

---

## 9. Useful links
- HST Notebooks index: https://spacetelescope.github.io/hst_notebooks/index.html
- HASP: https://archive.stsci.edu/missions-and-data/hst/hasp
- HSLA: https://archive.stsci.edu/missions-and-data/hst/hsla
- COS notebooks: https://spacetelescope.github.io/hst_notebooks/notebooks/COS/README.html
- STIS notebooks (see `extraction` / `1D_Extraction`): https://spacetelescope.github.io/hst_notebooks/notebooks/STIS/README.html
- HASP coadd tutorial: https://spacetelescope.github.io/hst_notebooks/notebooks/HASP/CoaddTutorial/CoaddTutorial.html
