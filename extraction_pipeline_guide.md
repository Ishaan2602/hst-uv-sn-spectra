# HST UV Supernova Extraction Pipeline - Master Guide

A single writeup of everything built for this project so far: the environment, the
reduction workflows (STIS + COS), defringing, coadds, the supernova discovery query,
the full end-to-end 2023ixf example, and the first automated pipeline. Written so you
can review the whole thing in one place and reproduce any piece. Normal day-to-day
status still lives in `project_state.md`, `logistics_notes.md`, `data_access.md`, and
`science_notes.md`; this file is the consolidated convenience copy.

---

## 0. What we're doing

Goal: a uniformly reduced repository of every supernova UV spectrum HST has taken, across
all the UV spectrographs (STIS and COS), reducing the data ourselves (our own trace,
background, cosmic-ray handling, defringing) instead of trusting MAST's default pipeline
products. The science driver is early-time/flash UV spectroscopy of core-collapse SNe.

This session got the full toolchain working end to end: discover targets -> download ->
reduce (STIS and COS, our own extraction) -> defringe -> coadd, plus a first automated
driver that chains it all.

---

## 1. Environment (the most important thing to understand)

The reduction tools are split across two places, and this is not optional.

**Native Windows, global Python 3.14.2** (`c:\Users\eluru\AppData\Local\Python\pythoncore-3.14-64\python.exe`)
- We do NOT use virtual environments. The global interpreter has astroquery, astropy,
  specutils, numpy, scipy, matplotlib, crds, stistools, nbformat.
- This runs everything that is pure Python: MAST queries, downloads, reading FITS,
  plotting, coadding, the discovery query.
- A `.venv` will silently break notebook kernels (it gets auto-selected and is empty of
  packages). If one appears, delete it. All notebooks use the global 3.14.2 kernel.

**WSL Debian, conda env `surf_uv`** (`~/miniforge3/envs/surf_uv`)
- Why: `stistools.x1d` / `ocrreject` / `defringe` and COS `calcos` are thin wrappers around
  the compiled HSTCAL C binaries (`cs0.e` .. `cs12.e`). HSTCAL has no native Windows build,
  only conda-forge for Linux/mac. So the actual reduction cannot run on native Windows.
- The env has: `hstcal` (conda-forge) + `calcos` (pip) + `stistools` + `crds` + the astropy stack.
- The project is visible from WSL at `/mnt/c/Users/eluru/CALTECH_SURF_2026/WORK`.

**Invoking WSL from here** (the pattern every reduction uses):
```
wsl.exe -d Debian -- bash -lc 'source ~/miniforge3/etc/profile.d/conda.sh && conda activate surf_uv && export oref=... && python /mnt/c/.../WORK/<script>.py <args>'
```

**WSL gotcha (cost me real time):** shell variables and loops (`$OBS`, `$e`, `for ... do`)
get expanded to EMPTY crossing Git Bash -> `wsl.exe`. So inside these one-liners, use
LITERAL paths only and do any looping in Python, not the shell. That is exactly why the
sweep drivers (`reduce_stis_sweep.py`, `reduce_cos_sweep.py`) loop in Python and take literal
args. `export VAR=literal` is fine (it is read internally by crds/calcos); just never
reference `$VAR` later in the same command text.

---

## 2. CRDS reference files

Calibration reference files come from CRDS and are resolved by env vars: STIS uses `oref`,
COS uses `lref`.

- Cache lives at `WORK/crds_cache`. STIS refs land in `crds_cache/references/hst/stis/`,
  COS in `crds_cache/references/hst/cos/`. Point `oref`/`lref` at those dirs (trailing slash).
- The cache is on the Windows side and reused fine from WSL over `/mnt/c`.
- Fresh-cache gotcha: CRDS errors writing files before it creates `mappings/hst`,
  `references/hst`, `config/hst`. Pre-create them (the scripts/`mkdir -p` do this).
- To sync refs for a dataset: `crds.assign_bestrefs([files], sync_references=True)` (Python),
  or `crds bestrefs --files ... --sync-references=1` (CLI).
- **Old-data gotcha:** for older observations the header ref pointers are stale. You must add
  `--update-bestrefs` (CLI) or use `assign_bestrefs` so the headers point at the synced
  current refs, or calcos/x1d will fail with "file not found". Hit this on the 2010 COS FUV data.

---

## 3. Data access and target discovery

**The landscape (full version in `data_access.md`):** MAST is the archive. The MAST Portal,
the HST Search Form, and `astroquery` are just different doors into it. HASP and HSLA are
higher-level *products*: HASP auto-coadds the pipeline `x1d`; HSLA combines across programs
and adds classification on top of HASP. Both start from the DEFAULT pipeline `x1d`, so they
inherit whatever the default trace/background got wrong, which is why we do our own.

**The naming problem:** `target_name` in MAST is whatever the proposer typed. SN2020fqv was
logged as `TESS-SN`; a name query misses it entirely. Use the resolver
(`Observations.query_criteria(objectname="SN2024iss", radius="0.05 deg", ...)`) for a single
known target.

**The discovery query (the whole point, in `uv_sn_query.ipynb`):** to find ALL UV SNe without
paper/abstract searches, query on `target_classification`, which holds the proposer's Phase II
object class (e.g. `EXT-STAR;SUPERNOVA TYPE IA`). Filtering HST spectra on
`target_classification='*upernova*'` across the STIS + COS instruments, then deduplicating by
coordinate (within ~5 arcsec), gives **140 unique HST UV supernovae** from ~1591 spectra,
written to `uv_sn_catalog.csv`. This catches SN2020fqv-as-TESS-SN that name-matching never
would. Top targets: SN1987A (203 spectra), SN2023ixf (133, COS+STIS), SN2010JL (72),
PTF11kly/SN2011fe (58), SN2014J (45). A future belt-and-suspenders improvement is a
SIMBAD/TNS coordinate cross-match.

**Product naming quirk:** STIS CCD data taken with CRSPLIT>1 gives `crj` (CR-rejected) and
`sx1` (1D); MAMA data gives `x1d`. Single CRSPLIT=1 exposures have no `crj` (see 5.).

---

## 4. STIS reduction

STIS CCD gratings: G230LB (NUV/blue), G430L (optical), G750L (red). The 2D image data flows
flt (flat-fielded) -> crj (CR-rejected, if multiple subexposures) -> x1d (1D extraction).

**Our workflow:**
1. `ocrreject(flt -> crj)` to combine CRSPLIT subexposures and reject cosmic rays. ONLY works
   if there are >=2 subexposures. Single CRSPLIT=1 exposures skip this and extract straight
   off the flt (x1d accepts flt, sfl, or crj).
2. `x1d(crj or flt -> x1d)` for the 1D extraction. Knobs: `extrsize` (aperture height),
   `a2center` (manual trace center), `maxsrch` (how far to search for the center; 0 = don't
   search, use a2center), `bk1offst`/`bk2offst` (background region offsets from the trace),
   `bk1size`/`bk2size` (background region heights).
3. Every output must be deleted before re-running; the pipeline will not overwrite.

**Helpers:** `reduce_stis.py` (single parametrized extraction), `reduce_stis_sweep.py` (loops
a set of params in Python and writes one x1d per setting).

**The extraction-region visualization** (the green-lines picture the PI shared): plot the 2D
flt with the extraction aperture (from `EXTRLOCY`/`EXTRSIZE`) and the two background regions
(from `BK1OFFST`/`BK2OFFST` + sizes). Helper `show_extraction_regions` in `stis_sandbox.ipynb`
and `1dspectrumreference.ipynb`.

**Findings from the interactive param play (SN2024iss G230LB, of8b02010):**
- The default extraction puts both background regions at offset -300/-320, i.e. both ~300px
  below the trace, because A2CENTER=894 is near the top detector edge (1024) and there is no
  room above. A tight zoom around the trace hides them. Moving them close (offset +/-14) made
  them visible straddling the trace, and barely changed the spectrum: for an isolated trace on
  a low uniform background, background placement is minor. It matters for traces on a
  structured/contaminated background (host-galaxy gradient, nearby source).
- `extrsize` sweep (1,3,5,7,9,11): extrsize=1 loses ~half the flux (only the trace peak),
  extrsize=3 already recovers nearly all of it (the trace is ~7px tall), and beyond ~7 you just
  add sky/noise. Sweet spot ext3-5.
- `a2center` centering: with maxsrch=0, miscentering by 4px drops the flux ~60-70% because the
  aperture barely overlaps the real trace. Getting the center right (or letting maxsrch find it)
  matters a lot.

---

## 5. COS reduction

COS is photon-counting TIME-TAG data (rawtag/corrtag event lists), processed by `calcos`
(pure Python, pip-installed into the env). It is NOT a simple 2D trace like STIS CCD.

**Run calcos:** on the `_asn.fits` to get a combined `x1dsum`, or on a single rawtag/corrtag to
get a per-exposure `x1d`. For FUV pass only the `rawtag_a`; calcos finds segment B. Output dir
must be empty (clear it first). Helper `reduce_cos.py`.

**Custom extraction = the XTRACTAB.** The `_1dx.fits` extraction table holds, per
(segment, aperture, grating, cenwave): `B_SPEC` (spectrum box center), `HEIGHT` (box height,
should be odd), `B_BKG1`/`B_BKG2` (background centers), and the bg heights. NUV has 3 stripes
(NUVA/NUVB/NUVC, shared background); FUV has FUVA/FUVB with per-segment background.
To customize: copy the XTRACTAB, edit the matching rows, write a new `_1dx.fits`, point the
rawtag/corrtag header `XTRACTAB` keyword at it, re-run calcos. Helpers `reduce_cos_custom.py`,
`reduce_cos_sweep.py` (NUV), `reduce_cos_fuv_sweep.py` (FUV).

**NUV vs FUV:** NUV uses BOXCAR extraction by default, so custom boxes just work. FUV defaults
to TWOZONE; to use boxcar custom boxes you must set the rawtag header `XTRCTALG=BOXCAR`,
`TRCECORR=OMIT`, `ALGNCORR=OMIT`.

**Two FUV gotchas:**
- Old FUV data (2010) needs `--update-bestrefs` so calcos finds the current refs (section 2).
- The XTRACTAB `B_SPEC` lives in the geometrically-corrected `YFULL` frame, not the raw `RAWY`.
  Plotting the boxes over `RAWY` puts them ~13-37px off the actual peaks; collapse the corrtag
  `YFULL` to align boxes with the trace.

**Findings from the interactive height sweeps:**
- NUV (SN2023ixf G230L, faint source): SNR peaks at box height ~9, then declines (the default
  57 is far too wide for a faint source). The narrow box roughly halves the flux scatter while
  preserving the continuum and the Mg II ~2800 emission. A ~25% SNR gain over default.
- FUV (SN2010jl G130M, bright source): SNR rises with height and plateaus by ~35-45 (the default
  41 is near-optimal). Opposite behavior.
- The lesson: optimal box height is source-brightness dependent. Faint sources are
  background-limited (narrow box to kill background); bright sources are source-limited (wide
  enough to catch the full profile). The default is tuned for bright point sources, so faint SNe
  benefit most from a custom narrow box.

**Extraction visualization:** collapse the corrtag counts onto the cross-dispersion axis to see
the stripes/segments, overlay the default vs custom boxes. In `cos_sandbox.ipynb`.

---

## 6. Defringing (G750L)

STIS G750L has fringing (fixed-pattern interference) at red wavelengths (worst past ~9000A). The
fix divides by a contemporaneous tungsten fringe flat (a CCDFLAT exposure taken in the same visit
right next to the science). For SN2024iss the science was of8b02040 and the flat of8b02050 (the
`_050` CCDFLAT in each visit is the fringe flat). For SN2023ixf: science oezt010e0, flat oezt010d0.

**Chain** (`reduce_defringe.py`, WSL):
1. `normspflat(flat_raw -> nsp, do_cal=True, wavecal=sci_wav)` - calibrate + normalize the flat.
2. `mkfringeflat(sci -> frr)` - shift and scale the normalized flat to match the science fringes.
3. `defringe(sci -> drj)` - divide it out. Output suffix is `_drj` for G750L regardless of whether
   the input is crj or flt.
4. `x1d(drj -> x1d)` - extract the defringed frame.

**Caveats:**
- Widen the shift search (`beg_shift=-2, end_shift=1, shift_step=0.1`); the default range hit its
  edge. After widening, the shift converged (best ~ -0.46px).
- The scale search can still hit the edge of its range (best scale 1.2 on the 2023ixf flat). The
  shift is the more important term; the scale-edge is a minor refinement worth revisiting if the
  red end needs to be pristine.
- Result: fringe oscillations past ~9000A are visibly reduced while the spectral shape is
  preserved. Residuals remain (normal for G750L).

---

## 7. Coadds

Combine spectra from different gratings/exposures onto one common wavelength axis, following the
PI's `combine_spectra_reference.ipynb` approach: resample each onto a common axis, then median.

- Common axis is finer in the UV, coarser to the red:
  `concatenate([arange(1650,3050,1.4), arange(3050,5600,2.7), arange(5600,10260,4.9)])`.
- Resampling: the PI uses specutils `FluxConservingResampler`. It is strict about a strictly
  increasing/decreasing axis and chokes on duplicate/non-monotonic wavelengths at grating splices.
  We fall back to `np.interp` (with `left=nan, right=nan`) after sorting + dropping non-increasing
  points (`clean()` helper). For a display coadd this is fine; flux-conserving rebinning is a later
  refinement.
- Combine with `nanmedian` across the stack: each component contributes where it has coverage,
  and overlaps get averaged.

**Caveats (real findings):**
- Do NOT median a high-SNR and a near-noise spectrum equally. In the full 2023ixf coadd the faint
  COS NUV fully overlaps STIS G230LB but is much noisier; a naive median dragged the combined NUV
  down to half the G230LB flux. Fix: combine the gratings that have good SNR (the 3 STIS gratings
  tile 1670-10250A) and keep the faint COS NUV as an overlay. A proper fix is SNR weighting.
- Grating flux offsets: the gratings don't always agree perfectly in their overlaps, so the median
  shows small steps at splices (~3000A). HASP/HSLA fix this with scaling; we don't yet.
- Our COS NUV coadd of the 6 SN2023ixf exposures matches the HASP coadd well except a ~2150-2350A
  region where a contaminated NUVA stripe leaks in. HASP's flux-checking rejects it; we don't yet.
  That discrepancy is a good illustration of what HASP's extra filtering buys.

STIS coadd in `stis_sandbox.ipynb`, COS coadd in `cos_sandbox.ipynb`.

---

## 8. Full 2023ixf reduction (end to end)

`full_2023ixf.ipynb` + `run_full_2023ixf.py`. One STIS epoch of SN2023ixf (prop 17205, the oezt01
visit): G230LB (oezt01040), G430L (oezt010h0), G750L (oezt010e0), and the contemporaneous CCDFLAT
(oezt010d0), combined with the COS NUV G230L spectrum.

- The notebook downloads flt/raw/wav, fires the WSL reduction, then loads the x1d files and builds
  the coadd.
- `run_full_2023ixf.py` (WSL) syncs the refs, runs x1d on the blue gratings, and does
  normspflat -> mkfringeflat -> defringe -> x1d for G750L.
- These are single CRSPLIT=1 exposures, so NO ocrreject; we extract straight off the flt.
- Result: a clean full HST UV-optical SED from 1670 to 10250A, our own extraction throughout, with
  the COS NUV overlaid (faint at this epoch, excluded from the median per section 7).

---

## 9. The automated pipeline (first pass)

`pipeline.py` (Windows orchestrator) + `reduce_epoch.py` (generalized WSL reducer). Ties the
session's pieces into one driver: discover -> pick an epoch -> download -> reduce -> coadd.

`pipeline.py` functions:
- `find_uv_sne()` - the classification query from section 3, deduped to unique SNe.
- `pick_stis_ccd_epoch(target)` - auto-finds a single STIS/CCD visit that has all three gratings
  (and the G750L CCDFLAT if present). Validated: it independently found the oezt03 visit for
  SN2023ixf. If a visit lacks a flat it degrades gracefully (G750L extracted without defringe).
- `fetch(ids, target)` - downloads flt/raw/wav for the chosen exposures.
- `reduce(ids, target)` - calls `reduce_epoch.py` in WSL via subprocess.
- `coadd(ids, target)` - combines the gratings and saves `output/<target>_coadd.csv` + `.png`.
- `run(target, ids=None)` - the whole chain; `ids` optional (auto-picks if omitted).

`reduce_epoch.py` (WSL): generalized epoch reducer. Takes a base dir and `grating=obsid` pairs
(any subset), syncs refs, x1d's each grating off the flt, and defringes G750L if a flat is given.

Run the demo: `python pipeline.py` (does SN2023ixf end to end, writes `output/SN2023ixf_coadd.*`).

**Limitations (it is a first pass):** the epoch auto-picker is simple (first visit with all 3
gratings; flat detection is by target name containing FLAT), the coadd is STIS-only and uses a
plain median (no SNR weighting or inter-grating flux scaling), and COS is not yet wired into the
auto pipeline (it is demonstrated separately in `cos_sandbox`). Natural next steps: SNR-weighted
coadd, inter-grating scaling, COS leg in the pipeline, SIMBAD/TNS cross-match in discovery.

---

## 10. File inventory

Notebooks:
- `uv_sn_query.ipynb` - the supernova discovery query -> `uv_sn_catalog.csv`.
- `stis_sandbox.ipynb` - STIS extraction, param sweeps, STIS coadd, G750L defringe comparison.
- `cos_sandbox.ipynb` - COS NUV + FUV custom extraction, height sweeps, extraction viz, COS coadd.
- `full_2023ixf.ipynb` - the end-to-end full reduction of one 2023ixf epoch.
- `1dspectrumreference.ipynb`, `combine_spectra_reference.ipynb` - the PI's reference notebooks.
- `spectra_sandbox.ipynb` - the original combined sandbox (superseded by the split ones; deletable).

Scripts (all the `reduce_*` and the pipeline run in WSL surf_uv, except pipeline.py which is the
Windows orchestrator):
- `reduce_stis.py` - single parametrized STIS extraction (ocrreject + x1d).
- `reduce_stis_sweep.py` - STIS param sweep (extrsize / centering), loops in Python.
- `reduce_cos.py` - default calcos run.
- `reduce_cos_custom.py` - COS custom XTRACTAB extraction.
- `reduce_cos_sweep.py` - COS NUV box-height sweep.
- `reduce_cos_fuv_sweep.py` - COS FUV box-height sweep.
- `reduce_defringe.py` - G750L defringe chain.
- `run_full_2023ixf.py` - the full 2023ixf epoch reduction.
- `reduce_epoch.py` - generalized epoch reducer used by the pipeline.
- `pipeline.py` - the automated orchestrator (Windows).

Docs: `project_context.md`, `project_state.md` (live status + todo log), `logistics_notes.md`,
`data_access.md`, `science_notes.md`, and this guide. Catalog/output: `uv_sn_catalog.csv`,
`output/*_coadd.{csv,png}`.

Not tracked in git (see `.gitignore`): `Data/` (downloads), `crds_cache/` (refs). Both refetchable.

---

## 11. All the gotchas in one place

- HSTCAL/calcos are Linux-only; reduction must run in WSL `surf_uv`, not native Windows.
- We use the global Python 3.14.2, no virtualenvs. A stray `.venv` breaks notebook kernels; delete it.
- WSL one-liners: shell `$vars`/loops expand to empty crossing Git Bash -> wsl.exe. Literal paths,
  loop in Python.
- CRDS fresh cache: pre-create `mappings/hst`, `references/hst`, `config/hst`.
- CRDS old data: needs `--update-bestrefs` / `assign_bestrefs` to repoint stale header refs.
- STIS single CRSPLIT=1 exposures: no crj, ocrreject errors ("needs more than one input"); extract
  off the flt directly.
- STIS default background can sit far off the trace (offset -300) and near a detector edge; it may
  be off-screen in a tight zoom.
- COS FUV defaults to TWOZONE; set XTRCTALG=BOXCAR to use custom XTRACTAB boxes.
- COS XTRACTAB B_SPEC is in YFULL, not RAWY.
- specutils resampler rejects non-monotonic axes (splice duplicates); sort+dedupe or use np.interp.
- Coadds: don't equally median high-SNR and near-noise spectra; watch inter-grating flux offsets.
- Defringe: widen the mkfringeflat shift range; the scale can still hit its edge.
- Every stistools output must be deleted before re-running (no overwrite).

---

## 12. Quick run reference

```
# discovery query (Windows)
python -c "import pipeline; print(len(pipeline.find_uv_sne()))"   # ~140

# one STIS extraction with custom params (WSL)
wsl.exe -d Debian -- bash -lc 'source ~/miniforge3/etc/profile.d/conda.sh && conda activate surf_uv && export oref=/mnt/c/.../crds_cache/references/hst/stis/ && python /mnt/c/.../reduce_stis.py <obsdir> <root> --out my.fits --extrsize 5 --a2center 894 --maxsrch 0'

# full automated pipeline on a target (Windows; auto-picks the epoch)
python -c "import pipeline; pipeline.run('SN2023ixf')"

# the bundled demo
python pipeline.py    # -> output/SN2023ixf_coadd.csv + .png
```
