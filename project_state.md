# Project State & Daily Log

## Current Status: Week 3 (June 25, 2026)
* **Core Goal:** STIS extraction working end to end. Now: vary STIS extraction params on purpose (PI says we kept everything default), do the full COS custom extraction, start our own coadds, and crack a general UV-supernova query. Slow, careful, interactive.

---

## Active Task List

### 1. Data Access & Environment Configuration
- [x] Install light packages on native Windows (`astroquery`, `astropy`, `specutils`, `numpy`, `scipy`, `matplotlib`). No conda on Windows, packages live on the machine (Python 3.14 kernel).
- [x] Map out the data access landscape (MAST vs Portal vs HST Search Form vs astroquery, and HASP vs HSLA). Written up in `data_access.md`.
- [x] Found hstcal blocker: `stistools`/`ocrreject`/`x1d`/`calcos` wrap compiled HSTCAL binaries with no native Windows build. Set up WSL Debian + Miniforge + `surf_uv` conda env (hstcal + stistools + crds). Details in `logistics_notes.md`.
- [x] Query + download a test target through `astroquery.mast` (working in `spectra_sandbox.ipynb`): SN2020fqv (Thursday, default x1d only) and SN2024iss (proper test).
- [x] Confirmed the SN query convention: use `objectname=` resolver + `radius`, not `target_name`.

### 2. Analytical Pipeline Logic (Bitmasking & Re-extraction)
- [x] DQ bitmask filtering in the notebook (bit 16 + bit 512): `(dq & 16 != 16) & (dq & 512 != 512)`.
- [x] Full STIS re-extraction on SN2024iss G230LB (of8b02010): our own `ocrreject` (flt -> crj) + `x1d` run in WSL via `reduce_stis.py`, output read back on Windows. Our extraction matches the pipeline sx1 (sanity check passed) and we visualized the 2D trace + extraction/background regions.
- [ ] Vary custom extraction params (`extrsize`, `bk1offst`, `bk2offst`, manual center, `maxsrch`) now that the toolchain works, to see effect on faint/contaminated traces.
- [x] COS test target + default products: SN2023ixf (M101, Type II, prop 17497, COS/NUV G230L). Pulled the x1dsum + HASP coadd; science-window plot shows a faint NUV continuum + Mg II 2800 emission. calcos inputs (corrtag/rawtag/flt) already downloaded under `./Data/2023ixf`.
- [x] COS full custom extraction (NUV + FUV done): installed `calcos` in the WSL env, synced COS refs (lref), ran the full custom extraction on SN2023ixf NUV G230L (lf8803d5q). `reduce_cos.py` = default, `reduce_cos_custom.py` copies the XTRACTAB, narrows the source box HEIGHT (57 -> 15), points a copy of the corrtag at it, re-runs calcos. Visualized in `cos_sandbox`. Finding: the narrow box roughly halves the flux scatter (NUVC std 4.84e-15 -> 2.49e-15) while preserving the continuum + Mg II 2800 emission. Interactive height sweep (5,9,15,25,41,57 via `reduce_cos_sweep.py`): on the NUVB science stripe, median flux rises with box height but per-pixel SNR PEAKS at h~9 (3.14) and falls to 2.50 at the default h=57, so the optimal box is ~9px (~25% better SNR than default). FUV done on SN2010jl G130M (2023ixf has no direct FUV) via `reduce_cos_fuv_sweep.py`: same XTRACTAB-edit workflow (already BOXCAR), two segments FUVA/FUVB. KEY CONTRAST: for this BRIGHT source SNR rises with box height and plateaus by h~35-45 (source-limited, want the full profile), whereas the FAINT NUV source peaked at h~9 (background-limited). So optimal box = source-brightness dependent.

### 3. Review & Theoretical Context
- [ ] Review Section 2 of the [HST CCD Spectra Reduction Paper](https://iopscience.iop.org/article/10.3847/2041-8213/ad7855/pdf) for realistic pipeline workflows.
- [x] Browse the interactive line velocity/shape models on the [LBL Supernova Spectra Page](https://supernova.lbl.gov/~dnkasen/tutorial/).

### 4. TODO from 6/22 meeting (active log)
- [x] **cell 9 viz**: orange bg lines were off-screen because the default bg sits at offset -300/-320 (both ~300px below the trace, since A2CENTER=894 is near the top edge). Fixed by re-extracting with close bg (a2center=894, extrsize=5, bk1offst=-8, bk2offst=+8, sizes 10) via `reduce_stis.py --out 230_x1d_closebg.fits`. Now in `stis_sandbox` the orange straddles the trace.
- [x] **STIS param play** (interactive, done this round): `reduce_stis.py` parametrized + `reduce_stis_sweep.py` driver. Findings on SN2024iss G230LB (of8b02010), logged for the PI:
    * **bg placement** barely changes the spectrum here (default offset -300/-320 vs close +/-14). This trace is isolated on a low uniform background, so sky region choice is minor. Expect it to matter for traces on a structured/contaminated background (host-galaxy gradient, nearby source).
    * **extrsize** (aperture height): extrsize=1 loses ~half the flux (only the trace peak). extrsize=3 already recovers nearly all of it (trace is ~7px tall). Beyond ~7 just adds sky/noise without flux. Sweet spot ext3-5.
    * **a2center** (centering): with maxsrch=0, miscentering by 4px drops the flux ~60-70% because the aperture barely overlaps the real trace. Getting the center right (or letting maxsrch find it) matters a lot.
- [x] **COS full custom extraction**: NUV G230L (SN2023ixf) + FUV G130M (SN2010jl) both done (see section 2). calcos in WSL + XTRACTAB box editing + interactive height sweeps. Key finding: optimal box height is source-brightness dependent (faint NUV peaks narrow ~h9; bright FUV keeps improving to ~h35-45).
- [x] **Our own coadds** for STIS and COS (instead of just using HASP). STIS done in `stis_sandbox`: SN2024iss G230LB(our extraction)+G430L+G750L resampled onto a common axis (specutils FluxConservingResampler, native dispersions per the PI reference) and median-combined into one 1670-10250 A spectrum. Gratings agree in the overlaps. Visible G750L fringing past ~8500 A (motivates the defringing todo). COS done in `cos_sandbox`: 6 SN2023ixf NUV G230L exposures x 3 stripes (18 spectra) resampled + median-combined, matches the HASP coadd well (continuum + Mg II 2800) except ~2150-2350 A where our naive median keeps a contaminated NUVA stripe that HASP's flux-checking rejects (motivates adding a sigma-clip / flux filter to our combiner).
- [x] **COS: visualize the extraction** (done during the COS extraction work): cross-dispersion count profiles with the default vs custom boxes overlaid, for both NUV (3 stripes) and FUV (FUVA/FUVB, in the geo-corrected YFULL frame). In `cos_sandbox`.
- [x] **Defringing** step for the G750L grating: done in `stis_sandbox` via `reduce_defringe.py` (WSL). SN2024iss G750L science of8b02040 + its contemporaneous fringe flat of8b02050 (50s CCDFLAT in the same visit). Chain: `stistools.defringe.normspflat` -> `mkfringeflat` (best shift -0.46px, scale 1.08 after widening the shift range off the edge) -> `defringe` -> `x1d`. Fringe oscillations past ~9000 A are clearly reduced (residuals remain, normal for G750L).
- [ ] **Full 2023ixf reduction in one go** (COS NUV + COS FUV).
- [x] **General UV-SN query** (the whole point): built `uv_sn_query.ipynb`. Instead of name-matching (misses SN2020fqv=`TESS-SN`) or paper searches, filter the HST obs table on `target_classification='*upernova*'` (the proposer Phase II class) across STIS+COS spectroscopy, then dedup by coordinates (5 arcsec). Result: 1591 SN-classified spectra -> 140 unique SNe (139 transients, 1 remnant, 102 with a UV grating), written to `uv_sn_catalog.csv`. Validated: SN2020fqv (TESS-SN), 2024iss, 2023ixf, 2010jl all caught. Caveat: classification is proposer-set so a SN observed under a non-SN program could slip through; a SIMBAD/TNS coordinate cross-match would be the belt-and-suspenders next step.
- [ ] **First automated pipeline pass** (where helper .py / OOP starts to make sense).

---

## Technical Reference & Code Snippets

### Standard Bitmasking Filter Template
```python
# Extracting array masks for data cleansing
dq512_mask = (x1d_data['DQ'] & 512 != 512)[0]
dq16_mask  = (x1d_data['DQ'] & 16 != 16)[0]
combined_clean_mask = dq16_mask & dq512_mask

clean_wavelength = wvl_array[combined_clean_mask]
clean_flux = flux_array[combined_clean_mask]
### 6/26 — full 2023ixf reduction (DONE)
- One STIS epoch of SN2023ixf (prop 17205, oezt01 visit): G230LB (oezt01040), G430L (oezt010h0), G750L (oezt010e0) + contemporaneous CCDFLAT (oezt010d0).
- These are single CRSPLIT=1 exposures, so NO ocrreject - extract straight from the FLT. `run_full_2023ixf.py` (WSL) syncs CRDS refs, x1d on the blue gratings, normspflat->mkfringeflat->defringe->x1d for G750L.
- Coadd in full_2023ixf.ipynb: load the 3 STIS x1d + COS NUV cspec, np.interp onto a common axis, nanmedian. Caveat: COS NUV is faint + fully overlapped by G230LB, so a naive median drags it down; combine the 3 STIS gratings (they tile 1670-10250) and keep COS NUV as overlay.
- Caveat: mkfringeflat scale hit the edge of its range (1.2); shift converged fine. Minor for the demo.
- Result: clean full HST UV-optical SED of SN2023ixf, our own extraction throughout.
