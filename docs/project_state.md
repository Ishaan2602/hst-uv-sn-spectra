# Project State & Daily Log

Log is reverse-chronological: most recent session at the top.

## Current Status: Week 3 (June 29, 2026)
* **Core Goal:** a uniformly reduced repository of all HST UV supernova spectra (STIS + COS), our own reduction throughout (trace, background, CR, defringe, coadd), not MAST defaults.
* **Latest:** full 2023ixf reduction + first automated pipeline done; workspace reorganized into `notebooks/ scripts/ docs/ output/`; reviewed the Bostroem et al. 2024 SN2023ixf paper to validate and refine our defringe + coadd; 6/29 PI meeting response received -- scope updated (see open items).
* Style: slow, careful, interactive.

---

## Open items / next steps (decisions to be made with PI / Wynn)

1. **Switch 2023ixf to main science epochs (days 14-66).** Wynn confirmed: the oezt01 early visit has saturation issues and should NOT be used as the primary product. We should do the full reduction for all epochs from day 14 onward. Wynn will send the early-epoch 1D files (days 3-11) separately to include in the repo without re-reducing. Action: identify the main-program obs IDs (prop likely 17315, E1 slit position, days 14/19/24/66), download, and run the full pipeline on each as a time series.
2. **2D aperture/background visualization saved per extraction.** Wynn explicitly asked for this: "output visualization plots of the 2-D image with aperture and backgrounds shown -- will help us identify any issues with the extractions." The `show_extraction_regions` function exists in `stis_sandbox.ipynb`; needs to be wired into `pipeline.py` / `reduce_epoch.py` so it saves a PNG for every extraction automatically.
3. **Inter-grating flux scaling in the coadd.** The Bostroem 2024 paper (Sec 2) aligns G430L and G750L to the G230LB flux by a constant percentage. Our coadd is a naive median. Tried it: G430L x0.595, G750L x0.913 (anchor G230LB). Wynn says the large G430L offset is likely a saturation artifact of the early epoch -- expect factors near 1 for the main science epochs. Revisit after switching to days 14-66.
4. **SNR weighting / faint-leg handling in the coadd.** A proper inverse-variance weighted combine would let the faint COS NUV leg contribute instead of being excluded.
5. **Sigma-clip / flux filter in the COS coadd.** Our 6-exposure NUV median diverges from HASP in ~2150-2350 A (a contaminated NUVA stripe). HASP drops exposures deviating >5% from the coadd median; adding that filter would match HASP.
6. **COS leg in the automated pipeline.** `pipeline.py` only handles STIS today; add a COS branch.
7. **SIMBAD/TNS cross-match for the UV-SN catalog.** Classification is proposer-set, so a SN under a non-SN program could be missed and a misclassified object included. A coordinate cross-match cleans both up.
8. **Systematic reduction of the full catalog.** Loop the pipeline over the 140 SNe in `output/uv_sn_catalog.csv`.
9. **CMFGEN model comparison.** The science goal: compare our reduced spectra to CMFGEN models for mass-loss rate / CSM radius. Not built yet.

---

## Daily Log (most recent first)

### 6/29 — PI meeting; Wynn's email response
Sent Wynn a summary email after the 6/26 session. His response (paraphrased; full text in `extraction_pipeline_guide.md` section 13):
- **Scaling offset is likely saturation artifact.** The G430L x0.595 factor is probably specific to the saturated oezt01 early epoch, not a real calibration issue. Expect near-1 factors for the main science epochs.
- **Clarification on COS NUV / SN2010jl:** to clarify in meeting -- the COS NUV discussion is about SN2023ixf; SN2010jl was a separate FUV test case.
- **Epoch scope update:** do NOT re-reduce the pre-14d data. Wynn will send the early-epoch 1D files directly. For us: do the full reduction on the main epochs (days 14, 19, 24, 66) from the primary program, at the E1 slit position.
- **2D aperture visualization:** Wynn explicitly asked for PNG outputs of the 2D image with aperture and background regions shown for every extraction. The helper exists (`show_extraction_regions`); needs to be wired into the pipeline.
- **Agreed next steps:** settle coadd scaling (revisit after switching epochs), add COS to the pipeline, robustify the discovery catalog, loop over all 140 SNe.

### 6/26 (PM) — Bostroem et al. 2024 paper review; defringe + coadd validation
Reviewed `docs/2023ixf_paper_review.md` (Circumstellar Interaction in the UV Spectra of SN 2023ixf, Bostroem et al. 2024), focus on Section 2 (Observations and Data Reduction). What it means for our pipeline:
- **Validates our DQ bitmask.** The paper masks bad pixels with DQ flags 16 (high dark rate) and 512 (bad reference pixel) -- exactly our `(dq & 16 != 16) & (dq & 512 != 512)`.
- **Validates our defringe approach.** They correct G750L fringing with a contemporaneous fringe flat via `stistools.defringe`, then extract with `stistools.x1d` -- our exact chain (normspflat -> mkfringeflat -> defringe -> x1d). Nuance: their fringe flat is taken through the smaller 0.3"x0.09" aperture (standard at the E1 slit position); ours is the visit's CCDFLAT exposure. Confirmed below: our flat (oezt010d0) is 52X0.1 while the science G750L (oezt010e0) is 52X0.2, so the flat IS the narrower slit -- same principle as the paper. (This epoch sits at the standard 52X0.2 position, not E1, consistent with it being the early visit.)
- **Refines our coadd.** They do inter-grating flux scaling (align G430L + G750L to G230LB by a constant percentage) and intra-visit exposure scaling (to the highest-flux exposure, first-degree poly), then median-combine excluding bad pixels. Our coadd is the first-order version (naive median, no scaling). This is the recipe for next-step #1.
- **Saturation caveat on our 2023ixf epoch.** Our full reduction used the GO-17205 (oezt01) early visit -- which the paper explicitly flags as the supplementary early data with pointing/acquisition/CCD-saturation problems (e.g. ~3200-5000 A saturated in visit 1). Those data are usable for feature ID, NOT absolute flux. So our pretty full-SED plot has flux-unreliable regions; label them. Confirmed below: our G430L is saturated across 3178-5022 A (668 px), matching the paper's stated 3200-5000 A visit-1 saturation.
- **Dark-rate artifact caveat.** The paper found several narrow "emission" features were localized high-dark-rate artifacts (verified against the 2D darks), not real. Relevant when we interpret narrow features.
- **CR strategy.** Main epochs took >=4 exposures per grating for auto CR rejection; our oezt01 epoch was single CRSPLIT=1 (consistent with it being the problematic early data, hence no ocrreject for us).
- Checked off the Section-2 paper-review task.

**Test results (ran the 3 checks on our oezt01 epoch, in `full_2023ixf.ipynb`):**
- **Fringe-flat aperture confirmed.** Science G750L (oezt010e0) APERTURE=52X0.2, fringe flat (oezt010d0) APERTURE=52X0.1 -- the flat is the narrower slit, same principle as the paper's small-aperture flat. Our defringe input is valid.
- **Saturation matches the paper on our own data.** DQ bit 256: G230LB 3 px and G750L 3 px (both minor), G430L **668 px spanning 3178-5022 A** -- essentially the paper's 3200-5000 A visit-1 saturation. So our G430L optical mid-band IS the saturated visit; flagged unreliable for absolute flux and shaded on the comparison plot.
- **Inter-grating scaling tried, kept as comparison.** Anchor G230LB: G430L x0.595, G750L x0.913. The 0.595 is a big (~40%) pull-down measured right at the grating edges (2950-3050 A, where G230LB throughput droops), so it is partly an edge artifact, not necessarily a true cal offset. Decision: naive median stays the primary coadd, scaled version shown as a comparison only, pending the PI's call on the absolute-flux anchoring. Overlap windows kept as-is with the caveat noted.
- Added three cells to `notebooks/full_2023ixf.ipynb` (paper-checks markdown, data-quality cell, naive-vs-scaled comparison) and saved `output/2023ixf_coadd_scaled.png`.
- **Bug fix:** the reorg find/replace had mangled the notebook's `../data` paths into `.../data` (three dots), which is why those cells had not re-run since the reorg. Fixed the download_dir + hst paths back to `../data`.

### 6/26 (AM) — full 2023ixf reduction, automated pipeline, workspace reorganization
- **Full 2023ixf reduction** (`notebooks/full_2023ixf.ipynb` + `scripts/run_full_2023ixf.py`): one STIS epoch (oezt01, prop 17205): G230LB (oezt01040), G430L (oezt010h0), G750L (oezt010e0) defringed with the contemporaneous CCDFLAT (oezt010d0). Single CRSPLIT=1 exposures so NO ocrreject -- extract straight off the FLT. Coadd loads the 3 STIS x1d + COS NUV cspec, np.interp onto a common axis, nanmedian. Clean full HST UV-optical SED, 1670-10250 A. COS NUV is faint + fully overlapped by G230LB so a naive median drags it down ~50%; combine the 3 STIS gratings and keep COS NUV as overlay.
- **First automated pipeline** (`scripts/pipeline.py` + `scripts/reduce_epoch.py`): discover (classification query) -> auto-pick a STIS CCD epoch -> download -> reduce in WSL -> coadd -> save `output/<sn>_coadd.{csv,png}`. Auto-picker validated (independently found a 3-grating visit for 2023ixf). Demo runs end to end.
- **Workspace reorganization**: flat root -> `notebooks/` (+ `notebooks/reference/`), `scripts/`, `docs/`, `output/`. Updated all internal paths (`../data/`, `../output/`, `../crds_cache/` in notebooks; `scripts/reduce_epoch.py` in pipeline). Deleted superseded `spectra_sandbox.ipynb`. PI reference notebooks left untouched.
- **Session fixes / issues**:
    * `.venv` kernel conflict: VS Code auto-selected a stray empty `.venv` for the new notebook, stalling kernel start. Deleted it; all notebooks use the global Python 3.14.2. No virtualenvs.
    * Specutils renamed `Spectrum1D` -> `Spectrum`; added a try/except import alias.
    * `FluxConservingResampler` rejects non-monotonic axes (splice duplicates); switched to `np.interp` after a `clean()` sort+dedupe. Flux-conserving rebinning is a later refinement.
    * mkfringeflat best-scale hit the top of its range (1.2) on the 2023ixf flat; shift converged fine (-0.46px). Minor.

### 6/25 — STIS param play, COS extraction, coadds, defringing, UV-SN query, first pipeline (the 6/22 meeting todos)
- **STIS param play** (`reduce_stis.py` parametrized + `reduce_stis_sweep.py`), SN2024iss G230LB (of8b02010):
    * **bg placement** barely changes the spectrum here (default offset -300/-320 vs close +/-14): isolated trace on low uniform background, so sky region choice is minor. Matters for structured/contaminated backgrounds.
    * **extrsize**: extrsize=1 loses ~half the flux (trace peak only); extrsize=3 recovers nearly all (trace ~7px tall); beyond ~7 just adds sky/noise. Sweet spot 3-5.
    * **a2center**: with maxsrch=0, miscentering by 4px drops the flux ~60-70%.
    * Earlier, fixed the cell-9 viz (default bg was off-screen at offset -300/-320 since A2CENTER=894 is near the top edge); re-extracting with close bg makes the orange straddle the trace.
- **COS custom extraction** (NUV SN2023ixf G230L + FUV SN2010jl G130M): `calcos` in WSL, XTRACTAB box editing, interactive height sweeps (`reduce_cos_custom.py`, `reduce_cos_sweep.py`, `reduce_cos_fuv_sweep.py`). KEY FINDING: optimal box height is source-brightness dependent. Faint NUV: SNR peaks at h~9 then falls (default 57 way too wide; ~25% SNR gain). Bright FUV: SNR keeps rising and plateaus by h~35-45 (default 41 near-optimal). FUV gotchas: old data needs `--update-bestrefs`; XTRACTAB B_SPEC lives in the geo-corrected YFULL frame, not RAWY; FUV defaults to TWOZONE (set XTRCTALG=BOXCAR for custom boxes).
- **Coadds** (STIS + COS): STIS SN2024iss G230LB(ours)+G430L+G750L resampled + median-combined into one 1670-10250 A spectrum; gratings agree in overlaps; visible G750L fringing past ~8500 A. COS 6 SN2023ixf NUV exposures x 3 stripes median-combined, matches HASP except ~2150-2350 A where our naive median keeps a contaminated NUVA stripe HASP rejects.
- **COS extraction viz**: cross-dispersion count profiles with default vs custom boxes, NUV (3 stripes) + FUV (FUVA/FUVB in YFULL).
- **Defringing** (`reduce_defringe.py`): SN2024iss G750L of8b02040 + contemporaneous fringe flat of8b02050. normspflat -> mkfringeflat (best shift -0.46px, scale 1.08 after widening the shift range) -> defringe -> x1d. Fringes past ~9000 A clearly reduced.
- **UV-SN query** (`uv_sn_query.ipynb`): filter HST obs on `target_classification='*upernova*'` across STIS+COS, dedup by coords (5 arcsec) -> 1591 spectra -> 140 unique SNe (139 transients, 1 remnant, 102 with a UV grating) -> `output/uv_sn_catalog.csv`. Catches SN2020fqv (logged `TESS-SN`) that name-matching misses.

### 6/22 — meeting: todo list set
- Out of the meeting: do the COS full custom extraction, do our own coadds (STIS + COS), visualize the COS extraction, defringe G750L, do the full 2023ixf in one go, build the general UV-SN query, and a first automated pipeline pass. Also vary STIS extraction params on purpose (PI: we had been keeping everything default).

### ~6/21 — STIS reduction toolchain on WSL
- Found the hstcal blocker: `stistools`/`ocrreject`/`x1d`/`calcos` wrap compiled HSTCAL binaries with no native Windows build. Set up WSL Debian + Miniforge + `surf_uv` conda env (hstcal + stistools + crds). Details in `logistics_notes.md`.
- First full STIS re-extraction on SN2024iss G230LB (of8b02010): our own `ocrreject` (flt -> crj) + `x1d` in WSL via `reduce_stis.py`; output matches the pipeline sx1 (sanity check passed). Visualized the 2D trace + extraction/background regions.
- WSL gotcha: shell vars/loops expand to empty crossing Git Bash -> wsl.exe; all scripts use literal paths and loop in Python.

### ~6/18 — environment + data access foundation
- Installed light packages on native Windows (astroquery, astropy, specutils, numpy, scipy, matplotlib); global Python 3.14, no conda on Windows.
- Mapped the data access landscape (MAST vs Portal vs HST Search Form vs astroquery; HASP vs HSLA) -> `data_access.md`.
- Query convention: use `objectname=` resolver + `radius`, not `target_name`. DQ bitmask bit 16 + bit 512.

---

## Task checklist (status)

### 1. Data access & environment
- [x] Light packages on native Windows; WSL `surf_uv` for hstcal/calcos.
- [x] Data access landscape mapped (`data_access.md`).
- [x] Query convention confirmed (`objectname=` + radius).

### 2. Extraction (STIS + COS)
- [x] DQ bitmask (16 + 512).
- [x] STIS re-extraction (SN2024iss G230LB) matches pipeline.
- [x] Vary STIS extraction params on purpose (extrsize / centering / bg) -- findings logged 6/25.
- [x] COS full custom extraction, NUV + FUV, with interactive height sweeps.

### 3. Review & theoretical context
- [x] Review Section 2 of the [HST CCD Spectra Reduction Paper](https://iopscience.iop.org/article/10.3847/2041-8213/ad7855/pdf) (Bostroem et al. 2024). Review in `docs/2023ixf_paper_review.md`; analysis logged 6/26.
- [x] Browse the [LBL Supernova Spectra Page](https://supernova.lbl.gov/~dnkasen/tutorial/).

### 4. Workflows
- [x] Coadds (STIS + COS).
- [x] COS extraction visualization.
- [x] G750L defringing.
- [x] Full 2023ixf reduction in one go.
- [x] General UV-SN query.
- [x] First automated pipeline pass.

---

## Technical Reference

### Standard bitmask filter template
```python
# Extracting array masks for data cleansing
dq512_mask = (x1d_data['DQ'] & 512 != 512)[0]
dq16_mask  = (x1d_data['DQ'] & 16 != 16)[0]
combined_clean_mask = dq16_mask & dq512_mask

clean_wavelength = wvl_array[combined_clean_mask]
clean_flux = flux_array[combined_clean_mask]
```
