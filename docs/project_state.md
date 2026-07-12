# Project State & Daily Log

Log is reverse-chronological: most recent session at the top.

---

## Workflow Principles (static - "review our workflow principles first")

Compiled from `copilot-instructions.md` (mirrored in VS Code `settings.json`) and the
process rules locked in over the weeks. These do not change session to session.

### Hard style rules (a PI reviews this code; AI involvement must not be detectable)
- Match the reference notebooks (`1dspectrumreference.ipynb`, `combine_spectra_reference.ipynb`):
  short lowercase abbreviated vars (wvl, flx, err, spec, tbdata, flist), terse casual inline
  comments ("u gotta delete it to run again"), no docstrings on trivial functions, no formal
  block headers or banner decorations.
- Comments explain the WHY, not the WHAT.
- NO emojis anywhere (code, comments, logs, chat). NO em-dashes.
- Banned words: robust, crucial, leverage, seamless, delve, dive in, furthermore, snippet.
- NEVER name-drop people/roles in CODE comments (no "Wynn", "PI", "the PI"). Use neutral refs
  ("the already-reduced early 1d spectra", not "Wynn's files"). Docs/.md files MAY name people;
  code comments may NOT.
- No unnecessary try/except; do not swallow errors, let them bubble. No redundant null checks
  when types already guarantee safety. No over-engineering beyond what's asked.

### Process rules (slow, careful, interactive)
- ALWAYS ask clarifying questions at real decision points (which dataset, instrument, extraction
  params, output path). Use the built-in question tool. Do not guess on consequential choices.
- DO NOT end the process until ALL outlined tasks and todos are complete.
- Maintain the active todo log in `project_state.md` and keep important todos in memory.
- Don't overwrite existing notebook cells much; prefer adding. May rewrite/rerun cells when a
  prior result was flawed (e.g. a biased sample) and note it.
- Keep commenting code (in the reference style) and documenting changes in `project_state.md`
  and the `pipeline_phaseX.md` docs as work proceeds.
- Helper `.py` / OOP only when actually useful or when building the full automated pipeline.
- USE MCPs for data access, library/code documentation, and up-to-date references.
- Aggressively refer to and learn from the existing docs compiled over prior weeks.

### Python environment (STRICT)
- NEVER create virtual environments (venv, uv venv, conda env) unless explicitly asked.
- Native Windows: global Python 3.14 for astroquery/astropy/specutils/plotting (NO hstcal).
- Reduction (hstcal/calcos/stistools/crds) runs in WSL Debian conda env `surf_uv`. Run pattern:
  `wsl.exe -d Debian -- bash -lc 'source ~/miniforge3/etc/profile.d/conda.sh && conda activate surf_uv && ...'`.
- Shell vars/loops expand to EMPTY crossing git-bash -> wsl.exe; use literal paths, loop in Python.

### Science scope
- A unified, uniformly reduced repository of ALL HST UV SN spectra (STIS + COS, not STIS-only).
- Our OWN reduction throughout (trace, background, CR rejection, defringe, coadd); MAST/HASP/HSLA
  defaults are not the final products.

---

## Current Status: Week 5 (July 9, 2026)
* **Core Goal:** a uniformly reduced repository of all HST UV supernova spectra (STIS + COS), our own reduction throughout (trace, background, CR, defringe, coadd), not MAST defaults.
* **Latest:** extended 2023ixf to the full 16-epoch time series (days 3-1094); conducted STIS automation investigation (automation_sandbox.ipynb); key findings: default bg hardcoded wrong, trace center reliable for bright sources but fails silently on faint ones. New docs: pipeline_phase3.md (today's work), pipeline_phase2.md P2.9 (defringe recap). See 7/9 log.
* Style: slow, careful, interactive.

---

## Open items / next steps (decisions to be made with PI / Wynn)

1. **[DONE] Switch 2023ixf to main science epochs (days 14-66).** Identified the main program as **prop 17313** (of43 visits, 52x0.2 E1), not 17315. Days 14/19/24/66 = visits of4301/02/03/05; of4304 skipped (the failed guide-star-acq visit, matches the paper's day-50 failure). Downloaded into per-epoch dirs `data/2023ixf/epochs/day{N}_of43xx/`, reduced with our own extraction (`reduce_epoch_ts.py`), and stacked into the time series (see 7/1-7/2 log). The pre-14d epochs use the supplied 1D files in `data/earlytime_2023ixf`.
2. **2D aperture/background visualization saved per extraction.** Wynn explicitly asked for this: "output visualization plots of the 2-D image with aperture and backgrounds shown -- will help us identify any issues with the extractions." The `show_extraction_regions` function exists in `stis_sandbox.ipynb`; needs to be wired into `pipeline.py` / `reduce_epoch.py` so it saves a PNG for every extraction automatically.
3. **[DONE] Inter-grating flux scaling in the coadd.** Added to the main-epoch coadd (anchor G230LB, constant-% in the overlaps, guarded against faint/edge artifacts). Unsaturated day 14/19/24 factors came out 0.99-1.08 (gratings agree), confirming the oezt01 x0.595 was a saturation artifact. Day 66's faint blue anchor is guarded. See the 7/1-7/2 log.
4. **SNR weighting / faint-leg handling in the coadd.** A proper inverse-variance weighted combine would let the faint COS NUV leg contribute instead of being excluded.
5. **Sigma-clip / flux filter in the COS coadd.** Our 6-exposure NUV median diverges from HASP in ~2150-2350 A (a contaminated NUVA stripe). HASP drops exposures deviating >5% from the coadd median; adding that filter would match HASP.
6. **COS leg in the automated pipeline.** `pipeline.py` only handles STIS today; add a COS branch.
7. **SIMBAD/TNS cross-match for the UV-SN catalog.** Classification is proposer-set, so a SN under a non-SN program could be missed and a misclassified object included. A coordinate cross-match cleans both up.
8. **Systematic reduction of the full catalog.** Loop the pipeline over the 140 SNe in `output/uv_sn_catalog.csv`.
9. **CMFGEN model comparison.** The science goal: compare our reduced spectra to CMFGEN models for mass-loss rate / CSM radius. Not built yet.

---

## Daily Log (most recent first)

### 7/9 — full 2023ixf time series (days 3-1094) + STIS automation investigation
Full details in `pipeline_phase3.md`.

**Full time series (16 epochs, days 3-1094):**
- Queried the full 2023ixf HST spectroscopy landscape: 28 visits, days 3-1094, most recent 2026-05-16.
- Validated COS epoch count: only 1 real COS epoch (lf8803, day 214, sep 0.1"); lf9256 is a different M101 target 413.7" away (cone-search false positive). PI was right.
- Downloaded all 16 late visits (days 183-1094) into `data/2023ixf/epochs/day{N}_{visit}/`.
- Reduced them our own way with new `scripts/reduce_stis_generic.py` (MAMA x1d off flt, G750M x1d off crj). 111 extractions, 0 failures on MAMA.
- G750M days 913/924 auto-trace failed ("Cannot extract") -- found trace at row ~898 off the 2D, re-extracted with a2center=898/maxsrch=0/extrsize=7.
- Day-577 (ofg001) is an aborted visit (zero exptime, calstis produces nothing) -- covered by ofg005 (day 619).
- Built full-series waterfall plot `output/2023ixf_fullseries.png` (16 epochs, rest-frame z=0.000804, broken at UV-Hα gaps). Shows smooth photospheric continuum -> UV-fading -> Mg II] 2800 + Hα emission (CSM interaction) at late times.
- Added `pipeline_phase2.md` P2.9 (consolidated defringe recap) and created `pipeline_phase3.md`.

**STIS automation investigation (automation_sandbox.ipynb):**
- Downloaded MAST sx1 headers for 5 SNe (SN2024iss, SN2023ixf, SN2010jl, SN2021yja, SN2009ip), all G230LB CCD 52X0.2.
- **A2CENTER**: 4 of 5 cluster at 893-894 (the standard E1 sub-position). SN2010jl lands at 912 -- verified by downloading the 2D crj: the trace genuinely IS at row 912 for that program (a different E1 sub-position). Pipeline auto-finder was correct when source is bright.
- **BK1OFFST/BK2OFFST**: **completely hardcoded at -300/-320 across every program, every SN.** For E1 sources at row ~893-912, this puts both bg regions off the bottom of the chip (~row 593-592). The bg subtraction is extracting from empty chip.
- **Bottom line**: trace center is reliable for bright sources; need profile-peak fallback for faint ones (can't hardcode E1 row -- it varies by program). Background must be placed relative to the found trace, not hardcoded.

### 7/1-7/2 — full 2023ixf time series (paper Figure 1 reproduction)
Built the day 3-66 NUV-to-NIR time series, our own reduction for the main epochs. Details also in `pipeline_phase2.md`.
- **Epoch/naming explainer** written and moved into a new `docs/pipeline_phase2.md` (phase-2 master doc); the phase-1 guide now points to it. Covers the rootname anatomy, the 3 gratings, the CCDFLAT fringe flat, and how a "day X" phase maps to an HST visit/program.
- **Figure 1 skeleton** added to `full_2023ixf.ipynb`: the 5 pre-14d epochs (already-reduced 1D spectra in `data/earlytime_2023ixf`, one per epoch, mapped 1st..5th = 3.25/4.25/5.25/8.25/11.25 d) plotted log-log, stacked with per-epoch offsets, the paper's color scheme, and the atomic-line guides. `output/2023ixf_timeseries.png`.
- **Main epochs identified + downloaded**: main program = **prop 17313** (of43, 52x0.2 E1), not 17315. Days 14/19/24/66 = of4301/02/03/05; of4304 is the failed guide-star-acq visit (skipped, matches the paper's day-50 failure). Each visit: 010 g230lb, 020 g430l, 030 g750l, 040 ccdflat. Per-epoch dirs `data/2023ixf/epochs/day{N}_of43xx/`.
- **Reduced the main epochs ourselves** (`scripts/reduce_epoch_ts.py`, a new crj-aware reducer; `reduce_epoch.py` left intact for backtracking). These are CRSPLIT>1, so we x1d off the pipeline cr-combined `crj` (g230lb/g430l) and defringe the g750l with the contemporaneous 0.3x0.09 E1 CCDFLAT (of43NN040) -- the paper's recipe. Defringe converged cleanly (shifts ~0.3-0.5 px, scales ~1.12, not at the range edges this time).
- **day-66 g230lb** blue flux is so low the auto trace-find fails ("Cannot extract"), exactly as the paper reports. Fixed by forcing the fixed E1 trace center (a2center=893.5, maxsrch=0) -- baked into the reducer as a retry on empty output. Also needed an extra CRDS bestrefs sync (some refs are date-specific and the loop's sync lagged).
- **Full 9-epoch figure**: coadd each main epoch's 3 gratings (dq clean + interp + nanmedian) and append to the same plot -> `output/2023ixf_timeseries_full.png`. Reproduces the paper Fig 1 evolution: hot blue UV + flash features early, UV fading through days 14-24, red P Cygni (Ha, Ca II) dominant by day 66 (its UV is noisy, as expected for the faint manually-extracted blue end).
- **Open decisions logged (pipeline_phase2.md P2.3)**: flux column FLUX vs FLUX_corr (using FLUX for now), narrow ISM-line clipping (queued), the rest-frame assumption, and the early(dereddened)-vs-main(raw) flux-system mismatch.
- **Comment-style pass**: stripped all PI/person name-drops and banner-style comments from the notebook; a formatter had also scrambled the time-series cells, so fixed the cell order and removed a duplicate plot cell.

**Handoff-note items (all cleared, same session):**
- **2010jl pipeline (STIS + FUV).** 2010jl has STIS/CCD G230LB + G430L (no G750L, i.e. 2-grating) plus COS FUV. Downloaded + reduced the obk002 epoch (prop 12242): extracted all 5 CRSPLIT exposures off their `crj`, coadded per grating, overlaid our COS FUV -> `output/2010jl_sed.png`. A bright IIn, so strong emission lines throughout; the two gratings agree at the ~3000 A splice. (pipeline.py still assumes 3 gratings; `reduce_epoch_ts.py` already handles any subset.)
- **STIS coadd -- inter-grating scaling.** Added the paper's constant-% grating alignment to the main-epoch coadd (anchor G230LB), guarding factors to 1 when the overlap sits on a faint/edge region. Unsaturated day 14/19/24 factors come out 0.99-1.08 -- the gratings genuinely agree, confirming the oezt01 0.595 was a saturation artifact (as Wynn predicted). Day 66's faint blue anchor gives an unreliable 1.88, so it's guarded. The time series now uses the scaled coadds (splice steps gone on days 14-24).
- **G750L defringe -- verified working.** Quantified the fringe residual (RMS of flux/median-smooth) on the day-14 G750L: defringe cuts it 30-48% over 8000-9700 A and visibly smooths the band. Residuals past ~9700 A remain (hardest reddest fringes + photon noise where the source is faint) -- normal for G750L. The earlier 'not working' was the oezt01 epoch (scale solver hit its boundary); the of43 epochs converge properly. `output/defringe_check_day14.png`.
- **COS NUVB box centering.** The apparent off-center box was a frame mismatch: `B_SPEC` lives in `YFULL` but the plot histogrammed `YCORR`. Redone in `YFULL`, the stripe peaks sit 0-3 px from `B_SPEC` (NUVB worst at 2.8 px), so the narrow custom box captures the flux fine. Fixed the plot to use `YFULL`.
- **COS FUV box-height / Ly-alpha.** The ~2x flux swing across the box-height sweep is a real point-source effect: the COS FUV cross-dispersion profile is broad, so a narrow box clips flux. From h11 to h45 the continuum grows 1.95x and the Ly-alpha peak 2.23x, both plateauing by h35-45. Use a box on the plateau (default 41 is fine) or aperture-correct.
- **Housekeeping.** Deleted three empty auto-generated `.venv*` dirs (uv spawns one per notebook-config call; they wrap the global 3.14.2 and are unused). Reconfirmed the no-venv rule.

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
