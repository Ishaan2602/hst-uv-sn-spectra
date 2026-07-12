# Logistics & Technical Notes

## Environment / Toolchain (set up 6/21)
* Native Windows can run the light stuff (astroquery, astropy, specutils, plotting) but NOT the reduction. `stistools.x1d`/`ocrreject` and COS `calcos` just wrap the compiled HSTCAL binaries (`cs0.e`..`cs12.e`), which have no native Windows build. See `data_access.md` section 7 for the full story.
* So the compiled reduction runs in WSL:
    * Distro: Debian (user `ishaang6`). Miniforge installed at `~/miniforge3`.
    * Conda env `surf_uv`: `hstcal` + astropy stack from conda-forge, then `pip install stistools crds` (those two arent on conda-forge). python 3.12.
    * hstcal binaries live in `~/miniforge3/envs/surf_uv/bin/` (cs2.e = ocrreject, cs6.e = x1d). They only resolve on PATH when the env is ACTIVATED, so always `conda activate surf_uv` first (stistools shells out to them by bare name).
* Run pattern from the Windows side:
    * `wsl.exe -d Debian -- bash -lc 'source ~/miniforge3/etc/profile.d/conda.sh && conda activate surf_uv && export oref="<refdir>/" && python <script> <args>'`
    * Project is visible in WSL at `/mnt/c/Users/eluru/CALTECH_SURF_2026/WORK`.
    * `reduce_stis.py` is the little driver: takes obsdir + rootname, runs ocrreject(flt->crj) then x1d(crj->x1d).
* CRDS reference files:
    * Cache at `WORK/crds_cache`. On a fresh Windows cache you must pre-make `mappings/hst`, `references/hst`, `config/hst` or CRDS errors writing into them.
    * `crds.assign_bestrefs([...], sync_references=True)` pulls bestrefs and stamps headers with `oref$...`. STIS refs land in `crds_cache/references/hst/stis/`. Point `oref` there (trailing slash). WSL reuses the same cache fine over `/mnt/c`.

## Query naming convention (astroquery.mast)
* `target_name` in MAST is whatever the proposer typed (SN2020fqv was logged `TESS-SN`), so name-matching misses. Use `Observations.query_criteria(objectname="SN2024iss", radius="0.05 deg", ...)` which runs the resolver (name -> coords) then cone searches. SN2024iss resolves fine.
* STIS CCD CRSPLIT data -> pipeline 1d product is `sx1` (from the CR-rejected `crj`), NOT `x1d`. MAMA data -> `x1d`. Filter products accordingly.

## Data & Extraction Methods
* Using `stistools` we put apertures along trace, put bg apertures at edges to subtract bg noise. And we use FLT and 2D.
* `FLT` = flat-fielded 2D images. `CRJ` = cosmic-ray rejected 2D images.
* In code: `./Data/2024iss/HST/of8b02010/230_1_crj.fits` — this `230_1` file will do more rigorous cosmic ray extraction.
* Right now we’re just getting the `x1d` file MAST gives us (we use the `.x1d.x1d` function to do our own trace).
* MAST's default pipeline traces for background regions might be totally wrong.
* NOTE: inside the `.x1d.x1d` function, you can manually enter the center (not currently shown in our example code).
* NOTE: Every time you make an `x1d` file, you must delete it to run the script again (the pipeline cannot override existing files).

## Instrument & Query Filters
* Check naming convention for `astroquery.mast` supernovae sources.
* Look at COS instrument, not just STIS.
* Focus on Instrument, not waveband.
* Relevant STIS detectors: STIS/NUV-MAMA, STIS/FUV-MAMA, STIS/CCD.
* COS has FUV filter G160M (but later we’d need to look if people use OTHER filters for older supernovae).
* MAST filtering bands: 1000-1600 is FUV, 1600-3200 is NUV.
* 230L: nUV. 430L: optical.

## General Observations
* Note: A single supernova might have multiple observations (epochs).

## PI Reference Image: 2D Trace + Extraction Regions
* PI shared a screenshot of a 2D spectroscopic viewer showing a raw/flat-fielded fits (FLT or CRJ).
* Speckled noise background, with a bright horizontal band across center = the 2D spectral trace of the SN.
* Guidelines overlaid mark the `stistools.x1d` extraction params:
    * dashed red line through the center of the bright trace
    * Two solid red lines above and below the bright trace show the main extraction aperture.
    * Tow pairs of orange lines above/below the trace = background extraction regions, placed off the trace to sample/subtract sky/detector background.
* This is the same idea as `show_extraction_regions` in `1dspectrumreference.ipynb` (there: red = extraction, orange = bg). The aperture follows `EXTRLOCY`/`EXTRSIZE`, bg follows `BK1OFFST`/`BK2OFFST` + `BK1SIZE`/`BK2SIZE`.