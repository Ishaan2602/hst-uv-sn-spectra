# Emission analysis skill (load via read_file when doing Mg II / Lya emission work)

## Current canonical state (phase 5, Sep 2026)

### Big finding: Bostroem+2026 Table 4 is OBSERVED flux
Table 4 fluxes are NOT extinction-corrected. Our "+15-25% offset" and the GGI window override were both
artifacts of assuming the opposite. Confirmed: Mg II observed/Table4 = 1.00 ± 0.03 across 6 epochs
and 2 SNe with different reddening factors.

### Measurement method
- Continuum: deg-1 polynomial fit to flanks at ±16000 km/s from line center (Mg II)
- Mg II window: **(-10000, +6000) km/s** default (confirmed safe per A3-A4 spill tests)
- Lya window: **(-10000, +5000) km/s** (blue edge extended from -8000 in phase 5; C III 1176 confirmed safe)
- Flux = model-independent direct integral over window (continuum-subtracted)
- No dereddening in flux measurement (Table 4 is observed; we store raw observed flux)

### Airglow (Lya only)
- `_airglow_bg`: reads sibling `*_x1d.fits` BACKGROUND column
- Subtracts geocoronal Lya if bg spike > 5 MAD; else fallback notch only if |v_geo| > 800 km/s
- Echelle: airglow in GROSS, not BACKGROUND (0.2" slit) -> always fallback path

### Shape models
Fit set to every clean-emission epoch; data picks the best by BIC:
- **gaussian** (3 params): symmetric peak, narrow cores
- **lorentzian** (3 params): IIn electron-scattering wings
- **skew-normal** (4 params): mild asymmetry
- **kwok** off-center-hole shell (5 params): broad blueshifted CSM shell (FWHM pegs 3000 km/s floor)

Selection: BIC, with dBIC >= 6 required for a k>3 model to beat the best k=3 model.

**Peg guard (phase 5 fix):** any model with a parameter within 2% of a bound is barred from winning.
If all models peg -> `best_model=null` + `shape_note`.

### P-Cygni / non-emission flags
- Photospheric epochs (Ia/Ibc): `pcygni=True, flux_reliable=False`
- Echelle excluded from emission thread (low-res pipeline not valid on echelle orders)
- Detection: `flux_reliable` requires `F / flux_err_syst >= 3`

### Lya core treatment
- Our recipe: subtract airglow, KEEP the ISM absorption trough (honest: what reached us minus airglow)
- Bostroem bridges the airglow cut region (interpolates over the core)
- The deficit vs Table 4 at Lya = core treatment difference, not flux calibration
- Optional: `F_ismcorr` = flux divided by H I transmission exp(-tau(N(HI))) -- intrinsic SN Lya before ISM absorbed it

### GGI specifics
- d41: G230L MAMA (PROPOSID 17614). Bostroem text says CCD G230LB -- their text is loose; our header is authoritative.
- d232 Mg II: default window = her Table 5 window. No narrow override (removed in phase 5).
- GGI E(B-V): total 0.154 (MW 0.070 + host 0.084).

### Key catalog files
- `catalog/emission_summary.csv`: per-epoch emission measurements
- `catalog/fuv_aod_summary.csv`: FUV AOD columns
- `reference/host_ebv.csv`: curated host E(B-V) per SN
- `linelists/csm_lines.csv`: emission line list

### Key scripts
- `scripts/emission_products.py`: full measurement pipeline (continuum, airglow, flux, shape)
- `scripts/build_products.py`: assembles catalog-level CSVs

### Do NOT reintroduce
- Narrow GGI window override (removed: it was tuned to hit 32.4 under wrong dereddening assumption)
- +15-25% offset claim (retracted: Table 4 is observed flux)
- GGI d232 Mg II "override" window (-6000, +4000): removed
- Lorentzian unfit bug (fixed: normalize all fits to O(1) before curve_fit)
- kwok-only peg guard: widened to all models in phase 5
- Width floors 500/300 km/s: replaced with 100 km/s
