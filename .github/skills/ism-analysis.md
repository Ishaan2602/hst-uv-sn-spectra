# ISM analysis skill (load via read_file when doing ISM absorption work)

## Current canonical state (phase 5, Sep 2026)

### Overview
Narrow foreground ISM lines (MW + host) printed on the SN UV continuum backlight.
Three analysis products: `ism.py` (automated CoG), `lya_nhi.py` (N(H I)), `absorption_products.py` (joins them).

### EW measurement
- Local deg-1 continuum to flanks, masking other known ISM cores
- **Joint deblend** for crowded groups (Mn II + Fe II picket fence 2577-2606; Cr/Zn clump 2056-2066):
  one shared continuum + sum of Gaussians, centers pinned at rest wavelengths, one shared instrumental width
- Two methods: (A) direct integral of (1-f/cont), (B) Gaussian fit area. Big A-B disagreement = saturation/blend.
- Error: Monte Carlo photon + continuum error; continuum dominates 10-25x for strong lines

### Detection gate (option A, phase 5)
- Keep a line if `EW / EW_err_tot >= 2`, where `EW_err_tot = hypot(photon, continuum)`
- Old gate (photon only) allowed continuum artifacts to pass -- fixed in phase 5

### Curve of growth
- Universal curve: `W/lam = (2b/c) * F(tau0)`, `tau0 = 1.497e-15 * N * f * lam / b`
- Fe II anchors b (most lines over widest f-range). At G230LB resolution b is instrument-dominated (~60 km/s).
- b-N degeneracy breaks at the flat part; Fe II pins it; other ions read off same b.
- b range plausibility: 30-120 km/s. logN(FeII) < 16 gate (> 16 = Ia photosphere / dense CSM / remnant).
- Need >= 4 Fe II lines to anchor reliably.

### Adopted flag (phase 5 fix)
- `adopted` = mean over all anchored epochs with the epoch scatter as the error (NOT single-epoch)
- Foreground ISM can't change between epochs; scatter proves photon error ~10x too small
- Reliability guards: drop Ia epochs, logN_err > 0.3, logN > 16, n_fe < 4

### Apparent optical depth (AOD, FUV only)
- Only for resolved FUV (G130M/G160M, R~15k-19k). G140L R~2000 is NOT resolved.
- `tau(v) = ln(I_cont/I_obs)`, `N = 3.768e14 / (f * lam) * integral tau dv` (Savage & Sembach 1991)
- Significance gate >= 3 sigma. Saturated lines -> lower limit. Non-detections -> upper limit.
- Same-ion cross-check: Fe II 1608 (AOD) vs Fe II CoG (NUV) -- the independent check.
- 2010jl result: NUV Fe II column (15.87) contradicted by AOD upper limit (< 14.4) -> CSM/photospheric, not ISM.

### N(H I) from damped Lya
- Voigt-Hjerting damped profile; b fixed (damping wings are b-independent at high N)
- Only for photospheric-backlight epochs with G140L coverage (not CSM-emission-contaminated)
- `F = cont(w) * exp(-tau)`. v_abs fixed at 0 (host frame).
- `lya_nhi.py` b = 25 km/s (damping regime). CHI2_GATE is an rms gate, not true chi2.

### Metallicity
- `[X/H] = log N(X) - log N(H I) - log(X/H)_sun`
- Zn tracks gas metallicity (undepleted); Fe reads low (depleted onto grains)
- `[Fe/Zn]` ~ -1.1 dex = dust depletion (warm-disk/halo range)
- Only for sightlines with BOTH a reliable NUV CoG epoch AND a Lya N(HI) detection

### Key catalog files
- `catalog/ism_cog_summary.csv`: per-epoch CoG results + adopted flags
- `catalog/lya_nhi_summary.csv`: N(HI) + v_abs per source
- `catalog/absorption_summary.csv`: joined columns + metallicity
- `reference/ism_columns.csv`: hand-curated reference values (4 well-measured sightlines)
- `linelists/ism_lines.csv`, `linelists/ism_lines_fuv.csv`: line lists
