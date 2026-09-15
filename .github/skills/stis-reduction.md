# STIS reduction skill (load via read_file when doing STIS pipeline work)

## Current canonical state (as of Aug 2026 full rerun)

### DQ bitmask (LOCKED)
- DROP: bits 512 (bad ref pixel), 256, 4
- KEEP: bit 16 (high dark rate = real data; removing it killed the Mg II doublet -- do NOT reintroduce)

### Extraction flow per detector type

**STIS CCD** (G230LB, G430L, G750L, G230M, G430M, G750M):
1. `ocrreject` (crj from multiple flt exposures) -- ONLY if CRSPLIT > 1
2. `x1d` off the crj (or off flt if CRSPLIT=1)
3. G750L only: `normspflat` -> `mkfringeflat` -> `defringe` (contemporaneous CCDFLAT) -> `x1d`
4. Adaptive extrsize: 98% enclosed flux; per-grating clamps; fallback to default if trace fails

**STIS MAMA** (G230L, G140L, G230M, FUV gratings):
- Photon-counting; no CR rejection step
- `x1d` directly off flt
- No defringe

**STIS echelle** (E140M, E230M):
- Cross-dispersed; trace finder does NOT apply
- `stis_extract.is_echelle()` detection -> skip finder, run x1d with SPTRCTAB defaults
- Multi-order 1D written under `STIS/ECHELLE/`

### Resampling (LOCKED -- no np.interp)
- Low-res (G230LB/G430L/G750L/G230L/G140L): FCR -> `COMMON_AXIS` (0.9-1.7x native, harmless)
- Medium-res (G750M/G430M/G140M/G230M) + echelle: FCR -> per-grating resel grid (~2px native)
- COS FUV: FCR -> resel grid (~6 px, ~0.06-0.07 A), preserving R~15k

### Known issues / do-not-reintroduce
- DQ 16+512 simultaneous drop: removed (killed 33-52% of NUV pixels, destroyed Mg II doublet)
- np.interp: removed (doesn't conserve flux, wrong error interpolation)
- Inter-grating scaling (_scale_leg): removed (calstis flux cal already agrees ~1%; scaling suppressed flux)
- Low-res FCR is a no-op vs np.interp -- confirmed on broad sample; medium-res/echelle is the win
- Faint trace (low-SNR G230LB): auto trace-find fails; retry with fixed E1 trace center (a2center=893.5, maxsrch=0)
- G750L defringe: use the NARROW slit flat (52x0.1 or 0.3x0.09 for E1), not the science slit

### Key scripts
- `scripts/stis_extract.py`: wrapper around calstis x1d, ocrreject, defringe
- `scripts/run_full_catalog.py`: orchestrator (--workers, --cos-only, --presync flags)
- `scripts/coadd.py`: FCR resampling + ivar_combine (sigma-clip active for >= 3 exposures)

### Environment
- WSL Debian, conda env `surf_uv`. Run via:
  `wsl.exe -d Debian -- bash -lc 'source ~/miniforge3/etc/profile.d/conda.sh && conda activate surf_uv && export oref=.../crds_cache/references/hst/stis/ && python /mnt/c/.../WORK/scripts/<script>.py'`
- CRDS cache: `crds_cache/`. STIS refs -> oref (trailing slash). COS refs -> lref (trailing slash).
- calcos: ALWAYS delete output dir before rerun (calcos needs an empty outdir).
- git.exe not plain git in WSL (git status stalls on NTFS).
