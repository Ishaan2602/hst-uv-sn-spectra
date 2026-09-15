# COS reduction skill (load via read_file when doing COS pipeline work)

## Current canonical state (as of Aug 2026 full rerun)

### Workflow
1. `calcos` on each raw corrtag file (standard correction chain)
2. Edit XTRACTAB reference file: adjust B_SPEC (trace center), HEIGHT (extraction box), B_BKG (background regions)
3. Point the header `XTRACTAB` keyword at the edited copy (`lref` dir)
4. Re-run calcos; inspect the per-stripe x1dsum outputs

### Key parameters
- NUV: boxcar extraction by default
- FUV: set `XTRCTALG=BOXCAR` for custom-height extraction
- `B_SPEC`: cross-dispersion center; compare SP_LOC to actual counts-image peak (SP_ERR ~30px on some targets)
- Height tuning: 98% enclosed flux (same philosophy as STIS adaptive extrsize)

### Resampling (LOCKED)
- Per-STRIPE FCR onto the resel grid (~6 FUV pixels, ~0.06-0.07 A), preserving R~15k
- DO NOT concatenate stripes first and then resample (the old gap-bridge approach)
- Gap-aware nan-fill via `resample_fcr` in `scripts/coadd.py` (GAP_K=5 for dead-edge guard)
- NUV NUVB stripe: each cenwave puts Mg II 2800 at a different detector spot; the per-stripe approach handles this

### Known issues / do-not-reintroduce
- Concatenate-then-FCR: removed (dead-edge pixels at gap boundaries bridge the gap -> flux loss; see COS day-214 2x flux bug)
- np.interp: removed (same as STIS)
- 82x downsampling onto 1.0A COMMON_AXIS for medium-res G130M/G160M: removed (buried lines in noise)
- TRCECORR=OMIT / ALGNCORR=OMIT run: trace sits at SP_LOC (reference nominal), not the actual peak

### Key scripts
- `scripts/run_full_catalog.py`: orchestrator (--cos-only flag to rerun COS only)
- `scripts/coadd.py`: resample_fcr, per-stripe gap-handling

### Environment
- Same WSL Debian surf_uv env as STIS. Replace `oref` with `lref` for COS refs.
- `lref` = `crds_cache/references/hst/cos/` (trailing slash).
- COS data: `data/<SN>/` and `data/cos_catalog/` for multi-epoch COS targets.
- LSF: non-Gaussian wings (mid-frequency wavefront errors). EW is LSF-invariant; profile fitting is not.
