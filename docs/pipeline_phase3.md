# HST UV Supernova Pipeline - Phase 3 (full 2023ixf time series + the automation push)

Continuation of `pipeline_phase2.md`. Phase 2 reproduced the paper's 9-epoch Figure 1 (days 3-66)
from our own reductions and cleared the handoff investigations. **Phase 3 (this file) is two things:**
(1) extend 2023ixf to the *full* HST time series that exists today -- out to ~1100 days -- reduced our
own way, and (2) begin "the automation push": understanding whether MAST default extraction can be
trusted so the pipeline can loop over the whole catalog without hand-holding each epoch.

Same spirit as before: read start to finish and reproduce any piece. Day-to-day status lives in
`project_state.md`; this is the consolidated phase-3 writeup.

---

## Quick-review summary (meeting prep)

Scan this in 5 minutes; sections below go deeper.

### What is "phase 3"?
The full 2023ixf spectral time series (days 3-1094, every HST spectrum that exists), each epoch
reduced our own way, plus the first concrete evidence for why extraction cannot be blindly automated.

### What each piece is (one line each)

**1. Full time series, days 3-1094 (P3.1).** Extended `full_2023ixf.ipynb` from 9 optical epochs to
**16 epochs**: the 9 optical (days 3-66) plus 7 UV/red late epochs (days 183-1094). Late epochs are
STIS MAMA (G140L FUV + G230L NUV) with two G750M H-alpha snippets at days 913/924. One waterfall
figure, `output/2023ixf_fullseries.png`, shows the whole campaign: hot blue photospheric continuum
early, fading into UV-only spectra dominated by **Mg II] 2800 and H-alpha emission** (CSM
interaction) at late times.

**2. Only one real COS epoch (P3.2).** The landscape query returned two COS visits; the PI said there
should be one. Verified with `SkyCoord.from_name` + separation: **lf8803** (day 214, NUV G230L) is the
SN (sep 0.1"); **lf9256** (day 632, FUV G160M) is a *different* M101 target, "M101-209+312", **413.7"
(~7') away**. The cone search picked it up by proximity; it is not the SN. PI was right. Lesson for
the catalog cross-match: cone-search membership != correct target; check the coordinate/name.

**3. Late-epoch reduction, our own way (P3.3).** New `scripts/reduce_stis_generic.py`: x1d off the
crj (CCD) or the flt (MAMA), any grating, no defringe/ocrreject. MAMA is photon-counting so there is
nothing to CR-reject or defringe -- extraction straight off the flt is the whole job. Ran it over all
16 late visits (111 extractions, 0 failures on the first pass).

**4. G750M auto-trace failure = the automation motivation (P3.4).** The two G750M H-alpha epochs
(days 913/924) both failed automatic extraction: "cross correlation to locate spectrum failed --
Cannot extract." At those phases the continuum is too faint for the pipeline to lock onto the trace.
Found the trace at CCD row ~898 straight off the 2D crj and re-extracted with a fixed center
(a2center=898, maxsrch=0, extrsize=7). **This is the concrete case for the automation investigation:**
on faint late traces the default extraction silently produces nothing, so a catalog-scale loop needs
its own trace-finding before it can trust x1d.

**5. Two data-quality catches.** Day-577 visit **ofg001** failed at the telescope (calstis reports
"Total exposure time = 0", an aborted visit; MAST only ever had raw) -- that phase is covered by
ofg005 (day 619), so nothing is lost. And the rest-wavelength axis now actually deredshifts the data
(divide by 1+z, z=0.000804 for M101/SN2023ixf) instead of just being labelled "rest".

**6. Science note.** COS has slightly higher UV resolution than STIS in general -- used to resolve hot
stars (subdwarfs/OB) in the Milky Way and narrow ISM/IGM absorption; STIS trades resolution for wider
grasp + long-slit spatial coverage. In `docs/science_notes.md`.

**7. STIS default-parameter investigation (P3.5).** Ran `automation_sandbox.ipynb` against 5 SNe
(SN2024iss, SN2023ixf, SN2010jl, SN2021yja, SN2009ip) using MAST sx1 headers. Key findings: trace
center is reliable for bright sources but varies by program (not a single fixed E1 row); background
offsets are hardcoded at -300/-320 everywhere and place both bg windows off the chip for E1 sources.
Automation loop needs a profile-peak trace fallback and relative bg placement. See P3.5.

**8. Full-series figure (P3.1, updated).** `output/2023ixf_fullseries.png` now includes the COS NUV
day-214 epoch (lf8803 HASP cspec, NUVC artifact zone 2150-2350 A masked), a log-phase colorbar, and
per-segment labels so the Hα stub labels don't collide with the UV labels.

### Still open / next
- COS NUV vs HASP param sweep (`cos_sandbox`): how close can we match/beat HASP on the NUV spectrum.
- Full "our-way" extraction pass for the STIS catalog once the trace question is settled.
- The late UV epochs currently use default-center MAMA extraction; revisit centers once the automation
  loop is built.

---

## P3.1 Full time series, days 3-1094

The phase-2 figure stopped at day 66 because that is where the *optical 3-grating* coverage stops.
Everything after is UV-only, from a string of later programs. The full landscape query (one row per
visit, phase from `t_min` - explosion MJD 60082.79) shows **28 visits, phase 3 to 1094 d**, most
recent 2026-05-16. Collapsing to the epochs we can actually build our way:

- **Days 3-66 (optical, 9 epochs):** unchanged from phase 2 -- the 5 pre-reduced early spectra plus
  our of43 reductions (days 14/19/24/66).
- **Days 183-1094 (UV/red, 7 epochs):** downloaded into `data/2023ixf/epochs/day{N}_{visit}/`.
  - Day 183: of8801, FUV-MAMA G140L only.
  - Day ~311: ofbp01/02/03 (FUV G140L) + ofbp04/05 (NUV G230L) -- one epoch.
  - Day 619: ofg005, FUV+NUV.
  - Day 723: ofg003 (NUV) + ofg004 (FUV, day 724) -- one epoch.
  - Day 913: ofrs01 (G750M H-alpha) + ofrs02 (FUV+NUV, day 914) -- one epoch.
  - Day 924: ofrs03 (G750M H-alpha) + ofrs04 (FUV+NUV) -- one epoch.
  - Day 1094: ofoz01 (FUV) + ofoz02 (NUV) -- one epoch.

The G750M cenwave is 6581, a ~570 A window centred on the **H-alpha region** (~6300-6900 A) -- not a
full red spectrum, a targeted late-time H-alpha snippet. Fringing is negligible below ~7500 A so no
defringe is needed there.

**The figure.** `output/2023ixf_fullseries.png` (cell near the bottom of `full_2023ixf.ipynb`): each
epoch normalized to its own median, `log10` + rank offset, colored by phase on a log-wavelength rest
axis, with the paper's atomic-line guides. Late epochs draw only shortward of ~3200 A (UV-only) plus
the two H-alpha stubs. The story reads cleanly: smooth blue photospheric continuum early -> UV fading
-> strong Mg II] 2800 and H-alpha emission at late times (the cool-dense-shell / CSM interaction the
paper flags at day 214, now traced all the way to day 1094).

Assembly detail: per visit we median-combine our x1d over the exposures within each grating
(`OPT_ELEM` read from the header, so grating grouping is automatic), then merge gratings onto a shared
axis; visits at the same phase are merged into one epoch.

## P3.2 How many COS epochs -- verifying the PI's claim

The landscape query's 3' cone returned two COS visits, but M101 is ~28' across and full of other HST
targets. `SkyCoord.from_name('SN2023ixf')` + separation settled it:

| visit | prop | grating | phase | target_name | sep to SN |
|-------|------|---------|-------|-------------|-----------|
| lf8803 | 17497 | NUV G230L | 214 d | SN2023IXF-COS | 0.1" -- **the SN** |
| lf9256 | 17494 | FUV G160M | 632 d | M101-209+312 | 413.7" (~7') -- **not the SN** |

So there is exactly **one** COS epoch of the SN (lf8803), as the PI said. lf9256 is a different M101
target that only shared the field. This is a small but important cross-match lesson for the eventual
catalog loop (P3.4 and the SIMBAD/TNS item): proximity in a cone search is not identity.

## P3.3 Reducing the late epochs our own way

`scripts/reduce_stis_generic.py` -- a deliberately small extractor for anything that is not the
optical 3-grating CCD case:
- Finds each visit's science exposures (raw files whose TARGNAME contains "IXF").
- x1d off the **crj** if present (CCD, already CR-combined by MAST) else the **flt** (MAMA).
- No ocrreject, no defringe. MAMA is photon-counting -- no CRSPLIT, no fringing -- so extraction
  straight off the flt is the entire reduction. G750M at H-alpha is below the fringing regime.
- One CRDS bestrefs sync up front, then loops.

First pass: 111 extractions, 0 failures for the MAMA epochs; the two G750M epochs needed the
fixed-center retry (P3.4). Sanity checks: FUV G140L covers ~1120-1720 A, NUV G230L ~1570-3160 A,
G750M ~6294-6863 A, all with sensible flux.

## P3.4 The automation motivation: when default extraction silently fails

The two G750M H-alpha epochs (days 913/924) are the first clean example of why extraction cannot be
blindly automated at catalog scale. Running x1d with defaults:

```
X1DCORR  PERFORM
Warning    Cross correlation to locate spectrum failed.
ERROR:    Cannot extract.
Warning    No rows were written; no table created.
```

At day 913+ the continuum is faint enough that the pipeline's cross-correlation cannot find the trace,
so it writes **no output** -- and exits 0, so a naive batch loop would record "success" with an empty
product. Collapsing the 2D crj along dispersion showed a clean trace peak at **row ~897** (the E1
nominal position; the row-1023 spike is the detector edge, not the source). Re-extracted with
a2center=898, maxsrch=0, extrsize=7 and both epochs came out fine.

**Why this matters for phase 3's second goal.** The plan is to loop the reduction over ~140 SNe. This
epoch proves that the loop cannot trust the automatic trace-finder on faint sources: it needs to
(a) detect the empty-output / "Cannot extract" case, and (b) fall back to finding the trace itself
from the 2D image before extracting. That is exactly the question `automation_sandbox.ipynb` (P3.5)
was built to answer.

## P3.5 STIS default-parameter investigation (`automation_sandbox.ipynb`)

Before building the catalog-scale loop, we need to know whether MAST's default extraction parameters
are trustworthy or whether every visit needs its own trace-finding first. The investigation queries
MAST sx1 headers for five SNe with G230LB CCD 52X0.2 exposures: SN2024iss (in-hand), SN2023ixf,
SN2010jl, SN2021yja, SN2009ip. No re-extraction -- just read the parameters the pipeline chose.

**Data collected.** One MAST sx1 per target, downloaded to `data/automation_test/{obs_id}/`. Key
columns read: `A2CENTER` (header), `EXTRLOCY` at column 512 (actual trace row), `EXTRSIZE`,
`BK1OFFST`, `BK2OFFST`, `BK1SIZE`, `BK2SIZE`.

### Finding 1: A2CENTER varies by program and is correct when the source is bright

| target | A2CENTER | EXTRLOCY | program notes |
|--------|----------|----------|---------------|
| SN2024iss | 894.0 | 894.0 | prop 17507, recent |
| SN2023ixf | 893.6 | 893.6 | prop 17313, recent |
| SN2021yja | 893.3 | 893.3 | prop 16178, recent |
| SN2009ip  | 893.4 | 893.4 | prop 13179, recent |
| SN2010jl  | 912.9 | 912.9 | prop 12242, older |

Four of five cluster at **893-894**, the standard E1 sub-position for recent programs. SN2010jl at
**912.9** looked like a pipeline miss, but we downloaded the 2D crj and plotted the spatial profile:
the trace genuinely peaks at row 912 for that program. Different programs can place the target at
slightly different E1 offsets; the pipeline's cross-correlation found it correctly.

**Implications for automation:**
1. The pipeline trace-finder works for bright sources -- no need to override A2CENTER on bright visits.
2. But it **silently fails on faint sources** (returns exit 0 with no output, as in the G750M case).
3. There is **no universal "E1 row" to hardcode** as a fallback. The row varies by program (~893-912
   seen so far). The only reliable fallback is to find the peak of the spatial profile from the 2D
   image itself.

### Finding 2: background offsets are hardcoded and wrong for E1 CCD sources

Every single program, every SN: `BK1OFFST = -300`, `BK2OFFST = -320`. This is a fixed XTRACTAB
default -- the pipeline never adapts it per visit.

For E1-position sources (trace at row ~893-912 on a 1024-row chip):
- bg1 center at row 893 - 300 = **593**, bg2 at **573** -- both well below the science trace and for
  a near-top-of-chip source these windows fall largely off the illuminated detector area.
- In practice the background is extracting from empty chip, giving near-zero sky. For UV spectra this
  is usually harmless (sky background is tiny compared to source flux), but it is wrong in principle
  and will become visible for faint late-time sources where the source and background fluxes are
  comparable.
- The 2024iss fix (both bg regions set to ±14 px from the trace, straddling it symmetrically) is the
  right approach: place bg relative to the found trace center, not from a hardcoded table entry.

For MAMA observations (trace at row ~400-524):
- BK1OFFST = -300 places bg1 around row 100-224 (near detector bottom).
- BK2OFFST = +300 places bg2 around row 700-824 (near detector top).
- These are marginally better (both windows exist on-chip), but still not adaptive.

### What to do in the automation loop

The required logic per visit:
1. Run default x1d. If output is empty (or "Cannot extract" in the log), proceed to step 2.
2. Collapse the 2D crj/flt along dispersion. Find the trace row as the peak of the spatial profile
   (search in the expected E1 region, ~800-980 for CCD; avoid the detector edge at row 1023).
3. Re-extract with a2center = found row, maxsrch=0, extrsize=7 (or sweep for optimal).
4. Set bg offsets relative to the found trace: e.g., BK1OFFST = -(extrsize//2 + 10),
   BK2OFFST = +(extrsize//2 + 10), size = 10. Never use the -300/-320 default.
5. Save a 2D visualization (aperture + bg overlay) per extraction for manual QC.

---

## P3.6 Full-series figure refinements

After the initial waterfall (`output/2023ixf_fullseries.png`), three improvements were applied to
cell 24 of `full_2023ixf.ipynb`:

1. **COS day-214 added.** The lf8803 HASP cspec (NUV G230L, already on disk) is loaded and injected
   into `series` between the day-183 and day-311 MAMA epochs. The NUVC artifact zone (2150-2350 A) is
   masked. The COS epoch shows broad Mg II] 2800 emission clearly at day 214 -- the first sign of the
   CSM interaction before the late-time MAMA epochs dominate.

2. **Per-segment labels.** Previously each epoch got one label at the rightmost point of the full
   array; for epochs with both a UV segment and a narrow Hα stub (~6300-6900 A), both labels landed
   at the Hα end and collided with each other. Now each disconnected segment gets its own label at its
   own rightmost point.

3. **Phase colorbar.** A `ScalarMappable` colorbar on the right edge shows the turbo colormap on a
   log phase axis (ticks at 3, 10, 30, 100, 300, 1000 days). Makes the temporal direction immediately
   readable without counting trace offsets.

---

## Session log

### 7/9 -- full time series + automation groundwork
- **Science note** added: COS vs STIS resolution (`docs/science_notes.md`).
- **Full landscape query**: 28 visits, days 3-1094, most recent 2026-05-16. Established that the
  optical 3-grating coverage is days 3-66 only; everything later is UV-only MAMA + two G750M H-alpha
  CCD epochs.
- **COS epoch validated** (P3.2): one real COS epoch (lf8803, day 214); lf9256 is a different M101
  target 7' away. Confirmed the PI's count.
- **Downloaded** all 16 late visits (days 183-1094) into per-epoch dirs.
- **Reduced them our own way** (P3.3) with new `scripts/reduce_stis_generic.py`: MAMA x1d off the flt,
  G750M x1d off the crj. 111 extractions, 0 failures on the MAMA.
- **G750M fixed-center retry** (P3.4): days 913/924 failed auto trace-finding; found row ~898 off the
  2D, re-extracted with a2center=898/maxsrch=0/extrsize=7. Flagged as the concrete automation case.
- **Data catches**: day-577 ofg001 is an aborted visit (zero exposure time; calstis produces nothing),
  covered by ofg005 (day 619); rest-wavelength axis now deredshifts by 1+z (z=0.000804).
- **STIS default-parameter investigation** (P3.5): `automation_sandbox.ipynb` -- downloaded MAST sx1
  headers for 5 SNe; A2CENTER is reliable for bright sources but program-dependent (893-912); bg at
  -300/-320 is hardcoded wrong everywhere. Checked SN2010jl 2D crj -- row 912 is correct for that
  program. Wrote up full conclusions.
- **Full figure refined** (P3.6): added COS day-214, per-segment labels, log-phase colorbar.
  Saved to `output/2023ixf_fullseries.png`.
- **Docs**: `pipeline_phase2.md` P2.9 defringe recap appended; `pipeline_phase3.md` created (this
  file). `automation_sandbox.ipynb` conclusion cells trimmed to pointers.

---

## TODOs

Priority order roughly: science-blocking first, then automation, then polish.

### Immediate / next session
- [ ] **COS NUV vs HASP param sweep** (`cos_sandbox`): run a parameter sweep on the lf8803 NUV
  extraction (box height, bg placement, NUVA masking) to see how well we can match or beat the HASP
  coadd, especially in the 2150-2350 A artifact zone.
- [ ] **Automation loop prototype**: write a small loop that takes a list of obs_ids, runs
  reduce_stis_generic.py, detects empty output ("Cannot extract"), falls back to 2D profile-peak trace
  finding, then extracts with relative bg. Test on the 5 SNe from the automation_sandbox investigation.
- [ ] **2D aperture visualization per extraction**: wire `show_extraction_regions` from stis_sandbox
  into the automation loop so every extraction saves a PNG of the 2D + aperture + bg overlay. The PI
  asked for this explicitly.

### Short-term
- [ ] **MAMA late-epoch trace centers**: the 16 late MAMA epochs (ofbp*, ofg*, etc.) were extracted
  with default centers (~400-524 as seen in the of43 params output). Verify these against the 2D
  profiles for a few visits; if they drift far from the trace peak, redo with fixed centers.
- [ ] **COS NUV for the full time series**: the lf8803 coadd in `full_2023ixf.ipynb` still uses the
  HASP cspec product. Once the param sweep picks an optimal extraction, switch to our own reduction.
- [ ] **Flux system consistency**: the early spectra (days 3-11) are E(B-V)-corrected and scaled;
  the of43 reductions (days 14-66) and all late epochs are raw flux. The waterfall uses normalized
  spectra so it's fine visually, but absolute-flux comparisons need all epochs on the same system.
- [ ] **ISM narrow-line masking**: the paper removed narrow ISM absorption features. Apply a
  line-list mask (sigma-clip or simple wavelength exclusion) to the early spectra and our reductions.

### Longer-term (catalog scale)
- [ ] **Full catalog extraction loop**: extend the automation loop prototype to all 138 STIS SNe
  in `output/uv_sn_catalog.csv`. Handle MAMA and CCD branches, detect empty outputs, save 2D viz.
- [ ] **COS branch in pipeline.py**: the main pipeline only handles STIS today; add COS calcos +
  custom XTRACTAB editing.
- [ ] **SIMBAD/TNS cross-match for the catalog**: flag any target whose cone-search position
  doesn't match the SIMBAD/TNS classification (the lf9256 / M101-209+312 case is the model).
- [ ] **SNR-weighted coadd**: inverse-variance weighted combine instead of nanmedian, so faint
  legs (COS NUV tails, G750M late epochs) contribute with appropriate down-weighting.
- [ ] **CMFGEN model comparison**: the end science goal -- compare our reduced spectra to CMFGEN
  model grids to extract mass-loss rates and CSM geometry.
