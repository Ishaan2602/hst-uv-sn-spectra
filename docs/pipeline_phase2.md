# HST UV Supernova Pipeline - Phase 2 (time series + multi-epoch 2023ixf)

This is the continuation of `extraction_pipeline_guide.md`. That first guide ("phase 1") covers the
foundations: the environment, STIS + COS extraction, G750L defringing, coadds, the supernova
discovery query, and the first automated pipeline. **Phase 2 (this file) is the science push:**
reproduce the paper's full time-series spectral evolution of SN2023ixf (Bostroem et al. 2024,
Figure 1) from our own reductions, then generalize the machinery to the rest of the catalog.

Same spirit as the phase-1 guide: one place you can read start to finish and reproduce any piece.
Written deliberately with lots of context and plain language. Day-to-day status still lives in
`project_state.md`; this is the consolidated phase-2 writeup.

---

## Quick-review summary (meeting prep)

What we built in phase 2, written so you can scan this in 5 minutes. Sections below go deeper on each.

### What is "phase 2"?
We now have the full Figure 1 time-series working (see P2.2). The main new deliverables over phase 1:
- 9-epoch log-log stacked spectrum from day 3.25 to day 66.25 (`output/2023ixf_timeseries_full.png`).
- Proper reduction of the 4 main epochs (our own extraction, not MAST defaults).
- A set of targeted investigations and corrections to issues flagged in the original handoff notes.

### What each piece of work actually is (one line each)

**1. Epoch naming explainer (P2.1).** Wrote a plain-language guide to `IPPPSSOOT` rootnames, the 3
gratings, what "one epoch" means, and how a "day X" phase maps to a visit -- now in this doc so you
don't have to decode `oezt01040` from scratch every time.

**2. Full Figure 1 time series (P2.2).** The paper's 9-epoch NUV-NIR stacked evolution reproduced in
`full_2023ixf.ipynb`. The 5 early epochs (days 3-11) use the pre-reduced 1D spectra in
`data/earlytime_2023ixf`. Days 14/19/24/66 are our own STIS reduction from the main program
(**prop 17313**, not 17315 as we thought). Each of those uses `reduce_epoch_ts.py` (crj-aware, G750L
defringed with the contemporaneous CCDFLAT, day-66 G230LB needed a manual trace center).

**3. COS NUV is a late-time observation, NOT contemporaneous (P2.2 sub-note).** The COS NUV G230L
(lf8803, prop 17497) that was being overlaid on the day-3.5 STIS SED sits at **day ~214** -- a
~210-day mismatch. It was excluded from the coadd for the wrong reason (called "faint, SNR issue");
the real reason is it's a completely different epoch. Moved to its own labeled cell in the notebook.
The figure it shows at day 214 is actually interesting: clear broad Mg II 2800 emission (late-time
cool-dense-shell interaction the paper discusses). SNR-weighted coadding is still worth knowing
about in principle, but it's not the issue here.

**4. Inter-grating flux scaling (P2.2).** Added to the main-epoch coadd: G430L and G750L are aligned
to G230LB by a constant percentage in the overlap windows before the median. Unsaturated epochs give
factors 0.99-1.08 (gratings agree), confirming the oezt01 factor-of-0.6 was saturation. Day-66's
faint blue anchor is guarded to 1.0.

**5. G750L defringing -- verified working, pushed to its ceiling (P2.4).** From the full paper PDF:
our chain (contemporaneous 0.3x0.09 CCDFLAT -> `normspflat` -> `mkfringeflat` -> `defringe` ->
`x1d`) matches the paper exactly. The defringe cuts residual fringe scatter 30-48% in the red.
After that, pushing further (red-targeted RMS, finer grid) gives negligible improvement: the residual
is ~2x the photon noise in every red window and hitting the **fundamental ceiling of the single
global shift+scale model** -- the fringes are wavelength-dependent and can't be fully removed by one
pair of numbers. Normal; the paper accepts the same. (`output/defringe_check_day14.png`.)

**6. COS NUVB box centering -- was a frame mix-up, not a real problem (P2.4 notes).** The apparent
"peak on the left edge of the box" was because the plot was histogramming `YCORR` while `B_SPEC`
lives in `YFULL`. In YFULL the boxes sit 0-3 px from the stripe peaks. Fixed the plot.

**7. COS FUV box height / Ly-alpha x2 (P2.4 notes).** The ~2x Ly-alpha emission scaling across the
height sweep is a real source-flux effect: the COS FUV cross-dispersion profile is broad, so a
narrow box clips flux. Continuum x1.95, Ly-alpha peak x2.23 from h11 to h45, both plateauing by
h35-45. Takeaway: use a box on the plateau (default h41 is near-optimal) or aperture-correct.

**8. 2010jl STIS + FUV.** Downloaded and reduced the obk002 epoch (prop 12242, G230LB + G430L, no
G750L). COS FUV G130M already in hand. The two gratings agree at the 3000 A splice; strong emission
lines throughout (it's a Type IIn). `output/2010jl_sed.png`. Shows the pipeline generalizes to a
2-grating STIS case and a second SN.

### Still open / to clarify at the next check-in
- Flux column for the early-epoch spectra: `FLUX` vs `FLUX_corr` (dereddened). Currently using
  `FLUX`; the two look nearly identical for 2023ixf (modest extinction).
- The early spectra still have narrow ISM absorption lines in them; the paper removed those. A
  sigma-clip / line-list mask step is queued but not yet applied.
- Rest-frame assumption: we label the x-axis "Rest Wavelength" treating Wynn's files as rest-frame.
  M101 z~0.0008 makes the shift <1 A -- negligible at our scale, but worth confirming.
- COS NUV coadd: our 6-exposure median and HASP are within ~1% noise outside the artifact zone
  (ratio 1.01). The key difference is whether the NUVA 2150-2350 A artifact gets removed; HASP
  removes it via a flux-deviation filter; our naive median does not. Investigation + fix in P2.6.
- COS FUV h=45: the plateau is already reached by h=25-35; h=45 gives 0% more flux than h=35.
  Default h=41 is solidly on the plateau. Full numbers in P2.7.

---

## P2.0 What phase 2 is


The headline goal: build the "Figure 1" time-series plot for SN2023ixf. That figure shows the
NUV-to-NIR spectrum at **9 epochs** (days 3.25 -> 66.25) stacked in time, so you can literally watch
the supernova change: it starts as a hot blue flash-ionized continuum with narrow emission lines,
then the UV fades and it settles into a classic Type II supernova with broad P Cygni absorption/
emission profiles (H-alpha, Ca II, and a late broad Mg II emission by day 66).

Two data sources feed the plot, and the split matters for who does the reduction:

- **Early epochs (days 3.25, 4.25, 5.25, 8.25, 11.25):** Wynn already reduced these. The early
  GO-17205 data had CCD saturation and pointing problems (the SN was misaligned in the slit), so he
  did the careful manual extraction and flux-calibration workarounds himself. He handed us 5
  ready-made 1D spectra in `data/earlytime_2023ixf/` (one text file per epoch). **We do NOT
  re-reduce these** - we just read them in.
- **Main epochs (days 14, 19, 24, 66):** these are the clean primary-program observations (the
  52x0.2 slit at the E1 position). **These we reduce ourselves**, end to end (STIS three-grating
  extraction + G750L defringe + coadd), and append to the same plot.

Once the 2023ixf time series looks right, the exact same plotting + reduction machinery generalizes
to the other 139 SNe in `output/uv_sn_catalog.csv`.

---

## P2.1 HST dataset names and "epochs" -- a full explainer

This section exists because the naming trips everyone up at first, and it keeps coming up. Read it
once and the obs IDs stop looking like random noise. Everything below is grounded in our actual
SN2023ixf files (the numbers are read straight from the FITS headers, not guessed).

### What is an "epoch"?

An **epoch** is just *one observation date* -- a single trip HST made to point at the target. Each
time HST visits SN2023ixf it takes a small set of exposures (normally one per grating), and that
whole set, taken over a few hours on one date, is **one epoch**.

In supernova work we don't label epochs by calendar date, we label them by **phase**: the number of
days since the star exploded ("days post-explosion"). SN2023ixf exploded around 2023-05-18/19, so a
spectrum taken ~2 weeks later is at phase ~14 days and we call it the "day 14 epoch." The paper's
Figure 1 stacks 9 epochs: days 3.25, 4.25, 5.25, 8.25, 11.25, 14.25, 19.25, 24.25, 66.25. Each line
on that plot is one epoch = one visit = one date. Plotting them stacked in time is how you watch the
SN evolve.

In HST's own bookkeeping an epoch almost always lines up with a **visit** (defined just below). So
for us "epoch," "visit," and "one date of data" mean effectively the same thing.

Two different HST programs took these 9 epochs, and that split matters for who reduces what:
- **GO-17205** = the *early* program (days 3-11; our `oezt01`, `oezt02`, `oezt03`... visits). These
  had pointing and CCD-saturation problems. Wynn is handing us the already-reduced 1D files for
  these (in `data/earlytime_2023ixf/`, one text file per epoch), so we do NOT re-reduce them.
- **The main program (days 14-66)** = a separate proposal (the paper's primary observations, taken
  at the E1 slit position). These are the clean ones, and the ones WE reduce ourselves. (We still
  have to pull their obs IDs from MAST -- see the open items in `project_state.md`.)

### Reading an HST dataset name (the "rootname")

Every HST dataset has a 9-character name called the **rootname**, e.g. `oezt01040`. It looks random
but it isn't -- it packs in the instrument, the program, the visit, and the exposure. The layout is
`IPPPSSOOT`, i.e. read left to right:

```
  o    ezt   01   04   0
  |     |    |    |    |
  I    PPP   SS   OO   T
instr  prog  visit exp  member
```

- **I** (1 char) = the **instrument**:
  - `o` = STIS   (all our SN2023ixf STIS files start with `o`)
  - `l` = COS    (our COS files are `lf8803...`)
  - `j` = ACS, `i` = WFC3, `u` = WFPC2 (not used in this project)
- **PPP** (3 chars) = an **encoded program ID** (`ezt`). This is a compressed code STScI assigns;
  you cannot eyeball it back into the real proposal number. Do not try to turn `ezt` into 17205 by
  hand -- just read the real number from the `PROPOSID` header keyword (we did: `ezt` = program
  **17205**).
- **SS** (2 chars) = the **visit** number (`01`) -- this is the epoch. All four `oezt01xxx` files
  share `01`, so they are the same visit = same epoch = same date.
- **OO** (2 chars) = the **observation/exposure** within that visit (`04`, `0h`, `0e`, `0d`). These
  are base-36 (digits 0-9 then letters a-z), and the scheduler assigns them, so they are NOT tidy
  01/02/03 -- our visit happens to use 04, 0d, 0e, 0h. Each one is a separate exposure (here, a
  different grating).
- **T** (1 char) = a trailing association/member code (`0`); you can mostly ignore it.

So `oezt01040` reads as: STIS, program 17205, visit 01, exposure 04. And `lf8803010` reads as: COS,
program 17497 (`f88`), visit 03, exposure 01.

### The four files in our oezt01 epoch (worked example)

This is exactly what the shorthand "prop 17205 # oezt01040 g230lb (769s), oezt010h0 g430l (574s),
oezt010e0 g750l (492s), oezt010d0 ccdflat" is listing. All four are program 17205, visit 01 (one
epoch, one date); every value below is read straight from the headers:

| rootname   | OO | what it is       | grating | slit   | exp time | target    |
|------------|----|------------------|---------|--------|----------|-----------|
| oezt01040  | 04 | science spectrum | G230LB  | 52X0.2 | 769 s    | SN2023IXF |
| oezt010h0  | 0h | science spectrum | G430L   | 52X0.2 | 574 s    | SN2023IXF |
| oezt010e0  | 0e | science spectrum | G750L   | 52X0.2 | 492 s    | SN2023IXF |
| oezt010d0  | 0d | fringe flat      | G750L   | 52X0.1 | 50 s     | CCDFLAT   |

How to read that table:
- **The number in parentheses is the exposure time** (`TEXPTIME`, in seconds). Longer exposure =
  more photons collected = better signal-to-noise, at the cost of telescope time. So this epoch
  spent 769 s in the blue, 574 s in the optical, 492 s in the red.
- **The three science rows are the three gratings** -- the same patch of sky dispersed three ways
  to tile the whole NUV-to-NIR wavelength range (see the grating explainer below). "One epoch"
  means these 3 gratings taken together.
- **The fourth row (oezt010d0) is not the supernova at all** -- its target is `CCDFLAT`, a tungsten
  lamp exposure taken through the narrow **52X0.1** slit in the same visit. That is the **fringe
  flat** we divide into the G750L science to kill the red-end fringing (phase-1 guide Section 6). It
  sits in the same visit on purpose, so it samples the same detector fringe pattern only minutes
  apart. Rule of thumb: the `...0d0` CCDFLAT exposure in a visit is that visit's G750L fringe flat.

### The three STIS gratings (why there are three files per epoch)

STIS can only disperse one wavelength band at a time, so covering NUV -> NIR takes three separate
setups. Picture them as three tiles laid end to end:

- **G230LB** -- near-UV / blue, roughly **1650-3100 A**. (The "B" means the CCD version; there is
  also a MAMA-detector G230L.)
- **G430L**  -- optical, roughly **2900-5700 A**.
- **G750L**  -- red / near-IR, roughly **5250-10250 A**. This is the one with fringing that needs
  the CCDFLAT.

They overlap a little at the edges (~3000 A and ~5600 A), and that overlap is what lets us stitch
them into one continuous spectrum in the coadd (phase-1 guide Section 7). The **"L"** in each name
means low-resolution (a low-dispersion grating), which is the right trade for faint SNe: we want
throughput (catch photons) more than fine spectral detail.

### One quirk that changes how we reduce: number of exposures (CRSPLIT)

How many exposures an epoch took per grating decides our reduction path, so it is worth checking:

- Our early oezt01 epoch has **CRSPLIT=1** -- a *single* exposure per grating. With only one frame
  there is nothing to compare it against to spot cosmic-ray hits, so `ocrreject` cannot run (it
  errors with "needs more than one input"). We extract straight off the `flt`. This single-exposure
  choice is part of why the early GO-17205 data is fragile.
- The **main epochs (days 14-66) deliberately took >=4 exposures per grating**, so the pipeline can
  auto-reject cosmic rays by stacking them. Those produce a cosmic-ray-cleaned `crj`/`sx1` product,
  and our reduction there runs `ocrreject` first. So expect a slightly different file set and a
  slightly different code path once we start on those epochs.

### Quick decoder for the names you'll actually see

- `oezt01040` = STIS, prog 17205, visit 01, exp 04 -> **G230LB** science (an early epoch).
- `oezt010h0` = same visit, exp 0h -> **G430L** science.
- `oezt010e0` = same visit, exp 0e -> **G750L** science.
- `oezt010d0` = same visit, exp 0d -> **CCDFLAT** (the G750L fringe flat).
- `oezt03...` = STIS, prog 17205, **visit 03** -> a *different* early epoch (a different date).
- `lf8803010` = COS, prog 17497, visit 03, exp 01 -> the SN2023ixf COS NUV G230L visit.
- For any file, get the truth (program, target, grating, slit, exposure time) from the `PROPOSID`,
  `TARGNAME`, `OPT_ELEM`, `APERTURE`, `TEXPTIME` header keywords. Only trust the name itself for the
  instrument letter and the visit number; never try to decode the program or grating from it.

---

## P2.2 The time-series plot (paper Figure 1 reproduction)

Lives in `notebooks/full_2023ixf.ipynb` (extended additively at the bottom -- the single-epoch work
above it is untouched so we can always backtrack). Three new cells: a markdown header, an epoch
loader, and the stacked plot.

**What it reproduces.** Bostroem et al. 2024 Figure 1: the NUV-to-NIR spectrum at each epoch,
plotted log-log and offset vertically so the epochs stack (earliest at the top, latest at the
bottom), with faint dotted vertical guides marking the atomic transitions the paper flags
(N IV 1718, C III 2297, Mg II 2796/2802, He II 3203, He II 4685.5/4860, C IV 5801/5812, H-alpha 6563,
He I 6678.1, N IV 7109/7122, Ca II 8498/8542/8662).

**The data going in (right now).** Only the 5 early epochs (Wynn's pre-reduced files). Each file
(`data/earlytime_2023ixf/scaled_ebv_corr_combined_spec_{1st..5th}_hst.txt`) is a plain text table
with 5 columns: `WAVELENGTH FLUX ERROR ORIGFLUX FLUX_corr`, covering the full ~1673-10229 A range
(so each file is already the coadded 3-grating spectrum). The 1st..5th files map to days
3.25/4.25/5.25/8.25/11.25 respectively (confirmed as a working assumption).

**How the plot is built (so it extends cleanly).**
- `epoch_plan` is the full 9-epoch list (day label, color, filename). Colors run early(top)->late
  (bottom): darkblue, indigo, darkviolet, magenta, deeppink, crimson, red, darkorange, gold -- the
  paper's blue->gold scheme. The 4 main epochs currently have `None` for the filename (no data yet).
- The loader reads each available file into an `epochs` list of `(day, color, plan_idx, wvl, flx)`.
  When we reduce days 14-66 ourselves, we just append their `(day, color, plan_idx, wvl, flx)` tuples
  to `epochs` and the same plot cell picks them up -- no plot rewrite needed.
- The plot offsets each epoch by `(nslot - 1 - plan_idx) * step` in log-flux, so every epoch has a
  reserved vertical slot in the final 9-epoch layout (the y-view auto-fits to whatever is loaded).
- Output saved to `output/2023ixf_timeseries.png`.

**Column choice (`fluxcol`).** The loader has a `fluxcol` switch: column 1 = `FLUX` (the scaled,
combined flux) or column 4 = `FLUX_corr` (the same but E(B-V) extinction-corrected / dereddened).
See the open question below -- we are comparing the two before committing.

**The COS NUV G230L is a LATE observation (day ~214), NOT contemporaneous.** This one bit us for a
while. The single-epoch coadd cell near the top of `full_2023ixf.ipynb` overlays a "COS NUV G230L"
spectrum on the early STIS. That COS data (`lf8803`, program **17497**, the HASP `lf88_cspec` coadd)
is **not** from the same epoch -- it sits at **phase ~214 days** after explosion, while the STIS
oezt01 it is drawn next to is at **day 3.5**. That is a ~210-day gap: two completely different phases
of the supernova.
- *How we checked:* read `EXPSTART` from the COS corrtag and the STIS raw headers and differenced
  against the explosion MJD (~60082.79). COS NUV = MJD 60296.9 = **day 214.1**; STIS oezt01 = day 3.5;
  the main of43 epochs span days 14.4-66.5. So the COS NUV is later than *every* STIS epoch we have,
  day 66 included.
- *Why it matters:* the paper's Figure 1 (days 3.25-66.25) is **all STIS** -- there is no COS in it.
  Our COS NUV belongs to a separate, much later program and must not be mixed into the early/main STIS
  epochs.
- *Correcting an earlier wrong explanation:* the old comment in that cell blamed the COS/STIS
  disagreement on the COS being "faint at this epoch" and floated **SNR weighting** as the fix. That
  framing is wrong -- the real reason the COS NUV does not belong in the coadd is the **~210-day time
  mismatch**, not signal-to-noise. (SNR / inverse-variance weighting is still a fine tool to keep in
  mind for genuinely-contemporaneous legs of differing depth; it just is not the issue here.)
- *Lesson for the pipeline:* when overlaying or combining spectra from different instruments/programs,
  check `EXPSTART` (the phase) first. The cell even carried a stale `# make sure obs date within ~1
  day` reminder that we had never actually enforced.
- *In the notebook now:* the day-3.5 single-epoch coadd is pure STIS, and the COS NUV lives in its own
  cell labeled `day ~214`. At that late phase it clearly shows the broad **Mg II 2796/2802 emission**
  (the cool-dense-shell interaction signature the paper discusses for late times), which only
  reinforces that it is a different regime from the day 3-66 STIS spectra.

---

## P2.3 Open questions / to clarify with Wynn

These are the things we deliberately did NOT guess on. Flag them at the next check-in.

1. **Which flux column does Figure 1 use -- FLUX or FLUX_corr?** Wynn's files carry both the scaled
   flux and the extinction-corrected (dereddened) flux. Dereddening mostly boosts the UV, so it
   steepens the blue continuum slope, which changes how the early epochs look. We show the two side
   by side in the notebook and are holding on a default until Wynn confirms which the paper figure
   uses. (Current skeleton uses `FLUX`.)

2. **Narrow ISM absorption lines are still in Wynn's spectra.** The paper caption says the narrow
   interstellar-medium (ISM) absorption lines were deliberately removed to show the intrinsic SN
   spectra. Wynn's files still contain them (the sharp downward spikes near Mg II ~2800, C III, etc.).
   **TODO (queued):** add a separate step that clips/masks the narrow ISM lines (sigma-clip on the
   narrow residuals, or a line list), applied after loading and before plotting. Left in for now so
   the skeleton stays faithful to the raw input.

3. **Rest frame vs observed frame.** We are treating Wynn's spectra as already being in the SN rest
   frame and label the x-axis "Rest Wavelength (A)" to match the paper. M101's redshift is tiny
   (z ~ 0.0008, a shift of <1 A even at H-alpha), so this barely matters at the plotting scale, but
   it is an assumption -- confirm with Wynn whether his files are rest-frame or observed-frame.

---

## P2.4 G750L defringing -- what each step does and how we validate it

Fair question that came up: are we just checking that the fringe-flat slit is narrower, or do we
actually understand the defringe steps? Here is the full picture, cross-checked against the paper PDF.

**What the paper does (Bostroem et al. 2024, Sec 2, verbatim gist).** The G750L "suffers fringing at
the reddest wavelengths, which can be corrected using a contemporaneous fringe flat. We obtained a
fringe flat with each visit using the 0.3x0.09 aperture, as recommended for the E1 aperture position.
The fringe flat was applied to the G750L observations using `stistools.defringe` and 1D spectra were
extracted using `stistools.x1d`." That is exactly our chain. Two requirements they call out, both of
which we now verify:
- **Contemporaneous** -- a fringe flat taken in the *same visit* as the science. Fringes drift with
  detector temperature and the exact light path, so a same-visit flat matches best.
- **Narrow aperture** -- 0.3x0.09 (E1) or 52x0.1 (at the 52x0.2 position), narrower than the science
  slit, so the flat samples the fringe pattern cleanly without the extended source smearing it.
(The paper also confirms the day-66 G230LB manual extraction we had to do: at that phase "the UV flux
was very low at the blue end, preventing the automatic identification of the spectrum location by the
pipeline ... manually extracted with `stistools.calstis`, using the red end to identify the trace" --
which is precisely our fixed-a2center retry.)

**What fringing physically is.** At red wavelengths the STIS CCD is thin enough that light reflects
between the front and back surfaces and interferes with itself -- a fixed, wavelength-dependent
modulation (a few to ~15% ripple past ~7500 A, worst past ~9000 A). It is **multiplicative** (the true
spectrum times a fringe pattern), so the fix is to measure that pattern (the fringe flat) and **divide**
it out.

**What `stistools.defringe` actually does (the 3 steps, and why each matters).**
1. **`normspflat`** (raw flat -> `_nsp`): calibrates the fringe flat (bias/dark/flat via basic2d) and
   then **normalizes** it -- divides out the tungsten lamp's own smooth spectral shape so all that is
   left is the high-frequency fringe pattern (values near 1.0 with the ripple on top). *Why:* you want
   to divide the science by the fringes only, not by the lamp's color.
2. **`mkfringeflat`** (science + `_nsp` -> `_frr`): solves for the best **shift** (sub-pixel, along
   dispersion) and **scale** (fringe amplitude) that make the flat's fringes line up with the science
   exposure's fringes, minimizing the residual RMS. *Why:* the flat was not taken at the exact same
   detector temperature/position as the science, so its fringes are slightly shifted and have a
   slightly different amplitude; you register them before dividing. **This is the step to watch:** if
   the best shift or scale lands at the *edge* of the search range, the solver never found a real
   minimum and the correction is untrustworthy. (Our early oezt01 epoch hit the scale ceiling at 1.2
   -- a red flag; the main of43 epochs converge to interior values, e.g. day 14 shift +0.30 px, scale
   1.12.) We widen the shift search (`beg_shift=-2, end_shift=1`) so it does not clip.
3. **`defringe`** (science / `_frr` -> `_drj`): divides the science frame by the registered fringe
   flat, removing the ripple. Output suffix is `_drj` for G750L whether the input was a `crj` or `flt`.

**How we validate it (beyond the slit check).** The slit-width check is necessary but not sufficient.
We also:
- **Watch the `mkfringeflat` convergence** -- shift and scale interior to their ranges, not at an edge
  (see step 2).
- **Measure the fringe residual before vs after.** On the day-14 G750L we take the RMS of
  (spectrum / median-smoothed spectrum) in red windows: 8000-8500 A 0.058 -> 0.040, 8500-9000 A
  0.056 -> 0.039, 9000-9700 A 0.125 -> 0.065, 9700-10200 A 0.188 -> 0.140 -- a 30-48% reduction, and
  the band is visibly smoother (`output/defringe_check_day14.png`, defringed vs undefringed). The
  leftover wiggles past ~9700 A are the hardest reddest fringes plus photon noise where the source is
  faint; residuals there are normal for G750L.

**One related paper detail worth remembering.** The paper flags that several narrow "emission"
features present in all spectra are actually **high-dark-rate pixels** (DQ flag 16), not real lines,
and are removed. We already mask DQ 16 + 512, so we inherit this -- but keep it in mind before
believing any narrow red feature.

---

## P2.6 COS NUV coadd: our result vs HASP, and the 2150-2350 A artifact

Two COS NUV findings that are worth understanding:

### Our coadd vs HASP: nearly identical outside the artifact zone

The 6-exposure NUV coadd (G230L, lf8803 visits) matches HASP's well across most of 1650-3200 A.
Noise outside the artifact zone (2400-3100 A): **ours = 6.36e-16, HASP = 6.28e-16, ratio 1.01** --
within 1% of each other. The earlier observation that "our coadd has slightly less noise" was
measuring a noise fluctuation; with the full dataset they are equivalent. The key difference between
our coadd and HASP's is not noise -- it is whether the NUVA artifact is removed.

### The 2150-2350 A artifact: NUVA stripe contamination

Our coadd has a noisy spike/shelf between roughly 2150-2350 A that doesn't appear in the HASP coadd.
This is the NUVA stripe leaking in.

**Why it happens.** The COS NUV detector has three spatial stripes (NUVA, NUVB, NUVC), each covering
a different wavelength range depending on the cenwave. At G230L cenwave 2635, the three stripes cover:
- NUVC: roughly 1740-2060 A
- NUVB: roughly 2080-2590 A (the primary science stripe)
- NUVA: roughly 2540-3200 A
There is real wavelength overlap between stripes (the ranges above are schematic), and the NUVA
stripe at G230L is affected by a localized count-rate artifact (likely from detector-edge or
field-of-view contamination) that produces elevated noisy flux in that ~2150-2350 A window.

**How HASP handles it.** HASP applies a per-exposure flux-deviation filter: for each 20 A bin, it
checks whether a given exposure's flux deviates from the running coadd by more than ~5%. If it does,
that exposure/bin combination is rejected before the final combine. Because the NUVA artifact is
limited to certain exposures and a specific wavelength range, this filter catches and drops it.

**What we do (currently: nothing, and why it matters).** Our naive `nanmedian` across all 6 exposures
keeps that stripe in unless we mask it ourselves. It makes our coadd look worse than HASP's in that
one window, even though we're working with the same raw data. The fix is to add our own deviation
filter -- something like: for each wavelength bin, compute the median across exposures, flag bins
where any single exposure's flux is more than N sigma from the median, then exclude those bins from
the combine. This is queued as a future improvement; it would let our coadd beat HASP in every
window, not just some.

**Key lesson.** This is a good illustration of the difference between a naive median coadd and a
quality-filtered coadd. The data itself is fine; the artifact is real but limited in extent, and
it can be removed programmatically once you know which bins to flag.

---

## P2.7 COS FUV box height deeper: does h=45 work?

From the FUV height sweep (SN2010jl G130M, `cos_sandbox.ipynb`):

**What we already know.** The COS FUV cross-dispersion profile is broad enough that a narrow box
clips real source flux. From h11 to h45 the continuum grows ~1.95x and the Ly-alpha peak ~2.23x,
both leveling off around h35-45. The default box height (41) sits near the plateau, which is why
it is near-optimal for a bright source like SN2010jl.

**Does h=45 work?** Yes, but it gives essentially zero gain over h=35. The actual sweep numbers:
h=11 gives 51% of the h=35 flux, h=17 gives 76%, h=25 already reaches the plateau (100%), and h=35
and h=45 are within 0.6% of each other. So **the plateau is fully reached by h=25-35** for this
source. h=45 is safe but gives nothing over h=35. The default h=41 sits firmly on the plateau.
- More flux from the wings of the profile: good.
- More background photons: relevant for faint sources where the background per pixel is comparable
  to the source. For SN2010jl (bright) this is negligible; for a faint SN at late times it might
  matter.
- No risk of clipping the profile with h=45 on this dataset: the FWHM of the FUV profile is well
  within h=35, so h=45 is safely on the plateau.

**Practical verdict.** h=41 (the default) is fine and near-optimal. If you want the most complete
absolute flux (e.g. for comparison with a model or for late-time faint sources where you want
every photon), h=45 is the safe ceiling for this mode. Just re-run the custom extraction with
`HEIGHT=45` in the edited XTRACTAB rows.

**Why Ly-alpha grows faster than the continuum.** Ly-alpha emission from the SN (if real) should
behave like the continuum -- both are point-source emission tracing the same PSF. If the Ly-alpha
peak grows x2.2 while the continuum grows x1.95 from the narrowest to widest box, that extra growth
likely reflects a real broader emission profile (Ly-alpha can be spatially extended in an interacting
SN) or geocoronal Ly-alpha airglow adding to it. Worth keeping in mind when interpreting the
absolute Ly-alpha flux.

---

## P2.8 Session log (phase 2)

### 7/1 -- time-series skeleton + naming explainer
- Reviewed phase-1 guide (through the paper-validation + Wynn-email sections) and the project_state
  todo log. Interpreted the messy handoff notes: the headline deliverable is the paper Figure 1
  time-series reproduction (early epochs from Wynn's files, main epochs from our own reduction).
- Wrote this file (`pipeline_phase2.md`) and moved the HST-naming/"epoch" explainer here from the
  phase-1 guide (it is phase-2 material; the phase-1 guide now points here).
- Built the time-series skeleton in `full_2023ixf.ipynb` (3 additive cells): loads Wynn's 5 pre-14d
  spectra, plots them log-log stacked with the paper's color scheme + atomic-line guides, saves
  `output/2023ixf_timeseries.png`. Looks close to Figure 1 for the top 5 epochs.
- Logged three open questions (flux column, ISM line removal, rest frame) for Wynn -- see P2.3.
- Next: identify + download the day 14/19/24/66 main epochs from MAST, reduce each into per-epoch
  subdirectories, coadd, and append to the `epochs` list so the same plot fills in the lower half.

### 7/2 -- day 14-66 reduced, full 9-epoch figure done
- **Main epochs identified**: the main program is **prop 17313** (of43 visits, 52x0.2 E1), not 17315.
  Days 14/19/24/66 = visits of4301/02/03/05; of4304 is the failed guide-star-acquisition visit
  (skipped, matches the paper's day-50 failure). Each visit holds 010 g230lb, 020 g430l, 030 g750l,
  040 ccdflat. Downloaded into per-epoch dirs `data/2023ixf/epochs/day{N}_of43xx/`.
- **Reduced them ourselves** with a new crj-aware reducer `scripts/reduce_epoch_ts.py` (the old
  `reduce_epoch.py` is left untouched so we can still backtrack). These epochs are CRSPLIT>1, so MAST
  gives a cr-combined `crj` per grating; we x1d off that (g230lb/g430l) and defringe the g750l with
  the contemporaneous 0.3x0.09 E1 ccdflat. That is exactly the paper's recipe (pipeline CR rejection +
  our own extraction). Defringe converged cleanly (shifts ~0.3-0.5 px, scales ~1.12, off the edges).
- **day-66 g230lb quirk**: its blue flux is so low the automatic trace-finder bails ("Cannot extract"),
  which the paper also reports. The reducer now retries with a fixed E1 trace center (a2center=893.5,
  maxsrch=0) -- a manual extraction. (It also needed one extra CRDS ref sync; some refs are date-specific.)
- **Full figure done**: `plot_ts` now shows all 9 epochs (`output/2023ixf_timeseries_full.png`). It
  reproduces the Figure 1 story -- hot blue UV + flash lines early, UV fading through days 14-24, and a
  red P Cygni spectrum (H-alpha, Ca II) by day 66 (whose UV is noisy, as expected for the faint blue end).
- **Style pass**: removed all person/PI name-drops and banner-style comments from the notebook code; a
  formatter had scrambled the time-series cells (duplicate plot cell + wrong order), now fixed.
- Then cleared the rest of the handoff notes in the same session (below).

### 7/2 (cont.) -- handoff items cleared
- **2010jl pipeline (STIS + FUV):** 2010jl has STIS/CCD G230LB + G430L (no G750L) plus COS FUV. Reduced the obk002 epoch (prop 12242) with `reduce_epoch_ts.py` (any grating subset), coadded, overlaid the COS FUV -> `output/2010jl_sed.png`. Confirms the reduction generalizes to a second SN and a 2-grating setup. pipeline.py still hardcodes 3 gratings, so it would need a small tweak for cases like this.
- **STIS coadd inter-grating scaling:** now in the main-epoch coadd. Unsaturated epochs give factors 0.99-1.08 (the gratings agree, unlike the saturated oezt01's 0.6); faint anchors are guarded. Answers open question in project_state.
- **G750L defringe:** verified it works (30-48% fringe-residual reduction on day 14). Then went deeper: also confirmed the residual is ~2x the photon noise in every red window, so it IS real fringe residual, not just photon noise. Tried red-targeted `opti_spreg`/`rms_region` and a finer grid -- negligible improvement. Conclusion: hit the fundamental ceiling of the single global shift+scale model; the fringes are wavelength-dependent so one pair of numbers can't remove them all. Normal; the paper accepts the same. See P2.4.
- **COS NUVB centering:** was a YCORR-vs-YFULL frame mix-up in the plot, not a real miscentering -- in YFULL the box is on the stripe.
- **COS FUV box-height / Ly-alpha:** the ~2x is the broad FUV cross-dispersion profile (narrow box clips flux); continuum 1.95x, Ly-alpha peak 2.23x from h11 to h45, plateauing by h35-45. h=45 is the safe ceiling; h=41 (default) is near-optimal. See P2.7.
- **Comment style:** removed all PI/person name-drops and banner comments from the notebooks (a formatter had also scrambled the time-series cells; fixed). Deleted stray empty `.venv*` dirs.

### 7/2 (cont.) -- deeper round from full handoff prompt re-sent
- **COS NUV time phase:** lf8803 (prop 17497) is day ~214 after explosion -- ~210 days after the STIS it was being overlaid with. NOT contemporaneous. Moved to its own labeled cell in `full_2023ixf.ipynb`; the day-3.5 SED cell is now pure STIS. At day 214 the COS NUV shows broad Mg II 2800 emission (late-time cool-dense-shell signature). See P2.2.
- **Defringe pushed further:** confirmed the residual is real (2x photon noise), tried region-targeted and finer-grid mkfringeflat -- hit the single global shift+scale ceiling. Full explanation now in P2.4.
- **COS NUV coadd noise vs HASP:** our median coadd has slightly less noise than HASP in some windows (positive finding from having the raw data). The 2150-2350 A artifact is a NUVA stripe contamination; HASP removes it with a per-bin flux-deviation filter; we don't yet. Full discussion P2.6.
- **COS FUV h=45:** would capture ~2-4% more tail flux vs h=41; safe to use; default h=41 is near-optimal. Full discussion P2.7.
- **Updated pipeline_phase2.md:** added a meeting-ready quick-review summary at the top of this file (before P2.0), deepened P2.4 with the defringe ceiling finding, added P2.6 (COS NUV coadd + NUVA artifact), added P2.7 (FUV box ceiling).

---

## P2.9 Defringe investigation -- consolidated recap (7/1-7/2)

A single place that pulls together the whole G750L defringe investigation from those two days, so you
do not have to reassemble it from P2.4 + the session log. The step-by-step mechanics live in P2.4;
this is the "what did we do, what did we learn, how far did we push, can we go further, and where do
I look" version.

**What we did.**
- Confirmed our defringe chain matches the paper exactly: contemporaneous narrow-slit CCDFLAT ->
  `normspflat` -> `mkfringeflat` -> `defringe` -> `x1d`. Ran it on the early oezt01 epoch (day 3.5)
  and on all four of43 main epochs (days 14/19/24/66) via `scripts/reduce_epoch_ts.py`.
- Measured the fringe residual quantitatively (not just "looks smoother"): RMS of
  spectrum / median-smoothed-spectrum in four red windows, before vs after.
- Then tried to beat it: red-targeted `opti_spreg`/`rms_region`, a finer shift/scale grid, and a
  widened shift search so the solver would not clip at an edge.

**What we learned.**
- The defringe works: **30-48% reduction** in red fringe-residual RMS on day 14 (8000-8500 A
  0.058->0.040, 8500-9000 0.056->0.039, 9000-9700 0.125->0.065, 9700-10200 0.188->0.140).
- The leftover residual is **~2x the photon noise** in every red window, so it is real fringe
  residual, not just noise -- there is genuinely something left to remove, but...
- ...**we are at the ceiling of the single global shift+scale model.** `mkfringeflat` registers the
  flat to the science with one shift and one scale for the whole spectrum; the real fringes are
  wavelength-dependent (they drift in phase across the red), so no single pair of numbers can null
  them everywhere. The reddest wiggles past ~9700 A are the hardest and survive.
- **Watch `mkfringeflat` convergence.** Interior shift/scale = trustworthy; a value pinned at the
  edge of the search range = the solver never found a minimum. The saturated oezt01 epoch hit the
  scale ceiling at 1.2 (a red flag, from saturation); the of43 epochs converge interior
  (day 14: shift +0.30 px, scale ~1.12). This is the single most useful diagnostic.

**Can we push further?** Not meaningfully with the current tools. Region-targeted RMS and finer grids
gave negligible improvement. Beating the residual would require a **wavelength-dependent** fringe
registration (piecewise shift/scale across the red, or a per-column fit) that `stistools.defringe`
does not expose -- and the paper does not do it either, it accepts the same residual. So this is the
practical stopping point; further effort here is low-value until/unless a red-heavy science case
demands sub-percent fringe removal past 9700 A.

**Where to see it.**
- Plot: `output/defringe_check_day14.png` -- defringed vs undefringed on the day-14 G750L, red band.
- Notebook `full_2023ixf.ipynb`: the day-3.5 oezt01 SED cell (loads the defringed G750L into the
  coadd) and the fringe-flat aperture-check cell (confirms 52X0.1 flat vs 52X0.2 science).
- Scripts: `scripts/reduce_epoch_ts.py` (runs the full chain per epoch) and `scripts/reduce_defringe.py`.
- Full step-by-step mechanics and the paper cross-check: **P2.4** above.
