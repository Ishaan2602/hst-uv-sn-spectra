# Interim Reports Notes

## Program Details
- SURF 10-week program
- Interim Report 1 due: ~July 3, 2026 (end of week 3)
- Interim Report 2 / Abstract due: ~July 31, 2026 (end of week 7)
- Submitted through SFP Online, must be mentor-approved (not co-mentor)

## First Interim Report Guidelines
1. Motivation, background, overview of ongoing group work. Include references.
2. Problem being worked on, how it fits into group work. Approach and methods.
3. Work done. What progress made?
4. Challenges so far + anticipated challenges going forward.

## Second Interim Report Guidelines
1. Work completed past month (experiments, data analysis). Exact technical specs,
   chronological order if possible.
2. Progress and observations. Are they in line with expectations?
3. Problems encountered: source and how addressing them.
4. Research goals for remainder. Have goals changed?

---

## Report 1 Content Notes (drafted July 6, 2026)

### Style
- Maintain roughly same language/style as the proposal
- No em-dashes, no AI-style wording
- Simple and direct, not boasty
- This is a science document by an undergrad 3 weeks in

### Figures used in Report 1
- Fig 1: `output/2023ixf_timeseries_full.png` -- 9-epoch NUV-NIR time series,
  our reproduction of Bostroem et al. 2024 Fig 1. Shows the main science product
  of the first 3 weeks.
  
### Key points not in the proposal that need to be added
- Both STIS and COS are used (proposal focused on STIS only for the pipeline)
- WSL requirement for HSTCAL binaries (no native Windows build)
- COS has a fundamentally different reduction path (photon-counting vs CCD)
- Discovery query result: 140 unique UV SNe from 1591 HST spectra
- Day-66 G230LB: faint blue end caused pipeline extraction failure, needed manual trace centering
- Early epoch saturation issue (GO-17205 data saturated 3178-5022 A in G430L)

### Proposal sections to copy largely verbatim (with fixes)
- Introduction/Background: mostly verbatim, fix em-dashes, add COS mention
- Objectives: mostly verbatim, fix objective 1 to say "stistools and calcos"
- Approach: update to mention both STIS and COS reduction paths

### Proposal language to NOT copy
- Work plan (weeks 1-10) -- replace with actual progress section
- References: keep the same ones; add Bostroem et al. 2024

---

## Report 2 Notes (to be filled in later, due ~July 31)

- Describe work completed weeks 4-7 in detail
- Focus on batch processing progress, CMFGEN comparison if started
- Exact technical specs, chronological order

