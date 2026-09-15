---
applyTo: '**/*.ipynb'
---

# Rules for every notebook cell

## Markdown cells
- **Title** (`# heading`): 1-2 line description max. Bullet list for multiple sections. No "This notebook does X, Y, Z" prose.
- **Section headers** (`##`/`###`): 2-4 terse bullets or 1-2 lines. Never a dense paragraph.
- **Result cells**: bullet list of findings. Not a run-on bold paragraph (`**Read-out.** The clip drops...`).
- **Banned constructs**: `**Goal:**`, `**Sections:**`, `**Overview:**`, numbered "Sections: 1. 2. 3." in title cells.
- **No internal doc references**: no "see analysis_phase5.md", no phase labels as citations.
- **No isolated prose paragraphs** starting mid-cell ("Across 33 SNe..." -> make it a bullet).

## Code cells
- Short lowercase abbreviated vars: `wvl`, `flx`, `err`, `spec`, `tbdata`, `flist`.
- Comments explain WHY, not WHAT. One short line max.
- No name-drops in comments. No emojis. No em-dashes.
- BANNED: `robust`, `crucial`, `leverage`, `seamless`, `delve`, `dive in`, `furthermore`, `snippet`.

## Figure rules
- `rcParams` block required per plot section (serif font, inward ticks, thick spines).
- Always `fig, ax = plt.subplots()`. Never `plt.plot()` directly.
- LaTeX axis labels. No sentence-style `ax.set_title()`. Short terse ID-style ok.
- Save via `save(fig, name)` helper to `PLOT_DIR`. Always display inline.
- Image budget: <= 15 plot cells per session turn.
