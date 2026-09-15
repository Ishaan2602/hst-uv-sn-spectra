---
applyTo: 'scripts/**/*.py'
---

# Script style rules

- Short lowercase abbreviated vars matching reference notebooks (`wvl`, `flx`, `err`, `spec`, `flist`).
- Comments explain WHY, not WHAT. One short line. No docstrings for trivial functions.
- No name-drops in comments (no "Wynn", "PI"). Refer to data neutrally.
- BANNED: `robust`, `crucial`, `leverage`, `seamless`, `delve`, `dive in`, `furthermore`, `snippet`.
- No emojis, no em-dashes.
- No error handling for cases that cannot happen. Let errors bubble up.
- No redundant null checks when types already guarantee safety.
- NEVER create virtual environments. Global Python 3.14 on Windows.
- No internal doc references in product output strings (no "phase-5", no "see analysis_phaseN.md").
- Public-facing product strings (JSON fields, CSV columns, README) must read as scientific text.
