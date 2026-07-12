# LaTeX Setup Explainer

This doc covers the actual software setup we have going, what all those extra files
in `interim_reports/` are, and other options you might want to know about. Does not
cover LaTeX syntax/macros since you already know that.

---

## What we have right now

### The editor: VS Code + LaTeX Workshop extension

You have the **LaTeX Workshop** extension installed. It adds two main things to VS Code:

1. **"Build LaTeX Project"** button (or `Ctrl+Alt+B`): runs `pdflatex` (or whichever
   recipe you configure) on the open `.tex` file and produces a `.pdf`.

2. **"View LaTeX PDF File"** button (or `Ctrl+Alt+V`): opens the compiled PDF in a
   split panel to the right of the editor. It stays synced -- when the PDF is
   recompiled, the viewer refreshes automatically.

3. **Auto-build on save**: by default, LaTeX Workshop recompiles every time you save
   the file (`Ctrl+S`). The PDF viewer updates a few seconds later once compilation
   finishes. This is the closest thing to Overleaf's auto-compile -- not live-typing,
   but save-triggers-compile. Compilation takes a couple of seconds for a short report
   like this one.

The extension also gives you autocomplete for `\commands`, go-to-definition for
`\ref{}` and `\cite{}` labels, and syntax highlighting.

### The compiler: pdflatex (via TeX Live in WSL)

LaTeX Workshop needs a TeX installation on the system. On your machine we're using
**TeX Live** inside WSL Debian (`/usr/share/texlive`), which is a full TeX distribution
(thousands of packages, fonts, etc). When you hit Build or save the `.tex` file,
LaTeX Workshop shells out to `pdflatex` inside WSL.

`pdflatex` compiles `.tex` -> `.pdf` in one pass, and usually needs to be run twice
to resolve cross-references (figure numbers, citations, page refs). LaTeX Workshop
handles this automatically -- it reruns until the references stabilize.

---

## The junk files in `interim_reports/`

When pdflatex compiles `interim_report1.tex`, it writes a bunch of auxiliary files
alongside the `.pdf`. Here is what each one is:

| File | What it is |
|------|-----------|
| `interim_report1.pdf` | The compiled PDF. This is the only output you actually care about. |
| `interim_report1.aux` | Auxiliary file. Stores cross-reference data (figure/table/section numbers, citation keys). Written on pass 1, read on pass 2 so `\ref{}` and `\cite{}` resolve correctly. Safe to delete; just gets regenerated. |
| `interim_report1.log` | Full compiler log. Every warning, error, font load, and package message from the last compile. If something looks wrong in the PDF, grep this file. Usually verbose (~several hundred lines). |
| `interim_report1.out` | Written by the `hyperref` package. Stores the PDF bookmarks/outlines (the clickable section list in PDF readers). |
| `interim_report1.fls` | File-list. Records every file read and written during compilation (useful for incremental build tools). LaTeX Workshop uses this. |
| `interim_report1.fdb_latexmk` | `latexmk` dependency database. `latexmk` is a smarter build driver than bare `pdflatex` -- it tracks which files changed and only reruns what is necessary. LaTeX Workshop uses `latexmk` under the hood by default. |

All of these (except the `.pdf`) can be deleted at any time without losing anything;
they all get regenerated on the next compile. You can add them to `.gitignore` if
you don't want to track them in git.

If you want to clean them manually:
```bash
cd interim_reports
rm -f *.aux *.log *.out *.fls *.fdb_latexmk
```

---

## How it compares to Overleaf

| | Our setup (VS Code + LaTeX Workshop + TeX Live in WSL) | Overleaf |
|---|---|---|
| **Where it runs** | Entirely local. No internet needed to compile. | Cloud. Requires internet. |
| **Compile trigger** | On save. | On save (or manually). Same speed. |
| **PDF viewer** | Side panel in VS Code. | Side panel in browser. |
| **Collaboration** | Edit the `.tex` file in git like any other code. | Built-in real-time co-editing. |
| **Package availability** | TeX Live has essentially everything. | Overleaf has essentially everything. Same. |
| **Version control** | Full git history since the `.tex` file lives in the repo. | Overleaf has its own internal history, not git (unless you use the git bridge on paid plans). |

In practice, our setup is equivalent to Overleaf for a solo author. The main
thing Overleaf adds is real-time co-editing, which is not relevant here since
you're the only author.

---

## Upgrading / other options

### latexmk vs pdflatex directly

LaTeX Workshop defaults to `latexmk`, which is smarter than bare `pdflatex`:
- Only reruns passes when needed (tracks file changes).
- Automatically runs BibTeX/Biber for bibliographies if you use `.bib` files.
- Handles the multi-pass reference resolution for you.

You probably already have this running via LaTeX Workshop without noticing.

### XeLaTeX / LuaLaTeX

Alternative compilers (instead of `pdflatex`) that support Unicode natively and
can use any system font. For this report, `pdflatex` is fine -- you're not doing
anything that needs Unicode or custom fonts.

### Using a `.bib` file (BibTeX)

Right now the bibliography is inline (`\begin{thebibliography}{99}...`), which is
fine for a short report. If the reference list grows or you want to reuse references
across documents (e.g., for the final paper), switch to a `.bib` file + BibTeX:

1. Create `refs.bib` in `interim_reports/`.
2. Replace the `\begin{thebibliography}` block with:
   ```latex
   \bibliographystyle{aasjournal}
   \bibliography{refs}
   ```
3. Add the `natbib` package (already included) and run the compile sequence:
   `pdflatex` -> `bibtex` -> `pdflatex` -> `pdflatex`.
   `latexmk` handles this sequence automatically.

### Overleaf sync

If you want to use Overleaf (e.g., to share with your mentor for comments), you
can upload the `.tex`, `.png` figures, and any `.bib` file directly to Overleaf.
The file structure is flat on Overleaf so adjust image paths from `../output/` to
just `stis_2d_extraction.png` etc. if you do this.

Alternatively, Overleaf's GitHub sync (paid feature) can stay in sync with the
repo automatically.

---

## Quick reference: what to do when things go wrong

**PDF is blank or has "??" for figure/ref numbers:**
Run build twice (LaTeX Workshop usually does this automatically with `latexmk`).

**"File not found" for an image:**
The path in `\includegraphics{}` is relative to the `.tex` file's location.
Our figures live in `../output/` and `../` (proposal_figure1.png) since the `.tex`
is in `interim_reports/`. Check that the file actually exists there.

**Weird spacing or layout changed:**
Delete the `.aux` and `.out` files and recompile clean.

**Build button does nothing / LaTeX Workshop not responding:**
Open the LaTeX Workshop output panel (View > Output, select "LaTeX Workshop") to
see what command it's trying to run and what error it hit.
