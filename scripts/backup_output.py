#!/usr/bin/env python3
# snapshot the current output tree into output_backups/ with a dated, phase-tagged name.
# run anytime: python scripts/backup_output.py   (optionally pass a source dir + a note)
# name format: output_<mon><dd>-<yy>_<function>_phase<N>   e.g. output_sep14-26_analysis_phase5

import os
import sys
import subprocess
from datetime import date

import paths

def backup_name(when=None):
    when = when or date.today()
    fn, ph = paths.project_phase()
    mon = when.strftime("%b").lower()
    return f"output_{mon}{when.strftime('%d')}-{when.strftime('%y')}_{fn}_phase{ph}"

def main():
    src = sys.argv[1] if len(sys.argv) > 1 else paths.OUT
    if not os.path.isdir(src):
        sys.exit(f"source tree not found: {src}")
    os.makedirs(paths.BACKUPS, exist_ok=True)
    name = backup_name()
    dst = os.path.join(paths.BACKUPS, name)
    if os.path.exists(dst):
        sys.exit(f"backup already exists (bump phase or delete it first): {dst}")
    print(f"backing up {src} -> {dst}")
    # robocopy is the fast native windows copy; /e = all subdirs incl empty, /nfl /ndl = quiet
    r = subprocess.run(["robocopy", src, dst, "/e", "/nfl", "/ndl", "/njh", "/njs"])
    # robocopy returns <8 on success (bit flags), >=8 is a real error
    if r.returncode >= 8:
        sys.exit(f"robocopy failed (code {r.returncode})")
    print(f"done: {name}")

if __name__ == "__main__":
    main()
