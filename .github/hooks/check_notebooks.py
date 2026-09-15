"""
scan all active notebooks for markdown cell anti-patterns (AI giveaways, banned words, name references).
run from WORK root: python .github/hooks/check_notebooks.py
"""
import json, os, re, sys

NOTEBOOKS_DIR = 'notebooks'
EXCLUDE = {'archive', 'reference'}

BANNED_WORDS = re.compile(r'\b(robust|crucial|leverage|seamless|delve|furthermore|snippet)\b', re.I)
AI_CONSTRUCTS = re.compile(r'\*\*(Goal|Sections|Overview|Background|Motivation):', re.I)
NAME_REFS = re.compile(r'\bWynn\b|\bthe PI\b')
INTERNAL_REFS = re.compile(r'analysis_phase\d+\.md|pipeline_phase\d+\.md|phase[- ]\d+ (item|point|task)')
DENSE_PARA = re.compile(r'^[A-Z][a-z].{100,}$', re.M)

issues = []

for fname in os.listdir(NOTEBOOKS_DIR):
    if not fname.endswith('.ipynb'): continue
    if any(ex in fname for ex in EXCLUDE): continue
    path = os.path.join(NOTEBOOKS_DIR, fname)
    nb = json.load(open(path, encoding='utf-8'))
    for i, cell in enumerate(nb['cells']):
        if cell['cell_type'] != 'markdown': continue
        src = ''.join(cell['source'])
        for pat, label in [
            (BANNED_WORDS, 'banned word'),
            (AI_CONSTRUCTS, 'AI construct'),
            (NAME_REFS, 'name reference'),
            (INTERNAL_REFS, 'internal doc reference'),
        ]:
            m = pat.search(src)
            if m:
                issues.append(f'{fname} cell {i+1} [{label}]: ...{src[max(0,m.start()-15):m.end()+30]}...')

if issues:
    print(f'ISSUES FOUND ({len(issues)}):')
    for s in issues: print(' ', s)
    sys.exit(1)
else:
    print(f'clean: no issues across {len(os.listdir(NOTEBOOKS_DIR))} notebooks checked')
