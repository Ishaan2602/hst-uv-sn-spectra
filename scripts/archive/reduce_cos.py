import os, sys, shutil
import calcos

# runs in wsl (surf_uv env). lref must be exported before calling (points at the cos refs).
# args: input (corrtag or asn), outdir. corrtag -> x1d, asn -> x1dsum.
inp = sys.argv[1]
outdir = sys.argv[2]

# calcos wont overwrite, so wipe the outdir first
if os.path.exists(outdir):
    shutil.rmtree(outdir)
os.makedirs(outdir, exist_ok=True)

calcos.calcos(inp, verbosity=0, outdir=outdir)
print("calcos done ->", os.path.basename(inp))
