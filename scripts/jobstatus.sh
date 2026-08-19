#!/bin/bash
# quick status of the running ism catalog job (writes per-epoch cog.csv into output5/<SN>/absorption/)
cd "$(dirname "$0")/.." || exit
echo "=== process ==="
ps -eo pid,etime,cmd | grep -E "ism\.py|run_full_catalog" | grep -v grep || echo "  no ism/pipeline process running"
echo
echo "=== progress (output5 absorption cog.csv) ==="
tot=$(find output5 -path "*/absorption/*_cog.csv" 2>/dev/null | wc -l)
fresh=$(find output5 -path "*/absorption/*_cog.csv" -newermt "-15 min" 2>/dev/null | wc -l)
echo "  total cog.csv: $tot    written in last 15min: $fresh"
echo "  newest 3:"
find output5 -path "*/absorption/*_cog.csv" -printf '%T+  %p\n' 2>/dev/null | sort | tail -3 | sed 's/^/    /'
echo
if [ -f output/ism_cog_summary.csv ]; then
  echo "=== summary exists (job wrote it at the end) : $(wc -l < output/ism_cog_summary.csv) rows ==="
else
  echo "=== summary output/ism_cog_summary.csv NOT yet written (job still running) ==="
fi
