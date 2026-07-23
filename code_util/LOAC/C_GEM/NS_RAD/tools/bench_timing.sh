#!/usr/bin/env bash
# Measure wall-clock runtime per river (full 2-year production config) for the schematic
# performance panel. Runs all four in parallel to a scratch dir (does not touch runs/),
# capturing each site's own wall time plus the overall parallel total. Serial total is
# the sum of the per-site times. Writes tools/../docs/performance_metrics.json.
set -euo pipefail
export PATH="$HOME/miniforge3/bin:$PATH"
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
CODE="$ROOT/code"
BENCH="/tmp/cgem_bench"
rm -rf "$BENCH"; mkdir -p "$BENCH"
: > "$BENCH/timings.txt"

overall_start=$(date +%s)
for s in colville kuparuk sagavanirktok canning; do
    ( mkdir -p "$BENCH/$s" && cd "$BENCH/$s"
      st=$(date +%s)
      CGEM_SITE="$s" CGEM_TS=48 CGEM_OUTPUT=nc PYTHONPATH="$CODE" \
          python3 "$CODE/main.py" > /dev/null 2>&1
      en=$(date +%s)
      echo "$s $((en - st))" >> "$BENCH/timings.txt" ) &
done
wait
overall_end=$(date +%s)
parallel_total=$((overall_end - overall_start))

# Emit a metrics JSON the schematic reads (per-site seconds, parallel + serial totals).
python3 - "$BENCH/timings.txt" "$parallel_total" "$ROOT/docs/performance_metrics.json" <<'PY'
import json, sys
times = {}
for line in open(sys.argv[1]):
    parts = line.split()
    if len(parts) == 2 and parts[0] != "PARALLEL_TOTAL":
        times[parts[0]] = int(parts[1])
parallel_total = int(sys.argv[2])
out = {
    "per_site_sec": times,
    "serial_total_sec": sum(times.values()),
    "parallel_total_sec": parallel_total,
    "config": "2-year run, DELTI=75 s, TS=48 (hourly), NetCDF, numba, 14 cores",
}
json.dump(out, open(sys.argv[3], "w"), indent=2)
print(json.dumps(out, indent=2))
PY
echo "BENCH DONE"
