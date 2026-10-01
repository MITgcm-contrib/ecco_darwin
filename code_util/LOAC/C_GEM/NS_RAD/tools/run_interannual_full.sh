#!/usr/bin/env bash
# Run C-GEM's FULL-FORCING INTERANNUAL site variants -- EVERY forcing category is
# genuinely multi-year (2005-2023, the shortest common window across every real
# source), not just discharge/DOC like tools/run_interannual.sh's variants. See
# sites/colville_interannual_full.py and CLAUDE.md -> "Interannual forcing" ->
# "Full-forcing interannual (2005-2023)". Mirrors run_interannual.sh exactly except
# for the record length (19 years here, not 44) and the site-variant suffix.
#
# Usage:
#   tools/run_interannual_full.sh                    # colville+kuparuk+sagavanirktok
#   tools/run_interannual_full.sh kuparuk             # just one
#   SERIAL=1 tools/run_interannual_full.sh            # one at a time
#   CGEM_MAXT_DAYS=30 tools/run_interannual_full.sh kuparuk    # quick smoke test
#
# Results land in runs/interannual_full/<site>/ ; stdout in
# runs/interannual_full/<site>/run.log. Requires the forcing files from
# tools/build_interannual_met.py, build_interannual_humidity.py,
# build_interannual_solar.py, build_interannual_surge.py, build_interannual_marine.py,
# and build_interannual_river_temp.py to already exist (the discharge/DOC slice is
# produced inline the first time this repo's forcing/ was built -- see CLAUDE.md).
#
# SAVE CADENCE IS DAILY (CGEM_TS=2880), same reasoning as run_interannual.sh: a 19-year
# 6-min-cadence run would still be multiple GB per site for no benefit this script's
# consumers (annual budgets, seasonal shape) need.

set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
CODE="$ROOT/code"
SITES=("${@:-}")
if [ -z "${SITES[*]}" ]; then
    SITES=(colville kuparuk sagavanirktok)
fi

# Full 2005-2023 record by default (19 years = 6935 days); override with CGEM_MAXT_DAYS.
MAXT_DAYS="${CGEM_MAXT_DAYS:-6935}"
WARMUP_DAYS="${CGEM_WARMUP_DAYS:-365}"
TS="${CGEM_TS:-2880}"

pids=()
for site in "${SITES[@]}"; do
    outdir="$ROOT/runs/interannual_full/$site"
    mkdir -p "$outdir"
    echo "-> ${site}_interannual_full  (output: runs/interannual_full/$site, $MAXT_DAYS days)"

    if [ -n "${SERIAL:-}" ]; then
        ( cd "$outdir" && CGEM_SITE="${site}_interannual_full" CGEM_MAXT_DAYS="$MAXT_DAYS" \
              CGEM_WARMUP_DAYS="$WARMUP_DAYS" CGEM_TS="$TS" PYTHONPATH="$CODE" \
              PYTHONWARNINGS="default" python3 "$CODE/main.py" 2>&1 | tee run.log )
    else
        ( cd "$outdir" && CGEM_SITE="${site}_interannual_full" CGEM_MAXT_DAYS="$MAXT_DAYS" \
              CGEM_WARMUP_DAYS="$WARMUP_DAYS" CGEM_TS="$TS" PYTHONPATH="$CODE" \
              PYTHONWARNINGS="default" python3 "$CODE/main.py" > run.log 2>&1 ) &
        pids+=($!)
    fi
done

if [ ${#pids[@]} -gt 0 ]; then
    echo "launched ${#pids[@]} run(s); waiting..."
    fail=0
    for pid in "${pids[@]}"; do wait "$pid" || fail=1; done
    [ $fail -eq 0 ] && echo "all runs finished" || { echo "a run FAILED -- check runs/interannual_full/*/run.log"; exit 1; }
fi
