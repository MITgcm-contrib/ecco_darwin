#!/usr/bin/env bash
# Run the FORTE campaign year site variants (sites/<river>_forte_campaign.py --
# CLIMATOLOGICAL forcing, tools/build_forte_campaign_forcing.py) at DEFINITIVE
# (6-min) cadence, long enough to cover the real FORTE monitored window: day-of-year
# 152 (Phase 1a, Jun 1) through 273 (mooring recovery, late September) -- see
# tools/build_osse_observations.py's PHASES/MOORING_WINDOW.
#
# MAXT_DAYS=638 = 365 (warmup, one climatological year) + 273 (day-of-year the
# campaign season ends, in a second climatological year). This is a SEPARATE
# framework from runs/osse_nature (2007/2008/2009 analog years, kept as-is for its
# own purposes) and runs/interannual_full (the real 19-year daily record, also kept
# as-is) -- see sites/<river>_forte_campaign.py's docstring for the full reasoning.
#
# Usage:
#   tools/run_forte_campaign.sh                    # colville+kuparuk+sagavanirktok
#   tools/run_forte_campaign.sh kuparuk             # just one
#   SERIAL=1 tools/run_forte_campaign.sh            # one at a time
#
# Results land in runs/forte_campaign/<site>/ ; stdout in
# runs/forte_campaign/<site>/run.log

set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
CODE="$ROOT/code"
SITES=("${@:-}")
if [ -z "${SITES[*]}" ]; then
    SITES=(colville kuparuk sagavanirktok)
fi

MAXT_DAYS="${CGEM_MAXT_DAYS:-638}"
WARMUP_DAYS="${CGEM_WARMUP_DAYS:-365}"
TS="${CGEM_TS:-12}"

pids=()
for site in "${SITES[@]}"; do
    outdir="$ROOT/runs/forte_campaign/$site"
    mkdir -p "$outdir"
    echo "-> ${site}_forte_campaign  (output: runs/forte_campaign/$site, $MAXT_DAYS days @ 6-min, climatological forcing)"

    if [ -n "${SERIAL:-}" ]; then
        ( cd "$outdir" && CGEM_SITE="${site}_forte_campaign" CGEM_MAXT_DAYS="$MAXT_DAYS" \
              CGEM_WARMUP_DAYS="$WARMUP_DAYS" CGEM_TS="$TS" PYTHONPATH="$CODE" \
              PYTHONWARNINGS="default" python3 "$CODE/main.py" 2>&1 | tee run.log )
    else
        ( cd "$outdir" && CGEM_SITE="${site}_forte_campaign" CGEM_MAXT_DAYS="$MAXT_DAYS" \
              CGEM_WARMUP_DAYS="$WARMUP_DAYS" CGEM_TS="$TS" PYTHONPATH="$CODE" \
              PYTHONWARNINGS="default" python3 "$CODE/main.py" > run.log 2>&1 ) &
        pids+=($!)
    fi
done

if [ ${#pids[@]} -gt 0 ]; then
    echo "launched ${#pids[@]} run(s); waiting..."
    fail=0
    for pid in "${pids[@]}"; do wait "$pid" || fail=1; done
    [ $fail -eq 0 ] && echo "all runs finished" || { echo "a run FAILED -- check runs/forte_campaign/*/run.log"; exit 1; }
fi
