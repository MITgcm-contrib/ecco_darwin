#!/usr/bin/env bash
# Run C-GEM's INTERANNUAL site variants (real 1980-2023 PWBM discharge + riverine TOC,
# everything else unchanged from the regular site) -- see sites/colville_interannual.py
# and CLAUDE.md -> "Interannual forcing". Mirrors run_sites.sh exactly; see that script's
# header for why each site needs its own directory and why runs are parallel by default.
#
# Usage:
#   tools/run_interannual.sh                    # colville+kuparuk+sagavanirktok, in parallel
#                                                # (canning excluded -- inherits EL=0, can't run)
#   tools/run_interannual.sh kuparuk             # just one
#   SERIAL=1 tools/run_interannual.sh            # one at a time
#   CGEM_MAXT_DAYS=30 tools/run_interannual.sh kuparuk    # quick smoke test instead of the
#                                                          # full 44-year record
#
# Results land in runs/interannual/<site>/ ; stdout in runs/interannual/<site>/run.log
# Requires the forcing files from tools/build_interannual_forcings.py to already exist.
#
# SAVE CADENCE IS DAILY BY DEFAULT (CGEM_TS=2880 -> TS*DELTI=86400s at these sites'
# DELTI=30s), not the definitive runs' default TS=12 (every 6 min). A 44-year run at the
# 6-min cadence would be ~22x runs/definitive's ~1-2 GB per site -- tens of GB each, and
# this machine only has ~60 GB free. Daily resolution is still plenty for annual budgets
# and seasonal shape (what tools/make_interannual_pdf.py reads); it is NOT enough for the
# sub-daily Hovmoller/profile detail runs/definitive's diagnostics PDF shows. Override
# with CGEM_TS if you need finer output and have the disk for it.

set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
CODE="$ROOT/code"
SITES=("${@:-}")
if [ -z "${SITES[*]}" ]; then
    SITES=(colville kuparuk sagavanirktok)
fi

# Full record by default (44 years); override with CGEM_MAXT_DAYS for a shorter run.
MAXT_DAYS="${CGEM_MAXT_DAYS:-16060}"
WARMUP_DAYS="${CGEM_WARMUP_DAYS:-365}"
TS="${CGEM_TS:-2880}"

pids=()
for site in "${SITES[@]}"; do
    outdir="$ROOT/runs/interannual/$site"
    mkdir -p "$outdir"
    echo "-> ${site}_interannual  (output: runs/interannual/$site, $MAXT_DAYS days)"

    if [ -n "${SERIAL:-}" ]; then
        ( cd "$outdir" && CGEM_SITE="${site}_interannual" CGEM_MAXT_DAYS="$MAXT_DAYS" \
              CGEM_WARMUP_DAYS="$WARMUP_DAYS" CGEM_TS="$TS" PYTHONPATH="$CODE" \
              PYTHONWARNINGS="default" python3 "$CODE/main.py" 2>&1 | tee run.log )
    else
        ( cd "$outdir" && CGEM_SITE="${site}_interannual" CGEM_MAXT_DAYS="$MAXT_DAYS" \
              CGEM_WARMUP_DAYS="$WARMUP_DAYS" CGEM_TS="$TS" PYTHONPATH="$CODE" \
              PYTHONWARNINGS="default" python3 "$CODE/main.py" > run.log 2>&1 ) &
        pids+=($!)
    fi
done

if [ ${#pids[@]} -gt 0 ]; then
    echo "launched ${#pids[@]} run(s); waiting..."
    fail=0
    for pid in "${pids[@]}"; do wait "$pid" || fail=1; done
    [ $fail -eq 0 ] && echo "all runs finished" || { echo "a run FAILED -- check runs/interannual/*/run.log"; exit 1; }
fi
