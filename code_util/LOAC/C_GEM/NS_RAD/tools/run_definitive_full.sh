#!/usr/bin/env bash
# Run the FULL-FORCING interannual site variants at DEFINITIVE (6-min) cadence for
# just the first 2 years of their real 2005-2023 record (2005-2006), instead of the
# usual 19-year daily-cadence run tools/run_interannual_full.sh does. This exists
# SOLELY to give Fig 3 of tools/make_publication_figures_full.py a tidal-resolved
# (6-min) reference to derive its daily-sampling aliasing correction from -- the
# existing runs/definitive/ serves that role for the ORIGINAL 1980-2023 publication
# figures, but it is a repeated-2022-climatology run, not drawn from the full-forcing
# period's real multi-year forcing at all, so it is not reused here.
#
# No new forcing files needed: the _interannual_full site variants' forcing files
# already span 2005-2023 at daily (met) and hourly (wind/air-temp/humidity/solar)
# resolution; this just reads the FIRST 730 days of them (real 2005-2006) at the
# model's native 6-min integration/output cadence instead of continuing to day 6935
# at daily output cadence. Same site variants, same forcing files, different
# CGEM_MAXT_DAYS/CGEM_TS -- not a new data pipeline.
#
# Usage:
#   tools/run_definitive_full.sh                    # colville+kuparuk+sagavanirktok
#   tools/run_definitive_full.sh kuparuk             # just one
#   SERIAL=1 tools/run_definitive_full.sh            # one at a time
#
# Results land in runs/definitive_full/<site>/ ; stdout in
# runs/definitive_full/<site>/run.log

set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
CODE="$ROOT/code"
SITES=("${@:-}")
if [ -z "${SITES[*]}" ]; then
    SITES=(colville kuparuk sagavanirktok)
fi

MAXT_DAYS="${CGEM_MAXT_DAYS:-730}"
WARMUP_DAYS="${CGEM_WARMUP_DAYS:-365}"
TS="${CGEM_TS:-12}"

pids=()
for site in "${SITES[@]}"; do
    outdir="$ROOT/runs/definitive_full/$site"
    mkdir -p "$outdir"
    echo "-> ${site}_interannual_full  (output: runs/definitive_full/$site, $MAXT_DAYS days @ 6-min, i.e. real 2005-2006)"

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
    [ $fail -eq 0 ] && echo "all runs finished" || { echo "a run FAILED -- check runs/definitive_full/*/run.log"; exit 1; }
fi
