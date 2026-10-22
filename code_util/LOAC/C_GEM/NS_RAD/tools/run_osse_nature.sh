#!/usr/bin/env bash
# Build the "nature run" (synthetic truth) for the FORTE campaign OSSE
# (tools/build_osse_observations.py): the *_interannual_full site variants at
# DEFINITIVE (6-min) cadence, run long enough to cover three analog flow-regime years
# -- 2007 (low), 2008 (average), 2009 (high), chosen because they are the cleanest
# consecutive low/average/high triplet in the full-forcing record for all three
# rivers simultaneously (checked directly against the discharge forcing: Colville
# 198/330/438 m3/s, Kuparuk 4.7/19.0/66.6, Sagavanirktok 31/77/172 mean annual
# discharge) -- and because they are EARLY in the 2005-2023 record, keeping this run
# affordable (NS-RAD has no restart/checkpoint capability, so reaching any later
# year means simulating every day before it at full 6-min cadence too).
#
# The real FORTE campaign is 2027-2028 (see the campaign plan this script's sibling
# tool reads); NS-RAD's forcing does not extend that far, so 2007/2008/2009 serve as
# flow-regime ANALOGS -- the campaign's real day-of-year phase windows are mapped
# onto each of these three years, not literally onto 2027/2028. See
# tools/build_osse_observations.py's docstring for the full reasoning and the
# satellite-availability resolution (SWOT/PACE are treated as a notional "if
# coordinated" channel, not gated by their real launch dates, since the whole
# analog-year substitution already breaks literal calendar correspondence).
#
# MAXT_DAYS=1825 = 365 (2005 warmup) + 365 (2006, unavoidable pass-through -- the
# model cannot skip ahead) + 365*3 (2007-2009, the actual study years). Saved at
# 6-min cadence for the WHOLE 1825 days (not just 2007-2009) rather than adding
# variable-cadence output logic to main.py/file_module.py for a few GB of savings --
# see run_definitive_full.sh's docstring for the identical reasoning on why the
# model's own integration timestep (DELTI) is unaffected by CGEM_TS (save cadence
# only); the extra disk is the 2006 pass-through's share, not a bigger run than
# needed.
#
# Usage:
#   tools/run_osse_nature.sh                    # colville+kuparuk+sagavanirktok
#   tools/run_osse_nature.sh kuparuk             # just one
#   SERIAL=1 tools/run_osse_nature.sh            # one at a time
#
# Results land in runs/osse_nature/<site>/ ; stdout in
# runs/osse_nature/<site>/run.log. Separate from runs/definitive_full/ (which stays
# at 2005-2006 for the publication-figure Fig-3 calibration, untouched by this) and
# from runs/interannual_full/ (the 19-year DAILY-cadence record, untouched).

set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
CODE="$ROOT/code"
SITES=("${@:-}")
if [ -z "${SITES[*]}" ]; then
    SITES=(colville kuparuk sagavanirktok)
fi

MAXT_DAYS="${CGEM_MAXT_DAYS:-1825}"
WARMUP_DAYS="${CGEM_WARMUP_DAYS:-365}"
TS="${CGEM_TS:-12}"

pids=()
for site in "${SITES[@]}"; do
    outdir="$ROOT/runs/osse_nature/$site"
    mkdir -p "$outdir"
    echo "-> ${site}_interannual_full  (output: runs/osse_nature/$site, $MAXT_DAYS days @ 6-min, i.e. real 2005-2009, study years 2007-2009)"

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
    [ $fail -eq 0 ] && echo "all runs finished" || { echo "a run FAILED -- check runs/osse_nature/*/run.log"; exit 1; }
fi
