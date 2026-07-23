#!/usr/bin/env bash
# Run C-GEM for one or more North Slope sites, each in its own output directory.
#
# main.py deletes every *.dat in the CURRENT working directory at startup and writes
# its results there, so every site must run from a separate directory or they will
# clobber each other. That is what this script arranges.
#
# The four runs are fully independent -- they share no state and write no common
# files -- so they are launched in parallel by default. Because variables.py
# allocates all model state at import time from config.M, separate processes are
# also the only way to run multiple sites without a deep refactor.
#
# Usage:
#   tools/run_sites.sh                          # all four, in parallel
#   tools/run_sites.sh kuparuk                  # just one
#   tools/run_sites.sh colville kuparuk         # a subset
#   SERIAL=1 tools/run_sites.sh                 # one at a time (easier to read logs)
#
# Results land in runs/definitive/<site>/ ; stdout in runs/definitive/<site>/run.log
# (runs/ is organized by run-type: definitive/ = these operational runs, regression_bnd/
#  = the independent-boundary check read by make_validation_pdf.py, idealized/ = the
#  verification fixture.)

set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
CODE="$ROOT/code"
SITES=("${@:-}")
if [ -z "${SITES[*]}" ]; then
    SITES=(colville kuparuk sagavanirktok canning)
fi

pids=()
for site in "${SITES[@]}"; do
    outdir="$ROOT/runs/definitive/$site"
    mkdir -p "$outdir"
    echo "-> $site  (output: runs/definitive/$site)"

    # PYTHONPATH so the model's modules resolve while the cwd is the output dir.
    # PYTHONWARNINGS so the placeholder-geometry warning is not swallowed.
    if [ -n "${SERIAL:-}" ]; then
        ( cd "$outdir" && CGEM_SITE="$site" PYTHONPATH="$CODE" \
              PYTHONWARNINGS="default" python3 "$CODE/main.py" 2>&1 | tee run.log )
    else
        ( cd "$outdir" && CGEM_SITE="$site" PYTHONPATH="$CODE" \
              PYTHONWARNINGS="default" python3 "$CODE/main.py" > run.log 2>&1 ) &
        pids+=($!)
    fi
done

if [ ${#pids[@]} -gt 0 ]; then
    echo "launched ${#pids[@]} run(s); waiting..."
    fail=0
    for pid in "${pids[@]}"; do wait "$pid" || fail=1; done
    [ $fail -eq 0 ] && echo "all runs finished" || { echo "a run FAILED -- check runs/definitive/*/run.log"; exit 1; }
fi
