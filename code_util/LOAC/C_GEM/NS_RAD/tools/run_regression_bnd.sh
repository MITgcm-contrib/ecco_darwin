#!/usr/bin/env bash
# Produce the REGRESSION-BOUNDARY runs that the validation document needs
# (runs/regression_bnd/<site>/), consumed only by tools/make_validation_pdf.py.
#
# WHY THIS RUN EXISTS. Kuparuk and Sagavanirktok use their own observed-blended water
# temperature as the upstream boundary (sites/<name>.py WATERTEMP_FILE). Validating their
# modelled temperature against the very USGS record that feeds that boundary would be
# near-tautological -- the model is partly being scored on its own input. Rerunning them
# with the pure air-temperature regression (river_watertemp_2022_degC.csv, which never
# sees the gauge) makes the comparison an INDEPENDENT skill test. Same physics, same
# geometry, different boundary source -- it is not an older or lesser model version.
#
# ONLY these two sites are built. Colville and Canning already use the regression file, so
# their regression-boundary run would be byte-identical to the definitive one, and
# make_validation_pdf only reads the two sites below anyway.
#
# Usage:
#   tools/run_regression_bnd.sh          # both sites, in parallel
#   tools/run_regression_bnd.sh kuparuk  # just one
#   SERIAL=1 tools/run_regression_bnd.sh # one at a time
#
# ~20 min. Results in runs/regression_bnd/<site>/ ; stdout in that dir's run.log.
# Rebuild the report afterwards with: tools/build_all.sh --figures-only

set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
CODE="$ROOT/code"
REGRESSION_FILE="river_watertemp_2022_degC.csv"

SITES=("${@:-}")
if [ -z "${SITES[*]}" ]; then
    SITES=(kuparuk sagavanirktok)
fi

pids=()
for site in "${SITES[@]}"; do
    outdir="$ROOT/runs/regression_bnd/$site"
    mkdir -p "$outdir"
    echo "-> $site  (regression T boundary: $REGRESSION_FILE)"
    (
        cd "$outdir"
        # Literal, separate VAR=value assignments: zsh does NOT word-split unquoted
        # expansions, so a packed "$EXTRA" string silently becomes one malformed
        # assignment and the override is lost. config echoes what it actually used.
        env CGEM_SITE="$site" \
            CGEM_WATERTEMP_FILE="$REGRESSION_FILE" \
            PYTHONPATH="$CODE" \
            PYTHONWARNINGS=default \
            python3 "$CODE/main.py" > run.log 2>&1
    ) &
    pids+=($!)
    if [ -n "${SERIAL:-}" ]; then wait "${pids[-1]}"; fi
done

fail=0
for pid in "${pids[@]}"; do wait "$pid" || fail=1; done

echo
for site in "${SITES[@]}"; do
    log="$ROOT/runs/regression_bnd/$site/run.log"
    # Verify the override actually took effect -- a run that silently fell back to the
    # site's observed boundary would produce plausible output that quietly defeats the
    # whole point of this run.
    if grep -q "\[watertemp\] upstream T boundary OVERRIDDEN" "$log" 2>/dev/null; then
        echo "  OK   $site  $(grep -m1 '\[watertemp\]' "$log")"
    elif grep -q "^t: " "$log" 2>/dev/null; then
        echo "  WARN $site ran but did NOT report the boundary override -- it may have"
        echo "       used its own observed series, which defeats the independence check."
        fail=1
    else
        echo "  FAIL $site -- see $log"
        fail=1
    fi
done

exit "$fail"
