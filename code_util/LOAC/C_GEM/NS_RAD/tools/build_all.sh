#!/usr/bin/env bash
# One command: run C-GEM for all four North Slope rivers, then regenerate every
# figure PDF and the movies from the fresh output.
#
#   tools/build_all.sh                   # run the 4 rivers, then all PDFs + movies
#   tools/build_all.sh --figures-only    # skip the model run; rebuild figures from existing runs
#   tools/build_all.sh --with-idealized  # ALSO run + verify + figure the idealized fixture
#   tools/build_all.sh --with-regression # ALSO build the regression-boundary runs the
#                                        #   validation PDF needs (+~20 min; see below)
#   SERIAL=1 tools/build_all.sh          # run the rivers one at a time (passed to run_sites.sh)
#
# THE VALIDATION PDF NEEDS A SECOND SET OF RUNS. runs/regression_bnd/ holds Kuparuk and
# Sagavanirktok rerun with the pure air-temperature regression as the upstream T boundary,
# instead of their own observed-blended series -- otherwise their modelled temperature
# would be scored against the very gauge record that feeds their boundary. Without it the
# report builds 5 of 6 parts and the model-vs-observation section is missing. It is not
# rebuilt by default because it doubles the run time and only changes when the temperature
# forcing does; pass --with-regression, or run tools/run_regression_bnd.sh once.
#
# The model run is the slow part (~15 min, all four in parallel). Figures take a
# couple of minutes. A single figure failing does not abort the rest — a PASS/FAIL
# summary is printed at the end and the exit code is non-zero if anything failed.
#
# Outputs:
#   runs/definitive/<site>/output.nc                  model state (per river)
#   docs/ns_rad_report.pdf                            THE deliverable: the six section
#                                                     PDFs below, stitched into one document
#   docs/movies/<site>_evolution.mp4 / _hovmoller.mp4 animations
#
# The six section PDFs (diagnostics, validation, geometry, model_schematic,
# model_summary, river_networks) are built as intermediates, stitched into
# ns_rad_report.pdf, then DELETED — only the combined report is kept.
set -uo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT"
PY="${PYTHON:-python3}"

FIGURES_ONLY=0
WITH_IDEALIZED=0
WITH_REGRESSION=0
for a in "$@"; do
    case "$a" in
        --figures-only) FIGURES_ONLY=1 ;;
        --with-idealized) WITH_IDEALIZED=1 ;;
        --with-regression) WITH_REGRESSION=1 ;;
        *) echo "unknown option: $a" >&2; exit 2 ;;
    esac
done

declare -a RESULTS
step() {                     # step "label" cmd...
    local label="$1"; shift
    echo; echo "==> $label"
    if "$@"; then RESULTS+=("PASS  $label"); else RESULTS+=("FAIL  $label"); fi
}

# --- 1. the model (all four rivers) ---
if [ "$FIGURES_ONLY" -eq 0 ]; then
    echo "### running C-GEM for all four rivers (this is the slow part) ###"
    if tools/run_sites.sh; then RESULTS+=("PASS  model run (4 rivers)")
    else RESULTS+=("FAIL  model run (4 rivers)"); echo "model run failed — see runs/definitive/*/run.log" >&2; fi
else
    echo "### --figures-only: skipping the model run ###"
fi

# --- 2. figures that DON'T need run output (always safe) ---
step "geometry PDF"           "$PY" tools/make_geometry_pdf.py
step "schematic PDF"          "$PY" tools/make_schematic_pdf.py
step "summary PDF"            "$PY" tools/make_summary_figures.py
step "river-network maps PDF" "$PY" tools/make_river_maps.py

# --- 3. figures that read runs/definitive ---
step "diagnostics PDF"        "$PY" tools/make_diagnostics_pdf.py
step "movies (4 rivers)"      "$PY" tools/make_movies.py

# --- 4. validation PDF needs the independent regression-boundary run ---
# Kuparuk + Sagavanirktok rerun with the air-temperature regression as the upstream T
# boundary instead of their own observed-blended series, so the temperature validation is
# an independent skill test rather than a comparison against the model's own input.
# Built on --with-regression, or once by hand with tools/run_regression_bnd.sh.
if [ "$WITH_REGRESSION" -eq 1 ] && [ "$FIGURES_ONLY" -eq 0 ]; then
    echo; echo "### building the regression-boundary runs (~20 min) ###"
    if tools/run_regression_bnd.sh; then RESULTS+=("PASS  regression-boundary runs")
    else RESULTS+=("FAIL  regression-boundary runs"); fi
elif [ "$WITH_REGRESSION" -eq 1 ]; then
    echo "### --figures-only: not rebuilding the regression-boundary runs ###"
fi

# (compgen -G is true if the glob matches ANYTHING; the regression run writes T.dat,
#  and the definitive run may be .nc or .dat — either present is enough.)
if compgen -G "runs/regression_bnd/*/T.dat" >/dev/null \
   || compgen -G "runs/regression_bnd/*/output.nc" >/dev/null; then
    step "validation PDF"     "$PY" tools/make_validation_pdf.py
else
    echo; echo "==> validation PDF: SKIPPED — runs/regression_bnd is empty, so the report"
    echo "    will build 5 of 6 parts and the model-vs-observation section will be MISSING."
    echo "    Fix: tools/run_regression_bnd.sh   (once, ~20 min), or --with-regression."
    RESULTS+=("SKIP  validation PDF (no runs/regression_bnd — report will be incomplete)")
fi

# --- 4b. combine the section PDFs into the single deliverable, then drop the parts ---
# (merges whichever of the six exist, in reading order -> docs/ns_rad_report.pdf)
step "combined report PDF"    "$PY" tools/make_report.py
# Keep only the stitched report: the six section PDFs are throwaway intermediates.
if [ -f docs/ns_rad_report.pdf ]; then
    rm -f docs/ns_rad_model_summary.pdf docs/ns_rad_model_schematic.pdf \
          docs/ns_rad_geometry.pdf docs/ns_rad_river_networks.pdf \
          docs/ns_rad_diagnostics.pdf docs/ns_rad_validation.pdf
    echo "    (removed the six section PDFs — ns_rad_report.pdf is the kept deliverable)"
fi

# --- 5. optional: the idealized verification fixture (separate ~15 min run) ---
if [ "$WITH_IDEALIZED" -eq 1 ]; then
    echo; echo "### --with-idealized: running + verifying the idealized fixture ###"
    step "idealized verify (full seasonal run)" "$PY" tools/verify_idealized.py --full
    step "idealized verification PDF"           "$PY" tools/make_idealized_verification_pdf.py
    step "idealized movies"                     "$PY" tools/make_movies.py --run runs/idealized idealized
fi

# --- summary ---
echo; echo "===================== build_all summary ====================="
fail=0
for r in "${RESULTS[@]}"; do
    printf '  %s\n' "$r"
    [[ "$r" == FAIL* ]] && fail=1
done
echo "============================================================="
[ "$fail" -eq 0 ] && echo "all steps OK" || echo "SOME STEPS FAILED (see above)"
exit "$fail"
