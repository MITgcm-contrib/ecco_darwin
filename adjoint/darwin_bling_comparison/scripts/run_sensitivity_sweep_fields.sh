#!/bin/bash
# Full 6-tracer field sweep at a single grid point (j=40,i=97), matched
# ICs (run scripts/match_ics_to_bling.py first). This eliminates a ~7%
# ALK / ~2% DIC provenance mismatch between the two models' starting
# states that was found to explain most of a ~2x discrepancy in
# d(pCO2)/dDIC sensitivity between Darwin's forward-FD and BLING's
# adjoint. All 6 tracers need a full rerun since every tracer's baseline
# state changes simultaneously when the ICs are matched (not just the one
# being perturbed), and this captures the FULL darwinPCO2Snap field (not
# just the point) for a proper cost-function FD comparison against
# BLING's adjoint. Physics/duration (30-day, nTimeSteps=2880) unchanged.
#
# Superseded by setup_and_launch_multipoint.sh for real use (relative
# Kay-style perturbation, many points at once, proper unperturbed
# baseline) -- kept here as the simpler single-point reference version.
#
# Prerequisite: run_plus/ and run_minus/ directories must already exist,
# each with input_ad/* symlinked in, a materialized `data` with
# nTimeSteps=2880 and checkpointing disabled, and mitgcmuv symlinked from
# ../build/mitgcmuv -- see the make_rundir() function in
# setup_and_launch_multipoint.sh for the exact pattern to follow.
#
# Configure via environment variable:
#   DARWIN_RUN_DIR -- Darwin experiment dir (default: current directory)
set -e
DARWIN="${DARWIN_RUN_DIR:-.}"
cd "$DARWIN"
mkdir -p pco2_fields

J=40; I=97

run_tracer () {
  local NAME="$1" ICFILE="$2" EPS="$3"
  echo "=========================================="
  echo "$(date) -- starting $NAME (eps=$EPS)"
  echo "=========================================="

  rm -f "run_plus/$ICFILE" "run_minus/$ICFILE"
  python3 scripts/make_perturbed_ic.py "input_ad/$ICFILE" "run_plus/$ICFILE" "run_minus/$ICFILE" "$J" "$I" "$EPS"

  ( cd run_plus && ./mitgcmuv > "output_field_${NAME}.txt" 2>&1 ) &
  local PID_PLUS=$!
  ( cd run_minus && ./mitgcmuv > "output_field_${NAME}.txt" 2>&1 ) &
  local PID_MINUS=$!

  local STATUS_PLUS=0 STATUS_MINUS=0
  wait $PID_PLUS || STATUS_PLUS=$?
  wait $PID_MINUS || STATUS_MINUS=$?
  echo "$(date) -- $NAME runs finished (plus=$STATUS_PLUS minus=$STATUS_MINUS)"

  if grep -q "STOP NORMAL END" "run_plus/output_field_${NAME}.txt" && grep -q "STOP NORMAL END" "run_minus/output_field_${NAME}.txt"; then
    for d in run_plus run_minus; do
      for f in "$d"/darwinPCO2Snap.0000001392.*.data "$d"/darwinPCO2Snap.0000001392.*.meta; do
        base=$(basename "$f")
        cp "$f" "pco2_fields/${NAME}_${d}_${base}"
      done
    done
    echo "$NAME: field saved to pco2_fields/"
  else
    echo "$NAME: WARNING -- one or both runs did not reach STOP NORMAL END, see output_field_${NAME}.txt"
  fi

  rm -f "run_plus/$ICFILE" "run_minus/$ICFILE"
  ln -sf "../input_ad/$ICFILE" "run_plus/$ICFILE"
  ln -sf "../input_ad/$ICFILE" "run_minus/$ICFILE"

  echo "$(date) -- $NAME done, disk: $(df -h / | tail -1)"
}

run_tracer DIC dic_darwin_init.bin 8.57
run_tracer ALK alk_darwin_init.bin 12.57
run_tracer O2  o2_darwin_init.bin  2.659
run_tracer NO3 no3_darwin_init.bin 0.3054
run_tracer PO4 po4_darwin_init.bin 0.02402
run_tracer FeT fet_darwin_init.bin 1.068e-5

echo "$(date) -- FIELD SWEEP DONE"
