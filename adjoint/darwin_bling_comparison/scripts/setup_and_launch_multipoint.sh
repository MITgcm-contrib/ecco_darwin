#!/bin/bash
# Set up and launch the multi-point FD-sensitivity extension: 6 tracers x
# (+eta/-eta) plus 1 unperturbed baseline = 13 concurrent runs, each with
# ALL 38 validation points perturbed simultaneously (well-separated, so no
# cross-contamination) using Kay-style relative (+/-5%) perturbations.
# M1QN3 campaign is done and disk has headroom, so these run fully
# concurrently (14 cores available) instead of the old sequential pattern.
#
# Configure via environment variable:
#   DARWIN_RUN_DIR -- Darwin experiment dir (default: current directory)
set -e
DARWIN="${DARWIN_RUN_DIR:-.}"
cd "$DARWIN"

POINTS_JSON='[[4,109],[6,63],[6,74],[6,87],[6,97],[7,10],[9,120],[10,25],[10,35],[13,53],[17,12],[17,63],[19,1],[19,101],[20,23],[24,37],[25,73],[25,91],[25,119],[27,55],[30,14],[30,25],[34,109],[38,77],[38,121],[39,46],[39,67],[39,89],[40,22],[41,56],[44,101],[44,111],[48,121],[50,80],[51,66],[54,108],[55,3],[40,97]]'
ETA=0.05

icfile_for () {
  case "$1" in
    DIC) echo dic_darwin_init.bin ;;
    ALK) echo alk_darwin_init.bin ;;
    O2)  echo o2_darwin_init.bin ;;
    NO3) echo no3_darwin_init.bin ;;
    PO4) echo po4_darwin_init.bin ;;
    FeT) echo fet_darwin_init.bin ;;
  esac
}

make_rundir () {
  local NAME="$1"
  local DIR="run_multi_${NAME}"
  mkdir -p "$DIR"
  for f in input_ad/*; do
    base=$(basename "$f")
    [ "$base" = "data" ] && continue
    ln -sf "../$f" "$DIR/$base"
  done
  ln -sf ../build/mitgcmuv "$DIR/mitgcmuv"
  sed -e 's/nTimeSteps = 34560,/nTimeSteps = 2880,/' \
      -e 's/pChkptFreq = 216000.,/pChkptFreq = 0.,/' \
      -e 's/chkptFreq  = 216000.,/chkptFreq  = 0.,/' \
      -e 's/dumpFreq   = 216000.,/dumpFreq   = 0.,/' \
      input_ad/data > "$DIR/data"
}

echo "$(date) -- setting up run directories"

make_rundir baseline

for T in DIC ALK O2 NO3 PO4 FeT; do
  for SIGN in plus minus; do
    make_rundir "${T}_${SIGN}"
  done
done

# perturb each tracer's file in its plus/minus dirs (replacing the symlink
# with a real perturbed file); baseline dir keeps all symlinks untouched
for T in DIC ALK O2 NO3 PO4 FeT; do
  IC=$(icfile_for "$T")
  rm -f "run_multi_${T}_plus/$IC" "run_multi_${T}_minus/$IC"
  python3 scripts/make_multipoint_perturbed_ic.py \
    "input_ad/$IC" "run_multi_${T}_plus/$IC" "run_multi_${T}_minus/$IC" \
    "$POINTS_JSON" "$ETA" "run_multi_${T}_plus/perturb_log.json"
  cp "run_multi_${T}_plus/perturb_log.json" "run_multi_${T}_minus/perturb_log.json"
done

echo "$(date) -- launching all 13 runs concurrently"

PIDS=()
NAMES=()
for RUNDIR in run_multi_baseline run_multi_DIC_plus run_multi_DIC_minus \
              run_multi_ALK_plus run_multi_ALK_minus run_multi_O2_plus run_multi_O2_minus \
              run_multi_NO3_plus run_multi_NO3_minus run_multi_PO4_plus run_multi_PO4_minus \
              run_multi_FeT_plus run_multi_FeT_minus; do
  ( cd "$RUNDIR" && ./mitgcmuv > output.txt 2>&1 ) &
  PIDS+=($!)
  NAMES+=("$RUNDIR")
  echo "  launched $RUNDIR (pid $!)"
done

echo "$(date) -- waiting for all 13 runs"
FAIL=0
for idx in "${!PIDS[@]}"; do
  if ! wait "${PIDS[$idx]}"; then
    echo "WARNING: ${NAMES[$idx]} exited nonzero"
    FAIL=1
  fi
done

echo "$(date) -- checking STOP NORMAL END for all runs"
for RUNDIR in run_multi_baseline run_multi_DIC_plus run_multi_DIC_minus \
              run_multi_ALK_plus run_multi_ALK_minus run_multi_O2_plus run_multi_O2_minus \
              run_multi_NO3_plus run_multi_NO3_minus run_multi_PO4_plus run_multi_PO4_minus \
              run_multi_FeT_plus run_multi_FeT_minus; do
  if grep -q "STOP NORMAL END" "$RUNDIR/output.txt"; then
    echo "  OK   $RUNDIR"
  else
    echo "  FAIL $RUNDIR -- see $RUNDIR/output.txt"
  fi
done

echo "$(date) -- MULTIPOINT SWEEP DONE"
