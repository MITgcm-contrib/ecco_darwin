#!/bin/bash
# job_production_chain.sh -- self-resubmitting multi-segment production run
# for the Kerguelen 3-year (2020-2023) regional integration.
#
# Task #13 (2026-08-09). Built after: Task #22's real 720-rank throughput
# validation (2.322 s/step), and the deltaT investigation that confirmed
# dt=90 is stable through the cold-start adjustment transient and cuts
# total steps 3x vs the original dt=30 (see DECISIONS.md). At dt=90:
# 1,052,160 total steps, ~28.3 days total compute, chkptFreq=7200 steps
# (~4.64 wallclock hours/interval), ~146 checkpoint intervals total ->
# ~6-7 segments at a 115h internal cap under the `long` queue's 120h max.
#
# No wallClockLimit support in this MITgcm build (confirmed via source
# grep), so segments rely on periodic chkptFreq boundaries + this script's
# own internal timeout, accepting a small amount of recompute waste
# (up to one chkptFreq interval, ~4.64h) on any walltime-triggered kill.
#
# Restart mechanics: MITgcm's pickup reader literally string-formats
# nIter0 into 'pickup.<iter>.data' -- it does NOT read rolling
# (ckptA/ckptB) checkpoints' .meta timeStepNumber to match nIter0
# (confirmed against read_pickup.F/ptracers_read_pickup.F/etc, 2026-08-09).
# find_restart_iter.py handles this: verifies the newest complete
# checkpoint (permanent or rolling, by exact expected byte size + a
# consistent timeStepNumber across all 5 filesets), and renames rolling
# ckptA/ckptB filesets to the numbered convention when needed.
#
# Usage:
#   RUN=/path/to/run_parent_dir qsub -v SEGMENT=1,RUN -W group_list=<group> \
#     scripts/job_production_chain.sh
#
# RUN is the directory CONTAINING the run/ subdirectory and this script's
# companion find_restart_iter.py. Set it, or edit the default below.
#PBS -l select=6:ncpus=120:mpiprocs=120:model=rom_ait
#PBS -q long
#PBS -l walltime=120:00:00
#PBS -N kerreg_prod
#PBS -j oe

RUN=${RUN:?set RUN to the directory containing run/ and find_restart_iter.py}
cd "$RUN/run"

module purge
module load comp-intel/2023.2.1 mpi-hpe/mpt.2.30 hdf4/4.2.12 hdf5/1.8.18_mpt netcdf/4.4.1.1_mpt python3/3.9.12

SCRIPTDIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
SEGMENT=${SEGMENT:-1}
TOTAL_STEPS=1052160
STALL_THRESHOLD_S=10800   # 3h with no new time_tsnumber => assume stalled
POLL_INTERVAL_S=300

LOGDIR="$RUN/run/chain_logs"
mkdir -p "$LOGDIR"

echo "PROD_SEGMENT_START segment=$SEGMENT $(date +%s)"

# Record the nIter0 this segment started from, to detect zero-progress
# segments later (would otherwise resubmit forever from the same point).
# 2026-09-29 FIX: the old pipeline emitted TWO lines -- the second grep
# also matched the '0' inside the NAME "nIter0" before reaching the value
# (" nIter0=158400" -> "0" and "158400"). START_ITER then held a two-line
# string, so the zero-progress guard below degraded to
# `[ 309600 -eq "0\n158400" ]`, which bash rejects as "integer expression
# expected": the test returned false and the guard NEVER fired, so a truly
# stuck chain would have resubmitted from the same iteration forever.
# Visible above in chain_logs/pbs_out.txt as the split line
# "RESOLVED_RESTART_ITER 309600 (started this segment at 0 / 158400)".
# Found while porting this script to fiji.
START_ITER=$(sed -n 's/^[[:space:]]*nIter0=[[:space:]]*\([0-9][0-9]*\).*/\1/p' data | head -1)

# Archive the previous segment's logs before this run overwrites them.
if ls STDOUT.0000 >/dev/null 2>&1; then
  PREV=$((SEGMENT - 1))
  mkdir -p "$LOGDIR/segment_${PREV}"
  mv STDOUT.* STDERR.* "$LOGDIR/segment_${PREV}/" 2>/dev/null
fi
rm -f STDOUT.* STDERR.*

# --- stall-detection watchdog: kill mpiexec early if truly stuck, rather
# --- than burning the full 115h cap on a hang. 3h threshold is well above
# --- the ~4.64h checkpoint interval's own natural monitor-print cadence
# --- (monitorFreq=86400s=960 steps -> prints roughly every ~37min at this
# --- throughput), with margin against a checkpoint-write pause.
(
  LASTSTEP=-1
  LASTCHANGE=$(date +%s)
  while true; do
    sleep "$POLL_INTERVAL_S"
    [ -f STDOUT.0000 ] || continue
    STEP=$(grep "time_tsnumber" STDOUT.0000 2>/dev/null | tail -1 | grep -oE '[0-9]+$')
    NOW=$(date +%s)
    if [ -n "$STEP" ] && [ "$STEP" != "$LASTSTEP" ]; then
      LASTSTEP=$STEP
      LASTCHANGE=$NOW
    elif [ $((NOW - LASTCHANGE)) -gt "$STALL_THRESHOLD_S" ]; then
      echo "STALL_DETECTED last_step=$LASTSTEP no_progress_for=$((NOW - LASTCHANGE))s -- killing mpiexec" >> "$LOGDIR/watchdog.log"
      pkill -9 -f mitgcmuv
      break
    fi
  done
) &
WATCHDOG_PID=$!

# 115h local cap, 5h buffer under the 120h PBS walltime for
# checkpoint-resolution + resubmission overhead.
timeout --kill-after=120 414000 mpiexec -np 720 ./mitgcmuv
RC=$?

kill "$WATCHDOG_PID" 2>/dev/null

echo "MPIEXEC_EXIT_CODE $RC"
echo "PROD_SEGMENT_END segment=$SEGMENT $(date +%s)"

# --- hard-failure check: never resubmit past a real numerical failure ---
if grep -qi "forrtl\|severe\|ABNORMAL\|floating" STDOUT.0000 2>/dev/null; then
  echo "PROD_CHAIN_HALTED_HARD_ERROR segment=$SEGMENT -- manual inspection required"
  exit 1
fi

# --- resolve the furthest verified-complete checkpoint ---
RESTART_ITER=$(python3 "$SCRIPTDIR/find_restart_iter.py" --run "$RUN/run" 2>>"$LOGDIR/watchdog.log")
if [ -z "$RESTART_ITER" ] || ! [[ "$RESTART_ITER" =~ ^[0-9]+$ ]]; then
  echo "PROD_CHAIN_HALTED_NO_VALID_CHECKPOINT segment=$SEGMENT -- manual inspection required"
  exit 1
fi
echo "RESOLVED_RESTART_ITER $RESTART_ITER (started this segment at $START_ITER)"

if [ "$RESTART_ITER" -eq "$START_ITER" ]; then
  echo "PROD_CHAIN_HALTED_ZERO_PROGRESS segment=$SEGMENT stuck_at=$RESTART_ITER -- manual inspection required"
  exit 1
fi

# --- completion check ---
if [ "$RESTART_ITER" -ge "$TOTAL_STEPS" ] || grep -q "PROGRAM MAIN: Execution ended Normally" STDOUT.0000 2>/dev/null; then
  echo "PROD_CHAIN_COMPLETE segment=$SEGMENT final_iter=$RESTART_ITER"
  exit 0
fi

# --- advance nIter0 and resubmit ---
sed -i "s/^ nIter0=.*/ nIter0=${RESTART_ITER},/" data

NEXT=$((SEGMENT + 1))
qsub -v SEGMENT=$NEXT,RUN -W group_list="${PBS_O_GROUP:-$(id -gn)}" "$SCRIPTDIR/job_production_chain.sh"
echo "RESUBMITTED_AS_SEGMENT $NEXT"
