#!/bin/bash
# Hybrid Darwin(forward)/pkg-dic(adjoint) loop, HAND-ROLLED gradient
# descent -- bypasses M1QN3/optim.x entirely, working around a
# reproducible bug in MITgcm's generic control-array (genarr3d) adjoint
# where whichever field lands in certain array slots comes back with
# extreme outlier gradient values, poisoning M1QN3's step sizing (see
# scripts/manual_gradient_step.py's docstring, or the conversation log,
# for the full diagnosis).
#
# Each cycle:
#   1. pkg/dic forward+adjoint with ONLY DIC registered (run_dic_A,
#      slot 1 -- empirically clean).
#   2. pkg/dic forward+adjoint with ONLY ALK registered (run_dic_B,
#      slot 1 -- empirically clean).
#   3. manual_gradient_step.py: read both adxx fields directly (never
#      touching M1QN3/optim.x), take a clipped steepest-descent step,
#      write the updated physical state into run_dic_A/run_dic_B's own
#      next-cycle ICs AND Darwin's IC.
#   4. Darwin forward run on the updated IC.
#   5. Compute Darwin's own SOCAT cost (logging/monitoring only -- no
#      formal line search/convergence check, this is a fixed-length
#      pilot).
set -e
HYBRID="$(cd "$(dirname "$0")/.." && pwd)"

NCYCLES=5
START_CYCLE="${1:-0}"

echo "$(date) -- hybrid Darwin/pkg-dic loop (manual gradient descent) starting at cycle ${START_CYCLE}"

for (( myiter=START_CYCLE; myiter<NCYCLES; myiter++ ))
do
    it=$(printf "%03d" $myiter)
    echo "=========================================="
    echo "$(date) -- starting hybrid cycle ${myiter}"
    echo "=========================================="

    # optimcycle stays 0 in run_dic_A/run_dic_B for EVERY cycle: we never use
    # M1QN3's xx/ecco_ctrl history, since our own physical state update is
    # already baked directly into dic_init.bin/alk_init.bin each cycle.
    # optimcycle>0 makes pkg/ctrl try to "unpack" a prior ecco_ctrl file that
    # M1QN3 would normally have written -- which doesn't exist here, and
    # crashed cycle 1 the first time this ran. Since every cycle now reuses
    # the SAME optimcycle=0 filenames, clean up the previous cycle's
    # ecco_cost/ecco_ctrl/xx/adxx files first so mitgcmuv_ad never trips
    # over a stale file it expects not to exist yet (the exact failure mode
    # that crashed the very first M1QN3-based attempt on this project).
    for D in run_dic_A run_dic_B; do
        rm -f "$HYBRID/$D"/ecco_cost_MIT_CE_000.opt* "$HYBRID/$D"/ecco_ctrl_MIT_CE_000.opt* "$HYBRID/$D"/OPWARM.opt*
        rm -f "$HYBRID/$D"/xx_ptr*.0000000000.* "$HYBRID/$D"/adxx_ptr*.0000000000.*
    done

    echo "$(date) -- [1/5] pkg/dic forward+adjoint, DIC only (run_dic_A)"
    ( cd "$HYBRID/run_dic_A" && ./mitgcmuv_ad > "output_optim_${it}.txt" 2>&1 )
    grep -q "global fc " "$HYBRID/run_dic_A/output_optim_${it}.txt" || { echo "FATAL: no 'global fc' in run_dic_A/output_optim_${it}.txt"; exit 1; }
    fc_a=$(grep "global fc " "$HYBRID/run_dic_A/output_optim_${it}.txt")
    echo "  pkg/dic (DIC-only) fc: ${fc_a}"

    echo "$(date) -- [2/5] pkg/dic forward+adjoint, ALK only (run_dic_B)"
    ( cd "$HYBRID/run_dic_B" && ./mitgcmuv_ad > "output_optim_${it}.txt" 2>&1 )
    grep -q "global fc " "$HYBRID/run_dic_B/output_optim_${it}.txt" || { echo "FATAL: no 'global fc' in run_dic_B/output_optim_${it}.txt"; exit 1; }
    fc_b=$(grep "global fc " "$HYBRID/run_dic_B/output_optim_${it}.txt")
    echo "  pkg/dic (ALK-only) fc: ${fc_b}"

    echo "$(date) -- [3/5] manual gradient-descent step (bypassing M1QN3)"
    python3 "$HYBRID/scripts/manual_gradient_step.py" "${myiter}"

    echo "$(date) -- [4/5] Darwin forward run (mitgcmuv)"
    ( cd "$HYBRID/run_darwin" && ./mitgcmuv > "output_darwin_${it}.txt" 2>&1 )
    grep -q "STOP NORMAL END" "$HYBRID/run_darwin/output_darwin_${it}.txt" || { echo "FATAL: Darwin did not reach STOP NORMAL END on cycle ${myiter}"; exit 1; }

    echo "$(date) -- [5/5] computing Darwin's own SOCAT cost"
    darwin_J=$(python3 "$HYBRID/scripts/compute_darwin_cost_dic.py" "${myiter}" | tail -1 | grep -oE '[0-9.eE+-]+$')
    echo "  Darwin's own: J = ${darwin_J}"

    echo "$(date) -- cycle ${myiter} done (pkg/dic DIC-only ${fc_a} | ALK-only ${fc_b} | Darwin J=${darwin_J})"

    echo "$(date) -- cleaning up cycle ${myiter}'s bulky adjoint tapes"
    rm -f "$HYBRID/run_dic_A"/tapes2.000.data "$HYBRID/run_dic_A"/tapes2.000.meta "$HYBRID/run_dic_A"/tapes3.000.data "$HYBRID/run_dic_A"/tapes3.000.meta
    rm -f "$HYBRID/run_dic_B"/tapes2.000.data "$HYBRID/run_dic_B"/tapes2.000.meta "$HYBRID/run_dic_B"/tapes3.000.data "$HYBRID/run_dic_B"/tapes3.000.meta
done
echo "$(date) -- hybrid loop (manual gradient descent) finished, ${NCYCLES} cycles"
