#!/bin/bash
# Full-year hybrid Darwin(forward)/BLING(adjoint) surrogate-gradient loop.
# Same mechanism as run_hybrid_loop.sh (see that file's header for the
# full 6-step description), scaled to nTimeSteps=34560 (1 year) and the
# full 12-month SOCAT cost (all of data.profiles' mult_profiles(m,13)=1,
# unmodified from the real BLING campaign) instead of the 1-month pilot's
# restricted single-month cost.
#
# numiter=1 in run_bling_1yr/data.optim -- mechanical shakeout of the
# full-year config (~2.4 days for this one cycle: ~7.3h BLING adjoint +
# ~49h Darwin forward) before committing to a multi-cycle campaign.
#
# After each cycle, BLING's own adjoint tapes (tapes2.000.data,
# tapes3.000.data, ~9.75GB/cycle) are deleted once that cycle's cost has
# been safely extracted and injected -- they're pure TAF checkpoint-
# recompute scratch space, not needed after the adjoint sweep completes,
# and disk (12GB free at the start of this pilot) can't hold more than
# about one cycle's worth.
set -e
HYBRID="$(cd "$(dirname "$0")/.." && pwd)"
cd "$HYBRID/run_bling_1yr"

myiter=0
nIter=1
nSim=5
nn=$(( myiter + nIter*nSim + 1 ))

echo "$(date) -- full-year hybrid loop starting"

while (( myiter < nn ))
do
    it=$(printf "%03d" $myiter)
    echo "=========================================="
    echo "$(date) -- starting full-year hybrid cycle ${myiter}"
    echo "=========================================="

    sed -i '' "s/^[[:space:]]*optimcycle[[:space:]]*=.*/ optimcycle=${myiter},/" data.optim
    grep -q "optimcycle=${myiter}," data.optim || { echo "FATAL: sed failed to set optimcycle=${myiter} in data.optim"; exit 1; }

    echo "$(date) -- [1/7] BLING full-year forward+adjoint (mitgcmuv_ad, ~7.3h)"
    ./mitgcmuv_ad > "output_optim_${it}.txt" 2>&1
    status=$?
    if (( status != 0 )); then
        echo "FATAL: mitgcmuv_ad exited with status ${status} on cycle ${myiter} -- see output_optim_${it}.txt"
        exit 1
    fi
    grep -q "global fc " "output_optim_${it}.txt" || { echo "FATAL: no 'global fc' in output_optim_${it}.txt"; exit 1; }
    bling_fc_line=$(grep "global fc " "output_optim_${it}.txt")
    echo "  BLING's own: ${bling_fc_line}"

    echo "$(date) -- [2/7] applying control guess to Darwin's physical IC"
    python3 "$HYBRID/scripts/apply_control_to_darwin_ic_1yr.py" "${myiter}"

    echo "$(date) -- [3/7] Darwin full-year forward run (mitgcmuv, ~49h)"
    ( cd "$HYBRID/run_darwin_1yr" && ./mitgcmuv > "output_darwin_${it}.txt" 2>&1 )
    darwin_status=$?
    if (( darwin_status != 0 )); then
        echo "FATAL: Darwin mitgcmuv exited with status ${darwin_status} on cycle ${myiter} -- see run_darwin_1yr/output_darwin_${it}.txt"
        exit 1
    fi
    grep -q "STOP NORMAL END" "$HYBRID/run_darwin_1yr/output_darwin_${it}.txt" || { echo "FATAL: Darwin did not reach STOP NORMAL END on cycle ${myiter}"; exit 1; }

    echo "$(date) -- [4/7] computing Darwin's own 12-month SOCAT cost"
    darwin_J=$(python3 "$HYBRID/scripts/compute_darwin_cost_1yr.py" "${myiter}" | tail -1 | grep -oE 'Darwin total J = [0-9.eE+-]+' | grep -oE '[0-9.eE+-]+$')
    echo "  Darwin's own: J = ${darwin_J}"

    echo "$(date) -- [5/7] patching ecco_cost with Darwin's J (BLING gradient untouched)"
    python3 "$HYBRID/scripts/patch_ecco_cost_1yr.py" "${myiter}" "${darwin_J}"

    echo "$(date) -- [6/7] optim.x (M1QN3 step)"
    prevlines=$(wc -l < m1qn3_output.txt 2>/dev/null || echo 0)
    ./optim.x > "opt_${it}.txt" 2>&1
    status=$?
    if (( status != 0 )); then
        echo "FATAL: optim.x exited with status ${status} on cycle ${myiter} -- see opt_${it}.txt"
        exit 1
    fi
    echo "$(date) -- cycle ${myiter} done (BLING ${bling_fc_line} | Darwin J=${darwin_J})"

    echo "$(date) -- [7/7] cleaning up cycle ${myiter}'s bulky BLING adjoint tapes"
    rm -f tapes2.000.data tapes2.000.meta tapes3.000.data tapes3.000.meta
    echo "  disk: $(df -h / | tail -1)"

    if tail -n +$((prevlines + 1)) m1qn3_output.txt 2>/dev/null | grep -q "m1qn3: output mode"; then
        echo "m1qn3 has finished (output mode reached)"
        break
    fi
    ((myiter++))
done
echo "$(date) -- full-year hybrid loop finished, last cycle=${myiter}"
