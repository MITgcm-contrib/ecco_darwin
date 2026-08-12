#!/bin/bash
# M1QN3 offline optimization driver for the 1-year (34560 timestep) SOCAT
# pCO2 adjoint assimilation. Each loop pass = 1 full forward+adjoint model
# evaluation (~7.3hr on a single core, see README "Performance") + 1 cheap
# optim.x call. Adapted from optim_m1qn3/src/run.sh.
#
# Run this from the experiment's run directory (with mitgcmuv_ad, optim.x,
# and all input files staged/symlinked in, per README.md "Building and
# running").
cd "$(dirname "$0")"

myiter=0
nIter=10
nSim=10
nn=$(( myiter + nIter*nSim + 1 ))

while (( myiter < nn ))
do
    it=$(printf "%03d" $myiter)
    echo "$(date) -- starting optimcycle ${myiter}"

    sed -i '' "s/^[[:space:]]*optimcycle[[:space:]]*=.*/ optimcycle=${myiter},/" data.optim
    grep -q "optimcycle=${myiter}," data.optim || { echo "FATAL: sed failed to set optimcycle=${myiter} in data.optim"; exit 1; }

    ./mitgcmuv_ad > output_optim_${it}.txt 2>&1
    status=$?
    if (( status != 0 )); then
        echo "FATAL: mitgcmuv_ad exited with status ${status} on optimcycle ${myiter} -- see output_optim_${it}.txt"
        exit 1
    fi
    echo "$(date) -- model run done for optimcycle ${myiter}"
    grep -q "global fc " output_optim_${it}.txt || { echo "FATAL: no 'global fc' found in output_optim_${it}.txt -- model likely did not complete"; exit 1; }
    grep "global fc " output_optim_${it}.txt

    # Only the lines optim.x appends THIS pass count toward output-mode
    # detection -- m1qn3_output.txt is append-only across warm restarts
    # (optim_sub.F), so a stale "output mode" line from an earlier run
    # must never be re-matched as if it were this cycle's result.
    prevlines=$(wc -l < m1qn3_output.txt 2>/dev/null || echo 0)

    ./optim.x > opt_${it}.txt 2>&1
    status=$?
    if (( status != 0 )); then
        echo "FATAL: optim.x exited with status ${status} on optimcycle ${myiter} -- see opt_${it}.txt"
        exit 1
    fi
    echo "$(date) -- optim.x done for optimcycle ${myiter}"

    if tail -n +$((prevlines + 1)) m1qn3_output.txt 2>/dev/null | grep -q "m1qn3: output mode"; then
        echo "m1qn3 has finished (output mode reached)"
        break
    fi
    ((myiter++))
done
echo "$(date) -- optimization loop finished, last optimcycle=${myiter}"
