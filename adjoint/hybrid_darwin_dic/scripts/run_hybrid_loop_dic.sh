#!/bin/bash
# Hybrid Darwin(forward)/pkg-dic(adjoint) surrogate-gradient M1QN3 loop,
# ALL 5 shared tracers (DIC, ALK, PO4, DOP, O2 -- pkg/dic's own tracer
# set), 1-month SOCAT cost, mechanical pilot (numiter=5 in
# run_dic/data.optim). Mirrors hybrid_darwin_bling/scripts/run_hybrid_loop.sh
# exactly; see that script's header for the full step-by-step rationale.
# Each pass:
#   1. pkg/dic forward+adjoint (mitgcmuv_ad) at the current control guess
#      -> writes ecco_cost_MIT_CE_000.opt<cycle> (pkg/dic's own fc + the
#      packed 5-tracer adjoint gradient).
#   2. Apply that SAME control guess to Darwin's own physical IC
#      (apply_control_to_darwin_ic_dic.py).
#   3. Darwin forward run (mitgcmuv) on that IC.
#   4. Compute Darwin's own cost J_Darwin against SOCAT month01
#      (compute_darwin_cost_dic.py).
#   5. Binary-patch ecco_cost_MIT_CE_000.opt<cycle>'s fc field with
#      J_Darwin (patch_ecco_cost_dic.py) -- pkg/dic's gradient payload is
#      left untouched.
#   6. optim.x reads (J_Darwin, pkg/dic gradient) and proposes the next
#      control guess.
set -e
HYBRID="$(cd "$(dirname "$0")/.." && pwd)"
cd "$HYBRID/run_dic"

myiter=0
nIter=5     # matches data.optim's numiter -- outer safety cap only,
nSim=5      # actual stopping is via m1qn3's own "output mode" signal
nn=$(( myiter + nIter*nSim + 1 ))

echo "$(date) -- hybrid Darwin/pkg-dic loop starting"

while (( myiter < nn ))
do
    it=$(printf "%03d" $myiter)
    echo "=========================================="
    echo "$(date) -- starting hybrid cycle ${myiter}"
    echo "=========================================="

    sed -i '' "s/^[[:space:]]*optimcycle[[:space:]]*=.*/ optimcycle=${myiter},/" data.optim
    grep -q "optimcycle=${myiter}," data.optim || { echo "FATAL: sed failed to set optimcycle=${myiter} in data.optim"; exit 1; }

    echo "$(date) -- [1/6] pkg/dic forward+adjoint (mitgcmuv_ad)"
    ./mitgcmuv_ad > "output_optim_${it}.txt" 2>&1
    status=$?
    if (( status != 0 )); then
        echo "FATAL: mitgcmuv_ad exited with status ${status} on cycle ${myiter} -- see output_optim_${it}.txt"
        exit 1
    fi
    grep -q "global fc " "output_optim_${it}.txt" || { echo "FATAL: no 'global fc' in output_optim_${it}.txt"; exit 1; }
    dic_fc_line=$(grep "global fc " "output_optim_${it}.txt")
    echo "  pkg/dic's own: ${dic_fc_line}"

    echo "$(date) -- [2/6] applying control guess to Darwin's physical IC"
    python3 "$HYBRID/scripts/apply_control_to_darwin_ic_dic.py" "${myiter}"

    echo "$(date) -- [3/6] Darwin forward run (mitgcmuv)"
    ( cd "$HYBRID/run_darwin" && ./mitgcmuv > "output_darwin_${it}.txt" 2>&1 )
    darwin_status=$?
    if (( darwin_status != 0 )); then
        echo "FATAL: Darwin mitgcmuv exited with status ${darwin_status} on cycle ${myiter} -- see run_darwin/output_darwin_${it}.txt"
        exit 1
    fi
    grep -q "STOP NORMAL END" "$HYBRID/run_darwin/output_darwin_${it}.txt" || { echo "FATAL: Darwin did not reach STOP NORMAL END on cycle ${myiter}"; exit 1; }

    echo "$(date) -- [4/6] computing Darwin's own SOCAT cost"
    darwin_J=$(python3 "$HYBRID/scripts/compute_darwin_cost_dic.py" "${myiter}" | tail -1 | grep -oE '[0-9.eE+-]+$')
    echo "  Darwin's own: J = ${darwin_J}"

    echo "$(date) -- [5/6] patching ecco_cost with Darwin's J (pkg/dic gradient untouched)"
    python3 "$HYBRID/scripts/patch_ecco_cost_dic.py" "${myiter}" "${darwin_J}"

    echo "$(date) -- [6/6] optim.x (M1QN3 step)"
    prevlines=$(wc -l < m1qn3_output.txt 2>/dev/null || echo 0)
    ./optim.x > "opt_${it}.txt" 2>&1
    status=$?
    if (( status != 0 )); then
        echo "FATAL: optim.x exited with status ${status} on cycle ${myiter} -- see opt_${it}.txt"
        exit 1
    fi
    echo "$(date) -- cycle ${myiter} done (pkg/dic ${dic_fc_line} | Darwin J=${darwin_J})"

    if tail -n +$((prevlines + 1)) m1qn3_output.txt 2>/dev/null | grep -q "m1qn3: output mode"; then
        echo "m1qn3 has finished (output mode reached)"
        break
    fi
    ((myiter++))
done
echo "$(date) -- hybrid loop finished, last cycle=${myiter}"
