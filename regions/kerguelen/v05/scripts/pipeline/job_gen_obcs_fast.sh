#!/bin/bash
# job_gen_obcs_fast.sh — run the optimized gen_obcs_fast.py (STEP3) on a
# COMPUTE node. Benchmarked at ~115s per 90-level (variable,boundary) task
# (job 24943932); 164 tasks total, parallelized across (boundary,variable)
# pairs -- see DECISIONS.md for the full extrapolation. Requesting a full
# sky_ele node's worth of parallelism.
#PBS -l select=1:ncpus=16:model=rom_ait
#PBS -q normal
#PBS -l walltime=6:00:00
#PBS -N kerobcsfast
#PBS -e /nobackup/dcarrol2/downscaling/kerguelen/obcs_fast_err.txt
#PBS -o /nobackup/dcarrol2/downscaling/kerguelen/obcs_fast_out.txt

D=/nobackup/dcarrol2/downscaling
K="$D/kerguelen"

export MAMBA_ROOT_PREFIX="$D/micromamba_root"
eval "$("$D/micromamba/bin/micromamba" shell hook --shell bash)"
micromamba activate downscaling

cd "$K"
python3 gen_obcs_fast.py -d "$K" -n kerguelen -bnd EWNS -seaice -bgc -p 14 -v

echo "GEN_OBCS_FAST_DONE"
ls -la "$K/forcings/OBCS/" | wc -l
ls -la "$K/forcings/OBCS/" | head -20
