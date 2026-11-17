#!/bin/bash
# Nares OBCS volume budget by strait end, for the original (east/north) and rotated OBCS.
#PBS -l select=1:ncpus=28:model=bro_ele
#PBS -q devel
#PBS -l walltime=0:30:00
#PBS -N nares_bud
#PBS -j oe
#PBS -o /nobackup/dcarrol2/downscaling/nares/prep/job_obcs_budget.out
set -u
D=/nobackup/dcarrol2/downscaling
export MAMBA_ROOT_PREFIX="$D/micromamba_root"
eval "$("$D/micromamba/bin/micromamba" shell hook --shell bash)"
micromamba activate downscaling
for s in OBCS OBCS_tf OBCS_rot; do python3 $D/nares/prep/obcs_boundary_budget.py $D/nares $s; done
