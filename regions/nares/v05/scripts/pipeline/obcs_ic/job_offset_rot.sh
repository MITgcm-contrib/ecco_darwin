#!/bin/bash
# Throughflow offset on the rotated OBCS (OBCS_rot -> OBCS_rot_tf), then the by-end budget of the result.
# Q = 0.85 - 0.57 = 0.28 Sv (0.57 = balance-split throughflow implied by OBCS_rot, job 25244231).
#PBS -l select=1:ncpus=28:model=bro_ele
#PBS -q devel
#PBS -l walltime=0:30:00
#PBS -N nares_off
#PBS -j oe
#PBS -o /nobackup/dcarrol2/downscaling/nares/prep/job_offset_rot.out
set -u
D=/nobackup/dcarrol2/downscaling
export MAMBA_ROOT_PREFIX="$D/micromamba_root"
eval "$("$D/micromamba/bin/micromamba" shell hook --shell bash)"
micromamba activate downscaling
python3 $D/nares/prep/nares_throughflow_offset.py $D/nares 0.28 OBCS_rot OBCS_rot_tf || exit 1
python3 $D/nares/prep/obcs_boundary_budget.py $D/nares OBCS_rot_tf
