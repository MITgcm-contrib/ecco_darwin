#!/bin/bash
# Rotate Nares OBCS velocities (east/north -> grid axes) into forcings/OBCS_rot; prints E/W transports.
#PBS -l select=1:ncpus=28:model=bro_ele
#PBS -q devel
#PBS -l walltime=0:30:00
#PBS -N nares_rot
#PBS -j oe
#PBS -o /nobackup/dcarrol2/downscaling/nares/prep/job_rotate_obcs.out
set -u
D=/nobackup/dcarrol2/downscaling
export MAMBA_ROOT_PREFIX="$D/micromamba_root"
eval "$("$D/micromamba/bin/micromamba" shell hook --shell bash)"
micromamba activate downscaling
python3 $D/nares/prep/rotate_obcs_velocity.py $D/nares
