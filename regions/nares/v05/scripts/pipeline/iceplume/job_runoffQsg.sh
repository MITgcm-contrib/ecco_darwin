#!/bin/bash
# Expand Humboldt daily discharge into yearly iceplume runoffQsg files (baseline).
# Runs on a compute node (no heavy work on pfe front ends).
#PBS -l select=1:ncpus=16:model=sky_ele
#PBS -q devel
#PBS -l walltime=0:30:00
#PBS -N nares_runoffQsg
#PBS -j oe
#PBS -o /nobackup/dcarrol2/downscaling/nares/forcings/iceplume/job_runoffQsg.out
N=/nobackup/dcarrol2/downscaling/nares
/nobackup/dcarrol2/downscaling/tools/cmvenv/bin/python $N/iceplume/gen_runoffQsg_files.py \
  $N/iceplume/Humboldt_Qsg_daily.npz $N/forcings/iceplume baseline
ls -la $N/forcings/iceplume
