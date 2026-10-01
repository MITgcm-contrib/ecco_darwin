#!/bin/csh
# STEP2 production run: diagnostics_vec-enabled global LLC270 Darwin BGC
# parent, restarted from the 2020-01-01 pickup (nIter0=736344), run 3
# years to 2023-01-01 (endtime=978307200.). Originally select=20:ncpus=40
# :model=sky_ele (np=767), switched to bro_ele 2026-07-31 to escape the
# sky_ele queue backlog (6751 nodes wanted vs 2196 total). Build's
# -axCORE-AVX2 optim path covers Broadwell fully, so the existing
# mitgcmuv binary runs unchanged -- no rebuild needed. 28 cores/node ->
# 28 nodes to keep >=767 slots for np=767.
#PBS -l select=28:ncpus=28:model=bro_ele
#PBS -q long
#PBS -l walltime=120:00:00
#PBS -N kerdv
#PBS -j oe
#PBS -m abe

module purge
module load comp-intel/2023.2.1 mpi-hpe/mpt.2.30 hdf4/4.2.12 hdf5/1.8.18_mpt netcdf/4.4.1.1_mpt
module list

umask 022
cd /nobackup/dcarrol2/downscaling/kerguelen/parent_run/run
limit stacksize unlimited
mpiexec -np 767 ./mitgcmuv
