#!/bin/bash
# socal 2-km nest run: 24 ranks on one sky_ele node (mitgcmuv is Intel-only, Kerguelen Bug I).
# qsub -l walltime=1:00:00 -v NSTEPS=720 job_run.sh   (1-day test at deltaT 120 s)
#PBS -l select=1:ncpus=40:mpiprocs=24:model=sky_ele
#PBS -q normal
#PBS -l walltime=1:00:00
#PBS -N socal_run
#PBS -j oe
module purge
module load comp-intel/2023.2.1 mpi-hpe/mpt.2.30 hdf4/4.2.12 hdf5/1.8.18_mpt netcdf/4.4.1.1_mpt
cd ${RUNDIR:-/nobackup/dcarrol2/downscaling/socal/regional_run/run_test1d}
mpiexec -np 24 ./mitgcmuv
grep -h "Execution ended Normally" STDOUT.0000 && echo RUN_OK
