#!/bin/bash
# Nares ice-arch test run (physics only), 192 ranks. Stage first: prep/stage_run_arch.sh <nTimeSteps>.
# Ask Dustin before qsub.
#PBS -l select=7:ncpus=28:mpiprocs=28:model=bro_ele
#PBS -q normal
#PBS -l walltime=1:00:00
#PBS -N nares_arch
#PBS -j oe
#PBS -o /nobackup/dcarrol2/downscaling/nares/regional_run/run_arch/job_run_arch.out
#PBS -W group_list=s1353
# 2026-10-21 (Dustin via Session watch): charge to s1353, not g26113
module purge
module load comp-intel/2023.2.1 mpi-hpe/mpt.2.30 hdf4/4.2.12 hdf5/1.8.18_mpt netcdf/4.4.1.1_mpt
cd ${RUN:-/nobackup/dcarrol2/downscaling/nares/regional_run/run_arch}
mpiexec -np 192 ./mitgcmuv
rc=$?
# MPT sometimes dies at launch ("Assertion failed at mmap.c:179", no STDOUT written): retry once
if [ ! -s STDOUT.0000 ]; then echo "no STDOUT after launch (rc=$rc): retrying mpiexec in 30 s"; sleep 30; mpiexec -np 192 ./mitgcmuv; rc=$?; fi
# MPT can return 0 after another rank aborts; trust MITgcm's own end-of-run line instead
if grep -q "Execution ended Normally" STDOUT.0000; then ok=0; else ok=1; fi
grep -h "ABNORMAL END\|STOP in" STDERR.0* 2>/dev/null | sort | uniq -c | head
echo "RUN_DONE mpiexec_rc=$rc normal_end=$((1-ok))"
exit $ok
