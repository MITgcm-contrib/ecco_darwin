#!/bin/bash
# Nares long arch run (2026-10-23): job_run_arch.sh plus a pre-flight check that every OBCS file holds enough hourly records
# for the run (nIter0..nIter0+nTimeSteps at deltaT 60 s from 2015-01-01 00Z). Ask Dustin before qsub.
# qsub -q long -l select=5:ncpus=40:mpiprocs=39:model=sky_ele -l walltime=40:00:00 -W group_list=s1353 -v RUN=<dir> <this>
#PBS -N nares_long
#PBS -j oe
#PBS -W group_list=s1353
module purge
module load comp-intel/2023.2.1 mpi-hpe/mpt.2.30 hdf4/4.2.12 hdf5/1.8.18_mpt netcdf/4.4.1.1_mpt
cd ${RUN:?set RUN}
exec > >(tee -a job_run.out) 2>&1
it0=$(grep -oP '^ *nIter0=\K\d+' data); nt=$(grep -oP '^ *nTimeSteps=\K\d+' data)
need=$(( (it0 + nt) * 60 / 3600 + 1 ))
for f in SIheff_east.bin SIheff_west.bin SIheff_north.bin SIheff_south.bin; do
  s=$(stat -L -c %s $f) || { echo "missing $f"; exit 1; }
  case $f in *east*|*west*) npt=360;; *) npt=1120;; esac
  rec=$(( s / (4 * npt) )); echo "$f: $rec hourly records (need >= $need)"; [ $rec -ge $need ] || { echo "PREFLIGHT FAIL"; exit 1; }
done
mpiexec -np 192 ./mitgcmuv
rc=$?
if [ ! -s STDOUT.0000 ]; then echo "no STDOUT after launch (rc=$rc): retrying mpiexec in 30 s"; sleep 30; mpiexec -np 192 ./mitgcmuv; rc=$?; fi
if grep -q "Execution ended Normally" STDOUT.0000; then ok=0; else ok=1; fi
grep -h "ABNORMAL END\|STOP in" STDERR.0* 2>/dev/null | sort | uniq -c | head
echo "RUN_DONE mpiexec_rc=$rc normal_end=$((1-ok))"
exit $ok
