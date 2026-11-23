#!/bin/bash
# Nares BGC scaling test (2026-11-16): job_run_arch_long.sh with the rank count from the environment (NP) instead of 192.
# qsub -q devel -l select=<NODES>:ncpus=40:mpiprocs=<PPN>:model=sky_ele -l walltime=2:00:00 -W group_list=s1353 -N nares_sc<TAG> -v RUN=<dir>,NP=<NRANK> <this>
# Ask Dustin before qsub.
#PBS -N nares_scaling
#PBS -j oe
#PBS -W group_list=s1353
module purge
module load comp-intel/2023.2.1 mpi-hpe/mpt.2.30 hdf4/4.2.12 hdf5/1.8.18_mpt netcdf/4.4.1.1_mpt
cd ${RUN:?set RUN}
NP=${NP:?set NP}
exec > >(tee -a job_run.out) 2>&1
echo "scaling run: NP=$NP nodes=$(sort -u $PBS_NODEFILE | wc -l) slots=$(wc -l < $PBS_NODEFILE) start $(date)"
grep -oE 'nPx *= *[0-9]+' $(dirname $(readlink -f mitgcmuv))/SIZE.h
it0=$(grep -oP '^ *nIter0=\K\d+' data); nt=$(grep -oP '^ *nTimeSteps=\K\d+' data)
need=$(( (it0 + nt) * 60 / 3600 + 1 ))
for f in SIheff_east.bin SIheff_west.bin SIheff_north.bin SIheff_south.bin; do
  s=$(stat -L -c %s $f) || { echo "missing $f"; exit 1; }
  case $f in *east*|*west*) npt=360;; *) npt=1120;; esac
  rec=$(( s / (4 * npt) )); echo "$f: $rec hourly records (need >= $need)"; [ $rec -ge $need ] || { echo "PREFLIGHT FAIL"; exit 1; }
done
mpiexec -np $NP ./mitgcmuv
rc=$?
if [ ! -s STDOUT.0000 ]; then echo "no STDOUT after launch (rc=$rc): retrying mpiexec in 30 s"; sleep 30; mpiexec -np $NP ./mitgcmuv; rc=$?; fi
if grep -q "Execution ended Normally" STDOUT.0000; then ok=0; else ok=1; fi
grep -h "ABNORMAL END\|STOP in" STDERR.0* 2>/dev/null | sort | uniq -c | head
echo "RUN_DONE mpiexec_rc=$rc normal_end=$((1-ok)) end $(date)"
exit $ok
