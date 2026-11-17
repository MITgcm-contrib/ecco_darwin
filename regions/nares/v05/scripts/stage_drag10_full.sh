#!/bin/bash
# Stage run_arch_D10 (2026-11-16): full Jan 1-Aug 31 2015 run of build_drag10 (ice-ocean drag on the top-10-m mean ocean
# velocity, SEAICE_OCNVEL_AVG). Namelists = run_arch_Y69m_sw exactly (nTimeSteps 349920 unchanged) plus the D10t daily list 12
# (UVEL, VVEL levels 1-7 = top 10.8 m, diags/daily/uv10). Physics change vs Y69m_sw: the drag velocity only (stop rule: one change per full run).
# Ask Dustin before qsub (same resources as Y69m_sw, job 25386502, which took 27.4 h):
# qsub -q long -l select=5:ncpus=40:mpiprocs=39:model=sky_ele -l walltime=40:00:00 -W group_list=s1353 -N nares_d10 -v RUN=<dir> prep/job_run_arch_long.sh
set -eu
R=/nobackup/dcarrol2/downscaling/nares/regional_run; Y=$R/run_arch_Y69m_sw; D=$R/run_arch_D10t; T=$R/run_arch_D10
[ -e $T ] && { echo "FATAL: $T exists"; exit 1; }
[ -x $R/build_drag10/mitgcmuv ] || { echo "FATAL: no build_drag10"; exit 1; }
mkdir $T; cd $T
for f in $Y/*; do b=$(basename $f); case $b in diags|STD*|job_run.out|*.ckpt*|*.bak_*|core*|crash*) continue;; esac
  if [ -L $f ]; then ln -s $(readlink $f) $b; elif [[ $b == pickup*.0000604872* ]]; then cp -p $f $b; elif [[ $b == pickup* ]]; then continue; else cp -p $f $b; fi; done
ln -sfn $R/build_drag10/mitgcmuv mitgcmuv
cp -p $D/data.diagnostics data.diagnostics   # = Y69m_sw + list 12 uv10 (verified by the diff below)
for d in $(cd $Y/diags && ls -d */); do mkdir -p diags/$d; done; mkdir -p diags/daily
echo "-- regular files differing from Y69m_sw (expect data.diagnostics only):"; for f in *; do [ -f $f ] && [ ! -L $f ] && { cmp -s $f $Y/$f || echo "  $f"; }; done
diff $Y/data data || true; diff $Y/data.diagnostics data.diagnostics || true
grep -E "^ *(nIter0|nTimeSteps|pChkptFreq)=" data; readlink mitgcmuv
echo "== dangling links"; find . -maxdepth 1 -xtype l | head; echo "(none above = ok)"
echo STAGED $T
