#!/bin/bash
# Stage run_arch_Y69m_sw (2026-11-09, Dustin "do 1: redo summer tuning on build_swfix; Y69m + ocean SW heating"): run_arch_Y69m namelists and inputs exactly,
# binary build_swfix (ocean SW heating restored for physics-only runs). Ask Dustin before qsub:
# qsub -q long -l select=5:ncpus=40:mpiprocs=39:model=sky_ele -l walltime=40:00:00 -W group_list=s1353 -N nares_y69msw -v RUN=<dir> prep/job_run_arch_long.sh
set -eu
R=/nobackup/dcarrol2/downscaling/nares/regional_run; Y=$R/run_arch_Y69m; T=$R/run_arch_Y69m_sw
[ -e $T ] && { echo "FATAL: $T exists"; exit 1; }
[ -x $R/build_swfix/mitgcmuv ] || { echo "FATAL: no build_swfix"; exit 1; }
mkdir $T; cd $T
for f in $Y/*; do b=$(basename $f); case $b in diags|STD*|job_run.out|*.ckpt*|*.bak_*|core*|crash*) continue;; esac
  if [ -L $f ]; then ln -s $(readlink $f) $b; elif [[ $b == pickup* ]]; then continue; else cp -p $f $b; fi; done
ln -sfn $R/build_swfix/mitgcmuv mitgcmuv
for d in $(cd $Y/diags && ls -d */); do mkdir -p diags/$d; done
echo "-- regular files differing from Y69m:"; for f in *; do [ -f $f ] && [ ! -L $f ] && { cmp -s $f $Y/$f || echo "  $f"; }; done
echo "-- links differing from Y69m:"; for f in *; do [ -L $f ] && [ "$(readlink $f)" != "$(readlink $Y/$f 2>/dev/null)" ] && echo "  $f -> $(readlink $f)"; done
grep -E "^ *(nIter0|nTimeSteps)=" data | tr -d " "; grep -i "SEAICE_strength\|presPow\|eccen\|AdvScheme\|areaLoss" data.seaice | grep -v "^ *#"
echo "pickup IC links: $(ls pickup*0000604872* 2>/dev/null | wc -l); diag dirs: $(ls diags | tr "\n" " ")"; echo STAGED $T
