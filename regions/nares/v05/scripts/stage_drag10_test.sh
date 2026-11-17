#!/bin/bash
# Stage run_arch_D10t (2026-11-15, Dustin fix (a)): 10-day cold-start test (Jan 1-10) of build_drag10 (ice-ocean drag on the
# top-10-m mean ocean velocity). Namelists = run_arch_Y69m_sw exactly, except nTimeSteps 14400 and one added daily list
# (UVEL, VVEL levels 1-7 = top 10.8 m) to compare the near-surface flow. Compare mouth/Robeson ice speed days 1-6 with Y69m_sw.
# Ask Dustin before qsub:
# qsub -q devel -l select=5:ncpus=40:mpiprocs=39:model=sky_ele -l walltime=2:00:00 -W group_list=s1353 -N nares_d10t -v RUN=<dir> prep/job_run_arch_long.sh
set -eu
R=/nobackup/dcarrol2/downscaling/nares/regional_run; Y=$R/run_arch_Y69m_sw; T=$R/run_arch_D10t
[ -e $T ] && { echo "FATAL: $T exists"; exit 1; }
[ -x $R/build_drag10/mitgcmuv ] || { echo "FATAL: no build_drag10"; exit 1; }
mkdir $T; cd $T
for f in $Y/*; do b=$(basename $f); case $b in diags|STD*|job_run.out|*.ckpt*|*.bak_*|core*|crash*) continue;; esac
  if [ -L $f ]; then ln -s $(readlink $f) $b; elif [[ $b == pickup*.0000604872* ]]; then cp -p $f $b; elif [[ $b == pickup* ]]; then continue; else cp -p $f $b; fi; done
ln -sfn $R/build_drag10/mitgcmuv mitgcmuv
sed -i 's/^ nTimeSteps=349920,/ nTimeSteps=14400,/' data
python3 - <<'PY'
s = open('data.diagnostics').read()
add = """# 2026-11-15 D10t (Dustin fix (a)): daily near-surface velocity, levels 1-7 (0-10.8 m)
  frequency(12) = 86400.,
   fields(1:2,12) = 'UVEL    ','VVEL    ',
   levels(1:7,12) = 1.,2.,3.,4.,5.,6.,7.,
   filename(12) = 'diags/daily/uv10',
"""
i = s.index('&DIAG_STATIS_PARMS'); j = s.rindex(' &', 0, i)
assert 'filename(12)' not in s
open('data.diagnostics', 'w').write(s[:j] + add + s[j:])
PY
for d in $(cd $Y/diags && ls -d */); do mkdir -p diags/$d; done
echo "-- regular files differing from Y69m_sw:"; for f in *; do [ -f $f ] && [ ! -L $f ] && { cmp -s $f $Y/$f || echo "  $f"; }; done
diff $Y/data data || true; diff $Y/data.diagnostics data.diagnostics || true
grep -E "^ *nIter0" data; readlink mitgcmuv; echo STAGED $T
