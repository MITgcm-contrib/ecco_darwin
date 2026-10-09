#!/bin/bash
# After the Pleiades GF ensemble: pull outputs (single rsync stream), then run
# global, hold-out-one-site and per-site Green's-function solves, then summarize.
set -e
cd ~/Documents/research/1-D
mkdir -p runs_pfe
rsync -a --include='*/' --include='mnc_out/***' --include='darwin_params.txt' --include='data*' \
      --include='gf_manifest.txt' --exclude='*' pfe:/nobackup/dcarrol2/1-D/runs/ runs_pfe/
O=gf_results
mkdir -p $O/global $O/cache
mkdir -p $O
python3 scripts/gf_solve.py runs_pfe/ctrl runs_pfe/gf obs $O/global/out > $O/global.log 2>&1
for s in HOT BATS HydroS PAPA PAP; do
  python3 scripts/gf_solve.py runs_pfe/ctrl runs_pfe/gf obs $O/holdout_$s/out --holdout $s > $O/holdout_$s.log 2>&1
  python3 scripts/gf_solve.py runs_pfe/ctrl runs_pfe/gf obs $O/only_$s/out --only $s > $O/only_$s.log 2>&1
done
python3 scripts/gf_summary.py $O > $O/summary.txt
cat $O/summary.txt
