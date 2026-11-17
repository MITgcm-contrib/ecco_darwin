#!/bin/bash
# Continue a Nares arch run from its Feb 1 pickup (iter 44640) for a long period (2026-10-23, arch evaluation idea 1).
# Same set-up as the parent (stage_arch_tc4.sh on input_baseline) with: nIter0 44640, nTimeSteps NT, pickups linked to the
# parent run, and every OBCS2_rot_tf link replaced by OBDIR (Jan-Sep set). Default: ICT -> run_arch_ICTs to Aug 31
# (NT = 305280 steps of 60 s = 212 d, Feb 1 00Z -> Sep 1 00Z). Ask Dustin before qsub.
# Usage: stage_arch_restart_long.sh [SRC=ICT] [NAME=ICTs] [NT=305280] [OBDIR=OBCS2_rot_tf_JanAug] [IT0=44640] [BUILD=build]
# 2026-10-24: IT0 (restart iteration, pickups taken from SRC) and BUILD (binary dir) added for run_arch_ICTs2.
set -eu
P=${1:-ICT}; X=${2:-ICTs}; NT=${3:-305280}; OBDIR=${4:-OBCS2_rot_tf_JanAug}; IT0=${5:-44640}; BUILD=${6:-build}
N=/nobackup/dcarrol2/downscaling/nares; R=$N/regional_run; SRC=$R/run_arch_$P; RUN=$R/run_arch_$X; IT=$(printf %010d $IT0)
[ -e $RUN ] && { echo "$RUN exists"; exit 1; }
[ -s $SRC/pickup.$IT.data ] && [ -s $SRC/pickup_seaice.$IT.data ] || { echo "no pickup/pickup_seaice $IT in $SRC"; exit 1; }
NAME=$X bash $N/prep/stage_arch_tc4.sh $NT >/dev/null
cd $RUN
python3 - $IT0 <<'PY'
import re, sys
s = open('data').read()
s, n = re.subn(r'^( *nIter0=)\d+,', r'\g<1>%s,' % sys.argv[1], s, count=1, flags=re.M); assert n == 1
open('data', 'w').write(s)
PY
ln -sfn $R/$BUILD/mitgcmuv mitgcmuv
for f in $(cd $SRC && ls pickup*.$IT.data pickup*.$IT.meta); do ln -sfn $SRC/$f $f; done
n=0
for f in *.bin; do
  t=$(readlink $f || true)
  case "$t" in */forcings/OBCS2_rot_tf/*) ln -sfn $N/forcings/$OBDIR/$f $f; n=$((n+1));; esac
done
echo "== $RUN: relinked $n OBCS files -> forcings/$OBDIR"; [ $n -eq 40 ] || { echo "expected 40 OBCS links"; exit 1; }
grep -n "nIter0=\|nTimeSteps=" data
echo "-- diff data vs $P:"; diff $SRC/data data || true
for f in data.obcs data.seaice data.diagnostics data.exf data.pkg; do cmp -s $f $SRC/$f && echo "$f = $P" || { echo "$f differs from $P:"; diff $SRC/$f $f | head -6; }; done
echo "pickups linked: $(ls pickup*.$IT.data | wc -l)"; readlink Nares_bathymetry.bin mitgcmuv
