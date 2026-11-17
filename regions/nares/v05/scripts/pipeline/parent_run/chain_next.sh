#!/bin/bash
# chain_next.sh -- after a 48-h segment of the shared parent dv extraction: stop if the run
# finished; otherwise restart from the newest COMPLETE pickup set and resubmit. Refuses to
# resubmit if no newer pickup exists (no progress) -> avoids an endless chain.
set -u
D=/nobackupp19/dcarrol2/downscaling/parent_V5r1_westAO_nares; R=$D/run
cd $R
if grep -q "Execution ended Normally" STDOUT.0000 2>/dev/null; then echo "CHAIN_DONE: run complete"; exit 0; fi
cur=$(grep -E "^ *nIter0=" data | sed -E "s/[^0-9]//g")
best=""
for f in $(ls pickup.[0-9]*.data 2>/dev/null | sort -r); do
  it=$(echo $f | sed -E "s/pickup\.([0-9]+)\.data/\1/")
  ok=1; for c in pickup pickup_darwin pickup_ggl90 pickup_ptracers pickup_seaice; do
    [ -s $c.$it.data ] && [ -s $c.$it.meta ] || { ok=0; continue; }
    # reject a pickup truncated by the walltime timeout: same byte size as the starting set
    ref=$(ls $c.[0-9]*.data 2>/dev/null | grep -v "$c.$it.data" | head -1)
    [ -n "$ref" ] && [ $(stat -L -c %s $c.$it.data) -ne $(stat -L -c %s $ref) ] && ok=0; done
  [ $ok -eq 1 ] && { best=$it; break; }
done
nb=$((10#${best:-0})); nc=$((10#$cur))
if [ -z "$best" ] || [ $nb -le $nc ]; then echo "CHAIN_STOP: no newer complete pickup than nIter0=$cur"; exit 1; fi
cp -p data data.seg_nIter0_$cur
sed -i -E "s/^ *nIter0=[0-9]+,/ nIter0=$nb,/" data
mkdir -p logs_seg_$cur; mv STDOUT.* STDERR.* logs_seg_$cur/ 2>/dev/null
echo "CHAIN_NEXT: nIter0 $cur -> $nb"; /PBS/bin/qsub $D/job_ECCO_darwin
