#!/bin/bash
# Stage the Nares BGC Jan 1 - Sep 1 2015 run in 3 segments (2026-10-30, Dustin "yes" to staging the Jan-Aug BGC run).
# Set-up = run_bgc_test2d (passed: = run_arch_Y16 physics + pan_AO/west_AO ecosystem, build_bgc, oasim_dTsolz 60) with:
#   segments at calendar-month pickups (Nares writes pickups at month starts):
#     s1 Jan 1 -> Apr 1  nIter0 0      NT 129600 (90 d, ~84 h at the measured 56 wall min per day)
#     s2 Apr 1 -> Jul 1  nIter0 129600 NT 131040 (91 d, ~85 h)
#     s3 Jul 1 -> Sep 1  nIter0 260640 NT  89280 (62 d, ~58 h)
#   (2 segments would need ~113 h each, too close to the 120-h long-queue limit.)
#   s2/s3 link pickup, pickup_seaice, pickup_ggl90, pickup_ptracers, pickup_darwin of the previous segment
#   (dangling until it ends; the chain uses afterok).
#   Diagnostics: Y16 physics lists without the hourly ones (Y16 has them) + daily surface BGC, daily CO2,
#   5-daily 3-D light/PP/Chl, monthly 3-D for all 33 tracers (west_AO file names).
# Ask Dustin before qsub.
set -eu
N=/nobackup/dcarrol2/downscaling/nares; R=$N/regional_run; T=$R/run_bgc_test2d; W=/nobackup/dcarrol2/downscaling/west_AO/regional_run/run_prod4_2015H1
SEGS=("s1 0 129600" "s2 129600 131040" "s3 260640 89280")
for s in "${SEGS[@]}"; do set -- $s; [ -e $R/run_bgc15_$1 ] && { echo "FATAL: $R/run_bgc15_$1 exists"; exit 1; }; done
prev=""
for s in "${SEGS[@]}"; do
  set -- $s; S=$1; IT0=$2; NT=$3; RUN=$R/run_bgc15_$S; IT=$(printf %010d $IT0)
  mkdir $RUN; cd $RUN
  for f in $T/*; do b=$(basename $f); case $b in diags|STD*|job_run.out|crash*|*.bak_*|*.0000002880.*|*.ckpt*|pickup*.ckpt*) continue;; esac
    if [ -L $f ]; then ln -s $(readlink $f) $b; else cp -p $f $b; fi; done
  ln -sfn $R/build_bgc/mitgcmuv mitgcmuv
  sed -i -E "s/^( *nTimeSteps=)[0-9]+,/\1$NT,/; s/^( *nIter0=)[0-9]+,/\1$IT0,/" data
  if [ -n "$prev" ]; then
    for p in pickup pickup_seaice pickup_ggl90 pickup_ptracers pickup_darwin; do
      for x in data meta; do ln -sfn $prev/$p.$IT.$x $p.$IT.$x; done; done
  fi
  mkdir -p diags/3hourly diags/6hourly diags/daily diags/5daily diags/monthly
  python3 - $W/data.diagnostics <<'PY'
import re, sys
west = open(sys.argv[1]).read()
# west_AO monthly 3-D tracer names: fields(1,k) = 'TRACnn' ... filename(k) = 'diags/monthly/<name>'
pairs = re.findall(r"fields\(1,(\d+)\) = '(TRAC\d\d) *',\s*\n\s*filename\(\1\) = 'diags/monthly/([A-Za-z0-9_]+)'", west)
names = {t: n for _, t, n in pairs}
assert len(names) == 33 and all(f'TRAC{i:02d}' in names for i in range(1, 34)), sorted(names)
L = ["# Nares BGC Jan-Aug 2015 (2026-10-30). Physics = run_arch_Y16 lists minus hourly. BGC names checked in code_bgc.",
     " &DIAGNOSTICS_LIST"]
def add(k, freq, fields, fname, lev1=False):
    L.append(f"  frequency({k}) = {freq},")
    L.append(f"   fields(1:{len(fields)},{k}) = " + ",".join(f"'{f:8s}'" for f in fields) + ",")
    if lev1: L.append(f"   levels(1,{k}) = 1.,")
    L.append(f"   filename({k}) = '{fname}',")
add(1, '10800.', ['SIuice', 'SIvice', 'SIarea', 'SIheff'], 'diags/3hourly/SIdyn')
add(2, '86400.', ['SIshear', 'SIdelta', 'SItensil', 'SIsig1', 'SIsig2', 'SIhsnow'], 'diags/daily/SIdeform')
add(3, '86400.', ['ETAN', 'MXLDEPTH', 'oceQnet'], 'diags/daily/surf')
add(4, '86400.', ['SHIfwFlx', 'SHIhtFlx', 'SHIuStar', 'SHIgammT'], 'diags/daily/SHI')
add(5, '2592000.', ['THETA', 'SALT', 'UVEL', 'VVEL'], 'diags/monthly/state')
add(6, '86400.', ['THETA', 'SALT'], 'diags/daily/surfTS', lev1=True)
add(7, '21600.', ['EXFuwind', 'EXFvwind', 'SItaux', 'SItauy', 'oceTAUX', 'oceTAUY', 'SIpress', 'SIzeta', 'SIuice', 'SIvice', 'SIheff', 'SIarea'], 'diags/6hourly/SIforce')
add(8, '86400.', ['Chl', 'PP', 'plankC', 'TRAC02', 'PARF', 'pH', 'pCO2'], 'diags/daily/bgcsurf', lev1=True)
add(9, '86400.', ['fluxCO2', 'gDICsurf'], 'diags/daily/co2')
add(10, '432000.', ['PARF', 'RTtrans', 'PP', 'Chl'], 'diags/5daily/bgc3d')
for i in range(1, 34):
    t = f'TRAC{i:02d}'; add(10 + i, '2592000.', [t], f'diags/monthly/{names[t]}')
L += [" &", " &DIAG_STATIS_PARMS", " /", ""]
open('data.diagnostics', 'w').write("\n".join(L))
print('diagnostics lists:', 10 + 33)
PY
  echo "== $RUN: $(grep -E '^ *(nIter0|nTimeSteps)=' data | tr -d ' \n')  pickups linked: $(ls pickup*.$IT.data 2>/dev/null | wc -l)  BGC IC links: $(ls pickup_pTr*.0000604872.data | wc -l)  OB tracer links: $(ls TRAC*_*.bin rDOC_*.bin CDOM_*.bin | wc -l)"
  prev=$RUN
done
cd $R/run_bgc15_s1; echo "-- s1 data/data.pkg/data.oasim vs run_bgc_test2d:"; for f in data data.pkg data.oasim data.darwin data.ptracers data.obcs data.seaice; do diff -q $T/$f $f >/dev/null && echo "$f same" || { echo "$f differs:"; diff $T/$f $f | head -4; }; done
echo STAGED
