#!/bin/bash
# run_arch_TC4 (2026-10-15, after TC3; Dustin "Yes"): = run_arch_TC2 set-up (8 d Jan, 8 constituents, C2 tides, OBCS sponge,
# baseline HO 2.0 + Nares_grid18) on the ORIGINAL 10 m minimum depth with the min-depth FILL bridges walled off
# (prep/wall_fill_bridges_nares.py -> grid/Nares_bathymetry_walled.bin: 571 cells -> land, incl. the SE-corner pocket and the
# i494 fjord; no new wet cells, so the baseline 10 m pickups are used unchanged). The SE-corner pocket held east-OB points
# j = 12..30 (1-based), now land -> OB_Ieast zeroed there. Ask Dustin before qsub.
set -eu
N=/nobackup/dcarrol2/downscaling/nares; R=$N/regional_run; NAME=${NAME:-TC4}; RUN=$R/run_arch_$NAME
BW=$N/grid/Nares_bathymetry_walled.bin
[ -s $BW ] || { echo "no $BW"; exit 1; }
[ -e $RUN ] && { echo "$RUN exists"; exit 1; }
bash $N/prep/stage_arch_tides.sh ${1:-11520} $NAME OBCS_tides_2015Jan_c2 | grep -E "^==|nTimeSteps|baseline"
cd $RUN
ln -sfn $BW Nares_bathymetry.bin
cp $R/run_arch_TC2/data.obcs data.obcs        # TC2 = C2 tides snippet + sponge (OBCS_PARM03); same tide dir
python3 - <<'P'
import re
s = open('data.obcs').read()
m = re.search(r'(OB_Ieast *=)([\s\d,]*?)(?=\n *[A-Za-z_#&])', s)
vals = [int(x) for x in re.findall(r'\d+', m.group(2))]
assert len(vals) == 360, len(vals)
assert all(vals[j] == 1120 for j in range(11, 30)), vals[8:32]
for j in range(11, 30):
    vals[j] = 0
body = '\n' + '\n'.join('   ' + ' '.join(f'{v:5d},' for v in vals[k:k + 15]) for k in range(0, 360, 15))
s = s[:m.start(2)] + body + s[m.end(2):]
s = re.sub(r'\n[ \t]*OB_Ieast =', '\n# TC4 2026-10-15: j=12..30 zeroed (SE-corner pocket walled off, wall_fill_bridges_nares.py)\n OB_Ieast =', s, count=1)
open('data.obcs', 'w').write(s)
print('OB_Ieast: zeroed j=12..30; nonzero now', sum(v > 0 for v in vals))
P
grep -n "useOBCSsponge\|spongeThickness\|tidalPeriod\|TC4" data.obcs | head -5
readlink Nares_bathymetry.bin pickup_SALT.0000604872.data Nares_grid18.face001.bin mitgcmuv
