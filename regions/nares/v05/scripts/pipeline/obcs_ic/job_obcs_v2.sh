#!/bin/bash
# Regenerate the Nares OBCS chain with the fixed gen_obcs_fast box-path masks (2026-10-09) in a SHADOW config
# (nares_obcsfix: links to everything in nares/ except forcings/), so running jobs that read nares/forcings/OBCS*
# are untouched. gen_obcs_fast -> rotate (child grid) -> budget -> offset Q = 0.85 - southward(OBCS_rot) -> OBCS_rot_tf.
#PBS -l select=1:ncpus=28:model=bro_ele
#PBS -q normal
#PBS -l walltime=2:00:00
#PBS -N nares_obcs2
#PBS -j oe
#PBS -o /nobackup/dcarrol2/downscaling/nares/prep/job_obcs_v2.out
set -u
D=/nobackup/dcarrol2/downscaling; N=$D/nares; X=$D/nares_obcsfix
export MAMBA_ROOT_PREFIX="$D/micromamba_root"
eval "$("$D/micromamba/bin/micromamba" shell hook --shell bash)"
micromamba activate downscaling
die() { echo "FATAL: $1"; exit 1; }
mkdir -p $X/forcings
for f in $N/*; do b=$(basename $f); [ "$b" = forcings ] || ln -sfn $f $X/$b; done
ln -sfn $N/forcings/pickups $X/forcings/pickups
cd $N
echo "===== gen_obcs_fast (fixed masks) $(date +%H:%M)"
python3 prep/gen_obcs_fast.py -d $X -n Nares -bnd EWNS -seaice -p 28 -v -disc off || die gen_obcs
python3 - <<P || die zerocheck
import numpy as np, xarray as xr
g = xr.open_dataset('$N/grid/Nares_ncgrid.nc'); hW = g['HFacW'].values; hS = g['HFacS'].values
NR, NY, NX = hW.shape[0], hW.shape[1], hS.shape[2]
u = np.fromfile('$X/forcings/OBCS/UVEL_west.bin', '>f4', count=NR*NY).reshape(NR, NY)
v = np.fromfile('$X/forcings/OBCS/VVEL_north.bin', '>f4', count=NR*NX).reshape(NR, NX)
w = hW[:, :, 1] > 0; n = hS[:, -2, :] > 0
print('west U: wet', w.sum(), 'zero', (u[w] == 0).sum(), '| north V: wet', n.sum(), 'zero', (v[n] == 0).sum())
assert (u[w] == 0).mean() < 0.01 and (v[n] == 0).mean() < 0.01, 'still zeros on wet interior faces'
P
echo "===== rotate $(date +%H:%M)"
python3 prep/rotate_obcs_velocity.py $X || die rotate
echo "===== budget OBCS_rot"
python3 prep/obcs_boundary_budget.py $X OBCS_rot | tee $X/budget_rot.txt || die budget
SOUTH=$(grep "implied southward" $X/budget_rot.txt | grep -oE '[-+][0-9.]+' | head -1)
Q=$(python3 -c "q=0.85-float('$SOUTH'); assert -0.5<q<1.5, q; print(f'{q:.3f}')") || die "offset out of range (south=$SOUTH)"
echo "southward(OBCS_rot) = $SOUTH Sv -> offset Q = $Q Sv"
python3 prep/nares_throughflow_offset.py $X $Q OBCS_rot OBCS_rot_tf || die offset
python3 prep/obcs_boundary_budget.py $X OBCS_rot_tf
echo "===== ALL DONE $(date +%H:%M)"
