#!/bin/bash
# Build the ecco_darwin downscaling config layout for Nares with SYMLINKS only (safe on pfe):
#   nares/parent/outputs/grid     -> LLC270 parent grid (same as kerguelen/pan_AO)
#   nares/parent/outputs/pickups  -> shared parent run (pickup*.0000604872 = 2015-01-01)
#   nares/parent/inputs/<bnd>_BC_mask.bin -> the Nares dv masks of the shared extraction
#   nares/Nares_ncgrid.nc         -> grid/Nares_ncgrid.nc
# parent/outputs/OBCS is created by job_nares_arch_inputs.sh (stable record snapshot).
set -eu
N=/nobackup/dcarrol2/downscaling/nares
PR=/nobackupp19/dcarrol2/downscaling/parent_V5r1_westAO_nares/run
mkdir -p $N/parent/outputs $N/parent/inputs
ln -sfn /nobackup/dcarrol2/downscaling/kerguelen/parent/outputs/grid $N/parent/outputs/grid
ln -sfn $PR $N/parent/outputs/pickups
for b in east west north south; do
  m=$(ls $PR/dv/nares_${b}_BC_mask.bin $PR/nares_${b}_BC_mask.bin 2>/dev/null | head -1)
  [ -n "$m" ] || { echo "missing mask for $b"; exit 1; }
  ln -sfn $m $N/parent/inputs/${b}_BC_mask.bin
done
ln -sfn $N/grid/Nares_ncgrid.nc $N/Nares_ncgrid.nc
ls -la $N/parent/outputs $N/parent/inputs $N/Nares_ncgrid.nc
