#!/bin/bash
# job_step1_stitch.sh — run stitch_ncgrid.py on a COMPUTE node. The stitched
# 3D fields (HFacC/HFacS/HFacW at Nr=90 x 1008 x 1170) are ~2.5-3GB combined,
# comparable in scale to gen_bathy.py's earlier login-node OOM (session cgroup
# ceiling ~1.5-2GB) -- running on a compute node per established STEP1 policy.
# Repad update (2026-08-08): grid is now 1170x1008 tiled at 39x42/30x24 (720
# ranks) -- see grid/mncs/mnc_0001..0720 (symlinked from gen_ncgrid/run/mnc_*).
#PBS -l select=1:ncpus=4:model=rom_ait
#PBS -q normal
#PBS -l walltime=1:00:00
#PBS -N kerstitch
#PBS -e /nobackup/dcarrol2/downscaling/kerguelen/stitch_err.txt
#PBS -o /nobackup/dcarrol2/downscaling/kerguelen/stitch_out.txt

D=/nobackup/dcarrol2/downscaling
K="$D/kerguelen"
UTILS="/nobackup/dcarrol2/v05_V5r1/ecco_darwin/regions/downscaling/utils"

export MAMBA_ROOT_PREFIX="$D/micromamba_root"
eval "$("$D/micromamba/bin/micromamba" shell hook --shell bash)"
micromamba activate downscaling

python3 "$UTILS/stitch_ncgrid.py" -d "$K/grid" -n kerguelen -z 90 -s 1170 1008 -p 39 42

echo "STITCH_DONE"
ls -la "$K/grid"/*_ncgrid.nc
