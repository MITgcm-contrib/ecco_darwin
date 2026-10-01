#!/bin/bash
# job_gen_pickups.sh — run gen_pickups.py (STEP3) on a COMPUTE node.
# Builds a Delaunay triangulation over the full parent LLC270 grid
# (947,700 points) plus several ~900MB regional-grid arrays -- same
# memory class that OOM-killed gen_bathy.py on the login node
# (GOTCHA 2/3 in DECISIONS.md), so this runs as a PBS job, not on login.
#PBS -l select=1:ncpus=8:model=sky_ele
#PBS -q normal
#PBS -l walltime=2:00:00
#PBS -N kerpickups
#PBS -e /nobackup/dcarrol2/downscaling/kerguelen/pickups_err.txt
#PBS -o /nobackup/dcarrol2/downscaling/kerguelen/pickups_out.txt

D=/nobackup/dcarrol2/downscaling
K="$D/kerguelen"
UTILS="/nobackup/dcarrol2/v05_V5r1/ecco_darwin/regions/downscaling/utils"

export MAMBA_ROOT_PREFIX="$D/micromamba_root"
eval "$("$D/micromamba/bin/micromamba" shell hook --shell bash)"
micromamba activate downscaling

cd "$K"
python3 "$UTILS/gen_pickups.py" \
  -d "$K" -n kerguelen -i 736344 -bgc -v

echo "GEN_PICKUPS_DONE"
ls -la "$K/forcings/pickups/"
