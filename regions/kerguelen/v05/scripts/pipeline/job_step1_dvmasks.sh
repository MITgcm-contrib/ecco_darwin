#!/bin/bash
# job_step1_dvmasks.sh — run gen_dvmasks.py on a COMPUTE node (loads the full
# parent LLC270 grid arrays + the 2.6GB kerguelen_ncgrid.nc; comparable in
# scale to stitch_ncgrid.py's 2.7GB usage, over the session-cgroup ceiling).
# -r 23: local LLC270 grid spacing near the Kerguelen domain (62-78E,54-45S),
# computed directly from the parent DXC/DYC dump (mean ~23-24km), not a
# guessed nominal value.
#PBS -l select=1:ncpus=4:model=rom_ait
#PBS -q normal
#PBS -l walltime=1:00:00
#PBS -N kerdvmasks
#PBS -e /nobackup/dcarrol2/downscaling/kerguelen/dvmasks_err.txt
#PBS -o /nobackup/dcarrol2/downscaling/kerguelen/dvmasks_out.txt

D=/nobackup/dcarrol2/downscaling
K="$D/kerguelen"
UTILS="/nobackup/dcarrol2/v05_V5r1/ecco_darwin/regions/downscaling/utils"

export MAMBA_ROOT_PREFIX="$D/micromamba_root"
eval "$("$D/micromamba/bin/micromamba" shell hook --shell bash)"
micromamba activate downscaling

python3 "$UTILS/gen_dvmasks.py" -d "$K" -n kerguelen \
  -bfl bathy270_filled_noCaspian_r4 -bnd EWNS -r 23 -v

echo "DVMASKS_DONE"
ls -la "$K/parent/inputs"
