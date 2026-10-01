#!/bin/bash
# job_step1_bathy.sh — run gen_bathy.py on a COMPUTE node. The login node's
# session memory cgroup OOM-killed it (scipy.griddata over ~8M GEBCO points
# onto a 1.16M-cell mitgrid). No internet needed (local mitgrid+GEBCO files),
# so this is a clean fit for PBS — unlike the conda env build.
#PBS -l select=1:ncpus=4:model=rom_ait
#PBS -q normal
#PBS -l walltime=2:00:00
#PBS -N kerbathy
#PBS -e /nobackup/dcarrol2/downscaling/kerguelen/bathy_err.txt
#PBS -o /nobackup/dcarrol2/downscaling/kerguelen/bathy_out.txt

D=/nobackup/dcarrol2/downscaling
K="$D/kerguelen"
UTILS="/nobackup/dcarrol2/v05_V5r1/ecco_darwin/regions/downscaling/utils"
GEBCO="$D/gebco/GEBCO_2024_CF.nc"

export MAMBA_ROOT_PREFIX="$D/micromamba_root"
eval "$("$D/micromamba/bin/micromamba" shell hook --shell bash)"
micromamba activate downscaling

read NROWS NCOLS < "$K/grid/mitgrid_dims.txt"
read WR WC < "$K/grid/wetpoint.txt"
echo "dims: $NROWS $NCOLS   wetpoint: $WR $WC"

python3 "$UTILS/gen_bathy.py" \
  -d "$K/grid" -g "$GEBCO" -n kerguelen \
  -s "$NROWS" "$NCOLS" \
  -cs 0.3 0.0 \
  -cw 1.00 1.14 \
  -wp "$WR" "$WC" \
  -v

echo "STEP1_BATHY_DONE"
touch "$K/grid/.bathy_done"
ls -la "$K/grid"
