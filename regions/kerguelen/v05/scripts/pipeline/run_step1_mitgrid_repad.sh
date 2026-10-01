#!/bin/bash
# run_step1_mitgrid_repad.sh -- repad variant of run_step1_mitgrid.sh (2026-08-08).
# Same bounding box as the original (62-78E, 54-45S) but a slightly finer
# resolution, chosen empirically (see DECISIONS.md) to target a highly
# composite Nx=1170,Ny=1008 instead of the original's large-prime-dominated
# 1157x1002 -- so a fine-grained MPI decomposition becomes possible (see
# Task #11-13 timing-benchmark scope discussion).
set -e
D=/nobackup/dcarrol2/downscaling
K="$D/kerguelen"
mkdir -p "$K/grid"

export MAMBA_ROOT_PREFIX="$D/micromamba_root"
eval "$("$D/micromamba/bin/micromamba" shell hook --shell bash)"
micromamba activate downscaling

UTILS="/nobackup/dcarrol2/v05_V5r1/ecco_darwin/regions/downscaling/utils"

echo "[step1-repad] running gen_mitgrid.py (target Nx=1170,Ny=1008)"
python3 "$UTILS/gen_mitgrid.py" \
  -d "$K/grid" -n kerguelen \
  -c 62 78 -54 -45 \
  -r 0.0136811 0.0089330 \
  -v 2>&1 | tee "$K/grid/mitgrid_run.log"

NROWS=$(grep -oE "Output shape: \([0-9]+, [0-9]+\)" "$K/grid/mitgrid_run.log" | grep -oE "[0-9]+" | sed -n 1p)
NCOLS=$(grep -oE "Output shape: \([0-9]+, [0-9]+\)" "$K/grid/mitgrid_run.log" | grep -oE "[0-9]+" | sed -n 2p)
echo "$NROWS $NCOLS" > "$K/grid/mitgrid_dims.txt"
echo "STEP1_MITGRID_DONE n_rows=$NROWS n_cols=$NCOLS"
touch "$K/grid/.mitgrid_done"
ls -la "$K/grid"
