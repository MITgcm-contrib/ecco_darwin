#!/bin/bash
# run_step1_mitgrid.sh — waits for the downscaling conda env to finish building,
# then runs gen_mitgrid.py for the Kerguelen Regional domain (STEP1.I).
# Domain: lon 62-78E, lat 54-45S; resolution tuned to ~1km at central lat 49.5S.
set -e
D=/nobackup/dcarrol2/downscaling
K="$D/kerguelen"
mkdir -p "$K/grid"

echo "[step1] waiting for micromamba env 'downscaling' to be ready..."
export MAMBA_ROOT_PREFIX="$D/micromamba_root"
MM="$D/micromamba/bin/micromamba"
while true; do
  if eval "$("$MM" shell hook --shell bash)" 2>/tmp/mitgrid_wait.err; then
    if micromamba activate downscaling 2>/dev/null; then
      if python3 -c "import simplegrid, pyproj, numpy" 2>/dev/null; then
        echo "[step1] env ready"
        break
      fi
    fi
  else
    echo "[step1] transient hiccup, retrying: $(cat /tmp/mitgrid_wait.err)"
  fi
  sleep 30
done

UTILS="$D/../v05_V5r1/ecco_darwin/regions/downscaling/utils"
# resolve absolute path robustly regardless of cwd
UTILS="/nobackup/dcarrol2/v05_V5r1/ecco_darwin/regions/downscaling/utils"

echo "[step1] running gen_mitgrid.py"
python3 "$UTILS/gen_mitgrid.py" \
  -d "$K/grid" -n kerguelen \
  -c 62 78 -54 -45 \
  -r 0.013831 0.0089833 \
  -v 2>&1 | tee "$K/grid/mitgrid_run.log"

# stable marker + parsed dims, so downstream steps never depend on a transient
# task-notification log path (bit us twice with grep-the-log-file coordination).
NROWS=$(grep -oE "Output shape: \([0-9]+, [0-9]+\)" "$K/grid/mitgrid_run.log" | grep -oE "[0-9]+" | sed -n 1p)
NCOLS=$(grep -oE "Output shape: \([0-9]+, [0-9]+\)" "$K/grid/mitgrid_run.log" | grep -oE "[0-9]+" | sed -n 2p)
echo "$NROWS $NCOLS" > "$K/grid/mitgrid_dims.txt"
echo "STEP1_MITGRID_DONE n_rows=$NROWS n_cols=$NCOLS"
touch "$K/grid/.mitgrid_done"
ls -la "$K/grid"
