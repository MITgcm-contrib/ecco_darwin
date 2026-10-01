#!/bin/bash
# job_gen_obcs_transport_correction.sh -- smoke test gen_obcs_transport_correction.py
# against the real STEP3 OBCS output (job 24944861, completed 2026-08-05).
# UVEL_*.bin/VVEL_*.bin backed up to forcings/OBCS_prebackup_transport_correction/
# before this runs, since the script overwrites them in place.
#PBS -l select=1:ncpus=4:model=rom_ait
#PBS -q normal
#PBS -l walltime=2:00:00
#PBS -N kertransport
#PBS -e /nobackup/dcarrol2/downscaling/kerguelen/transport_err.txt
#PBS -o /nobackup/dcarrol2/downscaling/kerguelen/transport_out.txt

D=/nobackup/dcarrol2/downscaling
K="$D/kerguelen"

export MAMBA_ROOT_PREFIX="$D/micromamba_root"
eval "$("$D/micromamba/bin/micromamba" shell hook --shell bash)"
micromamba activate downscaling

cd "$K"
python3 gen_obcs_transport_correction.py \
  -d "$K" -n kerguelen -bnd EWNS -i 736345 -v

echo "TRANSPORT_CORRECTION_DONE"
