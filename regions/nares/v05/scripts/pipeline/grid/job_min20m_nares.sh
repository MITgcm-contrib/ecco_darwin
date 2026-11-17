#!/bin/bash
# Nares 20 m minimum depth on COPIES (2026-10-15, Dustin D130). Same approach as WestAO job_min10m.sh: the 5 m/10 m-era
# inputs stay intact. grid/Nares_bathymetry.bin -> grid/Nares_bathymetry_min20m.bin (deepened in place, backup kept);
# forcings/pickups/* -> forcings/pickups_min20m/ (3-D pickups patched in place, npz backup). prep/deepen_min_depth_nares.py.
# hFacC of the current bathy from run_arch_TC2. Ask Dustin before qsub.
#PBS -l select=1:ncpus=40:model=sky_ele
#PBS -q devel
#PBS -l walltime=1:00:00
#PBS -N nares_min20m
#PBS -j oe
#PBS -o /nobackup/dcarrol2/downscaling/nares/prep/job_min20m_nares.out
set -eu
D=/nobackup/dcarrol2/downscaling
export MAMBA_ROOT_PREFIX="$D/micromamba_root"
eval "$("$D/micromamba/bin/micromamba" shell hook --shell bash)"
micromamba activate downscaling
N=$D/nares
B0=$N/grid/Nares_bathymetry.bin; B=$N/grid/Nares_bathymetry_min20m.bin; PK=$N/forcings/pickups_min20m
[ -e $B ] && { echo "FATAL: $B exists"; exit 1; }
[ -e $PK ] && { echo "FATAL: $PK exists"; exit 1; }
md5sum $B0                      # expect eb5e6532143745e7bab1d1e928cb3b82
cp -p $B0 $B
mkdir -p $PK; cp -p $N/forcings/pickups/* $PK/; ls $PK | wc -l
python3 $N/prep/deepen_min_depth_nares.py $B $N/grid/Nares_shelficeTopo.bin $PK $N/regional_run/run_arch_TC2 20
md5sum $B
echo MIN20M_DONE
