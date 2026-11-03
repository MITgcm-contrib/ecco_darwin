#!/bin/bash
# setup_run.sh <run dir> [nTimeSteps]: link binary, forcing, OBCS, pickups; copy namelists.
# INP=input_v2 (default input) picks the namelist dir.
set -e
S=/nobackup/dcarrol2/downscaling/socal; R=$S/regional_run; D=$1
mkdir -p $D/diags/daily $D/diags/monthly; cd $D
cp $R/${INP:-input}/* .
ln -sf $R/${BUILD:-build}/mitgcmuv .
ln -sf $R/socal_bathy_r8.bin .
for f in $S/forcings/OBCS/*.bin; do ln -sf $f .; done
for f in $S/forcings/pickups/pickup_*; do ln -sf $f .; done
for f in $S/forcings/${TIDES:-OBCS_tides}/*.bin; do [ -e "$f" ] && ln -sf $f .; done
mkdir -p diags/hourly
ln -sfn /nobackup/hzhang1/forcing/era5 ERA5
ln -sfn $S/forcings/era5_seafill era5_seafill
for y in 2024 2025 2026; do ln -sf /nobackup/dcarrol2/forcing/apCO2/NOAA_MBL/apCO2_$y .; done
ln -sf /nobackup/dcarrol2/v05_V5r1/ecco_darwin/v05/3deg/data_darwin/3deg_Mahowald_2009_soluble_iron_dust.bin .
if [ -n "$2" ]; then sed -i "s/^ endtime=.*/ nTimeSteps=$2,/" data; grep -n nTimeSteps data; fi
ls | wc -l
