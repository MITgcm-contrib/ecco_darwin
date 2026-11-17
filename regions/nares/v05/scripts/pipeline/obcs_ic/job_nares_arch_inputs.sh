#!/bin/bash
# Nares ice-arch test inputs (physics only: T, S, U, V, eta, sea ice), Jan-Jun 2015.
#   1. snapshot the first 4344 hourly boundary records (2015-01-01..06-30) of the live extraction
#   2. gen_pickups.py (ecco_darwin STEP3)  -> forcings/pickups (T, S, U, V, eta, sea ice, GGL90)
#   3. nares_seaice_ic.py                  -> pickup_seaice_si*.bin for data.seaice
#   4. nares_tref_sref.py                  -> tRef/sRef for input/data
#   5. gen_obcs_fast.py -seaice (no -bgc)  -> forcings/OBCS
#   6. nares_throughflow_offset.py         -> forcings/OBCS_tf (+0.46 Sv southward, net zero)
#   7. nares_shelfice_load.py              -> grid/Nares_shelficeLoadAnomaly.bin
# Run on a compute node (no analysis on pfe). Ask Dustin before qsub.
#PBS -l select=1:ncpus=28:model=bro_ele
#PBS -q normal
#PBS -l walltime=3:00:00
#PBS -N nares_arch_in
#PBS -j oe
#PBS -o /nobackup/dcarrol2/downscaling/nares/prep/job_nares_arch_inputs.out
set -u
D=/nobackup/dcarrol2/downscaling
N=$D/nares
PR=/nobackupp19/dcarrol2/downscaling/parent_V5r1_westAO_nares/run
UTILS=/nobackup/dcarrol2/v05_V5r1/ecco_darwin/regions/downscaling/utils
NREC=4344
export MAMBA_ROOT_PREFIX="$D/micromamba_root"
eval "$("$D/micromamba/bin/micromamba" shell hook --shell bash)"
micromamba activate downscaling
cd $N
step() { echo; echo "===== $1  $(date +%H:%M:%S)"; }
die() { echo "FATAL: $1"; exit 1; }

step "1 snapshot $NREC records"
rm -rf $N/parent/outputs/OBCS_2015H1
python3 prep/snapshot_dv.py $PR/dv $N/parent/inputs $N/parent/outputs/OBCS_2015H1 $NREC || die snapshot
ln -sfn $N/parent/outputs/OBCS_2015H1 $N/parent/outputs/OBCS

step "2 gen_pickups (physics + sea ice)"
python3 $UTILS/gen_pickups.py -d $N -n Nares -i 604872 -v -oh raise || die gen_pickups
ls -la $N/forcings/pickups

step "3 sea-ice IC split"
python3 prep/nares_seaice_ic.py $N/forcings/pickups 0000604872 || die seaice_ic

step "4 tRef/sRef"
python3 prep/nares_tref_sref.py $N/forcings/pickups $N/Nares_ncgrid.nc 0000604872 | tee $N/forcings/pickups/tref_sref_namelist.txt || die tref

step "5 OBCS (physics + sea ice, edges EWNS)"
python3 prep/gen_obcs_fast.py -d $N -n Nares -bnd EWNS -seaice -p 28 -v -disc off || die gen_obcs
ls -la $N/forcings/OBCS | head -50

step "6 throughflow offset"
python3 prep/nares_throughflow_offset.py $N 0.46 || die offset

step "7 Petermann load anomaly"
python3 prep/nares_shelfice_load.py $N $N/forcings/pickups/tref_sref.npz $N/grid/Nares_shelficeLoadAnomaly.bin || die load

step "ALL DONE"
exit 0
