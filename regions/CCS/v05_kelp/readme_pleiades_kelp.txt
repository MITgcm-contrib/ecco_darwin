# Instructions for building and running CCS regional simulation
# on pleiades (based on ecco_darwin/v05/llc270/readme2.txt).

==============
# 1. Get code
  git clone --depth 1 https://github.com/MITgcm-contrib/ecco_darwin.git
  git clone https://github.com/darwinproject/darwin3
  cd darwin3
  git checkout 24885b71
  mkdir build run

==============
# 2. Build executable
  cd build
  module purge
  module load comp-intel mpi-hpe python3
  module load hdf4/4.2.12 hdf5/1.8.18_mpt netcdf/4.4.1.1_mpt
  ../tools/genmake2 -of ../../ecco_darwin/regions/CCS/v05_kelp/code/linux_amd64_ifort+mpi_ice_nas \
  	-mo "../../ecco_darwin/regions/CCS/v05_kelp/code ../../ecco_darwin/regions/CCS/v05/code" -mpi
  make depend
  make -j

==============
# 3. Instructions for running simulation (1992-2023 period)
  cd ../run
  cp ../build/mitgcmuv .
  ln -sf /nobackup/dcarrol2/pub/regions/CCS/input/* .
  ln -sf /nobackup/hzhang1/pub/CCS_wave .
  ln -sf /nobackup/dcarrol2/forcing/apCO2/NOAA_MBL/* .
  mkdir diags diags/daily diags/monthly diags/monthly_kelp
  cp ../../ecco_darwin/regions/CCS/v05/input/* .
  cp ../../ecco_darwin/regions/CCS/v05_kelp/input/* .
  mpirun -np 20 ./mitgcmuv
#MM    or modify a job_ECCO_darwin as needed and then:
#    for long runs
#    qsub job_CCS_ED_MMlq_v2 
#    for 1-time step test run
#MM  qsub job_CCS_ED_MMdv_v2


=================================
use of wave forcing in MITgcm s/r:
=================================

#include "EXF_OPTIONS.h"
CBOP
#include "EXF_FIELDS.h"
CEOP

CBOB
#ifdef ALLOW_WAVE_FORCING
c assign   WvHeight
c assign   WvPeriod
#endif
CEOB


