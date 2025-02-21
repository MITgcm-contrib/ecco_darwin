# Instructions for building and running LR17 regional simulation
# on pleiades (based on ecco_darwin/v05/llc270/readme2.txt).

==============
# 1. Get code
  git clone git@github.com:MITgcm-contrib/ecco_darwin.git
  git clone git@github.com:darwinproject/darwin3
  cd darwin3
  git checkout 24885b71
  mkdir build run

==============
# 2. Build executable
  cd build
  module purge
  module load comp-intel mpi-hpe python3
  module load hdf4/4.2.12 hdf5/1.8.18_mpt netcdf/4.4.1.1_mpt
  ../tools/genmake2 -mo ../../ecco_darwin/regions/LR17/v05/code -mpi \
   -of ../../ecco_darwin/regions/LR17/v05/code/linux_amd64_ifort+mpi_ice_nas
  make depend
  make -j

==============
# 3. Instructions for running simulation (1992-2023 period)
  cd ../run
  ln -sf ../build/mitgcmuv .
  ln -sf /nobackup/dmenemen/ecco_darwin/LR17/run_template/* .
  ln -sf /nobackup/dcarrol2/forcing/apCO2/NOAA_MBL/* .
  mkdir diags diags/daily diags/monthly
  cp ../../ecco_darwin/regions/LR17/v05/input/* .
  mpirun -np 4 ./mitgcmuv
#    or modify a job_ECCO_darwin as needed and then:
#    qsub job_ECCO_darwin
