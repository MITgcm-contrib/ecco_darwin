# Instructions for building and running CCS regional simulation
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
  ../tools/genmake2 -mo ../../ecco_darwin/regions/CCS/v05_kelp/code -mpi \
   -of ../../ecco_darwin/regions/CCS/v05_kelp/code/linux_amd64_ifort+mpi_ice_nas
  make depend
  make -j

==============
# 3. Instructions for running simulation (1992-2023 period)
  cd ../run
  cp ../build/mitgcmuv .
#  ln -sf /nobackupp18/mmanizza/Kelp/CCS/run_test/* .
  ln -sf /nobackup/dmenemen/ecco_darwin/CCS_kelp/run_template/* .
  ln -s /nobackupp18/mmanizza/Kelp/CCS/run_test/init2 .
  ln -s /nobackupp18/mmanizza/Kelp/CCS/run_test/myobcs .
  ln -sf /nobackup/dcarrol2/forcing/apCO2/NOAA_MBL/* .
  mkdir diags diags/daily diags/monthly diags/monthly3
  cp ../../ecco_darwin/regions/CCS/v05_kelp/input/* .
  mpirun -np 20 ./mitgcmuv
#MM    or modify a job_ECCO_darwin as needed and then:
#    for long runs
#    qsub job_CCS_ED_MMlq_v2 
#    for 1-time step test run
#MM  qsub job_CCS_ED_MMdv_v2
