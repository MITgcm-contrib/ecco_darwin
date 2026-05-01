# Instructions for building and running a RedSea/kaust_v1s3
# ECCO-Darwin configuration on shaheen at KAUST

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
  rm *
  ../tools/genmake2 -mods ../../ecco_darwin/regions/RedSea/kaust_v1s3/code/  -of=../../ecco_darwin/regions/RedSea/kaust_v1s3/shaheen_build_options -mpi -make=gmake -rootdir=..
  make depend
  make -j 16

==============
# 3. Instructions for running simulation (1992-2023 period)
  cd ../run
  ln -sf ../build/mitgcmuv .
  ln -sf /nobackup/dmenemen/ecco_darwin/RedSea/run_template/* .
  ln -sf /nobackup/dcarrol2/forcing/apCO2/NOAA_MBL/* .
  mkdir diags diags/daily diags/monthly
  cp ../../ecco_darwin/regions/RedSea/v05/input/* .
  mpirun -np 8 ./mitgcmuv
