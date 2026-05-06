# Building RedSea/kaust_baby on arm64 macOS
# This is a development set-up, do not use for science

WORKDIR=~/mitgcm/

==============
# 1. Get code
  cd $WORKDIR
  git clone git@github.com:MITgcm-contrib/ecco_darwin.git
  git clone git@github.com:darwinproject/darwin3
  cd $WORKDIR/darwin3
  git checkout 24885b71
  mkdir build run

==============
# 2. Compile environment setup and darwin3 compatibility fixes
# See ecco_darwin/doc/MITgcm_on_Mac.txt for instructions on setting up
# compile environment (if needed) and fixing darwin3 compatibility issues.

==============
# 3. Build physics-only executable
  cd $WORKDIR/darwin3/build
  rm *
  ../tools/genmake2 -mo ../../ecco_darwin/regions/RedSea/kaust_baby/code -mpi \
   -of=../../ecco_darwin/regions/RedSea/kaust_baby/darwin_arm64_gfortran
  make depend
  make -j

==============
# 4. Run physics-only executable
  cd $WORKDIR/darwin3/run
  ln -sf ../build/mitgcmuv .

# Get forcing and configuration files from:
# --> https://nasa-ext.box.com/s/3d3qz47tvnhp2y8wbvd821rwdxk1m2un
# --> https://nasa-ext.box.com/s/mw2y1zu5z2ib81wqx5ywhb7bjs3lm4b9
# and deposit or link the contents of these directories
# inside the darwin3/run directory, for example,
  ln -sf <path_to_download_location>/NOAA_MBL/* .
  ln -sf <path_to_download_location>/RedSea/run_template/* .

# To save space, you can download only needed years for
# apCO2_* and era_xx_it42_v2

  mkdir diags diags/daily diags/monthly
  cp ../../ecco_darwin/regions/RedSea/kaust_baby/input/* .
  mpirun -np 8 ./mitgcmuv

==============
# 5. Build ECCO-Darwin executable
  cd $WORKDIR/darwin3/build
  rm *
  ../tools/genmake2 -of=../../ecco_darwin/regions/RedSea/kaust_baby/darwin_arm64_gfortran -mpi \
    -mo '../../ecco_darwin/regions/RedSea/kaust_baby/code_darwin ../../ecco_darwin/regions/RedSea/kaust_baby/code'
  make depend
  make -j

==============
# 6. Run ECCO-Darwin executable
  cd $WORKDIR/darwin3/run
  cp ../../ecco_darwin/regions/RedSea/kaust_baby/input_darwin/* .
  mpirun -np 8 ./mitgcmuv
