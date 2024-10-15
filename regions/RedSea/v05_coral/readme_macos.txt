# Instructions for building and running a Red Sea regional simulation
# on arm64 macOS (based on ecco_darwin/v05/llc270/readme2.txt).
# See ecco_darwin/doc/MITgcm_on_Mac.txt for additional instructions.

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
  ../tools/genmake2 -mo ../../ecco_darwin/regions/RedSea/v05/code -mpi \
   -of ../../ecco_darwin/regions/RedSea/v05/code/darwin_arm64_gfortran
  make depend
  make -j

==============
# 3. Instructions for running simulation (1992-2023 period)
  cd ../run
  ln -sf ../build/mitgcmuv .

# Get forcing and configuration files from
# https://nasa-ext.box.com/s/3d3qz47tvnhp2y8wbvd821rwdxk1m2un
# https://nasa-ext.box.com/s/mw2y1zu5z2ib81wqx5ywhb7bjs3lm4b9
# and deposit or link inside the darwin3/run directory, e.g.,
  ln -sf ~/NOAA_MBL/* .
  ln -sf ~/RedSea/run_template/* .

# To save space, you can download only needed years for
# apCO2_* and era_xx_it42_v2

  mkdir diags diags/daily diags/monthly
  cp ../../ecco_darwin/regions/RedSea/v05/input/* .
  mpirun -np 8 ./mitgcmuv
