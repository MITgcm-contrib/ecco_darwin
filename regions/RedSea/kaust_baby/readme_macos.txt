# Building physics-only RedSea/kaust_baby on arm64 macOS

==============
# 1. Get code
  git clone git@github.com:MITgcm-contrib/ecco_darwin.git
  git clone https://github.com/MITgcm/MITgcm
  cd MITgcm
  git checkout checkpoint67u
  mkdir build run

==============
# 2. Compile environment setup and darwin3 compatibility fixes
# See ecco_darwin/doc/MITgcm_on_Mac.txt for instructions on setting up
# compile environment (if needed) and fixing darwin3 compatibility issues.

==============
# 3. Build executable
  cd build
  rm *
  ../tools/genmake2 -mo ../../ecco_darwin/regions/RedSea/kaust_baby/code -mpi \
   -of=../../ecco_darwin/regions/RedSea/kaust_baby/darwin_arm64_gfortran
  make depend
  make -j

==============
# 4. Instructions for running simulation (1992-2023 period)
  cd ../run
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
