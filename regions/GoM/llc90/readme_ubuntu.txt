# Instructions for building and running GoML regional simulation
# on Ubuntu 24.04.01 Noble Numbat (based on ecco_darwin/v05/llc270/readme2.txt).
# See ecco_darwin/doc/Ubuntu_on_Windows.txt and MITgcm_on_Ubuntu.txt
# for additional instructions.

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
  export MPI_INC_DIR=/usr/lib/x86_64-linux-gnu/openmpi/include
  ../tools/genmake2 -mo ../../ecco_darwin/regions/LR17/v05/code -mpi \
   -of ../../ecco_darwin/regions/LR17/v05/code/darwin_arm64_gfortran
  make depend
  make -j

==============
# 3. Instructions for running simulation (1992-2023 period)
  cd ../run
  ln -sf ../build/mitgcmuv .

# Get forcing and configuraMITgcmtion files from:
# --> https://nasa-ext.box.com/s/3d3qz47tvnhp2y8wbvd821rwdxk1m2un
# --> https://nasa-ext.box.com/s/zionyzanq7h4jf4rdw7aieiuao515kmk
# and deposit or link the contents of these directories
# inside the darwin3/run directory, for example,
  ln -sf <path_to_download_location>/NOAA_MBL/* .
  ln -sf <path_to_download_location>/LR17/run_template/* .

# To save space, you can download only needed years for
# apCO2_* and era_xx_it42_v2

  mkdir diags diags/daily diags/monthly
  cp ../../ecco_darwin/regions/LR17/v05/input/* .
  mpirun -np 4 ./mitgcmuv
