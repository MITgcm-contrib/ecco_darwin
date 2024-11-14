# Instructions for building and running a Red Sea regional simulation
# on arm64 macOS (based on ecco_darwin/v05/llc270/readme2.txt).
# See ecco_darwin/doc/MITgcm_on_Mac.txt for additional instructions.

# This version replaces the single-column swfrac with swfrac2d to add
# capability to read space-time varying Jerlov water type.
# At present, swfrac2d is not consistent with pkg/layers and pkg/seaice

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
  ../tools/genmake2 -mo ../../ecco_darwin/regions/RedSea/v05_kpp/code -mpi \
   -of ../../ecco_darwin/regions/RedSea/v05/code/darwin_arm64_gfortran
  make depend
  make -j

==============
# 3. Instructions for running simulation (1992-2023 period)
  cd ../run
  ln -sf ../build/mitgcmuv .

# Get forcing and configuration files from:
# --> https://nasa-ext.box.com/s/3d3qz47tvnhp2y8wbvd821rwdxk1m2un
# --> https://nasa-ext.box.com/s/mw2y1zu5z2ib81wqx5ywhb7bjs3lm4b9
# --> https://nasa-ext.box.com/s/3tscnd3r3vbnyzcxdfh9avcbsvivj8cy
# and deposit or link the contents of these directories
# inside the darwin3/run directory, for example,
  ln -sf <path_to_download_location>/NOAA_MBL/* .
  ln -sf <path_to_download_location>/run_template/* .
  ln -sf <path_to_download_location>/swfrac/Jerlov_Kdclim12 * .

# To save space, you can download only needed years for
# apCO2_* and era_xx_it42_v2

  mkdir diags diags/daily diags/monthly
  cp ../../ecco_darwin/regions/RedSea/v05_kpp/input/* .
  mpirun -np 8 ./mitgcmuv
