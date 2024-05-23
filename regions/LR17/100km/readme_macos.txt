# Instructions for building and running LR17 regional simulation
# on macOS with 100-km horizontal grid spacing.
#
# See ecco_darwin/doc/MITgcm_on_Mac.txt for instructions on
# running MITgcm on Mac.
#
# Some LR17 input files needed are here:
# https://nasa-ext.box.com/s/6ubnsrzuvsx38hdz7rtowpkg8ra8nikw

==============
# 1. Get code
  git clone git@github.com:MITgcm/MITgcm.git
  git clone git@github.com:darwinproject/darwin3
  cd MITgcm
  mkdir build run

==============
# 2. Build executable
  cd build
  ../tools/genmake2 -mo ../../ecco_darwin/regions/LR17/100km/code -mpi
  make depend
  make -j

==============
# 3. Instructions for running simulation (1992-2023 period)
  cd ../run
  ln -sf ../build/mitgcmuv .
#    Get forcing and configuration files from
#    https://nasa-ext.box.com/s/3d3qz47tvnhp2y8wbvd821rwdxk1m2un
#    https://nasa-ext.box.com/s/zionyzanq7h4jf4rdw7aieiuao515kmk
#    and deposit or link inside the darwin3/run directory.
#    To save space, you can download only needed years for
#    apCO2_* and era_xx_it42_v2
  mkdir diags diags/daily diags/monthly
  cp ../../ecco_darwin/regions/LR17/v05/input/* .
  mpirun -np 4 ./mitgcmuv
