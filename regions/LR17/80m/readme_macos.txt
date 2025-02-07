# Instructions for building and running LR17 regional simulation
# on macOS with 100-km horizontal grid spacing.
#
# See ecco_darwin/doc/MITgcm_on_Mac.txt for instructions
# for running MITgcm on Mac.

==============
# 1 Download run_template from:
# https://nasa-ext.box.com/s/6ubnsrzuvsx38hdz7rtowpkg8ra8nikw
# Instructions below assume that run_template, MITgcm, and
# ecco_darwin are at same directory level

==============
# 2. Get code
  git clone git@github.com:MITgcm/MITgcm.git
  git clone git@github.com:MITgcm-contrib/ecco_darwin.git
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
  ln -sf ../../run_template/bathymetry/bathy_LR17_936x875 .
  ln -sf ../../run_template/era_xx_it42_v2 .
  cp ../../ecco_darwin/regions/LR17/100km/input/* .
  mpirun -np 4 ./mitgcmuv
