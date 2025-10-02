# Instructions for building and running LR17 regional simulation
# on macOS (based on ecco_darwin/v05/llc270/readme2.txt).
# See ecco_darwin/doc/MITgcm_on_Mac.txt for additional instructions.

==============
# 1. Get code
  git clone git@github.com:MITgcm-contrib/ecco_darwin.git
  git clone git@github.com:darwinproject/darwin3
  cd darwin3
  git checkout 24885b71
  cp ../ecco_darwin/doc/cogapp.py tools/darwin/cogapp
  mkdir build run

==============
# 2. Build executable
  cd build

# Set location of MPI_INC_DIR, for example,
  export MPI_INC_DIR=/opt/homebrew/include

  ../tools/genmake2 -mo ../../ecco_darwin/regions/LR17/v05/code -mpi \
   -of ../../ecco_darwin/regions/LR17/v05/code/darwin_arm64_gfortran
  make depend
  make -j

==============
# 3. Instructions for running simulation (1992-2023 period)
  cd ../run
  ln -sf ../build/mitgcmuv .

# Get contents of run_template from
# https://nasa-ext.box.com/s/6tjjkb01ibz8qnmtenf39fnnceof26mu
# for example,
  LR17_PATH=$HOME/Box/Public/LR17/run_template
  ln -sf $LR17_PATH/* .

# Get contents of NOAA_MBL from
# https://nasa-ext.box.com/s/ln8r6w9by0w0g6gadr8ibpe1dn9ddiae
# for example,
  ECCO_Darwin_PATH=$HOME/Box/Public/ECCO_Darwin
  ln -sf $ECCO_Darwin_PATH/NOAA_MBL/* .

  mkdir diags diags/daily diags/monthly
  cp ../../ecco_darwin/regions/LR17/v05/input/* .
  mpirun -np 4 ./mitgcmuv
