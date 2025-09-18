# Instructions for building and running Gulf of Mexico (GoM)
# regional simulation on macOS arm64 based on "version 2"
# of ecco_darwin/v06/1deg/readme_darwin_v4r5.txt

==============
# 1. Install gfortran, mpi, and python

# One way to do this is via the brew package manager
  /bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)"

# Install gfortran
  brew install gcc

# Install open-mpi
  brew install open-mpi

# Install python
  brew install python3
  cd /opt/homebrew/bin
  sudo ln -sf python3 python

==============
# 2. Get code
  git clone --depth 1 https://github.com/MITgcm-contrib/ecco_darwin.git
  git clone https://github.com/darwinproject/darwin3

# Use backport_ckpt68y for consistency with "version 2" of
# ecco_darwin/v06/1deg/readme_darwin_v4r5.txt
  cd darwin3
  git checkout backport_ckpt68y
  mkdir build run

==============
# 3. Build executable
  cd build
  ../tools/genmake2 -mo ../../ecco_darwin/regions/GoM/llc90/code
  make depend
  make -j

==============
# 4. Instructions for running simulation (1992-2023 period)
  cd ../run
  ln -sf ../build/mitgcmuv .

# Get contents of llc90/run_template from
# https://nasa-ext.box.com/s/a8iptfxasd95uxc4or9xmmgb68em8ext
# for example,
  GoM_PATH=$HOME/Box/Public/GoM/llc90/run_template
  ln -sf $GoM_PATH/* .

# Get contents of v06, oasim, and NOAA_MBL from
# https://nasa-ext.box.com/s/ln8r6w9by0w0g6gadr8ibpe1dn9ddiae
# for example,
  ECCO_Darwin_PATH=$HOME/Box/Public/ECCO_Darwin
  ln -sf $ECCO_Darwin_PATH/v06/* .
  ln -sf $ECCO_Darwin_PATH/oasim .
  ln -sf $ECCO_Darwin_PATH/NOAA_MBL/* .

  mkdir diags diags/3hourly diags/daily diags/monthly diags/budget
  mkdir diags/monthly/IOPS diags/monthly/PAR diags/monthly/RRS 
  cp -r ../../ecco_darwin/regions/GoM/llc90/input/* .
  ./mitgcmuv >& output.txt &
