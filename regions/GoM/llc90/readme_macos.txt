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

# Note that darwin3 includes a full copy of MITgcm;
# MITgcm head branch is downloaded to access latest tools/build_options
  git clone --depth 1 https://github.com/MITgcm/MITgcm.git

# Use backport_ckpt68y for consistency with "version 2" of
# ecco_darwin/v06/1deg/readme_darwin_v4r5.txt
  cd darwin3
  git checkout backport_ckpt68y
  mkdir build run

==============
# 3. Build executable
  cd build
  ../tools/genmake2 -mo ../../ecco_darwin/regions/GoM/llc90/code \
     -of ../../MITgcm/tools/build_options/darwin_arm64_gfortran	-mpi
  make depend
  make -j

==============
# 4. Instructions for running simulation (1992-2023 period)
  cd ../run
  ln -sf ../build/mitgcmuv .

# Get forcing and configuration files from:
# --> https://nasa-ext.box.com/s/3d3qz47tvnhp2y8wbvd821rwdxk1m2un
# --> https://nasa-ext.box.com/s/l7y193atfj5d4hxlwr8o1s3tvnbzecdp
# and deposit or link the contents of these directories
# inside the darwin3/run directory, for example,
  ln -sf <path_to_download_location>/NOAA_MBL/* .
  ln -sf <path_to_download_location>/GoM/llc90/run_template/* .

  mkdir diags diags/daily diags/monthly
  cp ../../ecco_darwin/regions/GoM/llc90/input_v4r5/* .
  cp -r ../../ecco_darwin/regions/GoM/llc90/input_darwin_v4r5/* .
  mpirun -np 4 ./mitgcmuv 
