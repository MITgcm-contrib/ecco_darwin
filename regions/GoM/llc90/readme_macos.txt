# Install gfortran
brew install gcc
# Extra package needed for parallel MPI compilation
brew install open-mpi
# If python is not available, one way to install is:
brew install python3
cd /opt/homebrew/bin
sudo ln -sf python3 python
# Make sure that /usr/local/bin is in $PATH, e.g.,
PATH=/usr/local/sbin:$PATH
export PATH

==============
# 1. Get code
git clone --branch backport_ckpt68g https://github.com/jahn/darwin3
git clone --depth 1 https://github.com/MITgcm-contrib/ecco_darwin.git
git clone --depth 1 https://github.com/MITgcm/MITgcm.git
cd darwin3
git checkout 24885b71
mkdir build run

# Note that, almost certainly, you will need to update the
# darwin3/tools/build_options. You can overwrite these files
# with those from MITgcm/tools/build_options downloaded above.

# Note that latest version of python (3.12.3) does not contain the
# imp module: https://docs.python.org/3.12/whatsnew/3.12.html#imp
# Need to replace darwin3/tools/darwin/cogapp/cogapp.py with
# ecco_darwin/doc/cogapp.py

==============
# 2. Build executable
  cd build
  ../tools/genmake2 -of ../../ecco_darwin/regions/GoM/llc90/code_v4r5/darwin_arm64_gfortran \
	-mo "../../ecco_darwin/regions/GoM/llc90/code_darwin_v4r5 ../../ecco_darwin/regions/GoM/llc90/code_v4r5" -mpi        
  make depend
  make -j

==============
# 3. Instructions for running simulation (1992-2023 period)
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
