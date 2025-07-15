# Instructions for building and running Gulf of Mexico (GoM)
# regional simulation on Windows 11 based on "version 2"
# of ecco_darwin/v06/1deg/readme_darwin_v4r5.txt
# These instructions assume WSL has already been downloaded and installed on your windows machine. 

==============
# 1. Install mpi and python

# One way to do this is via miniconda. Go to https://www.anaconda.com/download/success and download the Windows Miniconda Installer. Follow screen prompts to run
  
# Install gfortran
  brew install gcc

# Install open-mpi
  brew install open-mpi

# Install python
  brew install python3
  cd /opt/homebrew/bin
  sudo ln -sf python3 python

==============
# 2. Get code, make sure you are in the directory where you want the model to be. This will create a directory named "darwin3" where the rest of the instructions take place.
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
     -of ../../MITgcm/tools/build_options/darwin_arm64_gfortran
  make depend
  make -j

==============
# 4. Instructions for running simulation (1992-2023 period)
  cd ../run
  ln -sf ../build/mitgcmuv .

# To get the forcing and configuration files, navigate to the websites below and download the .zip files:
# --> https://nasa-ext.box.com/s/3d3qz47tvnhp2y8wbvd821rwdxk1m2un
# --> https://nasa-ext.box.com/s/l7y193atfj5d4hxlwr8o1s3tvnbzecdp
# --> https://nasa-ext.box.com/s/7m8wv9gj2cyf8nhpg2d8f2b4n92b8ftr
# --> https://nasa-ext.box.com/s/g2dchqvo0t70qnwk8s2d49owb2oorq33

# To unzip the downloads install unzip using WSL and then 
# navigate to where the downloaded zip files are and execute unzip.
   sudo apt-get install unzip
   unzip NOAA_MBL.zip
   unzip run_template.zip
   unzip v06.zip
   unzip oasim.zip


# Either move the newly unzipped directories inside 
# darwin3/run or link the contents of these directories
# inside the darwin3/run directory, by running the following lines of code
# while in the darwin3/run directory:
  ln -sf <path_to_download_location>/NOAA_MBL/* .
  ln -sf <path_to_download_location>/GoM/llc90/run_template/* .
  ln -sf <path_to_download_location>/ECCO_Darwin/v06/* .
  ln -sf <path_to_download_location>/ECCO_Darwin/oasim .

  mkdir diags diags/3hourly diags/daily diags/monthly diags/budget
  mkdir diags/monthly/IOPS diags/monthly/PAR diags/monthly/RRS 
  cp -r ../../ecco_darwin/regions/GoM/llc90/input/* .

# To start the model simulation, execute mitgcmuv. The following 
# code will execute the model run in the background of your terminal 
# and will put output messages in the file output.txt
  ./mitgcmuv >& output.txt &
