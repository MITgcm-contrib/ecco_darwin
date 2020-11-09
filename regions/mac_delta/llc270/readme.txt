# ========
#
# Mackenzie Delta regional setup based on LLC270
# WARNING: - Before starting make you have an Earthdata account (Or create it at: https://urs.earthdata.nasa.gov/users/new)
#          - If you run this setup on your laptop install a mpi library (as mpich or openmpi)
#
# ========

# ==============
# 1. Get code
# ==============
git clone https://github.com/MITgcm/MITgcm.git
svn checkout https://github.com/MITgcm-contrib/ecco_darwin/trunk/regions/mac_delta/llc270 Mac270

# Pleiades users skip to part 2.
# For the following requests you need your Earthdata username and WebDAV password (different from Earthdata password)
# Find it at :https://ecco.jpl.nasa.gov/drive
wget -r --no-parent --user=USERNAME --ask-password https://ecco.jpl.nasa.gov/drive/files/Version5/Alpha/era_xx
wget -r --no-parent --user=USERNAME --ask-password https://ecco.jpl.nasa.gov/drive/files/ECCO2/LLC270/Mac_Delta/run_template
mv ecco.jpl.nasa.gov/drive/files/Version5/Alpha/era_xx llc270_forcings/era_xx/
mv ecco.jpl.nasa.gov/drive/files/ECCO2/LLC270/Mac_Delta/run_template llc270_forcings/run_template_Mac270/
rm -r ecco.jpl.nasa.gov

# ================
# 2. Build executable for Mac Delta based on llc270 iteration 42 optimized solution
#    Prerequisite: 1. Get code
# ==============
cd MITgcm
mkdir build run
cd build

> On Pleidas follow intructions below:
   module purge
   module load comp-intel/2016.2.181 mpi-sgi/mpt.2.14r19 hdf4/4.2.12 hdf5/1.8.18_mpt netcdf/4.4.1.1_mpt
   ../tools/genmake2 -of ../tools/build_options/linux_amd64_ifort+mpi_ice_nas \
   -mo ../../llc270/code
   make depend
   make -j 16
 
 > On a laptop follow instructions below:
   export MPI_INC_DIR=PATH_TO_MPI_ENVIRONMENT_VARIABLE 
   # (path example on macintosh using homebrew: "/usr/local/opt/mpich/bin")
   ../tools/genmake2 -mpi -mo ../../Mac270/code
   make depend
   make -j 4

# ================
# 3. Run the setup
#    Prerequisite: 2. Build executable
#    Default running doesn't output diagnostic files
#    To enable diagnostics outputs follow instructions at the end of the readme
# ================
cd ../run
ln -sf ../build/mitgcmuv .
----------------
> On Pleiades:
  ln -sf /nobackup/hzhang1/forcing/era_xx .
  ln -sf /nobackup/hzhang1/pub/Mac_Delta270/run_template/* .
  cp ../../llc270/input/* .
  qsub job_Mac270_Bro
> On laptop:
  ln -sf ../../llc270_forcings/era_xx .
  ln -sf ../../llc270_forcings/run_template_Mac270/* .
  cp ../../Mac270/input/* .
  mpirun -np 4 ./mitgcmuv &
 ---------------
