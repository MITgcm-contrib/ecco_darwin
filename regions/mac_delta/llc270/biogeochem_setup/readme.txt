# ========
# Mackenzie Delta regional setup based on LLC270 with Darwin3
# WARNING: - Before starting make you have an Earthdata account (Or create it at: https://urs.earthdata.nasa.gov/users/new)
#          - If you run this setup on your laptop install a mpi library (as mpich or openmpi)
# ========

# ==============
# 1. Get code
# ==============

# For Carroll et al. (2020) Global Ocean Ecosystem
git clone https://github.com/darwinproject/darwin3
cd Darwin3
git checkout checkpoint67x

# For Arctic Ocean Ecosystem
git clone -b cdom-carbon https://github.com/jahn/darwin3

svn checkout https://github.com/MITgcm-contrib/ecco_darwin/trunk/regions/mac_delta/llc270/biogeochem_setup ECCO-darwin_Mac270

# ==============
# 2. Get forcings
# ==============
# Pleiades users skip to part 3.
# For the following requests you need your Earthdata username and WebDAV password (different from Earthdata password)
# Find it at :https://ecco.jpl.nasa.gov/drive

### Atmospheric forcings ###
wget -r --no-parent --user=USERNAME --ask-password https://ecco.jpl.nasa.gov/drive/files/Version5/Alpha/era_xx
mv ecco.jpl.nasa.gov/drive/files/Version5/Alpha/era_xx llc270_forcings/era_xx/

### llc270 ECCO-darwin forcings ###
# Should be revised
wget -r --no-parent --user=USERNAME --ask-password https://ecco.jpl.nasa.gov/drive/files/ECCO2/LLC270/Mac_Delta/run_template
mv ecco.jpl.nasa.gov/drive/files/ECCO2/LLC270/Mac_Delta/run_template llc270_forcings/run_template_Mac270/

### Freshwater tiver runoff ###
# (soon avalaible on Earth data)

### River Temperature forcings ###
# (soon avalaible on Earth data)

### Biogeochem runoffs ###
# (soon avalaible on Earth data)

# ================
# 3. Build executable for Mac Delta based on llc270 iteration 42 optimized solution
#    Prerequisite: 1. Get code and 2. Get forcings
# ==============

cd darwin3
mkdir build run
cd build

> On Pleidas follow intructions below:
  module purge
  module load comp-intel/2016.2.181 mpi-sgi/mpt.2.14r19 hdf4/4.2.12 hdf5/1.8.18_mpt netcdf/4.4.1.1_mpt
  mv  ../../ECCO-darwin_Mac270/code/packages.conf ../../ECCO-darwin_Mac270/code/packages.conf_org
  ln -s ../../ECCO-darwin_Mac270/code_darwin/packages.conf ../../ECCO-darwin_Mac270/code
  ../tools/genmake2 -of ../tools/build_options/linux_amd64_ifort+mpi_ice_nas \
  -mo '../../ECCO-darwin_Mac270/code_darwin ../../ECCO-darwin_Mac270/code'
  make depend
  make -j 16

> On a laptop follow instructions below:
  export MPI_INC_DIR=PATH_TO_MPI_ENVIRONMENT_VARIABLE
  mv  ../../ECCO-darwin_Mac270/code/packages.conf ../../ECCO-darwin_Mac270/code/packages.conf_org
  ln -s ../../ECCO-darwin_Mac270/code_darwin/packages.conf ../../ECCO-darwin_Mac270/code
  ../tools/genmake2 -mpi -mo '../../ECCO-darwin_Mac270/code_darwin ../../ECCO-darwin_Mac270/code'
  make depend
  make -j 16

# ================
# 4. Run the setup
#    Prerequisite: 3. Build executable
# ================
cd ../run
mkdir diags diags/daily diags/budget
ln -sf ../build/mitgcmuv .

> On Pleiades:
  ln -sf /nobackup/hzhang1/forcing/era_xx .
  ln -sf /nobackup/hzhang1/pub/Mac_Delta270/run_template/* .
  cp ../../Mac270/input_darwin/* .
  qsub job_Mac270_Bro

> On laptop:
  ln -sf ../../llc270_forcings/era_xx .
  ln -sf ../../llc270_forcings/run_template/* .
  ln -sf ../../llc270_forcings/ArcticGRO .
  ln -sf ../../llc270_forcings/river_temp .
  cp ../../mac270/input_darwin/* .
  mpirun -np 8 ./mitgcmuv
