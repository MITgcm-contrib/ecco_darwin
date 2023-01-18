# v05 1deg darwin3 simulation based on ECCOV4r4 set-up

# ========
# 1. Get code
git clone --branch backport-c66g https://github.com/jahn/darwin3
git clone --depth 1 https://github.com/MITgcm-contrib/ecco_darwin.git
cd darwin3

mkdir -p ECCOV4/release4
cd ECCOV4/release4
git clone https://github.com/ECCO-GROUP/ECCO-v4-Configurations.git
mv ECCO-v4-Configurations/ECCOv4\ Release\ 4/code .
rm -rf ECCO-v4-Configurations

# ================
# 2. Build executable
# Prerequisite: 1. Get code
mkdir build run
cd build
rm *
module load comp-intel/2020.4.304 mpi-hpe/mpt.2.25 hdf4/4.2.12 hdf5/1.8.18_mpt netcdf/4.4.1.1_mpt
../../../tools/genmake2 -mods=../code -optfile=../../../tools/build_options/linux_amd64_ifort+mpi_ice_nas -mpi
make depend
make all


../../../../ecco_darwin/v05/1deg/code_darwin

==============
# 3. Instructions for running simulation (1992-2017 period)

cd ../run
rm -rf *
ln -sf ../build/mitgcmuv .

# YOURUSERNAME = add your Earthdata username here
# for new users, visit https://urs.earthdata.nasa.gov/users/new

wget -r --no-parent --user ${YOURUSERNAME} --ask-password \
https://ecco.jpl.nasa.gov/drive/files/Version4/Release4/input_forcing
wget -r --no-parent --user ${YOURUSERNAME} --ask-password \
https://ecco.jpl.nasa.gov/drive/files/Version4/Release4/input_init
wget -r --no-parent --user ${YOURUSERNAME} --ask-password
https://ecco.jpl.nasa.gov/drive/files/Version4/Release4/input_ecco

module load comp-intel/2020.4.304 mpi-hpe/mpt.2.25 hdf4/4.2.12 hdf5/1.8.18_mpt netcdf/4.4.1.1_mpt
../../../tools/genmake2 -mods=../code -optfile=../../../tools/build_options/linux_amd64_ifort+mpi_ice_nas -mpi
make depend
make all
