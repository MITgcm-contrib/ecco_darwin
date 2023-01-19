# v05 1deg Darwin3 simulation based on ECCOV4r4 set-up
https://www.ecco-group.org/products-ECCO-V4r4.htm
https://ecco-group.org/docs/v4r4_reproduction_howto.pdf

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
../../../tools/genmake2 -mods='../../../../ecco_darwin/v05/1deg/code_darwin ../code' -optfile=../../../tools/build_options/linux_amd64_ifort+mpi_ice_nas -mpi
make depend
make -j 16

==============
# 3. Instructions for running simulation (1992-2017 period)

cd ../run
rm -rf *
ln -sf ../build/mitgcmuv .

# USERNAME='add your Earthdata username here'
# for new users, visit https://urs.earthdata.nasa.gov/users/new
# wget -r --no-parent --user ${USERNAME} --ask-password \
# https://ecco.jpl.nasa.gov/drive/files/Version4/Release4/input_forcing
# wget -r --no-parent --user ${USERNAME} --ask-password \
# https://ecco.jpl.nasa.gov/drive/files/Version4/Release4/input_init
# wget -r --no-parent --user ${USERNAME} --ask-password
# https://ecco.jpl.nasa.gov/drive/files/Version4/Release4/input_ecco

# INPUTDIR='./ecco.jpl.nasa.gov/drive/files/Version4/Release4/'
INPUTDIR='/nobackup/hzhang1/pub/Release4/'

ln -s ${INPUTDIR}/input_init/NAMELIST/* .
ln -s ${INPUTDIR}/input_init/error_weight/ctrl_weight/* .
ln -s ${INPUTDIR}/input_init/error_weight/data_error/* .
ln -s ${INPUTDIR}/input_init/* .
ln -s ${INPUTDIR}/input_init/tools/* .
ln -s ${INPUTDIR}/input_ecco/*/* .
ln -s ${INPUTDIR}/input_forcing/eccov4r4* .
python mkdir_subdir_diags.py

rm data.exf data.ptracers data.pkg data.diagnostics
cp ../../../../ecco_darwin/v05/1deg/input_darwin/* .
ln -sf /nobackup/dcarrol2/forcing/apCO2/NOAA_MBL/* .
ln -sf /nobackup/dcarrol2/v05_1deg/forcing/iron_dust/* .
ln -sf /nobackup/dcarrol2/v05_1deg/pickup/* .
# qsub job_ECCOV4r4
