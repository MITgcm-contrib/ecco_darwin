# v05 1deg Darwin3 simulation based on ECCOV4r5 set-up with RADIv2 metamodel
# include daily point source river freshwater runoff from JRA55-do and 
# biogechemical runoff capabilities from Savelli et al
# default set to set of equations for global (sed_depth_threshold = 0 in data.darwin)
# change sed_depth_threshold to desired depth for deep/coastal sets of equations
# for runs with RADI undef, please use data.diagnostics_noRADI and data.darwin_noRADI
# by default comes with data.diagnostics and data.darwin with RADI
#code base: c68g
# Version similar to v4r4 w/o ctrl/smooth, similar to v4r4

# "Release5/" folder is also available at ECCO Drive
# To Download, one needs to have an Earthdata account
# (Or create it at https://urs.earthdata.nasa.gov/users/new)
# For using wget, one needs an Earthdata username and WebDAV password (different from Earthdata password)
# Find it at https://ecco.jpl.nasa.gov/drive
# and https://ecco-group.org/docs/wget_download_multiple_files_and_directories.pdf for more detail
wget -r --no-parent --user=USERNAME --ask-password https://ecco.jpl.nasa.gov/drive/files/ECCO2/LLC90/Release5


# 1. Get code
git clone --branch backport_ckpt68g https://github.com/jahn/darwin3
git clone --depth 1 https://github.com/MITgcm-contrib/ecco_darwin.git
cd darwin3

# ================
# 2. Build executable
# Prerequisite: 1. Get code
mkdir build run
cd build
rm *
module load comp-intel/2020.4.304 mpi-hpe/mpt hdf4/4.2.12 hdf5/1.8.18_mpt netcdf/4.4.1.1_mpt python3/3.9.12
MOD="../../ecco_darwin/v05/1deg_RADIv2"
../tools/genmake2 -of ../tools/build_options/linux_amd64_ifort+mpi_ice_nas \
	-mo "${MOD}/code_darwin_v4r5_v2 ../../ecco_darwin/v05/1deg_runoff/code_darwin_v4r5_v2 ../../ecco_darwin/v05/1deg/code_darwin_v4r5_v2 ../../ecco_darwin/v05/1deg_runoff/code_v4r5_v2 ../../ecco_darwin/v05/1deg/code_v4r5_v2" -mpi
make depend
make -j 16

==============
# 3. Instructions for running simulation (1992-2024 period)

cd ../run
rm -rf *
mkdir -p diags
ln -sf ../build/mitgcmuv .

INPUTDIR='/nobackup/hzhang1/pub/Release5'
ln -s ${INPUTDIR}/input_bin/* .
ln -s ${INPUTDIR}/TBADJ .
cp ../../ecco_darwin/v05/1deg/input_v4r5_v2/* .


rm data data.pkg data.diagnostics
cp ../../ecco_darwin/v05/1deg/input_darwin_v4r5_v2/* .
cp ../../ecco_darwin/v05/1deg_runoff/input_v4r5_v2/* .
cp ${MOD}/input_darwin_v4r5_v2/* .
ln -sf /nobackup/dcarrol2/forcing/apCO2/NOAA_MBL/* .
ln -sf /nobackup/dcarrol2/pub/1deg/v05/V4r5/* .
ln -sf /nobackup/rsavelli/LOAC/ECCO_V4r5/freshwater_runoff/* .
ln -sf /nobackup/rsavelli/LOAC/ECCO_V4r5/bgc_runoff/* .
mkdir diags/3hourly diags/daily diags/monthly diags/budget


mv pickup_ptracers.0000000001.data pickup_ptracers.0000000002.data
mv pickup_ptracers.0000000001.meta pickup_ptracers.0000000002.meta
# qsub job_ECCOv4r5
