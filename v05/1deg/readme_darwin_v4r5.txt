# v05 1deg Darwin3 simulation based on ECCOV4r5 set-up
#code base: c68g
#two versions:
#1 code_v4r5 + input_v4r5: 	including ctrl/smooth, similar to llc270
#2 code_v4r5_v2 input_v4r5_v2:	w/o ctrl/smooth, similar to v4r4
#				but for place holder now

# "Release5/" folder is also available at ECCO Drive
# To Download, one needs to have an Earthdata account
# (Or create it at https://urs.earthdata.nasa.gov/users/new)
# For using wget, one needs an Earthdata username and WebDAV password (different from Earthdata password)
# Find it at https://ecco.jpl.nasa.gov/drive
# and https://ecco-group.org/docs/wget_download_multiple_files_and_directories.pdf for more detail
wget -r --no-parent --user=USERNAME --ask-password https://ecco.jpl.nasa.gov/drive/files/ECCO2/LLC90/Release5

# ========
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
module load comp-intel/2020.4.304 mpi-hpe/mpt.2.25 hdf4/4.2.12 hdf5/1.8.18_mpt netcdf/4.4.1.1_mpt python3/3.9.12
MOD="../../ecco_darwin/v05/1deg"
../tools/genmake2 -of ../tools/build_options/linux_amd64_ifort+mpi_ice_nas \
	-mo '../../ecco_darwin/v05/1deg/code_darwin_v4r5 ../../ecco_darwin/v05/1deg/code_v4r5' -mpi
make depend
make -j 16

==============
# 3. Instructions for running simulation (1992-2017 period)

cd ../run
rm -rf *
mkdir -p diags
ln -sf ../build/mitgcmuv .

INPUTDIR='/nobackup/hzhang1/pub/Release5'
ln -s ${INPUTDIR}/input_bin/* .
ln -s ${INPUTDIR}/input_forcing/* .
cp ${MOD}/input_v4r5/* .

rm data data.pkg data.diagnostics
cp ${MOD}/input_darwin_v4r5/* .
ln -sf /nobackup/dcarrol2/forcing/apCO2/NOAA_MBL/* .
ln -sf /nobackup/dcarrol2/pub/1deg/v05/V4r5/* .
mkdir diags/3hourly diags/daily diags/monthly diags/budget
mv pickup_ptracers.0000000001.data pickup_ptracers.0000000002.data
mv pickup_ptracers.0000000001.meta pickup_ptracers.0000000002.meta
# qsub job_ECCOV4
