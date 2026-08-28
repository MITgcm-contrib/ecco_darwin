# OAEMIP v05 LLC270 Darwin3 simulation based on ECCOV5r1 set-up with OAE experiments

# Code uses DAWRIN_OPTIONS.h flag DARWIN_OAE and leverages darwin_add_surfforc.Fto add alkalinity at surface from forcing files built from https://github.com/MITgcm-contrib/ecco_darwin/blob/master/code_util/OAEMIP/build_OAE_forcing_masks.m.

# Code includes a passive alkalinity tracer that tracks added alkalinity advection and diffusion (ALKtr). Runs can be started/initialized without modifying initial pickup.ptracers file. Last expected record in pickup.ptracers is skipped when nIter0 lower or equal PTRACERS_OAE_nopickupIter in data and data.ptracers, respectively. After PTRACERS_OAE_nopickupIter, added alkalinity tracer is being read from pickup files.

# Version 5 Release1  folder is also available at ECCO Drive
# To Download, one needs to have an Earthdata account
# (Or create it at https://urs.earthdata.nasa.gov/users/new)
# For using wget, one needs an Earthdata username and WebDAV password (different from Earthdata password)
# Find it at https://ecco.jpl.nasa.gov/drive
# and https://ecco-group.org/docs/wget_download_multiple_files_and_directories.pdf for more detail

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
MOD="../../ecco_darwin/v05/llc270_oaemip"
cp ../../ecco_darwin/v05/llc270/code_darwin/packages.conf .
sed -i '/ctrl/d;/smooth/d' packages.conf
sed -i '$ashelfice'        packages.conf
../tools/genmake2 -of ../tools/build_options/linux_amd64_ifort+mpi_ice_nas \
-mo "$MOD/code_darwin ../../ecco_darwin/v05/1deg_oaemip/code_darwin_v4r5 ../../ecco_darwin/v05/1deg/code_darwin_v4r5_v2 ../../ecco_darwin/v05/llc270/code_darwin ../../ecco_darwin/v05/llc270/code_v5r1" -mpi
make depend
make -j 16

==============
# 3. Instructions for running simulation (1992-2024 period)

cd ../run
rm -rf *
mkdir -p diags
ln -sf ../build/mitgcmuv .

==============
# 3. Instructions for running simulation (1992-2025 period)
INPUT="/nobackup/hzhang1/pub/llc270_FWD/v5r1"
ln -sf $INPUT/input_bin/* .
ln -sf $INPUT/input_darwin_bin/* .
ln -sf /nobackup/dcarrol2/forcing/apCO2/NOAA_MBL/* .
cp ../../ecco_darwin/v05/llc270/input_v5r1/* .
cp ../../ecco_darwin/v05/llc270/input2/data.{darwin,diagnostics,gchem,pkg,ptracers,traits}* .
sed -i "/nIter0=2/a \ pickupSuff='0000000002',"  data
cp $MOD/input_darwin/* .
mkdir diags/3hourly diags/daily diags/monthly diags/budget

==============
# copy AOE forcing files
ln -sf /nobackup/rsavelli/OAEMIP/forcings/test/LLC_270/zero/* .
ln -sf /nobackup/rsavelli/OAEMIP/forcings/test/LLC_270/nonzero/Ainjection_LLC_270_2003 .
mv Ainjection_LLC_270_2003 Ainjection0_LLC_270_2003
# qsub job_v5r1
