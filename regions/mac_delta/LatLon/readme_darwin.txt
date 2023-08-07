#=========
#
# Mackenzie Delta regional setup based on LatLon w/ Darwin
#
#=========

#=========
#1 Get code

git clone --depth 1 https://github.com/MITgcm-contrib/ecco_darwin.git
git clone https://github.com/darwinproject/darwin3
cd darwin3
git checkout 24885b71
mkdir build run

#=========
#2 Build executable

cd build
sed "s|tidalComponents = 10|tidalComponents = 8|" \
	../pkg/obcs/OBCS_PARAMS.h>OBCS_PARAMS.h

diff ../pkg/obcs/OBCS_PARAMS.h OBCS_PARAMS.h

module purge
module load comp-intel/2020.4.304 mpi-hpe/mpt hdf4/4.2.12 hdf5/1.8.18_mpt netcdf/4.4.1.1_mpt python3/3.9.12
MOD="../../ecco_darwin/regions/mac_delta/LatLon"
../tools/genmake2 -of ../tools/build_options/linux_amd64_ifort+mpi_ice_nas \
 -mo '../../ecco_darwin/regions/mac_delta/LatLon/code_darwin ../../ecco_darwin/regions/mac_delta/LatLon/code' -mpi

make depend
make -j 16

#=========
#3 Run model

cd ../run
ln -sf ../build/mitgcmuv .
ln -sf /nobackup/hzhang1/pub/Mac_Delta_LatLon/run_template/* .
ln -sf /nobackup/hzhang1/forcing/era5 ERA5
ln -sf /nobackupp12/mwood7/Darwin/darwin3/configurations/downscaled_greenland/L1/L1_mac_delta/exf . 
ln -sf /nobackupp12/mwood7/Darwin/darwin3/configurations/downscaled_greenland/L1/L1_mac_delta/obcs . 
cp ${MOD}/input/* .
cp ${MOD}/input_darwin/* .
ln -sf /nobackup/dcarrol2/pub/regions/mac_delta/* .

mkdir diags
mkdir diags/EtaN_day_snap
mkdir diags/EtaN_mon_mean
mkdir diags/SI_daily_snap
mkdir diags/TS_surf_daily_snap
mkdir diags/vel_surf_daily_snap
mkdir diags/TS_AW_daily_snap
mkdir diags/state_3D_mon_snap
mkdir diags/vel_3D_mon_snap
mkdir diags/vel_3D_mon_mean
mkdir diags/BGC_daily_consts
mkdir diags/BGC_daily_DO
mkdir diags/BGC_daily_PO
mkdir diags/BGC_daily_misc
mkdir diags/BGC_daily_cx
mkdir diags/BGC_daily_Chl
 
#qsub job_Mac_Ivy

