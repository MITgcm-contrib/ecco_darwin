#=========
#
# Mackenzie Delta regional setup based on LatLon
#
#=========

#=========
#1 Get code
git clone https://github.com/MITgcm/MITgcm.git
git clone https://github.com/MITgcm-contrib/ecco_darwin

 cd MITgcm
 mkdir build run

#=========
#2 Build executable
 cd build
 sed "s|tidalComponents = 10|tidalComponents = 8|" \
 	../pkg/obcs/OBCS_PARAMS.h>OBCS_PARAMS.h
   diff ../pkg/obcs/OBCS_PARAMS.h OBCS_PARAMS.h

 module purge
 module load comp-intel mpi-hpe hdf4/4.2.12 hdf5/1.8.18_mpt netcdf/4.4.1.1_mpt
 MOD="../../ecco_darwin/regions/mac_delta/LatLon"
 ../tools/genmake2 -of ../tools/build_options/linux_amd64_ifort+mpi_ice_nas \
                   -mo ${MOD}/code -mpi
 make depend
 make -j 16

#=========
#3 Run model
 cd ../run
 ln -sf ../build/mitgcmuv .
 ln -sf /nobackup/hzhang1/pub/Mac_Delta_LatLon/run_template/* .
 ln -sf /nobackup/hzhang1/forcing/era5 ERA5
 cp ${MOD}/input/* .
#change data, ... 
 qsub job_Mac_Ivy


#=========
# regarding pickup generation
pickup.0000000001.data
	corresponding to 2006/1/1
	1.6 m water is added to EtaN to make area mean SSH approximately zero
pickup.0015147648.data_1.1 (in case of "deltaT = 25.")
	corresponding to 2018/1/1
	1.1 m water is added to EtaN to make area mean SSH approximately zero
pickup.0015147648.data (for "deltaT = 25.")
	link from pickup.0015147648.data_1.1
pickup.0007573824.data (for "deltaT = 50.")
	link from pickup.0015147648.data
	link from pickup.0015147648.data_1.1
pickup.0007573824.data
	



