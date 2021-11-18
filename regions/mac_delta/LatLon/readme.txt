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
#revert to 68c b/c https://github.com/MITgcm/MITgcm/pull/545 broken
# git checkout checkpoint68c
#no revert b/c https://github.com/MITgcm/MITgcm/pull/562 fixed
 mkdir build run

#=========
#2 Build executable
 cd build
 sed "s|tidalComponents = 10|tidalComponents = 8|" \
 	../pkg/obcs/OBCS_PARAMS.h>OBCS_PARAMS.h
   diff ../pkg/obcs/OBCS_PARAMS.h OBCS_PARAMS.h

 sed "s|heffTooHeavy = dzSurf \* 0.2 _d 0|heffTooHeavy = 2.0 _d 0|" \
        ../pkg/seaice/seaice_growth.F>seaice_growth.F
   diff ../pkg/seaice/seaice_growth.F seaice_growth.F 

 module purge
 module load comp-intel/2016.2.181 mpi-sgi/mpt.2.14r19 hdf4/4.2.12 hdf5/1.8.18_mpt netcdf/4.4.1.1_mpt
 MOD="../../ecco_darwin/regions/mac_delta/LatLon"
 ../tools/genmake2 -of ../tools/build_options/linux_amd64_ifort+mpi_ice_nas \
                   -mo ${MOD}/code
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

