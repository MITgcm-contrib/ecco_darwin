#=========
#
# Guld of Mexico regional setup based on LLC4320 grid + LLC270 IC/BC
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
MOD="../../ecco_darwin/regions/GoM"
../tools/genmake2 -of ../tools/build_options/linux_amd64_ifort+mpi_ice_nas \
                  -mo ${MOD}/code -mpi
make depend
make -j 16

#=========
#3 Run model
 
cd ../run
mkdir -p diags
ln -sf ../build/mitgcmuv .
ln -sf /nobackup/hzhang1/pub/GoM_960/new/run_template/* .
ln -sf /nobackup/hzhang1/forcing/era5 ERA5
cp ${MOD}/input/* .
#change data, ... 
qsub job_gom960_bro

#=========
#4 Result
1995-1998 runs
https://ecco.jpl.nasa.gov/drive/files/ECCO2/LLC4320/GoM960/RUNOFF_exps/readme

