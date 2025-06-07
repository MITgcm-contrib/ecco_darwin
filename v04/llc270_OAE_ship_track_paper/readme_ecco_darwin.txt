# Instructions for llc270 ECCO-Darwin simulation
# Instructions are specific to pleiades and need to be adjusted for other platforms
# iter42/input is available at https://data.nas.nasa.gov/ecco/data.php?dir=/eccodata/llc_270/iter42/input
# ecco_darwin_v4/input is available at https://data.nas.nasa.gov/ecco/data.php?dir=/eccodata/llc_270/ecco_darwin_v4/input
# forcing/era_xx is available at https://ecco.jpl.nasa.gov/drive/files/Version5/Alpha/era_xx

#This solution is documented in:
#Carroll, D., Menemenlis, D., Adkins, J. F., Bowman, K. W., Brix, H., & Dutkiewicz, S., et al. (2020). 
#The ECCO-Darwin data-assimilative global ocean biogeochemistry model: Estimates of seasonal to multidecadal 
#surface ocean pCO2 and air-sea CO2 flux. Journal of Advances in Modeling Earth Systems, 12, e2019MS001888. 
#https://doi.org/10.1029/2019MS001888

#and modified for ship track OAE simulations in Sijia et al. (2025), Science Advances

#Note: you will need to generate initial condition for calcium as the last ptracer, 
#where CA = 1.028 _d -2*salinity/35. _d 0 * 1. _d 3
     
==============
# 1. Get code

git clone --depth 1 https://github.com/MITgcm-contrib/ecco_darwin.git
mkdir build run
cd build

==============
# 2. Build executable

module purge
module load comp-intel mpi-hpe/mpt hdf4/4.2.12 hdf5/1.8.18_mpt netcdf/4.4.1.1_mpt python3/3.9.12
../tools/genmake2 -of ../../ecco_darwin/v04/llc270_OAE_ship_track_paper/code/linux_amd64_ifort+mpi_ice_nas \
  -mo '../../ecco_darwin/v04/llc270_OAE_ship_track_paper/code_darwin ../../ecco_darwin/v04/llc270_OAE_ship_track_paper/code' -mpi
make depend
make -j 16

==============
# 3. Instructions for running simulation (1992-2023 period)

cd ../run
ln -sf ../build/mitgcmuv .
ln -sf /nobackupp19/dmenemen/public/llc_270/iter42/input/* .
ln -sf /nobackupp19/dmenemen/public/llc_270/ecco_darwin_v4/input/darwin_forcing/* .
ln -sf /nobackupp19/dmenemen/public/llc_270/ecco_darwin_v4/input/darwin_initial_conditions/pickup_ptracers_experiment_18.data pickup_ptracers.0000000001.data
ln -sf /nobackupp19/dmenemen/public/llc_270/ecco_darwin_v4/input/darwin_initial_conditions/pickup_ptracers.0000000001.meta .
ln -sf /nobackup/hzhang1/forcing/era_xx .
cp ../../ecco_darwin/v04/llc270_JAMES_paper/input/* .
cp ../../ecco_darwin/v04/llc270_JAMES_paper/input_darwin/* .
qsub job_ECCO_darwin

==============