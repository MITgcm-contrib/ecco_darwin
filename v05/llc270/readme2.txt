# iter42/input is available at https://data.nas.nasa.gov/ecco/data.php?dir=/eccodata/llc_270/iter42/input
# ecco_darwin_v5/input is available at https://data.nas.nasa.gov/ecco/data.php?dir=/eccodata/llc_270/ecco_darwin_v5/input
# forcing/era_xx is available at https://ecco.jpl.nasa.gov/drive/files/Version5/Alpha/era_xx

# Instructions for building and running ECCO-Darwin v05 with Darwin 3

#This solution is documented in:
#Carroll, D., Menemenlis, D., Dutkiewicz, S., Lauderdale, J. M., Adkins, J. F., Bowman, K. W., et al. (2022). 
#Attribution of space-time variability in global-ocean dissolved inorganic carbon. Global Biogeochemical Cycles, 
#36, e2021GB007162. https://doi.org/10.1029/2021GB007162

==============
# 1. Get code

git clone --depth 1 https://github.com/MITgcm-contrib/ecco_darwin.git
git clone https://github.com/darwinproject/darwin3
cd darwin3
git checkout 24885b71
mkdir build2 run2
cd build2

==============
# 2. Build executable

module purge
module load comp-intel mpi-hpe hdf4/4.2.12 hdf5/1.8.18_mpt netcdf/4.4.1.1_mpt
../tools/genmake2 -of ../../ecco_darwin/v05/llc270/code/linux_amd64_ifort+mpi_ice_nas \
  -mo ../../ecco_darwin/v05/code2 -mpi
make depend
make -j 16

==============
# 3. Instructions for running simulation (1992-2023 period)

cd ../run2
ln -sf ../build2/mitgcmuv .
ln -sf /nobackupp19/dmenemen/public/llc_270/iter42/input/* .
ln -sf /nobackupp19/dmenemen/public/llc_270/ecco_darwin_v5/input/darwin_initial_conditions/* .
ln -sf /nobackupp19/dmenemen/public/llc_270/ecco_darwin_v5/input/darwin_forcing/* .
ln -sf /nobackup/hzhang1/forcing/era_xx_it42_v2 .
cp ../../ecco_darwin/v05/llc270/input2/* .
mkdir diags diags/3hourly diags/daily diags/monthly diags/budget
# modify job_ECCO_darwin as needed
qsub job_ECCO_darwin

