# iter42/input is available at https://data.nas.nasa.gov/ecco/data.php?dir=/eccodata/llc_270/iter42/input
# ecco_darwin_v4/input is available at https://data.nas.nasa.gov/ecco/data.php?dir=/eccodata/llc_270/ecco_darwin_v4/input
# forcing/era_xx is available at https://ecco.jpl.nasa.gov/drive/files/Version5/Alpha/era_xx
# LOACv1.4.0_HJ available at https://ecco.jpl.nasa.gov/drive/files/ECCO2/LLC270/LOAC/LOACv1.4.0_HJ
# LOACriver_temp available at https://ecco.jpl.nasa.gov/drive/files/ECCO2/LLC270/LOAC/LOACriver_temp

# Instructions for building and running ECCO-Darwin v4 with Darwin 3,
# JRA55-do, and globalNEWS nutrient runoff

==============
# 1. Get code

git clone --depth 1 https://github.com/darwinproject/darwin3
git clone --depth 1 https://github.com/MITgcm-contrib/ecco_darwin.git
cd darwin3
mkdir build run
cd build

==============
# 2. Build executable for ECCO-Darwin v4 with Darwin 3, JRA55-do, and globalNEWS nutrient runoff

module purge
module load comp-intel/2016.2.181 mpi-sgi/mpt.2.14r19 hdf4/4.2.12 hdf5/1.8.18_mpt netcdf/4.4.1.1_mpt
../tools/genmake2 -of ../../ecco_darwin/v05/v5_llc270_jra55do_nutrients/code/linux_amd64_ifort+mpi_ice_nas \
  -mo '../../ecco_darwin/v05/v5_llc270_jra55do_nutrients/code_darwin ../../ecco_darwin/v05/v5_llc270_jra55do_nutrients/code'
make depend
make -j 16

==============
# 3. Instructions for running ECCO-Darwin v4 with Darwin 3, JRA55-do, and globalNEWS nutrient runoff
# for 1992-2018 period

cd ../run
ln -sf ../build/mitgcmuv .
ln -sf /nobackupp2/dmenemen/public/llc_270/iter42/input/* .
ln -sf /nobackupp2/dmenemen/public/llc_270/ecco_darwin_v4/input/* .
ln -sf /nobackup/hzhang1/forcing/era_xx .
ln -sf /nobackup/hzhang1/forcing/jra55_do/LOACv1.4.0_HJ .
ln -sf /nobackup/hzhang1/forcing/jra55_do/LOACriver_temp .
Ln -sf /nobackup/dcarrol2/LOAC/write_bin/globalNEWS/* .
cp ../../ecco_darwin/v05/v5_llc270_jra55do/input/* .
cp ../../ecco_darwin/v05/v5_llc270_jra55do/input_darwin/* .
ln -sf /nobackup/dcarrol2/pickup/v05/* .
# modify job_llc270_fdH as needed
qsub job_llc270_fdH
