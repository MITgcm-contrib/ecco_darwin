# iter42/input is available at https://data.nas.nasa.gov/ecco/data.php?dir=/eccodata/llc_270/iter42/input
# ecco_darwin_v5/input is available at https://data.nas.nasa.gov/ecco/data.php?dir=/eccodata/llc_270/ecco_darwin_v5/input
# forcing/era_xx is available at https://ecco.jpl.nasa.gov/drive/files/Version5/Alpha/era_xx

# Instructions for building and running ECCO-Darwin v05 with Darwin 3

==============
# 1. Get code

git clone https://github.com/darwinproject/darwin3
cd darwin3/pkg/darwin
git checkout 24885b71
cd ../../../
git clone --depth 1 https://github.com/MITgcm-contrib/ecco_darwin.git
cd darwin3
mkdir build run
cd build

==============
# 2. Build executable

module purge
module load comp-intel/2020.4.304 mpi-hpe/mpt.2.25 hdf4/4.2.12 hdf5/1.8.18_mpt netcdf/4.4.1.1_mpt
../tools/genmake2 -of ../../ecco_darwin/v05/llc270/code/linux_amd64_ifort+mpi_ice_nas \
  -mo '../../ecco_darwin/v05/llc270/code_darwin ../../ecco_darwin/v05/llc270/code'
make depend
make -j 16

==============
# 2. Instructions for running simulation (1992-2020 period)

cd ../run
ln -sf ../build/mitgcmuv .
ln -sf /nobackupp2/dmenemen/public/llc_270/iter42/input/* .
ln -sf /nobackupp2/dmenemen/public/llc_270/ecco_darwin_v5/input/darwin_initial_conditions/* .
ln -sf /nobackupp2/dmenemen/public/llc_270/ecco_darwin_v5/input/darwin_forcing/* .
ln -sf /nobackup/hzhang1/forcing/era_xx .
ln -sf /nobackup/hzhang1/pub/llc270_FWD/input/19920101/to2021/xx*42.data .
cp ../../ecco_darwin/v05/llc270/input/* .
mkdir diags diags/3hourly diags/daily diags/monthly diags/budget
# modify job_ECCO_darwin as needed
qsub job_ECCO_darwin