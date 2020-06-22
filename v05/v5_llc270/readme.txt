# iter42/input is available at https://data.nas.nasa.gov/ecco/data.php?dir=/eccodata/llc_270/iter42/input
# ecco_darwin_v4/input is available at https://data.nas.nasa.gov/ecco/data.php?dir=/eccodata/llc_270/ecco_darwin_v4/input
# forcing/era_xx is available at https://ecco.jpl.nasa.gov/drive/files/Version5/Alpha/era_xx

# Instructions for building and running ECCO-Darwin v4 with Darwin 3

==============
# 1. Get code

git clone --depth 1 https://github.com/darwinproject/darwin3
git clone --depth 1 https://github.com/MITgcm-contrib/ecco_darwin.git
cd darwin3
mkdir build run
cd build

==============
# 2. Build executable for ECCO-Darwin v4 with Darwin 3

module purge
module load comp-intel/2016.2.181 mpi-sgi/mpt.2.14r19 hdf4/4.2.12 hdf5/1.8.18_mpt netcdf/4.4.1.1_mpt
../tools/genmake2 -of ../../ecco_darwin/v05/v5_llc270/code/linux_amd64_ifort+mpi_ice_nas \
  -mo '../../ecco_darwin/v05/v5_llc270/code_darwin ../../ecco_darwin/v05/v5_llc270/code'
make depend
make -j 16

==============
# 2. Instructions for running ECCO-Darwin v4 with Darwin 3 for 1992-2018 period

cd ../run
ln -sf ../build/mitgcmuv .
ln -sf /nobackupp2/dmenemen/public/llc_270/iter42/input/* .
ln -sf /nobackupp2/dmenemen/public/llc_270/ecco_darwin_v4/input/* .
ln -sf /nobackup/hzhang1/forcing/era_xx .
cp ../../ecco_darwin/v05/v5_llc270/input/* .
ln -sf /nobackup/ojahn/ecco_darwin/v4_llc270/darwin3/pickup_ptracers_optimized.0000000001.data pickup_ptracers.0000000001.data
ln -sf /nobackup/ojahn/ecco_darwin/v4_llc270/darwin3/pickup_ptracers_optimized.0000000001.meta pickup_ptracers.0000000001.meta
# modify job_llc270_fdH as needed
qsub job_llc270_fdH




