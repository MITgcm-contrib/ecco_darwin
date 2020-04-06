# Instructions for building and running ECCO-Darwin v4 with Darwin 3 and JRA55-do 

==============
# 1. Get code

git clone --depth 1 https://github.com/darwinproject/darwin3
git clone --depth 1 https://github.com/MITgcm-contrib/ecco_darwin.git
cd darwin3
mkdir build run
cd build

==============
# 2. Build executable for ECCO-Darwin v4 with Darwin 3 and JRA55-do

module purge
module load comp-intel/2016.2.181 mpi-sgi/mpt.2.14r19 hdf4/4.2.12 hdf5/1.8.18_mpt netcdf/4.4.1.1_mpt
../tools/genmake2 -of ../../ecco_darwin/v05/v5_llc270_jra55do/code/linux_amd64_ifort+mpi_ice_nas \
  -mo '../../ecco_darwin/v05/v5_llc270_jra55do/code_darwin ../../ecco_darwin/v05/v5_llc270_jra55do/code'
make depend
make -j 16

==============
# 3. Instructions for running ECCO-Darwin v4 with Darwin 3 and JRA55-do for 1992-2018 period

cd ../run
ln -sf ../build/mitgcmuv .
ln -sf /nobackupp2/dmenemen/public/llc_270/iter42/input/* .
ln -sf /nobackupp2/dmenemen/public/llc_270/ecco_darwin_v4/input/* .
ln -sf /nobackup/hzhang1/forcing/era_xx .
ln -sf /nobackup/hzhang1/forcing/jra55_do/LOACv1.4.0_HJ .
ln -sf /nobackup/hzhang1/forcing/jra55_do/LOACriver_temp .
cp ../../ecco_darwin/v05/v5_llc270_jra55do/input/* .
cp ../../ecco_darwin/v05/v5_llc270_jra55do/input_darwin/* .
ln -sf /nobackup/ojahn/ecco_darwin/v4_llc270/darwin3/pickup_ptracers_optimized.0000000001.data pickup_ptracers.0000000001.data
ln -sf /nobackup/ojahn/ecco_darwin/v4_llc270/darwin3/pickup_ptracers_optimized.0000000001.meta pickup_ptracers.0000000001.meta
# modify job_llc270_fdH as needed
qsub job_llc270_fdH