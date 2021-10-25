# Instructions for building and running v06 ECCO-Darwin

==============
# 1. Get code

git clone https://github.com/darwinproject/darwin3
git clone --depth 1 https://github.com/MITgcm-contrib/ecco_darwin.git
cd darwin3
mkdir build run
cd build

==============
# 2. Build executable for v06 ECCO-Darwin

module purge
module load comp-intel/2016.2.181 mpi-sgi/mpt.2.14r19 hdf4/4.2.12 hdf5/1.8.18_mpt netcdf/4.4.1.1_mpt
../tools/genmake2 -of ../../ecco_darwin/v06/llc270/code/linux_amd64_ifort+mpi_ice_nas \
  -mo '../../ecco_darwin/v06/llc270/code_darwin ../../ecco_darwin/v06/llc270/code'
make depend
make -j 16

==============
# 2. Instructions for running v06 ECCO-Darwin for 1992-2020 period

cd ../run
ln -sf ../build/mitgcmuv .
ln -sf /nobackupp2/dmenemen/public/llc_270/iter42/input/* .
ln -sf /nobackupp2/dmenemen/public/llc_270/ecco_darwin_v5/input/darwin_initial_conditions/* .
ln -sf /nobackupp2/dmenemen/public/llc_270/ecco_darwin_v5/input/darwin_forcing/* .
ln -sf /nobackup/hzhang1/forcing/era_xx .
ln -sf /nobackup/hzhang1/pub/llc270_FWD/input/19920101/to2021/xx*42.data .
cp ../../ecco_darwin/v06/llc270/input/* .
ln -sf /nobackup/ojahn/ecco_darwin/v06/llc270/data_darwin/* .
ln -sf /nobackup/dcarrol2/forcing/apCO2/NOAA_MBL/* .
mkdir diags diags/3hourly diags/daily diags/monthly diags/budget diags
mkdir diags/monthly/IOPS diags/monthly/PAR diags/monthly/RRS 

# modify job_ECCO_darwin as needed
qsub job_ecco_darlin