# Instructions for building and running v06 ECCO-Darwin
# latest MITgcm code (c69e as of 2025/05/20)

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
module load comp-intel mpi-hpe/mpt hdf4/4.2.12 hdf5/1.8.18_mpt netcdf/4.4.1.1_mpt python3/3.9.12
../tools/genmake2 -of ../../ecco_darwin/v06/llc270/code_physics/linux_amd64_ifort+mpi_ice_nas \
  -mo '../../ecco_darwin/v06/llc270/code_darwin ../../ecco_darwin/v06/llc270/code_physics'
make depend
make -j 16

==============
# 2. Instructions for running v06 ECCO-Darwin for 1992-2024 period

cd ../run
ln -sf ../build/mitgcmuv .
ln -sf /nobackup/hzhang1/pub/llc270_FWD/nbp19_dmenemen_public_llc270/* .
ln -sf /nobackup/hzhang1/forcing/era_xx_it42_v2 .
cp ../../ecco_darwin/v06/llc270/input_physics/* .
cp ../../ecco_darwin/v06/llc270/input_darwin/* .
ln -sf /nobackup/ojahn/ecco_darwin/v06/llc270/data_darwin/* .
ln -sf /nobackup/ojahn/forcing/oasim .
ln -sf /nobackup/dcarrol2/forcing/apCO2/NOAA_MBL/* .
ln -sf /nobackup/rsavelli/LOAC/write_bin/jra55_do/v1.4.0/LLC_270/* .
ln -sf /nobackup/dcarrol2/pub/LLC_270/v06/* .
mkdir diags diags/3hourly diags/daily diags/monthly diags/budget
mkdir diags/monthly/IOPS diags/monthly/PAR diags/monthly/RRS

# modify job_ECCO_darwin as needed
qsub job_ECCO_darwin
