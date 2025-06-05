# Instructions for building and running ECCO-Darwin v05 with Darwin 3 
# experiments back in time from 1985/01/01
# adapted from readme.txt

==============
# 1. Get code

git clone --depth 1 https://github.com/MITgcm-contrib/ecco_darwin.git
git clone https://github.com/darwinproject/darwin3
cd darwin3
git checkout 24885b71
mkdir build run
cd build

==============
# 2. Build executable

module purge
module load comp-intel mpi-hpe/mpt hdf4/4.2.12 hdf5/1.8.18_mpt netcdf/4.4.1.1_mpt python3/3.9.12
../tools/genmake2 -of ../../ecco_darwin/v05/llc270/code/linux_amd64_ifort+mpi_ice_nas \
  -mo '../../ecco_darwin/v05/llc270/code_darwin ../../ecco_darwin/v05/llc270/code' -mpi
make depend
make -j 16

==============
# 3. Instructions for running simulation (1985-2022 period)
mkdir run_1985 run_1985_clim95 run_1985_linearCO2 run_1985_clim95_linearCO2 run_1985_clim99 run_1985_clim99_linearCO2

cd run_1985
ln -sf ../build/mitgcmuv .
ln -sf /nobackup/hzhang1/pub/llc270_FWD/nbp19_dmenemen_public_llc270/* .
ln -sf /nobackup/dcarrol2/forcing/apCO2/1984_NOAA_MBL/apCO2* .
ln -sf /nobackup/dcarrol2/pub/LLC_270/v05_1985_on/pickup_ptracers_1985_on.data \
	pickup_ptracers.0000000001.data
ln -sf /nobackup/hzhang1/forcing/exf_1985 .
ln -sf /nobackup/hzhang1/pub/llc270_FWD/input/19850101/pickup* .
cp ../../ecco_darwin/v05/llc270/input_1985/* .
mkdir diags diags/3hourly diags/daily diags/monthly diags/budget
# modify job_ECCO_darwin as needed
qsub job_ECCO_darwin

==============
# 4. Instructions for running simulation (1985-2022 period w/ repetitive year 1995)

cd run_1985_clim95
ln -sf ../build/mitgcmuv .
ln -sf /nobackup/hzhang1/pub/llc270_FWD/nbp19_dmenemen_public_llc270/* .
ln -sf /nobackup/dcarrol2/forcing/apCO2/1984_NOAA_MBL/apCO2* .
ln -sf /nobackup/dcarrol2/pub/LLC_270/v05_1985_on/pickup_ptracers_1985_on.data \
	pickup_ptracers.0000000001.data
ln -sf /nobackup/hzhang1/forcing/exf_1995 .
ln -sf /nobackup/hzhang1/pub/llc270_FWD/input/19850101/pickup* .
cp ../../ecco_darwin/v05/llc270/input_1995/* .
mkdir diags diags/3hourly diags/daily diags/monthly diags/budget
# modify job_ECCO_darwin as needed
qsub job_ECCO_darwin

==============
# 5. Instructions for running simulation (1985-2022 period w/ linear apCO2)

cd run_1985_linearCO2
ln -sf ../build/mitgcmuv .
ln -sf /nobackup/hzhang1/pub/llc270_FWD/nbp19_dmenemen_public_llc270/* .
ln -sf /nobackup/hzhang1/pub/linear_apCO2_forcing/apCO2* .
ln -sf /nobackup/dcarrol2/pub/LLC_270/v05_1985_on/pickup_ptracers_1985_on.data \
	pickup_ptracers.0000000001.data
ln -sf /nobackup/hzhang1/forcing/exf_1985 .
ln -sf /nobackup/hzhang1/pub/llc270_FWD/input/19850101/pickup* .
cp ../../ecco_darwin/v05/llc270/input_1985/* .
mkdir diags diags/3hourly diags/daily diags/monthly diags/budget
# modify job_ECCO_darwin as needed
qsub job_ECCO_darwin

==============
# 6. Instructions for running simulation (1985-2022 period w/ repetitive year 1995+linear apCO2)

cd run_1985_clim95_linearCO2
ln -sf ../build/mitgcmuv .
ln -sf /nobackup/hzhang1/pub/llc270_FWD/nbp19_dmenemen_public_llc270/* .
ln -sf /nobackup/hzhang1/pub/linear_apCO2_forcing/apCO2* .
ln -sf /nobackup/dcarrol2/pub/LLC_270/v05_1985_on/pickup_ptracers_1985_on.data \
	pickup_ptracers.0000000001.data
ln -sf /nobackup/hzhang1/forcing/exf_1995 .
ln -sf /nobackup/hzhang1/pub/llc270_FWD/input/19850101/pickup* .
cp ../../ecco_darwin/v05/llc270/input_1995/* .
mkdir diags diags/3hourly diags/daily diags/monthly diags/budget
# modify job_ECCO_darwin as needed
qsub job_ECCO_darwin

==============
# 7. Instructions for running simulation (1985-2022 period w/ repetitive year 1999)

cd run_1985_clim99
ln -sf ../build/mitgcmuv .
ln -sf /nobackup/hzhang1/pub/llc270_FWD/nbp19_dmenemen_public_llc270/* .
ln -sf /nobackup/dcarrol2/forcing/apCO2/1984_NOAA_MBL/apCO2* .
ln -sf /nobackup/dcarrol2/pub/LLC_270/v05_1985_on/pickup_ptracers_1985_on.data \
	pickup_ptracers.0000000001.data
ln -sf /nobackup/hzhang1/forcing/exf_1999 .
ln -sf /nobackup/hzhang1/pub/llc270_FWD/input/19850101/pickup* .
cp ../../ecco_darwin/v05/llc270/input_1999/* .
mkdir diags diags/3hourly diags/daily diags/monthly diags/budget
# modify job_ECCO_darwin as needed
qsub job_ECCO_darwin

==============
# 8. Instructions for running simulation (1985-2022 period w/ repetitive year 1999+linear apCO2)

cd run_1985_clim99_linearCO2
ln -sf ../build/mitgcmuv .
ln -sf /nobackup/hzhang1/pub/llc270_FWD/nbp19_dmenemen_public_llc270/* .
ln -sf /nobackup/hzhang1/pub/linear_apCO2_forcing/apCO2* .
ln -sf /nobackup/dcarrol2/pub/LLC_270/v05_1985_on/pickup_ptracers_1985_on.data \
	pickup_ptracers.0000000001.data
ln -sf /nobackup/hzhang1/forcing/exf_1999 .
ln -sf /nobackup/hzhang1/pub/llc270_FWD/input/19850101/pickup* .
cp ../../ecco_darwin/v05/llc270/input_1999/* .
mkdir diags diags/3hourly diags/daily diags/monthly diags/budget
# modify job_ECCO_darwin as needed
qsub job_ECCO_darwin

