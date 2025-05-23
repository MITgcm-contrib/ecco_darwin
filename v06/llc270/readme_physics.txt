# Instructions for building and running v06 ECCO-Darwin llc270 physics
# latest MITgcm code (c69e as of 2025/05/20)

==============
# 1. Get code

git clone https://github.com/darwinproject/darwin3
git clone --depth 1 https://github.com/MITgcm-contrib/ecco_darwin.git
cd darwin3
mkdir build run
cd build

==============
# 2. Build executable for v06 ECCO-Darwin llc270 physics

module purge
module load comp-intel mpi-hpe/mpt hdf4/4.2.12 hdf5/1.8.18_mpt netcdf/4.4.1.1_mpt python3/3.9.12
../tools/genmake2 -of ../../ecco_darwin/v06/llc270/code/linux_amd64_ifort+mpi_ice_nas \
  -mo ../../ecco_darwin/v06/llc270/code_physics
make depend
make -j 16

==============
# 3. Instructions for running v06 ECCO-Darwin llc270 physics for 1992-2024 period

cd ../run
ln -sf ../build/mitgcmuv .
ln -sf /nobackupp19/dmenemen/public/llc_270/iter42/input/* .
ln -sf /nobackup/hzhang1/forcing/era_xx_it42_v2 .
cp ../../ecco_darwin/v06/llc270/input_physics/* .
mkdir diags diags/3hourly diags/daily diags/monthly diags/budget
mkdir diags/monthly/IOPS diags/monthly/PAR diags/monthly/RRS
# modify job_ECCO_darwin as needed
qsub job_ECCO_darwin

