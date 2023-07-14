# Instructions for building and running CCS regional simulation
# Based on ecco_darwin/v05/llc270/readme.txt

==============
# 1. Get code

git clone git@github.com:MITgcm-contrib/ecco_darwin.git
git clone git@github.com:darwinproject/darwin3
cd darwin3
git checkout 24885b71
mkdir build run
cd build

==============
# 2. Build executable

module purge
module load comp-intel mpi-hpe python3
module load hdf4/4.2.12 hdf5/1.8.18_mpt netcdf/4.4.1.1_mpt
../tools/genmake2 -mo ../../ecco_darwin/regions/CCS/code -mpi \
         -of ../../ecco_darwin/regions/CCS/code/linux_amd64_ifort+mpi_ice_nas
make depend
make -j

==============
# 3. Instructions for running simulation (1992-2023 period)

cd ../run
ln -sf ../build/mitgcmuv .
ln -sf /nobackup/dmenemen/CCS_kelp/run_template/* .
#ln -sf /nobackupp19/dmenemen/public/llc_270/iter42/input/* .
#ln -sf /nobackupp19/dmenemen/public/llc_270/ecco_darwin_v5/input/darwin_initial_conditions/* .
ln -sf /nobackupp19/dmenemen/public/llc_270/ecco_darwin_v5/input/darwin_forcing/* .
ln -sf /nobackup/hzhang1/forcing/era_xx .
#ln -sf /nobackup/hzhang1/pub/llc270_FWD/input/19920101/to2023/xx*42.data .
mkdir diags diags/3hourly diags/daily diags/monthly diags/budget
cp ../../ecco_darwin/regions/CCS/input/* .
mpirun -np 16 ./mitgcmuv
# or modify job_ECCO_darwin as needed and then:
# qsub job_ECCO_darwin
