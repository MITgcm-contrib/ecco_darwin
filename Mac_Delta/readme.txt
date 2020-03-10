# Build executable for Mackenzie Delta regional model based on LLC4320
git clone https://github.com/MITgcm/MITgcm.git
git clone https://github.com/MITgcm-contrib/ecco-darwin.git


cd MITgcm
mkdir build
cd build
module purge
module load comp-intel/2016.2.181 mpi-sgi/mpt.2.14r19 hdf4/4.2.12 hdf5/1.8.18_mpt netcdf/4.4.1.1_mpt
../tools/genmake2 -of \
    ../tools/linux_amd64_ifort+mpi_ice_nas \
    -mo ../../MITgcm_contrib/ecco_darwin/Mac_Delta/code
make depend
make -j 16


==============
# Instructions for running Mackenzie Delta regional model based on LLC4320
cd ..
mkdir run
cd run
mkdir diags
ln -sf ../build/mitgcmuv .
ln -sf /nobackup/dmenemen/forcing/ECMWF_operational/EOG* .
ln -sf /nobackup/hzhang1/pub/Mac_Delta/run_template/* .
cp ../../MITgcm_contrib/ecco_darwin/Mac_Delta/input/* .
# modify job_mac_Bro_15 as needed
qsub job_mac_Bro_15

