==============
# Build executable for forward-only llc270 iteration 42 optimized solution
cvs -d :pserver:cvsanon:cvsanon@mitgcm.org:/u/gcmpack co -D "11/28/17" MITgcm_code
cvs -d :pserver:cvsanon:cvsanon@mitgcm.org:/u/gcmpack co MITgcm_contrib/ecco_darwin/v5_llc270
cd MITgcm
mkdir build
cd build
module purge
module load comp-intel/2016.2.181 mpi-sgi/mpt.2.14r19 hdf4/4.2.12 hdf5/1.8.18_mpt netcdf/4.4.1.1_mpt
../tools/genmake2 -of \
    ../../MITgcm_contrib/ecco_darwin/v5_llc270/code/linux_amd64_ifort+mpi_ice_nas \
    -mo ../../MITgcm_contrib/ecco_darwin/v5_llc270/code
make depend
make -j 16


==============
# Instructions for running forward-only llc270 optimized solution (1992-2015)
cd ..
mkdir run
cd run
mkdir diags
ln -sf ../build/mitgcmuv .
ln -sf /nobackupp2/dmenemen/public/llc_270/iter42/input/* .
ln -sf /nobackup/hzhang1/forcing/era_xx .
ln -sf /nobackup/hzhang1/forcing/jra55_do/LOACv1.4.0_HJ .
ln -sf /nobackup/hzhang1/forcing/jra55_do/LOACriver_temp .
cp ../../MITgcm_contrib/ecco_darwin/v5_llc270/input/* .
# modify job_llc270_fdH as needed
qsub job_llc270_fdH
