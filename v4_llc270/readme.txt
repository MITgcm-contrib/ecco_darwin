# Build executable for Hong Zhang's forward-only llc270 optimized solution

cvs -d :pserver:cvsanon:cvsanon@mitgcm.org:/u/gcmpack co MITgcm_code
cd MITgcm
mkdir build
cd build
module purge
module load comp-intel/2016.2.181 mpi-sgi/mpt.2.14r19 hdf4/4.2.12 hdf5/1.8.18_mpt netcdf/4.4.1.1_mpt
../tools/genmake2 -of \
    ../../MITgcm_contrib/ecco_darwin/v4_llc270/code/linux_amd64_ifort+mpi_ice_nas \
    -mo ../../MITgcm_contrib/ecco_darwin/v4_llc270/code
make depend
make -j 16

==============
# Instructions for running forward-only llc270 optimized solution (2001-2015)

cd MITgcm
mkdir run
cd run
ln -sf ../build/mitgcmuv .
cp ../../MITgcm_contrib/ecco_darwin/v4_llc270/input/* .
ln -sf /nobackup/hzhang1/obs/input/bathy270_filled_noCaspian_r4 .
ln -sf /nobackup/hzhang1/forcing/era_xx .
ln -sf /nobackup/hzhang1/obs/input/pickup* .
ln -sf /nobackup/hzhang1/obs/input/runoff-2d-Fekete-1deg-mon-V4-SMOOTH.bin .
ln -sf /nobackup/hzhang1/obs/pri_err/smooth* .
ln -sf /nobackup/hzhang1/obs/input/tile* .
ln -sf /nobackup/hzhang1/pub/llc270_FWD/input/w*data .
ln -sf /nobackup/hzhang1/obs/optim33/xx_* .

# modify job_llc270_fdH as needed
qsub job_llc270_fdH
