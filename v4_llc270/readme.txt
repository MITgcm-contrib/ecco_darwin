==============
# Build executable for forward-only llc270 iteration 42 optimized solution
cvs -d :pserver:cvsanon:cvsanon@mitgcm.org:/u/gcmpack co -D "11/28/17" MITgcm_code
cvs -d :pserver:cvsanon:cvsanon@mitgcm.org:/u/gcmpack co MITgcm_contrib/ecco_darwin/v4_llc270
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
# Instructions for running forward-only llc270 optimized solution (1992-2015)
cd ..
mkdir run
cd run
mkdir diags
ln -sf ../build/mitgcmuv .
ln -sf /nobackupp2/dmenemen/public/llc_270/iter42/input/* .
ln -sf /nobackup/hzhang1/forcing/era_xx .
cp ../../MITgcm_contrib/ecco_darwin/v4_llc270/input/* .
# modify job_llc270_fdH as needed
qsub job_llc270_fdH


==============
# Creating input directory for llc270 iteration 42
mkdir /nobackupp2/dmenemen/public/llc_270/iter42/input
cd /nobackupp2/dmenemen/public/llc_270/iter42/input
cp /nobackup/hzhang1/obs/input/tile* .
cp /nobackup/hzhang1/obs/input/bathy270_filled_noCaspian_r4 .
cp /nobackup/hzhang1/obs/pri_err/smooth* .
cp /nobackup/hzhang1/pub/llc270_FWD/input/19920101/pickup*meta .
cp /nobackup/hzhang1/pub/llc270_FWD/input/19920101/pickup.0000000001.data_it42_v4 pickup.0000000001.data
cp /nobackup/hzhang1/pub/llc270_FWD/input/19920101/pickup_seaice.0000000001.data .
cp /nobackup/hzhang1/pub/llc270_FWD/input/19920101/pickup_ggl90.0000000001.data .
cp /nobackup/hzhang1/pub/llc270_FWD/input/wkapgmFld.data .
cp /nobackup/hzhang1/pub/llc270_FWD/input/wkaprediFld.data .
cp /nobackup/hzhang1/pub/llc270_FWD/input/wdiffkrFld.data .
cp /nobackup/hzhang1/pub/llc270_FWD/input/wprecip.data .
cp /nobackup/hzhang1/pub/llc270_FWD/input/wlwdown.data .
cp /nobackup/hzhang1/pub/llc270_FWD/input/wswdown.data .
cp /nobackup/hzhang1/pub/llc270_FWD/input/waqh.data .
cp /nobackup/hzhang1/pub/llc270_FWD/input/watemp.data .
cp /nobackup/hzhang1/pub/llc270_FWD/input/wuwind.data .
cp /nobackup/hzhang1/pub/llc270_FWD/input/wvwind.data .
cp /nobackup/hzhang1/obs/optim42/xx_kapgm.0000000042.data .
cp /nobackup/hzhang1/obs/optim42/xx_kapredi.0000000042.data .
cp /nobackup/hzhang1/obs/optim42/xx_diffkr.0000000042.data .
cp /nobackup/hzhang1/pub/llc270_FWD/input/19920101/to2018/xx_precip.0000000042.data .
cp /nobackup/hzhang1/pub/llc270_FWD/input/19920101/to2018/xx_lwdown.0000000042.data .
cp /nobackup/hzhang1/pub/llc270_FWD/input/19920101/to2018/xx_swdown.0000000042.data .
cp /nobackup/hzhang1/pub/llc270_FWD/input/19920101/to2018/xx_aqh.0000000042.data .
cp /nobackup/hzhang1/pub/llc270_FWD/input/19920101/to2018/xx_atemp.0000000042.data .
cp /nobackup/hzhang1/pub/llc270_FWD/input/19920101/to2018/xx_uwind.0000000042.data .
cp /nobackup/hzhang1/pub/llc270_FWD/input/19920101/to2018/xx_vwind.0000000042.data .
cp /nobackup/hzhang1/obs/input/runoff-2d-Fekete-1deg-mon-V4-SMOOTH.bin .
