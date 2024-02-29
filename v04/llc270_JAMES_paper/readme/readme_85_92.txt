# Instructions for llc270 physical simulation without Darwin,
# including modifications for 1985-1992 back-in-time extension.

# Instructions are specific to pleiades and need to be adjusted for other platforms
# iter42/input is available at https://data.nas.nasa.gov/ecco/data.php?dir=/eccodata/llc_270/iter42/input
# ecco_darwin_v4/input is available at https://data.nas.nasa.gov/ecco/data.php?dir=/eccodata/llc_270/ecco_darwin_v4/input
# forcing/era_xx is available at https://ecco.jpl.nasa.gov/drive/files/Version5/Alpha/era_xx

==============
# Build executable for forward-only llc270 iteration 42 optimized solution
 cvs -d :pserver:cvsanon:cvsanon@mitgcm.org:/u/gcmpack co -D "11/28/17" MITgcm_code
 git clone https://github.com/MITgcm-contrib/ecco_darwin
 cd MITgcm
 mkdir build run
 cd build
 module purge
 module load comp-intel/2020.4.304 mpi-hpe hdf4/4.2.12 hdf5/1.8.18_mpt netcdf/4.4.1.1_mpt
 MOD="../../ecco_darwin/v04/llc270_JAMES_paper"
 ../tools/genmake2 -of ${MOD}/code/linux_amd64_ifort+mpi_ice_nas \
                   -mo "${MOD}/code ${MOD}/code_85_92 "
 make depend
 make -j 16

==============
# Instructions for running forward-only llc270 optimized solution (1992-2017)
 cd ../run
 mkdir -p diags
 ln -sf ../build/mitgcmuv .
 ln -sf /nobackupp19/dmenemen/public/llc_270/iter42/input/* .
 ln -sf /nobackup/hzhang1/forcing/era_xx .
 cp ${MOD}/input/* .
 cp ${MOD}/input_85_92/* .
 qsub job_ECCO_darwin

==============
# Creating input directory for llc270 iteration 42
 mkdir /nobackupp19/dmenemen/public/llc_270/iter42/input
 cd /nobackupp19/dmenemen/public/llc_270/iter42/input
 cp /nobackup/hzhang1/obs/input/tile* .
 cp /nobackup/hzhang1/obs/input/bathy270_filled_noCaspian_r4 .
 cp /nobackup/hzhang1/obs/pri_err/smooth* .
 cp /nobackup/hzhang1/pub/llc270_FWD/input/19920101/pickup*meta .
 cp /nobackup/hzhang1/pub/llc270_FWD/input/19920101/pickup.0000000001.data_it42_xxr8 pickup.0000000001.data
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
