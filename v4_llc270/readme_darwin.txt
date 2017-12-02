# Build executable for ECCO-Darwin version 4
cvs -d :pserver:cvsanon:cvsanon@mitgcm.org:/u/gcmpack co MITgcm_code
cvs -d :pserver:cvsanon:cvsanon@mitgcm.org:/u/gcmpack co MITgcm_contrib/darwin/pkg/darwin
cvs -d :pserver:cvsanon:cvsanon@mitgcm.org:/u/gcmpack co MITgcm_contrib/ecco_darwin/v4_llc270
cd MITgcm/pkg
ln -sf ../../MITgcm_contrib/darwin/pkg/darwin .
cd ..
mkdir build
cd build
module purge
module load comp-intel/2016.2.181 mpi-sgi/mpt.2.14r19 hdf4/4.2.12 hdf5/1.8.18_mpt netcdf/4.4.1.1_mpt
../tools/genmake2 -of \
  ../../MITgcm_contrib/ecco_darwin/v4_llc270/code/linux_amd64_ifort+mpi_ice_nas -mo \
  '../../MITgcm_contrib/ecco_darwin/v4_llc270/code_darwin ../../MITgcm_contrib/ecco_darwin/v4_llc270/code'
make depend
make -j 16


==============
# Instructions for running ECCO-Darwin Version 4 for 2009-2015 period
cd ..
mkdir run
cd run
ln -sf ../build/mitgcmuv .
cp ../../MITgcm_contrib/ecco_darwin/v4_llc270/input/* .
cp ../../MITgcm_contrib/ecco_darwin/v4_llc270/input_darwin/* .
cp data_2009 data
cp data.ctrl_2009 data.ctrl
ln -sf /nobackup/hzhang1/obs/input/bathy270_filled_noCaspian_r4 .
ln -sf /nobackup/hzhang1/forcing/era_xx .
ln -sf /nobackup/hzhang1/pub/llc270_FWD/input/20090101/pickup* .
ln -sf /nobackup/hzhang1/obs/input/runoff-2d-Fekete-1deg-mon-V4-SMOOTH.bin .
ln -sf /nobackup/hzhang1/obs/pri_err/smooth* .
ln -sf /nobackup/hzhang1/obs/input/tile* .
ln -sf /nobackup/hzhang1/pub/llc270_FWD/input/w*data .
ln -sf /nobackup/hzhang1/obs/optim33/xx_* .
ln -sf /nobackup/hbrix/ICBC/pco2a_blended_* .
ln -sf pco2a_blended_2012 pco2a_blended_2013
ln -sf pco2a_blended_2012 pco2a_blended_2014
ln -sf pco2a_blended_2012 pco2a_blended_2015
ln -sf pco2a_blended_2012 pco2a_blended_2016
ln -sf ~dmenemen/CMS/run_template_cg1/darwin* .
ln -sf /nobackup/dcarrol2/temp/pickup_llc270_dic.gf_cg1.0000078912.data pickup_dic.0000078912.data
ln -sf /nobackup/dcarrol2/temp/pickup_llc270_ptracers.gf_cg1.0000078912.data pickup_ptracers.0000078912.data
ln -sf /nobackup/dcarrol2/temp/pickup*.meta .
ln -sf /nobackup/hzhang1/cs510/ICBC_2009/pickup*.0000000001.??ta .
ln -sf /nobackup/dmenemen/forcing/runoff/runoff-360x180x12.bin .
ln -sf /nobackup/hzhang1/cs510/run_template/tile00* .
# modify job_llc270_fdH as needed
qsub job_llc270_fdH
