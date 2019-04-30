==============
# Build executable for ECCO-Darwin version 4
cvs -d :pserver:cvsanon:cvsanon@mitgcm.org:/u/gcmpack co -D "11/28/17" MITgcm_code
cvs -d :pserver:cvsanon:cvsanon@mitgcm.org:/u/gcmpack co MITgcm_contrib/ecco_darwin/v4_llc270
cvs -d :pserver:cvsanon:cvsanon@mitgcm.org:/u/gcmpack co -D "03/22/18" MITgcm_contrib/darwin/pkg/darwin
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
# Instructions for running ECCO-Darwin Version 4 for 1992-2018 period
cd ..
mkdir run
cd run
ln -sf ../build/mitgcmuv .
ln -sf /nobackupp2/dmenemen/public/llc_270/iter42/input/* .
ln -sf /nobackupp2/dmenemen/public/llc_270/ecco_darwin_v4/input/* .
ln -sf /nobackup/hzhang1/forcing/era_xx .
cp ../../MITgcm_contrib/ecco_darwin/v4_llc270/input/* .
cp ../../MITgcm_contrib/ecco_darwin/v4_llc270/input_darwin/* .
# modify job_llc270_fdH as needed
qsub job_llc270_fdH


==============
# Creating input directory for ECCO-Darwin version 4
mkdir /nobackupp2/dmenemen/public/llc_270/ecco_darwin_v4/input
cd /nobackupp2/dmenemen/public/llc_270/ecco_darwin_v4/input
cp /nobackup/dcarrol2/temp/pickup_ptracers_optimized.0000000001.data pickup_ptracers.0000000001.data
cp /nobackup/dcarrol2/temp/pickup_dic_optimized.0000000001.data pickup_dic.0000000001.data
cp /nobackup/dcarrol2/temp/pickup_ptracers.0000210384.meta pickup_ptracers.0000000001.meta 
cp /nobackup/dcarrol2/temp/pickup_dic.0000210384.meta pickup_dic.0000000001.meta
cp /nobackup/dcarrol2/forcing/apCO2/* .
cp /nobackup/dcarrol2/temp/llc270_Mahowald_2009_soluble_iron_dust.bin .
