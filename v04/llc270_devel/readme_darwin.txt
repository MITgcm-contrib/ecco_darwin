# development version of ECCO-Darwin with nonlinear water-column dissolution and
# DIC/ALK fluxes from bottom sediments

# iter42/input is available at https://data.nas.nasa.gov/ecco/data.php?dir=/eccodata/llc_270/iter42/input
# ecco_darwin_v4/input is available at https://data.nas.nasa.gov/ecco/data.php?dir=/eccodata/llc_270/ecco_darwin_v4/input
# forcing/era_xx is available at https://ecco.jpl.nasa.gov/drive/files/Version5/Alpha/era_xx

==============
# Build executable for ECCO-Darwin version 4
cvs -d :pserver:cvsanon:cvsanon@mitgcm.org:/u/gcmpack co -D "11/28/17" MITgcm_code
cvs -d :pserver:cvsanon:cvsanon@mitgcm.org:/u/gcmpack co -D "03/22/18" MITgcm_contrib/darwin/pkg/darwin
git clone --depth 1 https://github.com/MITgcm-contrib/ecco_darwin.git
cd MITgcm/pkg
ln -sf ../../MITgcm_contrib/darwin/pkg/darwin .
cd ..
mkdir build
cd build
module purge
module load comp-intel/2016.2.181 mpi-sgi/mpt.2.14r19 hdf4/4.2.12 hdf5/1.8.18_mpt netcdf/4.4.1.1_mpt
../tools/genmake2 -of \
  ../../ecco_darwin/v04/llc270_devel/code/linux_amd64_ifort+mpi_ice_nas -mo \
  '../../ecco_darwin/v04/llc270_devel/code_darwin ../../ecco_darwin/v04/llc270_devel/code'
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
cp ../../ecco_darwin/v04/llc270_devel/input/* .
cp ../../ecco_darwin/v04/llc270_devel/input_darwin/* .
# modify job_llc270_fdH as needed
qsub job_llc270_fdH
