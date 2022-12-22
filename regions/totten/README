%Contents will no longer updated here. Please check 
https://github.com/MITgcm-contrib/ecco_polar_regional/tree/main/East_Antarctica_Totten_EccoDarwin/Nakayama_in_prep

==============
%How to build executable for Totten regional ECCO-Darwin run

v05 (20220303)
(build)
mkdir build5
cd build5
../../tools/genmake2 -of ../code/linux_amd64_ifort+mpi_ice_nas -mo '../code_darwin5 ../code'
make depend
make -j 16

(Note)
note that in darwin_sinking.F we need to modify this line to avoid problems in ice shelf code.

     IF (dzdn .GT. 0 _d 0) THEN

     IF (dzdn .GT. 0 d 0 .AND. dzup .GT. 0 d 0) THEN

(build_old) is created using darwin_old for debug purpose.

cp ../build/mitgcmuv .
cp ../input*s*/* .
ln -sf /nobackup/hzhang1/forcing/era_xx .
qsub run12_sandy_tracer_init.pbs

Note that the content of input and code_darwin are updated following changes in v05. 


%OLD (darwin3-dev)
v03 (This code is from Oliver Darwin v03 (according to Dustin))
# Build executable for ECCO-Darwin version 4 with Darwin3
git clone git://gud.mit.edu/darwin3-dev darwin3
cd darwin3
git checkout ecco_darwin_v4_llc270_darwin3
cd ecco_darwin/regions/totten/
mkdir build
cd build
module purge
module load comp-intel/2016.2.181 mpi-sgi/mpt.2.14r19 hdf4/4.2.12 hdf5/1.8.18_mpt netcdf/4.4.1.1_mpt
../../../tools/genmake2 -of \
  ../code/linux_amd64_ifort+mpi_ice_nas -mo '../code_darwin ../code'
make depend
make -j 16

==============
# Instructions for running ECCO-Darwin Version 4 for 1992-2018 period
cd ..
mkdir run
cd run
cp ../build/mitgcmuv .
cp ../input*s*/* .
ln -sf /nobackup/hzhang1/forcing/era_xx .
qsub run12_sandy_tracer_init.pbs
