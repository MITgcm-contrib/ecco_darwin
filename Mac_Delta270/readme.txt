# forcing/era_xx is available at https://ecco.jpl.nasa.gov/drive/files/Version5

==============
# Build executable for Mac Delta based on llc270 iteration 42 optimized solution
 cvs -d :pserver:cvsanon:cvsanon@mitgcm.org:/u/gcmpack co -D "11/28/17" MITgcm_code
 git clone https://github.com/MITgcm-contrib/ecco_darwin.git

 cd MITgcm
 mkdir build run
 cd build
 module purge
 module load comp-intel/2016.2.181 mpi-sgi/mpt.2.14r19 hdf4/4.2.12 hdf5/1.8.18_mpt netcdf/4.4.1.1_mpt
 ../tools/genmake2 -of \
     ../tools/build_options/linux_amd64_ifort+mpi_ice_nas \
    -mo ../../ecco_darwin/Mac_Delta270/code
 make depend
 make -j 16

==============
# Instructions for running Mackenzie Delta regional model based on LLC4320
cd ..
mkdir run
cd run
mkdir diags
ln -sf ../build/mitgcmuv .
ln -sf /nobackup/hzhang1/forcing/era_xx .
ln -sf /nobackup/hzhang1/pub/Mac_Delta270/run_template/* .
cp ../../ecco_darwin/Mac_Delta270/input/* .
qsub job_Mmac270_Bro

