# Instructions for llc270 verification experiment
# Instructions are specific to pleiades and need to be adjusted for other platforms
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
 mkdir build run
 cd build
 module purge
 module load comp-intel/2020.4.304 mpi-hpe/mpt.2.25 hdf4/4.2.12 hdf5/1.8.18_mpt netcdf/4.4.1.1_mpt
 ../tools/genmake2 -ieee -of \
  ../../ecco_darwin/v04/llc270_JAMES_paper/code/linux_amd64_ifort+mpi_ice_nas -mo \
  '../../ecco_darwin/v04/llc270_JAMES_paper/code_darwin ../../ecco_darwin/v04/llc270_JAMES_paper/code'
 make depend
 make -j 16

==============
# Instructions for running ECCO-Darwin Version 4 for 1992-2018 period
 cd ../run
 ln -sf ../build/mitgcmuv .
 ln -sf /nobackupp2/dmenemen/public/llc_270/iter42/input/* .
 ln -sf /nobackupp2/dmenemen/public/llc_270/ecco_darwin_v4/input/darwin_forcing/* .
 ln -sf /nobackupp2/dmenemen/public/llc_270/ecco_darwin_v4/input/darwin_initial_conditions/pickup_ptracers_experiment_18.data pickup_ptracers.0000000001.data
 ln -sf /nobackupp2/dmenemen/public/llc_270/ecco_darwin_v4/input/darwin_initial_conditions/pickup_ptracers.0000000001.meta .
 ln -sf /nobackup/hzhang1/forcing/era_xx .
 cp ../../ecco_darwin/v04/llc270_JAMES_paper/input/* .
 cp ../../ecco_darwin/v04/llc270_JAMES_paper/input_darwin/* .
 mv data_test data
 qsub job_ECCO_darwin_test
# Compare to verification output
 diff <(grep %MON STDOUT.0000) <(grep %MON ../../ecco_darwin/v04/llc270_JAMES_paper/results/STDOUT.0000)
