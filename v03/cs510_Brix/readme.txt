# Build executable for ECCO-Darwin version 3
# with circa-2011 MITgcm code, equivalent to
# /nobackup/hbrix/MITgcm_110502/build_darwin_p6/mitgcmuv
# Note: this version has different ICs from V2 and
# uses a linear piston velocity formulation. The improved
# V3 air-sea CO2 flux uptake, both globally and in the Southern Ocean,
# results primarily from changing the piston velocity formulation
# from quadratic to linear. By implementing various bug fixes and adjusted
# ICs, we were able to return to the quadratic formulation in V4. 

cvs -d :pserver:cvsanon:cvsanon@mitgcm.org:/u/gcmpack co -D "2011-05-12 10:30:47" MITgcm_code
cvs -d :pserver:cvsanon:cvsanon@mitgcm.org:/u/gcmpack co -D "2011-05-12 10:30:47" MITgcm_contrib/darwin/pkg/darwin
git clone --depth 1 https://github.com/MITgcm-contrib/ecco_darwin.git
cd MITgcm/pkg
ln -sf ../../MITgcm_contrib/darwin/pkg/darwin .
cd ..
mkdir build
cd build
module purge
module load comp-intel/2020.4.304 mpi-hpe hdf4/4.2.12 hdf5/1.8.18_mpt netcdf/4.4.1.1_mpt
../tools/genmake2 -of \
    ../../ecco_darwin/v03/cs510_brix/code/linux_amd64_ifort+mpi_ice_nas \
    -mo ../../ecco_darwin/v03/cs510_brix/code
make depend
make -j 16

==============
# Instructions for running Version 3 (cg1) to 2009-2013.

cd ..
mkdir run
cd run
ln -sf ../build/mitgcmuv .
cp ../../ecco_darwin/v03/cs510_brix/input/* .
ln -sf ~dmenemen/CMS/run_template_cg1/darwin* .
ln -sf ~dmenemen/CMS/run_template_cg1/DIFFKR_2_20_1_lat6070_cube81 .
ln -sf /nobackup/hzhang1/cs510/run_template/GEBCO_510x6x510_ver06_dig.bin .
ln -sf /nobackup/hzhang1/forcing/jra25/jra25_dlw* .
ln -sf /nobackup/hzhang1/cs510/ICBC_2009_iter26/jra25_xx_* .
ln -sf /nobackup/hbrix/ICBC/pco2a_blended_* .
ln -sf pco2a_blended_2012 pco2a_blended_2013
ln -sf pco2a_blended_2012 pco2a_blended_2014
ln -sf /nobackup/hzhang1/cs510/ICBC_2009_iter26/pickup.0000000001.data_xx pickup.0000078912.data
ln -sf /nobackup/hzhang1/cs510/ICBC_2009/pickup.0000000001.meta pickup.0000078912.meta
ln -sf /nobackup/hzhang1/cs510/ICBC_2009/pickup_seaice.0000000001.data pickup_seaice.0000078912.data
ln -sf /nobackup/hzhang1/cs510/ICBC_2009/pickup_seaice.0000000001.meta pickup_seaice.0000078912.meta
ln -sf /nobackup/hbrix/ICBC/pickup_greensf/pickup_dic.gf_cg1.0000078912.data pickup_dic.0000078912.data
ln -sf /nobackup/hbrix/ICBC/pickup_greensf/pickup_dic.gf_cg1.0000078912.meta pickup_dic.0000078912.meta
ln -sf /nobackup/hbrix/ICBC/pickup_greensf/pickup_ptracers.gf_cg1.0000078912.data pickup_ptracers.0000078912.data
ln -sf /nobackup/hbrix/ICBC/pickup_greensf/pickup_ptracers.gf_cg1.0000078912.meta pickup_ptracers.0000078912.meta
ln -sf /nobackup/dmenemen/forcing/runoff/runoff-360x180x12.bin .
ln -sf /nobackup/hzhang1/cs510/run_template/tile00* .

# modify run_darwin_450 as needed
qsub run_darwin_450
