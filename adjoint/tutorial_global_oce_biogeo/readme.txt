# Instructions for building and running adjoint of tutorial_global_oce_biogeo
# adapted from MITgcm/verification/tutorial_global_oce_biogeo 
# using TAF
# 2 cases: 1) total CO2 uptake over all months
#	   2) just  CO2 uptake in last month

==============
# 1. Get code

git clone https://github.com/MITgcm-contrib/ecco_darwin.git
git clone https://github.com/MITgcm/MITgcm.git

==============
# 2. Build executable
cd MITgcm
mkdir build run
cd build

MOD='../../ecco_darwin/dic/tutorial_global_oce_biogeo'
module purge
module load comp-intel mpi-hpe hdf4/4.2.12 hdf5/1.8.18_mpt netcdf/4.4.1.1_mpt

#to build case 1)  
../tools/genmake2 -of ../tools/build_options/linux_amd64_ifort+mpi_ice_nas \
  -mo $MOD/code_ad -mpi
#to build case 2)  
../tools/genmake2 -of ../tools/build_options/linux_amd64_ifort+mpi_ice_nas \
  -mo "$MOD/code_ad2 $MOD/code_ad" -mpi

make depend
make -j16 adall

==============
# 3. Instructions for running simulation (12 months)

cd ../run
ln -sf ../build/mitgcmuv_ad .
ln -sf ../verification/tutorial_global_oce_biogeo/input/*.bin .
ln -sf ../verification/tutorial_global_oce_biogeo/input/pick* .
ln -sf ../verification/isomip/input_ad/ones_64b.bin .
cp $MOD/input_ad/* .
qsub job_dic_adj

