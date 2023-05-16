# v05 1deg Darwin3 simulation based on ECCOV4r5 set-up
#code base: c68g
#two versions:
#1 code_v4r5 + input_v4r5: 	including ctrl/smooth, similar to llc270
#2 code_v4r5_v2 input_v4r5_v2:	w/o ctrl/smooth, similar to v4r4

#version 1
# ========
# 1. Get code
git clone https://github.com/MITgcm/MITgcm.git -b checkpoint68g
git clone --depth 1 https://github.com/MITgcm-contrib/ecco_darwin.git
cd MITgcm

# ================
# 2. Build executable
# Prerequisite: 1. Get code
mkdir build run
cd build
rm *
module load comp-intel/2020.4.304 mpi-hpe/mpt.2.25 hdf4/4.2.12 hdf5/1.8.18_mpt netcdf/4.4.1.1_mpt
MOD="../../ecco_darwin/v05/1deg"
../tools/genmake2 -of ../tools/build_options/linux_amd64_ifort+mpi_ice_nas \
		  -mo ${MOD}/code_v4r5 -mpi
make depend
make -j 16

==============
# 3. Instructions for running simulation (1992-2019 period)

cd ../run
rm -rf *
mkdir -p diags
ln -sf ../build/mitgcmuv .

INPUTDIR='/nobackup/hzhang1/pub/Release5'

ln -s ${INPUTDIR}/input_bin/* .
ln -s ${INPUTDIR}/input_forcing/* .
cp ${MOD}/input_v4r5/* .

qsub job_v4r5

#version 2
# ========
# 1. Get code
git clone https://github.com/MITgcm/MITgcm.git -b checkpoint68g
git clone --depth 1 https://github.com/MITgcm-contrib/ecco_darwin.git
cd MITgcm

# ================
# 2. Build executable
# Prerequisite: 1. Get code
mkdir build2 run2
cd build2
rm *
module load comp-intel/2020.4.304 mpi-hpe/mpt.2.25 hdf4/4.2.12 hdf5/1.8.18_mpt netcdf/4.4.1.1_mpt
MOD="../../ecco_darwin/v05/1deg"
../tools/genmake2 -of ../tools/build_options/linux_amd64_ifort+mpi_ice_nas \
		  -mo ${MOD}/code_v4r5_v2 -mpi
make depend
make -j 16

==============
# 3. Instructions for running simulation (1992-2019 period)

cd ../run2
rm -rf *
mkdir -p diags
ln -sf ../build2/mitgcmuv .

INPUTDIR='/nobackup/hzhang1/pub/Release5'

ln -s ${INPUTDIR}/input_bin/* .
ln -s ${INPUTDIR}/input_forcing/TBADJ .
cp ${MOD}/input_v4r5_v2/* .

qsub job_v4r5
