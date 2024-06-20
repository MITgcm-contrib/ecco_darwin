# ECCOV4r5 set-up
#based on https://github.com/MITgcm-contrib/llc_hires/blob/master/llc_90/ecco_v4r5/readme_v4r5_68y.txt
#code base: c68y


# ========
# 1. Get code
git clone https://github.com/MITgcm-contrib/llc_hires.git
git clone https://github.com/MITgcm/MITgcm.git
cd MITgcm
git checkout checkpoint68y

# ================
# 2. Build executable
# Prerequisite: 1. Get code
mkdir build run
cd build
rm *
module load comp-intel mpi-hpe hdf4 hdf5/1.8.18_mpt netcdf/4.4.1.1_mpt
MOD="../../ecco_darwin/v06/1deg"
../tools/genmake2 -of ../tools/build_options/linux_amd64_ifort+mpi_ice_nas \
		  -mo ${MOD}/code_v4r5_v2 -mpi
make depend
make -j 16

==============
# 3. Instructions for running simulation (1992-2024 period)

cd ../run
rm -rf *
mkdir -p diags
ln -sf ../build/mitgcmuv .

INPUTDIR='/nobackup/hzhang1/pub/Release5'

ln -s ${INPUTDIR}/input_bin/* .
ln -s ${INPUTDIR}/TBADJ .
cp ${MOD}/input_v4r5_v2/* .

qsub job_v4r5
