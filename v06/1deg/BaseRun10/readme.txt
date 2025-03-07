# v06 1deg Darwin3 simulation based on ECCOV4r5 set-up
#based on https://github.com/MITgcm-contrib/llc_hires/blob/master/llc_90/ecco_v4r5/readme_v4r5_68y.txt
#https://github.com/jahn/ecco_darwin/tree/v06-1deg-oasim-v2
#code base: c68y


# ========
# 1. Get code
git clone https://github.com/darwinproject/darwin3
git clone --depth 1 https://github.com/MITgcm-contrib/ecco_darwin.git
cd darwin3/pkg/darwin
git checkout darwin_ckpt68y
cd ../../

# ================
# 2. Build executable
# Prerequisite: 1. Get code
mkdir build run
cd build
rm *
module load comp-intel mpi-hpe hdf4/4.2.12 hdf5/1.8.18_mpt netcdf/4.4.1.1_mpt python3/3.9.12
MOD="../../ecco_darwin/v06/1deg/BaseRun10"
../tools/genmake2 -of ../tools/build_options/linux_amd64_ifort+mpi_ice_nas \
	-mo ${MOD}/code -mpi
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
cp ${MOD}/input/* .
ln -sf /nobackup/dcarrol2/forcing/apCO2/NOAA_MBL/* .
ln -sf /nobackup/ojahn/ecco_darwin/v06/llc270/data_darwin/* .
ln -sf /nobackup/ojahn/forcing/oasim .
ln -sf /nobackup/dcarrol2/pub/1deg/v06/V4r5/* .
ln -sf /nobackup/rsavelli/LOAC/ECCO_V4r5/bgc_runoff/* .
mv pickup_ptracers.0000000001.data pickup_ptracers.0000000002.data
mv pickup_ptracers.0000000001.meta pickup_ptracers.0000000002.meta
mkdir diags diags/3hourly diags/daily diags/monthly diags/budget
mkdir diags/monthly/IOPS diags/monthly/PAR diags/monthly/RRS 

# qsub job_ECCOv4r5
