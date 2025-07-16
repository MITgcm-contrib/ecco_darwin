# v06 1deg Darwin3 simulation based on ECCOV4r5 set-up
# 
#Arctic Ocean ecosystem:
# ========
# 1. Get code
git clone https://github.com/darwinproject/darwin3
git clone --depth 1 https://github.com/MITgcm-contrib/ecco_darwin.git
cd darwin3/pkg/darwin
git checkout backport_ckpt68y
cd ../../

# ================
# 2. Build executable
# Prerequisite: 1. Get code
mkdir build run
cd build
rm *
module load comp-intel mpi-hpe hdf4/4.2.12 hdf5/1.8.18_mpt netcdf/4.4.1.1_mpt python3/3.9.12
MOD="../../ecco_darwin/v06/1deg"
../tools/genmake2 -of ../tools/build_options/linux_amd64_ifort+mpi_ice_nas \
	-mo '../../ecco_darwin/v06/1deg/code_darwin_v4r5_v2_AO ../../ecco_darwin/v06/1deg/code_v4r5_v2' -mpi
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

rm data data.pkg data.diagnostics
cp ${MOD}/input_darwin_v4r5_v2_AO/* .
ln -sf /nobackup/dcarrol2/forcing/apCO2/NOAA_MBL/* .
ln -sf /nobackup/ojahn/ecco_darwin/v06/llc270/data_darwin/* .
ln -sf /nobackup/ojahn/forcing/oasim .
ln -sf /nobackup/dcarrol2/pub/1deg/v06/V4r5/* .
ln -sf /nobackup/rsavelli/LOAC/ECCO_V4r5/bgc_runoff/* .
######### TO CHANGE
mv pickup_ptracers_AO.0000000001.data pickup_ptracers.0000000002.data
mv pickup_ptracers_AO.0000000001.meta pickup_ptracers.0000000002.meta
######### TO CHANGE
mkdir diags diags/3hourly diags/daily diags/monthly diags/budget
mkdir diags/monthly/IOPS diags/monthly/PAR diags/monthly/RRS 

# qsub job_ECCOv4r5
