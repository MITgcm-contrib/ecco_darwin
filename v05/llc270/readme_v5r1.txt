# Instructions for building and running ECCO-Darwin v05 with Darwin 3
# based on ECCO-v5r1

#This solution is documented in:
#Carroll, D., Menemenlis, D., Dutkiewicz, S., Lauderdale, J. M., Adkins, J. F., Bowman, K. W., et al. (2022). 
#Attribution of space-time variability in global-ocean dissolved inorganic carbon. Global Biogeochemical Cycles, 
#36, e2021GB007162. https://doi.org/10.1029/2021GB007162

==============
# 1. Get code

git clone --depth 1 https://github.com/MITgcm-contrib/ecco_darwin.git
git clone https://github.com/darwinproject/darwin3
cd darwin3
git checkout darwin_ckpt68g
mkdir build2 run2
cd build2

==============
# 2. Build executable

module purge
module load comp-intel mpi-hpe/mpt hdf4/4.2.12 hdf5/1.8.18_mpt netcdf/4.4.1.1_mpt python3/3.9.12
MOD="../../ecco_darwin/v05/llc270"
cp $MOD/code_darwin/packages.conf .
sed -i '/ctrl/d;/smooth/d' packages.conf
sed -i '$ashelfice'        packages.conf
../tools/genmake2 -of ../tools/build_options/linux_amd64_ifort+mpi_ice_nas \
  -mo "$MOD/code_darwin $MOD/code_v5r1" -mpi
make depend
make -j 16

==============
# 3. Instructions for running simulation (1992-2025 period)
INPUT="/nobackup/hzhang1/pub/llc270_FWD/v5r1"

cd ../run2
ln -sf ../build2/mitgcmuv .
ln -sf $INPUT/input_bin/* .
ln -sf $INPUT/input_darwin_bin/* .
ln -sf /nobackup/dcarrol2/forcing/apCO2/NOAA_MBL/* .
cp $MOD/input_v5r1/* .
cp $MOD/input2/data.{darwin,diagnostics,gchem,pkg,ptracers,traits}* .
sed -i "/nIter0=2/a \ pickupSuff='0000000002',"  data
mkdir diags diags/3hourly diags/daily diags/monthly diags/budget
# modify job_v5r1 as needed between "devel" and "long"
qsub job_v5r1

