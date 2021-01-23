# ========
# Mackenzie Delta regional setup based on LLC270 with Darwin3
# ========

# ==============
# 1. Get code
# ==============
git clone --depth 1 https://github.com/darwinproject/darwin3
svn checkout https://github.com/MITgcm-contrib/ecco_darwin/trunk/regions/mac_delta/llc270 Mac270

# ================
# 2. Build executable for Mac Delta based on llc270 iteration 42 optimized solution
#    Prerequisite: 1. Get code
# ==============
cd darwin3
mkdir build run
cd build
module purge
module load comp-intel/2016.2.181 mpi-sgi/mpt.2.14r19 hdf4/4.2.12 hdf5/1.8.18_mpt netcdf/4.4.1.1_mpt
mv  ../../Mac270/code/packages.conf ../../Mac270/code/packages.conf_org
ln -s ../../Mac270/code_darwin/packages.conf ../../Mac270/code
../tools/genmake2 -of ../tools/build_options/linux_amd64_ifort+mpi_ice_nas \
-mo '../../Mac270/code ../../Mac270/code_darwin'
make depend
make -j 16

# ================
# 3. Run the setup
#    Prerequisite: 2. Build executable
#    Default running doesn't output diagnostic files
#    To enable diagnostics outputs follow instructions at the end of the readme
# ================
cd ../run
ln -sf ../build/mitgcmuv .
ln -sf /nobackup/hzhang1/forcing/era_xx .
ln -sf /nobackup/hzhang1/pub/Mac_Delta270/run_template/* .
cp ../../Mac270/input_darwin/* .
mkdir diags diags/3hourly diags/daily diags/monthly diags/budget
qsub job_Mac270_Bro
