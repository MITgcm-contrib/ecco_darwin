# v05 3deg darwin3 verification experiment with volume, salt, salinity, DIC, and FeT budget
# diagnostics, initially based on MITgcm/verification/tutorial_global_oce_biogeo

# ========
# 1. Get code
 git clone https://github.com/darwinproject/darwin3
 cd darwin3/pkg/darwin
 git checkout 24885b71
 cd ../../../
 git clone --depth 1 https://github.com/MITgcm-contrib/ecco_darwin.git
 cd darwin3
 mkdir build run
 cd build

# ================
# 2. Build executable
#    Prerequisite: 1. Get code
 ../tools/genmake2 -ieee -mo \
 '../../ecco_darwin/v05/3deg/code ../../ecco_darwin/v05/llc270_sediment/code_darwin ../../ecco_darwin/v05/llc270/code'
 make depend
 make -j 8

# ======================
# 3. Run verification setup
#    Prerequisite: 2. Build executable
 cd ../run
 ln -sf ../build/mitgcmuv .
 cp ../../ecco_darwin/v05/llc270_sediment/input/data.diagnostics .
 cp ../../ecco_darwin/v05/3deg/input/* .
 ln -sf ../../ecco_darwin/v05/3deg/data_darwin/* .
 rm data.ctrl data.exch2 data.smooth
 mkdir diags diags/3hourly diags/daily diags/monthly diags/budget
 ./mitgcmuv > output.txt
# Compare to verification output
 diff <(grep %MON output.txt) <(grep %MON ../../ecco_darwin/v05/3deg/results/output_3deg.txt)

# ============================
# 4. Build and run MPI executable
#    Prerequisite: 1. Get code
 cd ../build
 rm *
 module load comp-intel/2020.4.304 mpi-hpe/mpt.2.25 hdf4/4.2.12 hdf5/1.8.18_mpt netcdf/4.4.1.1_mpt
 cp ../../ecco_darwin/v05/3deg/code/SIZE.h_mpi SIZE.h
 ../tools/genmake2 -of  ../../ecco_darwin/v05/llc270/code/linux_amd64_ifort+mpi_ice_nas  \
 -mo '../../ecco_darwin/v05/3deg/code ../../ecco_darwin/v05/llc270/code_darwin ../../ecco_darwin/v05/llc270/code' -mpi
  make depend
 make -j 8

 cd ../run
 rm -rf *
 ln -sf ../build/mitgcmuv .
 cp ../../ecco_darwin/v05/llc270_sediment/input/data* .
 cp ../../ecco_darwin/v05/3deg/input/* .
 ln -sf ../../ecco_darwin/v05/3deg/data_darwin/* .
 mkdir diags diags/3hourly diags/daily diags/monthly diags/budget
 mv data_mpi data
 rm data.ctrl data.exch2 data.smooth
# modify job_3deg_darwin3 as needed
 qsub job_3deg_darwin3
 
# ============================
# 5. MATLAB code for computing volume, salt, salinity, DIC, and Fe budgets
#    Prerequisite: 4. Build and run MPI executable
#    Can be executed as soon as 3 or more days of output are available
 cd ../../ecco_darwin/v05/3deg/matlab
# start MATLAB
# with gcmfaces use: *budget_with_gcmfaces*.m 
# without gcmfaces use: *budget_without_gcmfaces*.m
