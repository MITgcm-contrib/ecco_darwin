# v05 3deg Darwin3 setup for Carbon Dioxide Removal (CDR) simulations  

# ========
# 1. Get code
 git clone https://github.com/darwinproject/darwin3
 cd darwin3/pkg/darwin
 git checkout 24885b71
 cd ../../../
 git clone --depth 1 https://github.com/MITgcm-contrib/ecco_darwin.git
 cd darwin3
 mkdir build_cdr run_cdr
 cd build_cdr

# ================
# 2. Build executable
#    Prerequisite: 1. Get code


../tools/genmake2 -ieee -mo '../../ecco_darwin/v05/llc270_CDR/code_darwin ../../ecco_darwin/v05/3deg/code ../../ecco_darwin/v05/llc270/code_darwin ../../ecco_darwin/v05/llc270/code'
 make depend
 make -j 8

# ======================
# 3. Run experiment (verification setup)
#    Prerequisite: 2. Build executable
 cd ../run_cdr
 ln -sf ../build_cdr/mitgcmuv .
 cp ../../ecco_darwin/v05/llc270/input/data* .
 cp ../../ecco_darwin/v05/3deg/input/* .
 ln -sf ../../ecco_darwin/v05/3deg/data_darwin/* .
 rm data.ctrl data.exch2 data.smooth
 cp ../../ecco_darwin/v05/llc270_CDR/input/* .     

# Modify "data.darwin": set the runofffiles, start dates and periods for ALK, Si and Fe enhancments
# Modify data: specify whether you want to start from initial conditions or from the pickup files 
#             currently simulation starts from "00000859080" pickups which are also needed
 
mkdir diags diags/3hourly diags/daily diags/monthly diags/budget
 ./mitgcmuv > output.txt
# Compare to verification output (this would make sense only if yu start with initial files)
# diff <(grep %MON output.txt) <(grep %MON ../../ecco_darwin/v05/3deg/results/output_3deg.txt)

# ============================
# 4. Build and run MPI executable
#    Prerequisite: 1. Get code
 cd ../build_cdr
 rm *
 module load comp-intel mpi-hpe hdf4/4.2.12 hdf5/1.8.18_mpt netcdf/4.4.1.1_mpt
 cp ../../ecco_darwin/v05/llc270_CDR/input/SIZE_16.h SIZE.h
 ../tools/genmake2 -of ../../ecco_darwin/v05/llc270/code/linux_amd64_ifort+mpi_ice_nas -mo '../../ecco_darwin/v05/llc270_CDR/code_darwin ../../ecco_darwin/v05/3deg/code ../../ecco_darwin/v05/llc270/code_darwin ../../ecco_darwin/v05/llc270/code' -mpi

 make depend
 make -j 16

 cd ../run_cdr
 rm -rf *
 ln -sf ../build_cdr/mitgcmuv .
 cp ../../ecco_darwin/v05/llc270/input/data* .
 cp ../../ecco_darwin/v05/3deg/input/* .
 ln -sf ../../ecco_darwin/v05/3deg/data_darwin/* .
 cp ../../ecco_darwin/v05/llc270_CDR/input/* .  
# mv data_mpi data (not really necessary) rather modify "monitorFreq = 43200.0," 

# Modify "data.darwin": set the runofffiles, start dates and periods for ALK, Si and Fe enhancments
# Modify data: specify whether you want to start from initial conditions or from the pickup files 
#             currently simulation starts from "00000859080" pickups which are also needed
 
 mkdir diags diags/3hourly diags/daily diags/monthly diags/budget
 rm data.ctrl data.exch2 data.smooth


# modify job_3deg_darwin3 as needed

 qsub job_3deg_darwin3
 
