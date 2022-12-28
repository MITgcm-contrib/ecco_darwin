# v05 3deg Darwin3 setup for Carbon Dioxide Removal (CDR) simulations  

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

 cp ../../ecco_darwin/v05/llc270_CDR/input/DIAGNOSTICS_SIZE.h .
../tools/genmake2 -ieee -mo '../../ecco_darwin/v05/llc270_CDR/code_darwin ../../ecco_darwin/v05/3deg/code ../../ecco_darwin/v05/llc270/code_darwin ../../ecco_darwin/v05/llc270/code'
 make depend
 make -j 8

# ======================
# 3. Run verification setup
#    Prerequisite: 2. Build executable
 cd ../run
 ln -sf ../build/mitgcmuv .
 cp ../../ecco_darwin/v05/llc270/input/data* .
 cp ../../ecco_darwin/v05/3deg/input/* .
 ln -sf ../../ecco_darwin/v05/3deg/data_darwin/* .
 rm data.ctrl data.exch2 data.smooth
 cp ../../ecco_darwin/v05/llc270_CDR/input/* .     (modify data.darwin)
 mkdir diags diags/3hourly diags/daily diags/monthly diags/budget
 ./mitgcmuv > output.txt
# Compare to verification output
 diff <(grep %MON output.txt) <(grep %MON ../../ecco_darwin/v05/3deg/results/output_3deg.txt)

# ============================
# 4. Build and run MPI executable
#    Prerequisite: 1. Get code
 cd ../build
 rm *
 module load comp-intel mpi-hpe hdf4/4.2.12 hdf5/1.8.18_mpt netcdf/4.4.1.1_mpt
 cp ../../ecco_darwin/v05/llc270_CDR/input/SIZE_16.h SIZE.h
 cp ../../ecco_darwin/v05/llc270_CDR/input/DIAGNOSTICS_SIZE.h .
 ../tools/genmake2 -of ../../ecco_darwin/v05/llc270/code/linux_amd64_ifort+mpi_ice_nas -mo '../../ecco_darwin/v05/llc270_CDR/code_darwin ../../ecco_darwin/v05/3deg/code ../../ecco_darwin/v05/llc270/code_darwin ../../ecco_darwin/v05/llc270/code' -mpi

 make depend
 make -j 16

 cd ../run
 rm -rf *
 ln -sf ../build/mitgcmuv .
 cp ../../ecco_darwin/v05/llc270/input/data* .
 cp ../../ecco_darwin/v05/3deg/input/* .
 ln -sf ../../ecco_darwin/v05/3deg/data_darwin/* .
 cp ../../ecco_darwin/v05/llc270_CDR/input/* .   (modify data.darwin)

 mkdir diags diags/3hourly diags/daily diags/monthly diags/budget
# mv data_mpi data      ... modify data (time step + nITER0 + nTimeSteps)
 rm data.ctrl data.exch2 data.smooth
# from model_input copy restarts and directory to the right experiment

# modify job_3deg_darwin3 as needed
 qsub job_3deg_darwin3
 
