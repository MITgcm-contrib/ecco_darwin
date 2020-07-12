# 3deg darwin3 verification experiment with volume, salt, salinity, DIC, and FeT budget
# diagnostics, initially based on MITgcm/verification/tutorial_global_oce_biogeo

# ========
# 1. Get code
 git clone --depth 1 https://github.com/darwinproject/darwin3
 git clone --depth 1 https://github.com/MITgcm-contrib/ecco_darwin.git
 cd darwin3/
 mkdir build run
 cd build

# ================
# 2. Build executable
#    Prerequisite: 1. Get code
 ../tools/genmake2 -ieee -mo \
 '../../ecco_darwin/v05/v5_3deg/code ../../ecco_darwin/v05/v5_llc270/code_darwin ../../ecco_darwin/v05/v5_llc270/code'
 make depend
 make -j 8

# ======================
# 3. Run verification setup
#    Prerequisite: 2. Build executable
 cd ../run
 ln -sf ../build/mitgcmuv .
 cp ../../ecco_darwin/v05/v5_3deg/input/data* .
 cp ../../ecco_darwin/v05/v5_3deg/input_3deg/* .
 cp ../../ecco_darwin/v05/v5_3deg/input_darwin/data* .
 ln -sf ../../ecco_darwin/v05/v5_3deg/data_darwin/* .
 rm data.exch2
 mkdir diags diags/3hourly diags/daily diags/monthly diags/budget
 ./mitgcmuv > output.txt
# Compare to verification output
 diff <(grep %MON output.txt) <(grep %MON ../../ecco_darwin/v05/v5_3deg/results/output_3deg.txt)

# ============================
# 4. Build and run MPI executable
#    Prerequisite: 1. Get code
 cd ../build
 rm *
 cp ../../ecco_darwin/v05/v5_3deg/code/SIZE.h_mpi SIZE.h
 ../tools/genmake2 -mpi -ieee -mo \
 '../../ecco_darwin/v05/v5_3deg/code ../../ecco_darwin/v05/v5_3deg/code_darwin ../../ecco_darwin/v05/v5_3deg/code'
 make depend
 make -j 8
 cd ../run
 ln -sf ../build/mitgcmuv .
 cp ../../ecco_darwin/v05/v5_3deg/input/data* .
 cp ../../ecco_darwin/v05/v5_3deg/input_3deg/* .
 cp ../../ecco_darwin/v05/v5_3deg/input_darwin/data* .
 Ln -sf ../../ecco_darwin/v05/v5_3deg/data_darwin/* .
 mkdir diags diags/3hourly diags/daily diags/monthly diags/budget
 mv data_mpi data
 rm data.exch2
 mpirun -np 8 ./mitgcmuv &
# Monitor run
 tail -f STDOUT.0000 | grep advcfl_W

# ============================
# 5. MATLAB code for computing volume, salt, salinity, DIC, and Fe budgets
#    Prerequisite: 4. Build and run MPI executable
#    Can be executed as soon as 3 or more months of output are available
 cd ../../ecco_darwin/v05/v5_3deg/matlab
# start MATLAB
# if using gcmfaces: *budget_v4_3deg_with_gcmfaces.m 
