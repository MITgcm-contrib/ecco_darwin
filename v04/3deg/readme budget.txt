# Verification experiment, initially based on
# MITgcm/verification/tutorial_global_oce_biogeo 
# with volume, salt, salinity, and DIC budget

# ========
# 1. Get code
 cvs -d :pserver:cvsanon:cvsanon@mitgcm.org:/u/gcmpack co -D "11/28/17" MITgcm_code
 cvs -d :pserver:cvsanon:cvsanon@mitgcm.org:/u/gcmpack co -D "03/22/18" MITgcm_contrib/darwin/pkg/darwin
 git clone --depth 1 https://github.com/MITgcm-contrib/ecco_darwin.git
 cd MITgcm/pkg
 ln -sf ../../MITgcm_contrib/darwin/pkg/darwin .
 cd ..
 mkdir build run

# ================
# 2. Build executable
#    Prerequisite: 1. Get code
 cd build
 ../tools/genmake2 -ieee -mo \
  '../../ecco_darwin/v04/3deg/code ../../ecco_darwin/v04/llc270_JAMES_budget/code_darwin ../../ecco_darwin/v04/llc270_JAMES_budget/code'
 make depend
 make -j 8

# ======================
# 3. Run verification setup
#    Prerequisite: 2. Build executable
 cd ../run
 ln -sf ../build/mitgcmuv .
 cp ../../ecco_darwin/v04/llc270_JAMES_paper/input/data* .
 cp ../../ecco_darwin/v04/3deg/input_budget/*data* .
 ln -sf ../../ecco_darwin/v04/3deg/data/* .
 mkdir diags
 ./mitgcmuv > output.txt

# ============================
# 4. Build and run MPI executable
#    Prerequisite: 1. Get code
 cd build
 rm *
 cp ../../ecco_darwin/v04/3deg/code/SIZE.h_mpi SIZE.h
 ../tools/genmake2 -mpi -mo \
  '../../ecco_darwin/v04/3deg/code ../../ecco_darwin/v04/llc270_JAMES_budget/code_darwin ../../ecco_darwin/v04/llc270_JAMES_budget/code'
 make depend
 make -j 8
 cd ../run
 ln -sf ../build/mitgcmuv .
 cp ../../ecco_darwin/v04/llc270_JAMES_paper/input/data* .
 cp ../../ecco_darwin/v04/3deg/input_budget/*data* .
 mv data_mpi data
 ln -sf ../../ecco_darwin/v04/3deg/data/* .
 mkdir diags
 mpirun -np 8 ./mitgcmuv &
# Monitor run
 tail -f STDOUT.0000 | grep advcfl_W

# ============================
# 5. MATLAB code for computing volume, salt, salinity, and DIC budgets
#    Prerequisite: 4. Build and run MPI executable
#    Can be executed as soon as 3 or more months of output are available
 cd ../../MITgcm_contrib/ecco_darwin/llc270_JAMES_budget/matlab
# start matlab
# if using gcmfaces: budget_3deg_with_gcmfaces
# if not using gcmfaces: budget_3deg_without_gcmfaces
