# Verification experiment with volume, salt, salinity, and DIC budget
# diagnostics, initially based on MITgcm/verification/tutorial_global_oce_biogeo

# ========
# 1. Get code
 cvs -d :pserver:cvsanon:cvsanon@mitgcm.org:/u/gcmpack co -D "11/28/17" MITgcm_code
 cvs co MITgcm_contrib/ecco_darwin/v4_3deg/data
 cvs co MITgcm_contrib/ecco_darwin/v4_llc270_JAMES_budget
 cvs -d :pserver:cvsanon:cvsanon@mitgcm.org:/u/gcmpack co -D "03/22/18" MITgcm_contrib/darwin/pkg/darwin
 cd MITgcm/pkg
 ln -sf ../../MITgcm_contrib/darwin/pkg/darwin .
 cd ..
 mkdir build run

# ================
# 2. Build executable
#    Prerequisite: 1. Get code
 cd build
 ../tools/genmake2 -ieee -mo \
  '../../MITgcm_contrib/ecco_darwin/v4_llc270_JAMES_budget/code_3deg ../../MITgcm_contrib/ecco_darwin/v4_llc270_JAMES_budget/code_darwin ../../MITgcm_contrib/ecco_darwin/v4_llc270_JAMES_budget/code'
 make depend
 make -j 8

# ======================
# 3. Run verification setup
#    Prerequisite: 2. Build executable
 cd ../run
 ln -sf ../build/mitgcmuv .
 cp ../../MITgcm_contrib/ecco_darwin/v4_llc270_JAMES_budget/input/data* .
 cp ../../MITgcm_contrib/ecco_darwin/v4_llc270_JAMES_budget/input_darwin/data* .
 cp ../../MITgcm_contrib/ecco_darwin/v4_llc270_JAMES_budget/input_3deg/*data* .
 ln -sf ../../MITgcm_contrib/ecco_darwin/v4_3deg/data/* .
 mkdir diags
 ./mitgcmuv > output.txt
# Compare to verification output
 diff output.txt ../../MITgcm_contrib/ecco_darwin/v4_llc270_JAMES_budget/results/output_3deg.txt

# ============================
# 4. Build and run MPI executable
#    Prerequisite: 1. Get code
 cd build
 rm *
 cp ../../MITgcm_contrib/ecco_darwin/v4_llc270_JAMES_budget/code/SIZE.h_mpi SIZE.h
 ../tools/genmake2 -mpi -mo \
  '../../MITgcm_contrib/ecco_darwin/v4_llc270_JAMES_budget/code ../../MITgcm_contrib/ecco_darwin/v4_llc270_JAMES_budget/code_darwin ../../MITgcm_contrib/ecco_darwin/v4_llc270_JAMES_budget/code'
 make depend
 make -j 8
 cd ../run
 mkdir diags
 ln -sf ../build/mitgcmuv .
 cp ../../MITgcm_contrib/ecco_darwin/v4_llc270_JAMES_budget/input/data* .
 cp ../../MITgcm_contrib/ecco_darwin/v4_llc270_JAMES_budget/input_darwin/data* .
 cp ../../MITgcm_contrib/ecco_darwin/v4_llc270_JAMES_budget/input/*data* .
 mv data_mpi data
 ln -sf ../../MITgcm_contrib/ecco_darwin/v4_3deg/data/* .
 mpirun -np 8 ./mitgcmuv &
# Monitor run
 tail -f STDOUT.0000 | grep advcfl_W

# ============================
# 5. MATLAB code for computing volume, salt, salinity, and DIC budgets
#    Prerequisite: 4. Build and run MPI executable
#    Can be executed as soon as 3 or more months of output are available
 cd ../../MITgcm_contrib/ecco_darwin/v4_llc270_JAMES_budget/matlab
# start matlab
# if using gcmfaces: budget_v4_3deg_with_gcmfaces
# if not using gcmfaces: budget_v4_3deg_without_gcmfaces
