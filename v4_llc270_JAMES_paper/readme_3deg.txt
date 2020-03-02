# Verification experiment, initially based on
# MITgcm/verification/tutorial_global_oce_biogeo

# ========
# 1. Get code
 cvs -d :pserver:cvsanon:cvsanon@mitgcm.org:/u/gcmpack co -D "11/28/17" MITgcm_code
 cvs co MITgcm_contrib/ecco_darwin/v4_3deg/data
 cvs co MITgcm_contrib/ecco_darwin/v4_llc270_JAMES_paper
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
  '../../MITgcm_contrib/ecco_darwin/v4_llc270_JAMES_paper/code_3deg ../../MITgcm_contrib/ecco_darwin/v4_llc270_JAMES_paper/code_darwin ../../MITgcm_contrib/ecco_darwin/v4_llc270_JAMES_paper/code'
 make depend
 make -j 8

# ======================
# 3. Run verification setup
#    Prerequisite: 2. Build executable
 cd ../run
 ln -sf ../build/mitgcmuv .
 cp ../../MITgcm_contrib/ecco_darwin/v4_llc270_JAMES_paper/input/data* .
 cp ../../MITgcm_contrib/ecco_darwin/v4_llc270_JAMES_paper/input_darwin/data* .
 cp ../../MITgcm_contrib/ecco_darwin/v4_llc270_JAMES_paper/input_3deg/*data* .
 ln -sf ../../MITgcm_contrib/ecco_darwin/v4_3deg/data/* .
 ./mitgcmuv > output.txt
# Compare to verification output
 diff <(grep %MON output.txt) <(grep %MON ../../MITgcm_contrib/ecco_darwin/v4_llc270_JAMES_paper/results/output_3deg.txt)
