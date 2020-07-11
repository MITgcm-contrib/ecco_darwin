# Verification experiment, initially based on
# MITgcm/verification/tutorial_global_oce_biogeo

# ========
# 1. Get code
 cvs -d :pserver:cvsanon:cvsanon@mitgcm.org:/u/gcmpack co -D "11/28/17" MITgcm_code
 git clone --depth 1 https://github.com/MITgcm-contrib/ecco_darwin.git
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
  '../../ecco_darwin/v04/v4_3deg/code ../../ecco_darwin/v04/v4_llc270_JAMES_paper/code_darwin ../../ecco_darwin/v04/v4_llc270_JAMES_paper/code'
 make depend
 make -j 8

# ======================
# 3. Run verification setup
#    Prerequisite: 2. Build executable
 cd ../run
 ln -sf ../build/mitgcmuv .
 cp ../../ecco_darwin/v04/v4_llc270_JAMES_paper/input/data* .
 cp ../../ecco_darwin/v04/v4_llc270_JAMES_paper/input_darwin/data* .
 cp ../../ecco_darwin/v04/v4_3deg/input/*data* .
 ln -sf ../../ecco_darwin/v04/v4_3deg/data/* .
 rm data.exch2
 ./mitgcmuv > output.txt
# Compare to verification output
 diff <(grep %MON output.txt) <(grep %MON ../../ecco_darwin/v04/v4_llc270_JAMES_paper/results/output_3deg.txt)
