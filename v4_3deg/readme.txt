# Verification experiment, initially based on
# MITgcm/verification/tutorial_global_oce_biogeo

# ========
# Get code
 cvs -d :pserver:cvsanon:cvsanon@mitgcm.org:/u/gcmpack co -D "11/28/17" MITgcm_code
 cvs -d :pserver:cvsanon:cvsanon@mitgcm.org:/u/gcmpack co MITgcm_contrib/ecco_darwin/v4_3deg
 cvs -d :pserver:cvsanon:cvsanon@mitgcm.org:/u/gcmpack co -D "03/22/18" MITgcm_contrib/darwin/pkg/darwin
 cvs -d :pserver:cvsanon:cvsanon@mitgcm.org:/u/gcmpack co -D "11/28/17" MITgcm/verification/tutorial_global_oce_biogeo
 cd MITgcm/pkg
 ln -sf ../../MITgcm_contrib/darwin/pkg/darwin .
 cd ..

# ================
# Build executable
 mkdir build
 cd build
 ../tools/genmake2 -mo ../../MITgcm_contrib/ecco_darwin/v4_3deg/code
 make depend
 make -j 8
 cd ..

# ======================
# Run verification setup
 mkdir run
 cd run
 ln -sf ../build/mitgcmuv .
 cp ../../MITgcm_contrib/ecco_darwin/v4_3deg/input/*data* .
 ln -sf ../verification/tutorial_global_oce_biogeo/input/*bin .
 rm lev_* shi_* tren_t*
 ln -sf ../../MITgcm_contrib/ecco_darwin/v4_3deg/data/runof* .
 ln -sf ../../MITgcm_contrib/ecco_darwin/v4_3deg/data/*_2000 .
 ln -sf ../../MITgcm_contrib/ecco_darwin/v4_3deg/data/*.0005184000 .
 ln -sf ../../MITgcm_contrib/ecco_darwin/v4_3deg/data/ptracers* .
 ln -sf ../../MITgcm_contrib/ecco_darwin/v4_3deg/data/3deg* .
 ./mitgcmuv > output.txt

# ==============================
# Compare to verification output
 diff output.txt ../../MITgcm_contrib/ecco_darwin/v4_3deg/results/output.txt
