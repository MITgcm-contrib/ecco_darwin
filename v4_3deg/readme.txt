# Verification experiment, initially based on
# MITgcm/verification/tutorial_global_oce_biogeo

# ========
# Get code
 cvs co MITgcm_contrib/ecco_darwin/v4_3deg
 git clone git@github.com:MITgcm/MITgcm.git

# ================
# Build executable
 cd MITgcm
 mkdir build
 cd build
 ../tools/genmake2 -mo ../../MITgcm_contrib/ecco_darwin/v4_3deg/code
 make depend
 make -j 8

# ======================
# Run verification setup
 cd ..
 mkdir run
 cd run
 ln -sf ../build/mitgcmuv .
 cp ../../MITgcm_contrib/ecco_darwin/v4_3deg/input/* .
 ln -sf ../verification/tutorial_global_oce_biogeo/input/*.bin .
 ln -sf ../verification/tutorial_global_oce_biogeo/input/pickup* .
 ./mitgcmuv > output.txt

# ==============================
# Compare to verification output
 diff output.txt ../../MITgcm_contrib/ecco_darwin/v4_3deg/results/output.txt
