# Verification experiment, initially based on
# MITgcm/verification/tutorial_global_oce_biogeo

# ========
# Get code
 cvs -d :pserver:cvsanon:cvsanon@mitgcm.org:/u/gcmpack co -D "11/28/17" MITgcm_code
 cvs -d :pserver:cvsanon:cvsanon@mitgcm.org:/u/gcmpack co -D "11/28/17" MITgcm/verification/tutorial_global_oce_biogeo
 cvs co MITgcm_contrib/ecco_darwin/v4_3deg
 cvs co MITgcm_contrib/ecco_darwin/v4_llc270
 cvs -d :pserver:cvsanon:cvsanon@mitgcm.org:/u/gcmpack co -D "03/22/18" MITgcm_contrib/darwin/pkg/darwin
 cd MITgcm/pkg
 ln -sf ../../MITgcm_contrib/darwin/pkg/darwin .
 cd ..

# ================
# Build executable
 mkdir build
 cd build
 ../tools/genmake2 -ieee -mo \
  '../../MITgcm_contrib/ecco_darwin/v4_3deg/code ../../MITgcm_contrib/ecco_darwin/v4_llc270/code_darwin ../../MITgcm_contrib/ecco_darwin/v4_llc270/code'
 make depend
 make -j 8
 cd ..

# ======================
# Run verification setup
 mkdir run
 cd run
 ln -sf ../build/mitgcmuv .
 cp ../../MITgcm_contrib/ecco_darwin/v4_llc270/input/data* .
 cp ../../MITgcm_contrib/ecco_darwin/v4_llc270/input_darwin/data* .
 cp ../../MITgcm_contrib/ecco_darwin/v4_3deg/input/*data* .
 ln -sf ../verification/tutorial_global_oce_biogeo/input/bathy.bin .
 ln -sf ../../MITgcm_contrib/ecco_darwin/v4_3deg/data/* .
 ./mitgcmuv > output.txt

# ==============================
# Compare to verification output
 diff output.txt ../../MITgcm_contrib/ecco_darwin/v4_3deg/results/output.txt


# ============================
# Build and run MPI executable
 cd build
 rm *
 cp ../../MITgcm_contrib/ecco_darwin/v4_3deg/code/SIZE.h_mpi SIZE.h
 ../tools/genmake2 -mpi -mo \
  '../../MITgcm_contrib/ecco_darwin/v4_3deg/code ../../MITgcm_contrib/ecco_darwin/v4_llc270/code_darwin ../../MITgcm_contrib/ecco_darwin/v4_llc270/code'
 make depend
 make -j 8
 cd ../run
 mpirun -np 8 ./mitgcmuv
