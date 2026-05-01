# Instructions for building and running a physics-only
# RedSea/kaust_v1s3 configuration on shaheen at KAUST

==============
# 1. Get code
  git clone git@github.com:MITgcm-contrib/ecco_darwin.git
  git clone https://github.com/MITgcm/MITgcm
  cd MITgcm
  git checkout checkpoint67u
  mkdir build run

==============
# 2. Build executable
  cd build
  rm *
  ../tools/genmake2 -mods ../../ecco_darwin/regions/RedSea/kaust_v1s3/code/  -of=../../ecco_darwin/regions/RedSea/kaust_v1s3/shaheen_build_options -mpi -make=gmake -rootdir=..
  make depend
  make -j 16

==============
# 3. Instructions for running simulation
  cd ../run
  ln -sf ../build/mitgcmuv .
  ln -sf [INPUT_FILE_DIR]/* .
  mkdir diag2D diag3D
  cp ../../ecco_darwin/regions/RedSea/kaust_v1s3/input/* .
  mpirun -np 1649 ./mitgcmuv
