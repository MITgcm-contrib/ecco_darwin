# Building physics-only RedSea/kaust_v1s3 on shaheen at KAUST

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
