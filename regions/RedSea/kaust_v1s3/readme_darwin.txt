# Instructions for building and running a RedSea/kaust_v1s3
# ECCO-Darwin configuration on shaheen at KAUST

WORKDIR=/scratch/wangy0m/ECCO_Darwin/

==============
# 1. Get code
  cd $WORKDIR
  git clone git@github.com:MITgcm-contrib/ecco_darwin.git
  git clone git@github.com:darwinproject/darwin3
  cd $WORKDIR/darwin3
  git checkout 24885b71
  mkdir build run

==============
# 2. Build executable
  cd $WORKDIR/darwin3/build
  rm *
  module load python
  ../tools/genmake2 -of=../../ecco_darwin/regions/RedSea/kaust_v1s3/shaheen_build_options -mpi -make=gmake \
   -mo '../../ecco_darwin/regions/RedSea/kaust_v1s3/code_darwin ../../ecco_darwin/regions/RedSea/kaust_v1s3/code'
  make depend
  make -j 16
  % if make -j 16 not working, simply: make

==============
# 3. Instructions for running simulation (currently only for 1996 to resolve seasonal cycle)
  cd $WORKDIR/darwin3/run
  ln -sf ../build/mitgcmuv .
  cp ../../ecco_darwin/regions/RedSea/kaust_v1s3/input/* .
  sbatch run_sbatch.sh
