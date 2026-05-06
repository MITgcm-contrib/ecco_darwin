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
  ../tools/genmake2 -of=../../ecco_darwin/regions/RedSea/kaust_v1s3/shaheen_build_options -mpi -make=gmake \
   -mo '../../ecco_darwin/regions/RedSea/kaust_v1s3/code_darwin ../../ecco_darwin/regions/RedSea/kaust_v1s3/code'
  make depend
  make -j 16

==============
# 3. Instructions for running simulation (1992-2023 period)
  cd $WORKDIR/darwin3/run
  ln -sf ../build/mitgcmuv .
  ln -sf $WORKDIR/Data_for_testing/* .
  mkdir diags diags/daily diags/monthly
  cp ../../ecco_darwin/regions/RedSea/kaust_v1s3/input/* .
  cp ../../ecco_darwin/regions/RedSea/kaust_v1s3/input_darwin/* .
  sbatch run_sbatch.sh
