# Instructions for building and running a physics-only
# RedSea/kaust_v1s3 configuration on shaheen at KAUST

==============
# 1. Get code
  git clone git@github.com:MITgcm-contrib/ecco_darwin.git
  git clone git@github.com:MITgcm/MITgcm
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
# 3. Prepare input data
 # if you do not have a shaheen account, please use this 2 links to download the bathymetry and init+exf+obcs files:
  # links for bathymetry data: https://wetransfer.com/downloads/c0996a454998dee5fe0babf9516d293c20260502153153/5acf7cbffbe0c39e30849c3c9ce5d6b120260502153202/a9075d
  # link for init+exf+obcs data: https://wetransfer.com/downloads/8e4f4a407d86c8234673b57e294a22ef20260502145256/6be66a43e7022c22eb3d375b77316bc520260502151736/ed1c57
 # if you have shaheen account, please go /scratch/wangy0m/ECCO_Darwin/Data_for_testing for the bathymetry+init+exf+obcs files
  cp -r /scratch/wangy0m/ECCO_Darwin/Data_for_testing [YOUR_INPUTDATA_FILE_DIR] 

==============
# 4. Instructions for running simulation (only physics part now for testing)
  cd ../run
  ln -sf ../build/mitgcmuv .
  ln -sf [YOUR_INPUTDATA_FILE_DIR]/* .
  cp ../../ecco_darwin/regions/RedSea/kaust_v1s3/input/* .

 # if you have shaheen account, please use the following commands to run (you may need to change project name, currently it is using k10023):
  rm -rf data data.exch2.mpi
  sbatch run_sbatch.sh
 # if don't have shaheen account, please use the following commands to run:
  mkdir diag2D diag3D
  rm -rf data.temp2 run_sbatch.sh tmp.start
  mpirun -np 1649 ./mitgcmuv
