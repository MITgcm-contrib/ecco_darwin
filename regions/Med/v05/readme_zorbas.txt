# Instructions for building and running Mediterranean regional simulation
# on Ubuntu Crete HCMR's zorbas cluster (https://hpc.hcmr.gr/)

==============
# 1. Log in to zorbas login node
  (see https://hpc.hcmr.gr/infrastructure/ for details)

==============
# 2. Get code
  git clone git@github.com:MITgcm-contrib/ecco_darwin.git
  git clone git@github.com:darwinproject/darwin3
  cd darwin3
  git checkout 24885b71
  mkdir build run

==============
# 3. Run a single-CPU verification experiment
  cd darwin3/verification
  ./testreport -t lab_sea

==============
# 4. Log in to hydrogen (the intermediate node)
#    and run an mpi verification experiment
  ssh hydrogen
  cd darwin3/verification
  ./testreport -mpi -t lab_sea

==============
# 5. Build executable
  cd build
  export MPI_INC_DIR=/usr/lib/x86_64-linux-gnu/openmpi/include
  ../tools/genmake2 -mo ../../ecco_darwin/regions/Med/v05/code -mpi
  make depend
  make -j

==============
# 6. Instructions for running simulation (1992-2023 period)
  cd ../run
  ln -sf ../build/mitgcmuv .

# Get forcing and configuration files from:
# --> https://nasa-ext.box.com/s/3d3qz47tvnhp2y8wbvd821rwdxk1m2un
# --> https://nasa-ext.box.com/s/khtbuge4wvt5yleyigcjbqvdjebh3xcn
# and deposit or link the contents of these directories
# inside the darwin3/run directory, for example,
  ln -sf <path_to_download_location>/NOAA_MBL/* .
  ln -sf <path_to_download_location>/Med/run_template/* .

  mkdir diags diags/daily diags/monthly
  cp ../../ecco_darwin/regions/Med/v05/input/* .
  mpirun -np 8 ./mitgcmuv

==============
# 7. Run simulation with sbatch
  sbatch run_on_zorbas.sh

# a few more sbatch commands:
  squeue  # to look at running jobs
  sinfo   # for information on available queues
  scancel # for canceling a running job

# notes and things to try
- 8 cpus (2 nodes * 4 cpus, 36x36 tiles): ~20 minutes per month

- 108 cpus (3 nodes * 36 cpus, 12x8 tiles): ~6 minutes per month 
