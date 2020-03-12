#
# Mackenzie Delta regional setup based on LLC4320
# WARNING: Before starting make you have an Earthdata account (Or Create it at: https://urs.earthdata.nasa.gov/users/new)

# ========
# 1. Get code
git clone https://github.com/MITgcm/MITgcm.git
svn checkout https://github.com/MITgcm-contrib/ecco_darwin/trunk/Mac_Delta
# For the following requests you need your Earthdata username and WebDAV password (different from Earthdata password)
# Find it at :https://ecco.jpl.nasa.gov/drive
wget -r --no-parent --user=USERNAME --ask-password https://ecco.jpl.nasa.gov/drive/files/ECCO2/LLC4320/Mac_Delta/EOG/
wget -r --no-parent --user=USERNAME --ask-password https://ecco.jpl.nasa.gov/drive/files/ECCO2/LLC4320/Mac_Delta/run_template 
mv ecco.jpl.nasa.gov/drive/files/ECCO2/LLC4320/Mac_Delta/run_template Mac_Delta/
mv ecco.jpl.nasa.gov/drive/files/ECCO2/LLC4320/Mac_Delta/EOG/* Mac_Delta/input/
rm -r ecco.jpl.nasa.gov
cd MITgcm
mkdir build run

# ================
# 2. Build MPI executable
#    Prerequisite: 1. Get code
cd build
cp ../../Mac_Delta/code/* .

    # On Pleiades follow the instructions below:
    module purge
    module load comp-intel/2016.2.181 mpi-sgi/mpt.2.14r19 hdf4/4.2.12 hdf5/1.8.18_mpt netcdf/4.4.1.1_mpt
    ../tools/genmake2 -of ../tools/build_options/linux_amd64_ifort+mpi_ice_nas -mo ../../Mac_Delta/code
    
    # On other machine, use the following command with build option file ("-of ..") compatible with your machine
    # Following example works on Computecanada:
    ../tools/genmake2 -mpi -of ../tools/build_options/linux_amd64_ifort+gcc -mo ../../Mac_Delta/code
    
make depend
make -j 16

# ================
# 3. Run the setup
#    Prerequisite: 2. Build executable
cd ../run
mkdir diags
ln -sf ../build/mitgcmuv .
cp ../../Mac_Delta/input/*data* .
ln -sf .../../Mac_Delta/input/EOG* .
ln -sf .../../Mac_Delta/run_template/* .
# Run the job (Running on supercomputer might request sbatch submission)
mpirun -np 2227 ./mitgcmuv

# ================
# 4. sbacth request file example (running on computecanada)
------------------------
#!/bin/bash
#SBATCH --time=01:00:00
#SBATCH --account=def-account
#SBATCH --nodes=56
#SBATCH --ntasks-per-node=40
#SBATCH --mem=0
#SBATCH --output=Job_%j.out
#SBATCH --mail-user=clement.bertin1@univ-lr.fr
#SBATCH --mail-type=ALL

mpirun -np 2227 ./mitgcmuv
------------------------

