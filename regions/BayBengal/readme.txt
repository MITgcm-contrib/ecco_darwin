Regional Bay of Bengal ECCO-Darwin cutout
Collaboration with Smitha Ratheesh
Indian Space Research Organization
smitha@sac.isro.gov.in

Binary files available in NASA Box:
https://nasa-ext.box.com/s/lp61pcgv3jmusn3gprl6pvr2tirt9sn6

# Instructions for building and running v06 Bay of Bengal
# latest MITgcm code (c69e as of 2025/05/20)

==============
# 1. Get code

git clone https://github.com/darwinproject/darwin3
git clone --depth 1 https://github.com/MITgcm-contrib/ecco_darwin.git
cd darwin3
mkdir build run
cd build

==============
# 2. Build executable for v06 Bay of Bengal

module purge
module load comp-intel mpi-hpe/mpt hdf4/4.2.12 hdf5/1.8.18_mpt netcdf/4.4.1.1_mpt python3/3.9.12

######## To change once cut-out available ########
../tools/genmake2 -of ../../ecco_darwin/[add cutout model path here]/code_physics/linux_amd64_ifort+mpi_ice_nas \
  -mo '../../ecco_darwin/[add cutout model path here]/code_darwin ../../ecco_darwin/[add cutout model path here]/code_physics'
######## 

make depend
make -j 16

==============
# 2. Instructions for running v06 ECCO-Darwin Bay of Bengal for 1992-2024 period

cd ../run
ln -sf ../build/mitgcmuv .
######## Copy inputs in run directory     ########
######## To change once cut-out available ########
# copy initial and boundary conditions
# copy forcings
# copy river forcings
mkdir diags diags/3hourly diags/daily diags/monthly diags/budget
mkdir diags/monthly/IOPS diags/monthly/PAR diags/monthly/RRS

# modify job_ECCO_darwin as needed
qsub job_ECCO_darwin


