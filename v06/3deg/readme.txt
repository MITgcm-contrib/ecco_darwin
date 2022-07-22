# v06 3deg darwin3 verification experiment
# initially based on MITgcm/verification/tutorial_global_oce_biogeo

# ========
# 1. Get code
git clone https://github.com/darwinproject/darwin3
git clone --depth 1 https://github.com/MITgcm-contrib/ecco_darwin.git
cd darwin3
mkdir build run
cd build


# ================
# 2. Build executable
#    Prerequisite: 1. Get code
 ../tools/genmake2 -ieee -mo \
 '../../ecco_darwin/v06/3deg/code ../../ecco_darwin/v06/llc270/code_darwin ../../ecco_darwin/v06/llc270/code'
 make depend
 make -j 8

# ======================
# 3. Run verification setup
#    Prerequisite: 2. Build executable
 cd ../run
 ln -sf ../build/mitgcmuv .
 cp ../../ecco_darwin/v06/llc270/input/data* .
 cp ../../ecco_darwin/v06/3deg/input/* .
 ln -sf ../../ecco_darwin/v06/3deg/data_darwin/* .
 rm data.ctrl data.exch2 data.smooth
 mkdir diags diags/3hourly diags/daily diags/monthly diags/budget
 mkdir diags/monthly/IOPS diags/monthly/PAR diags/monthly/RRS 
 ./mitgcmuv > output.txt
# Compare to verification output
 diff <(grep %MON output.txt) <(grep %MON ../../ecco_darwin/v06/3deg/results/output_3deg.txt)

# ============================
# 4. Build and run MPI executable
#    Prerequisite: 1. Get code
 cd ../build
 rm *
 module load comp-intel/2020.4.304 mpi-hpe/mpt.2.25 hdf4/4.2.12 hdf5/1.8.18_mpt netcdf/4.4.1.1_mpt
 cp ../../ecco_darwin/v06/3deg/code/SIZE.h_mpi SIZE.h
 ../tools/genmake2 -of  ../../ecco_darwin/v06/llc270/code/linux_amd64_ifort+mpi_ice_nas  \
 -mo '../../ecco_darwin/v06/3deg/code ../../ecco_darwin/v06/llc270/code_darwin ../../ecco_darwin/v06/llc270/code'
 make depend
 make -j 16

 cd ../run
 rm -rf *
 ln -sf ../build/mitgcmuv .
 cp ../../ecco_darwin/v06/llc270/input/data* .
 cp ../../ecco_darwin/v06/3deg/input/* .
 ln -sf ../../ecco_darwin/v06/3deg/data_darwin/* .
 mkdir diags diags/3hourly diags/daily diags/monthly diags/budget
 mkdir diags/monthly/IOPS diags/monthly/PAR diags/monthly/RRS 
 mv data_mpi data
 rm data.ctrl data.exch2 data.smooth
 mpirun -np 8 ./mitgcmuv &
# Monitor run
 tail -f STDOUT.0000 | grep advcfl_W
