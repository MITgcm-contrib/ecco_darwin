# v06 3deg darwin3 verification experiment
# 128x64x15 global lat-lon grid, ported from v05/3deg

# ========
# 1. Get code

 git clone https://github.com/darwinproject/darwin3
 git clone --depth 1 https://github.com/MITgcm-contrib/ecco_darwin.git
 cd darwin3/pkg/darwin
 git checkout backport_ckpt68y
 cd ../../
 mkdir build run
 cd build

# ================
# 2. Build executable (serial, for the verification run)
#    Prerequisite: 1. Get code

 MOD="../../ecco_darwin/v06/3deg"
 ../tools/genmake2 -ieee -mo "${MOD}/code ../../ecco_darwin/v06/llc270/code_physics"
 make depend
 make -j 8

# ======================
# 3. Build initial conditions
#    Prerequisite: 1. Get code

# v05/3deg ships 31 ptracer IC files; v06 needs 36, and the plankton indices
# do NOT correspond to the same organisms between versions. Build them with:

 cd ../../ecco_darwin/v06/3deg/matlab
# start MATLAB
 make_v06_3deg_initial_conditions

# This writes ptracers_v06_3deg_01..30.0000000001 into v06/3deg/data_darwin/.
# Tracers 31-36 (Chl01..Chl06) are intentionally not written -- data.darwin
# sets darwin_chlInitBalanced=T so Darwin derives Chl from biomass.
# Read the header of that script before trusting the plankton mapping.

# ======================
# 4. Run verification setup
#    Prerequisite: 2. Build executable, 3. Build initial conditions

 cd ../../../../darwin3/run
 ln -sf ../build/mitgcmuv .
 cp ../../ecco_darwin/v06/3deg/input/* .
 ln -sf ../../ecco_darwin/v05/3deg/data_darwin/* .        # grid, bathy, EXF forcing
 ln -sf ../../ecco_darwin/v06/3deg/data_darwin/* .        # v06 ptracer ICs
 mkdir -p diags/3hourly diags/daily diags/monthly diags/budget
 mkdir -p diags/monthly/IOPS diags/monthly/PAR diags/monthly/RRS
 ./mitgcmuv > output.txt

# ======================
# 5. Required v06 input that does NOT ship with v05/3deg
#    These are the files you must supply before the run will start.

# On NAS all of these are already available; link them in:

 ln -sf /nobackup/ojahn/forcing/oasim .
 ln -sf /nobackup/ojahn/ecco_darwin/v06/llc270/data_darwin/OPTICS_COEFF2 .
 ln -sf /nobackup/ojahn/ecco_darwin/v06/llc270/data_darwin/mitgcm_3he_flux.bin .
 ln -sf /nobackup/ojahn/ecco_darwin/v06/llc270/data_darwin/CESM-MIMI_*.bin .
 ln -sf /nobackup/dcarrol2/forcing/apCO2/NOAA_MBL/* .

# ======================
# 6. Build and run MPI executable
#    Prerequisite: 1. Get code

 cd ../build
 rm *
 module load comp-intel mpi-hpe hdf4/4.2.12 hdf5/1.8.18_mpt netcdf/4.4.1.1_mpt python3/3.9.12
 cp ../../ecco_darwin/v06/3deg/code/SIZE.h_mpi SIZE.h
 MOD="../../ecco_darwin/v06/3deg"
 ../tools/genmake2 -of ../../ecco_darwin/v06/llc270/code_physics/linux_amd64_ifort+mpi_ice_nas \
   -mo "${MOD}/code ../../ecco_darwin/v06/llc270/code_physics" -mpi
 make depend
 make -j 8

 cd ../run
 rm -rf *
 ln -sf ../build/mitgcmuv .
 cp ../../ecco_darwin/v06/3deg/input/* .
 ln -sf ../../ecco_darwin/v05/3deg/data_darwin/* .
 ln -sf ../../ecco_darwin/v06/3deg/data_darwin/* .
# plus the section 5 links
 mkdir -p diags/3hourly diags/daily diags/monthly diags/budget
 mkdir -p diags/monthly/IOPS diags/monthly/PAR diags/monthly/RRS
 mv data_mpi data
# modify job_3deg_darwin3 as needed
 qsub job_3deg_darwin3
