# v06 1D_ocean_ice_column_darwin idealized experiment
# initially based on MITgcm/verification/1D_ocean_ice_column

# ========
# 1. Get code
cd ${DARWIN_PATH}/darwin3-darwin
 mkdir project
 cd project
    DOWNLOAD AND UNZIP HERE THE EXPERIMENTS
 cd 1D_ocean_ice_column_darwin06
 mkdir build
 cd build

# ================
# 2. Build executable
#    Prerequisite: 1. Get code
 ../../../tools/genmake2 -mods ../code -optfile «/PATH/TO/OPTFILE»
 make depend
 make

# ======================
# 3. Run verification setup
#    Prerequisite: 2. Build executable
 cd ../run
 ln -s ../input/* .
 cp ../build/mitgcmuv .
 ./mitgcmuv
