#!/bin/bash
# socal 2-km nest, tides build (code_tides/obcs_add_tides.F ramp + code/): build mitgcmuv (v05_V5r1 darwin3, code/ = Kerguelen regional code + code_darwin
# merged, seaice/salt_plume dropped, SIZE.h 24 x 1 tiles of 7 x 101 x 60).
# Same recipe as kerguelen/regional_run/job_build.sh (-mcmodel=medium, no -ax dispatch).
#PBS -l select=1:ncpus=16:model=sky_ele
#PBS -q normal
#PBS -l walltime=1:00:00
#PBS -N socal_btide
#PBS -j oe
#PBS -o /nobackup/dcarrol2/downscaling/socal/regional_run/build_tides_out.txt

module purge
module load comp-intel/2023.2.1 mpi-hpe/mpt.2.30 hdf4/4.2.12 hdf5/1.8.18_mpt netcdf/4.4.1.1_mpt python3/3.9.12

DARWIN3=/nobackup/dcarrol2/v05_V5r1/darwin3
R=/nobackup/dcarrol2/downscaling/socal/regional_run
mkdir -p "$R/build_tides"; cd "$R/build_tides"
rm -rf ./*   # fresh build dir (make CLEAN equivalent)
"$DARWIN3/tools/genmake2" -rootdir "$DARWIN3" \
  -of "$DARWIN3/tools/build_options/linux_amd64_ifort+mpi_ice_nas" \
  -mpi -mo "$R/code_tides $R/code"
sed -i '/^LINK = /{/-mcmodel=medium/!s/$/ -mcmodel=medium -Wl,--no-relax/}' Makefile
sed -i '/^FFLAGS = /{/-mcmodel=medium/!s/$/ -mcmodel=medium/}' Makefile
sed -i '/^F90FLAGS = /{/-mcmodel=medium/!s/$/ -mcmodel=medium/}' Makefile
sed -i 's/ -axCORE-AVX2 -xSSE4\.2//' Makefile
make depend > /dev/null
make -j 16
echo "BUILD_DONE"; ls -la mitgcmuv
