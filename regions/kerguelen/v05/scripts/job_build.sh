#!/bin/bash
# job_build.sh -- build the Kerguelen ~1 km regional executable on Pleiades.
#
# Usage:
#   DARWIN3=/path/to/darwin3 ECCO_DARWIN=/path/to/ecco_darwin \
#     qsub -v DARWIN3,ECCO_DARWIN scripts/job_build.sh
# or edit the two defaults below.
#
#PBS -l select=1:ncpus=16:model=sky_ele
#PBS -q normal
#PBS -l walltime=1:00:00
#PBS -N kerg_build
#PBS -j oe

set -e

DARWIN3=${DARWIN3:?set DARWIN3 to your darwin3 checkout}
ECCO_DARWIN=${ECCO_DARWIN:?set ECCO_DARWIN to your ecco_darwin checkout}
BUILD=${BUILD:-$DARWIN3/build}

CFG="$ECCO_DARWIN/regions/kerguelen/v05"

module purge
module load comp-intel/2023.2.1 mpi-hpe/mpt.2.30 hdf4/4.2.12 hdf5/1.8.18_mpt \
            netcdf/4.4.1.1_mpt python3/3.9.12

mkdir -p "$BUILD"
cd "$BUILD"

# Full object wipe. SIZE.h / *_SIZE.h carry PARAMETERs that are compiled into
# every routine dimensioning an array against them, so a partial rebuild that
# leaves one stale .o gives a binary whose routines disagree about an array
# bound -- a silent memory bug, not a compile error.
rm -f ./*.o ./*.mod mitgcmuv

# ORDER MATTERS: code_darwin FIRST. genmake2 resolves -mo left to right and
# the first match wins. Both directories contain OBCS_OPTIONS.h and
# DIAGNOSTICS_SIZE.h, so code_darwin's copies are the ones compiled -- notably
# ALLOW_OBCS_SPONGE, which is #define'd there. Reversing this silently changes
# the build.
"$DARWIN3/tools/genmake2" -rootdir "$DARWIN3" \
  -of "$DARWIN3/tools/build_options/linux_amd64_ifort+mpi_ice_nas" \
  -mpi -mo "$CFG/code_darwin $CFG/code"

# AMD (rom_ait) portability. Required there, harmless elsewhere:
#  - the stock options file emits Intel-vendor-locked CPU dispatch
#    (-axCORE-AVX2) that traps on AMD;
#  - this configuration's static arrays exceed the small code model, hence
#    -mcmodel=medium, and --no-relax so the linker does not undo it.
sed -i 's/ -axCORE-AVX2 -xSSE4\.2//' Makefile
sed -i '/^LINK = /{/-mcmodel=medium/!s/$/ -mcmodel=medium -Wl,--no-relax/}' Makefile
sed -i '/^FFLAGS = /{/-mcmodel=medium/!s/$/ -mcmodel=medium/}' Makefile
sed -i '/^F90FLAGS = /{/-mcmodel=medium/!s/$/ -mcmodel=medium/}' Makefile

grep -nE '^FOPTIM|^F90OPTIM|mcmodel|no-relax' Makefile

make depend
make -j 16

echo "BUILD_DONE"
ls -la mitgcmuv
