#!/bin/bash
# job_build.sh -- compile the diagnostics_vec-enabled parent LLC270/Darwin
# config for pan_AO's parent-side boundary extraction (mirrors fiji's/
# Kerguelen's STEP2 job_build.sh exactly). BUILD ONLY (genmake2+make), not
# the production run -- cheap/fast.
#PBS -l select=1:ncpus=16:model=sky_ele
#PBS -q normal
#PBS -l walltime=1:00:00
#PBS -N parentWANbuild
#PBS -e /nobackupp19/dcarrol2/downscaling/parent_V5r1_westAO_nares/build_err.txt
#PBS -o /nobackupp19/dcarrol2/downscaling/parent_V5r1_westAO_nares/build_out.txt

module purge
module load comp-intel/2023.2.1 mpi-hpe/mpt.2.30 hdf4/4.2.12 hdf5/1.8.18_mpt netcdf/4.4.1.1_mpt python3/3.9.12

DARWIN3=/nobackup/dcarrol2/v05_V5r1/darwin3
ECCODARWIN_LLC270=/nobackup/dcarrol2/v05_V5r1/ecco_darwin/v05/llc270
F=/nobackupp19/dcarrol2/downscaling/parent_V5r1_westAO_nares

cd "$F/build"

# 2026-10-01: force a full object rebuild. DIAGNOSTICS_VEC_SIZE.h changed
# (nVEC_mask 15 -> 20, west boundary restored on the disc rim), and that
# PARAMETER is compiled into every routine that dimensions a vec array. A
# partial rebuild that left even one stale .o would give a binary whose
# routines disagree about the array bound -- a silent memory bug, not a
# compile error. The objects are cheap; the ambiguity is not.
# NOTE: deliberately does NOT touch the source/header symlinks. genmake2 will
# not replace symlinks that already exist, but the MITgcm tree is unchanged
# here (same $DARWIN3 as the 09-29 build), so they are all still correct. The
# 09-30 gen_ncgrid failure -- a build dir made of two trees -- only arises
# when $DARWIN3 itself changes, and then the fix is to wipe build/ entirely.
rm -f ./*.o ./*.mod mitgcmuv

"$DARWIN3/tools/genmake2" -rootdir "$DARWIN3" \
  -of "$DARWIN3/tools/build_options/linux_amd64_ifort+mpi_ice_nas" \
  -mpi -mo "$F/code_darwin $ECCODARWIN_LLC270/code_v5r1"
make depend
make -j 16

echo "BUILD_DONE"
ls -la mitgcmuv
