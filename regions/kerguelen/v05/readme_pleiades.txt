# Instructions for building and running the Kerguelen Plateau ~1 km regional
# simulation on Pleiades. Downscaled from ECCO-Darwin LLC270 v05 (v5r1) using
# the ecco_darwin/regions/downscaling toolkit (STEP1-STEP4).

==============
# 0. Domain summary
  Grid        1170 x 1008 x 90, lat/lon, uniform spacing
  Extent      62.000 - 78.007 E, 54.000 - 44.996 S
  Resolution  0.0136811 lon x 0.0089330 lat  (~0.99 km x ~0.99 km)
  Period      2020-01-01 to 2023-01-01 (3 years), deltaT = 90 s
              endtime = 94694400 -> 1,052,160 timesteps
  Decomp      sNx=39 sNy=42, nPx=30 nPy=24 -> 720 MPI ranks
  Parent      ECCO-Darwin LLC270 v05 (v5r1), 31-tracer Darwin
  Boundaries  all four open (N/S/E/W), prescribed hourly from the parent
              via pkg/diagnostics_vec extraction + OBCS generation

==============
# 1. Get code
  git clone git@github.com:MITgcm-contrib/ecco_darwin.git
  git clone git@github.com:darwinproject/darwin3
  cd darwin3
  mkdir build run

==============
# 2. Build executable
  cd build
  module purge
  module load comp-intel/2023.2.1 mpi-hpe/mpt.2.30 python3/3.9.12
  module load hdf4/4.2.12 hdf5/1.8.18_mpt netcdf/4.4.1.1_mpt

  ../tools/genmake2 -rootdir .. -mpi \
    -of ../tools/build_options/linux_amd64_ifort+mpi_ice_nas \
    -mo "../../ecco_darwin/regions/kerguelen/v05/code_darwin ../../ecco_darwin/regions/kerguelen/v05/code"

# NOTE THE ORDER: code_darwin comes FIRST. genmake2 resolves -mo directories
# left to right and the FIRST match wins. Both directories contain
# OBCS_OPTIONS.h and DIAGNOSTICS_SIZE.h, so code_darwin's copies are the ones
# that actually get compiled. In particular ALLOW_OBCS_SPONGE is governed by
# code_darwin/OBCS_OPTIONS.h, not by code/OBCS_OPTIONS.h. Reversing the order
# silently changes the build.

# AMD (rom_ait) portability -- required if running on rom_ait, harmless
# otherwise. The stock options file emits Intel-vendor-locked CPU dispatch
# that traps on AMD, and this configuration's static arrays exceed the small
# code model:
  sed -i 's/ -axCORE-AVX2 -xSSE4\.2//' Makefile
  sed -i '/^LINK = /{/-mcmodel=medium/!s/$/ -mcmodel=medium -Wl,--no-relax/}' Makefile
  sed -i '/^FFLAGS = /{/-mcmodel=medium/!s/$/ -mcmodel=medium/}' Makefile
  sed -i '/^F90FLAGS = /{/-mcmodel=medium/!s/$/ -mcmodel=medium/}' Makefile

  make depend
  make -j 16

# Or use the packaged job script, which does all of the above:
#   qsub scripts/job_build.sh

==============
# 3. Stage the run directory
  cd ../run
  ln -sf ../build/mitgcmuv .
  cp ../../ecco_darwin/regions/kerguelen/v05/input/* .
  cp ../../ecco_darwin/regions/kerguelen/v05/input_darwin/* .
# input_darwin is copied SECOND and deliberately overwrites EXACTLY THREE
# files from input/: data.diagnostics, data.obcs and data.pkg. The input/
# versions are the physics-only configuration; the input_darwin/ versions add
# the BGC tracers, their boundary conditions and their diagnostics. Copying in
# the other order gives a physics-only run that still starts normally and
# reports no error -- check data.pkg has usePTRACERS/useGCHEM = .TRUE. if in
# any doubt.

# Atmospheric / surface forcing (shared, not in this repo). data.exf refers to
# the ERA5 files through the directory name, so link the DIRECTORY, not its
# contents:
  ln -sf /nobackup/hzhang1/forcing/era5 ERA5
  ln -sf /nobackup/dcarrol2/forcing/apCO2/NOAA_MBL/* .
  ln -sf /nobackup/dcarrol2/v05_V5r1/ecco_darwin/v05/3deg/data_darwin/3deg_Mahowald_2009_soluble_iron_dust.bin .

# Grid, bathymetry, boundary conditions and initial conditions. These are
# GENERATED, not stored in the repo -- see PIPELINE.md. They are large
# (~100 GB of OBCS + ~50 GB of pickups):
  ln -sf <config>/grid/kerguelen_bathymetry.bin .
  ln -sf <config>/forcings/OBCS/*.bin .            # 164 files
  ln -sf <config>/forcings/pickups/* .             #  78 files

  mkdir -p diags chain_logs

==============
# 4. Run
# Interactively, or from your own single-segment PBS script:
  mpiexec -np 720 ./mitgcmuv
# This will not finish a 3-year integration in one queue window -- see below.

# Production: the integration is ~28 days of compute at 720 ranks, which
# exceeds the `long` queue's 120 h ceiling, so it must be chained across
# segments. scripts/job_production_chain.sh self-resubmits:
  qsub -v SEGMENT=1 -W group_list=<group> scripts/job_production_chain.sh

# It caps itself at 115 h (5 h under the PBS walltime) to leave room for
# checkpoint resolution and resubmission, resolves the furthest
# verified-complete checkpoint with find_restart_iter.py, rewrites nIter0,
# and resubmits. It halts rather than resubmitting on: a hard numerical
# error in STDOUT, no valid checkpoint, or zero progress since the segment
# started.

# Why find_restart_iter.py exists: MITgcm's pickup reader string-formats
# nIter0 into 'pickup.<iter>.data'. It does NOT inspect the rolling
# ckptA/ckptB .meta timeStepNumber to find a matching restart. The script
# verifies the newest complete checkpoint (permanent or rolling) by exact
# expected byte size AND a consistent timeStepNumber across all five
# filesets, then renames rolling checkpoints to the numbered convention.

==============
# 5. Notes on this configuration

# OBCS sponge. useOBCSsponge = .TRUE. in data.obcs relaxes U/V toward the
# prescribed boundary profile over a short sponge zone. This damps a
# near-boundary velocity flicker: OBCS_u1_adv_Tr's per-timestep
# outflow/inflow test keys off the sign of the local normal velocity, which
# can oscillate at the edge in an eddying flow and produced a tracer
# reflection artifact at the east and north boundaries. Note that stock
# pkg/obcs has no ptracer sponge (obcs_sponge.F defines only
# OBCS_SPONGE_U/V/T/S), so tracers are not relaxed directly.

# Barotropic transport correction. The OBCS velocity files have a uniform
# per-timestep correction applied so the regional boundary carries the same
# net volume transport as the parent. Without it, interpolation error across
# the ~23x jump in resolution accumulates into ETAN drift under
# nonlinFreeSurf/exactConserv. See PIPELINE.md STEP3.

# Cold start. input/data ships nIter0=0 with hydrogThetaFile /
# hydrogSaltFile pointing at the STEP3-generated pickup fields at parent
# iteration 736344 (2020-01-01). The production chain rewrites nIter0 in
# place as it advances, so a run directory mid-chain will not match the
# repo copy -- that is expected.
