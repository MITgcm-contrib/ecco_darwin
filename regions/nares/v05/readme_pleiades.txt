# Instructions for building and running the Nares Strait 1 km regional
# simulation on Pleiades. Downscaled from ECCO-Darwin LLC270 v05 (V5r1) with
# the ecco_darwin/regions/downscaling toolkit plus Nares-specific scripts.
#
# Physics frozen 2026-11-20 at run_arch_D10 (Dustin). code/ and input/
# reproduce D10, with the 20x20 production tiling (678 ranks). Every 2015
# tuning run used 40x40 tiles (192 ranks). See ../OPEN_QUESTIONS.md.

==============
# 0. Domain summary
  Grid        1120 x 360 x 90, rotated lat-lon, one exch2 facet; i runs S -> N
              along the strait (i=1 North Water / N Baffin Bay, i=1120 Lincoln Sea)
  Rotation    rotated equator = great circle through A (83.6 N, 55 W) and
              B (76.0 N, 73.5 W); +-180 km across, 60 km past A, 120 km past B
  Extent      80.80 W - 40.52 W, 74.37 N - 84.94 N (corner values)
  Resolution  999.6 - 1000.0 m in both directions (analytic metrics)
  Vertical    90 levels, 1.00 m at the surface (same delR as pan_AO / west_AO)
  Bathymetry  GEBCO 2026 cell means, 10 m minimum depth (deepen only), 571
              min-depth fill-bridge cells walled off (Nares_bathymetry_walled.bin);
              Petermann floating tongue = static pkg/shelfice cavity from
              BedMachine v6 (1200 cells, draft 9-586 m)
  Period      2015-01-01 00:00 to 2015-09-01 00:00 (243 d), deltaT = 60 s,
              nTimeSteps = 349920, cold start (nIter0 = 0).
  Decomp      sNx = sNy = 20, 56 x 18 = 1008 tiles, 330 blank (data.exch2),
              678 active -> SIZE.h nPx = 678, nPy = 1 -> 678 MPI ranks
              (the 2015 tuning runs: 40x40, 252 tiles, 60 blank, 192 ranks)
  Parent      ECCO-Darwin LLC270 v05 (V5r1), 31-tracer Darwin; IC = parent
              pickup iteration 604872 (2015-01-01); boundaries from the shared
              west_AO + Nares pkg/diagnostics_vec extraction (hourly, 2015-2025)
  Boundaries  all four edges open (W i=1, E i=1120, S j=1, N j=360), prescribed
              hourly (T, S, U, V, sea ice: area, heff, hsnow, uice, vice; 33
              ptracers in BGC runs), useOBCSbalance (OBCSbalanceSurf), sponge
              10 cells (6 h at the edge -> 10 d inside), barotropic tides with
              8 constituents (M2 S2 N2 K2 K1 O1 P1 Q1, FES2014 currents scaled
              to the FES2022 prism budget, phase zero at 2015-01-01).
              Velocities rotated to the grid axes; +0.331 Sv southward
              throughflow offset added (net zero over the two strait ends).
  Forcing     ERA5 hourly (shared /nobackupp17/hzhang1/forcing/era5); no EXF
              runoff. Humboldt subglacial discharge goes through pkg/iceplume
              (yearly runoffQsg files, daily), which is compiled but switched
              OFF (useICEPLUME = .FALSE.) in every run so far.
  Sea ice     LSR (JFNK compiled in, off), TEM, SEAICE_strength 6.875e4,
              eccen 1.6, advScheme 77, HO 2.0, area_max 0.98, no tensile /
              side / basal drag, ice load capped at 6 m of water (code_capload),
              thinned-mouth IC (siHEFF_ICT), ice-ocean drag on the 0-10 m
              mean ocean current (code_drag10). = frozen run_arch_D10
  Mixing      GGL90 with the parent set (alpha 30, TKEmin 1e-7, TKEbottom
              1e-6, mxlMaxFlag 2, mxlSurfFlag T).
  BGC         33 ptracers (31 parent + rDOC + CDOM), pan_AO/west_AO v5
              ecosystem (RADIv2 emulator + limiter), radtrans + OASIM
              (oasim_dTsolz 60), under-ice light (RT_iceTransmit).

==============
# 1. Get code
  git clone https://github.com/MITgcm-contrib/ecco_darwin.git
  git clone https://github.com/darwinproject/darwin3
  cd darwin3
  git checkout darwin_ckpt69k          # = 1925cd1e1, the tag the Nares tree is built on
  cp -r ../ecco_darwin/regions/nares/v05/pkg/iceplume pkg/
  mkdir build run

# The Nares model tree on Pleiades (/nobackup/dcarrol2/downscaling/nares/darwin3,
# branch `nares`, HEAD 8eb545d38) = darwin_ckpt69k + exactly two additions,
# both reproduced by this repository:
#   1. pkg/iceplume (github.com/mhwood/iceplume e4d4ac9 + 13 fixes): the
#      15 files in v05/pkg/iceplume/, copied above. See its README.md.
#   2. a SEAICE_SWFrac backport from darwin3 master (pkg/seaice/SEAICE.h,
#      seaice_init_fixed.F; 51 added lines): shipped as overrides in
#      v05/code/, so genmake2 -mo picks them up. The darwin branch has the
#      Darwin extensions but an older pkg/seaice without SEAICE_SWFrac;
#      master has SEAICE_SWFrac but no Darwin extensions. seaice_growth.F
#      (code/) needs both.
# The core files iceplume patches are in code/ as well (PROVENANCE.md).
# Verified 2026-11-16 against the Pleiades tree: a build from these two
# directories (build_v05test) compiles the same sources as build_drag10 plus
# the 18 code_darwin files. Not verified: the same build from a fresh public
# darwin_ckpt69k clone (the SWFrac backport files were only checked as
# overrides on the Pleiades tree, where pkg/seaice already carries them).

==============
# 2. Build executable
  cd build
  module purge
  module load comp-intel/2023.2.1 mpi-hpe/mpt.2.30 hdf4/4.2.12 hdf5/1.8.18_mpt netcdf/4.4.1.1_mpt python3/3.9.12

# Physics-only binary (= build_drag10 sources, with the 20x20 SIZE.h):
  ../tools/genmake2 -rootdir .. -mpi \
    -of ../tools/build_options/linux_amd64_ifort+mpi_ice_nas \
    -mo "../../ecco_darwin/regions/nares/v05/code"

# BGC binary (code_darwin FIRST, then code):
  ../tools/genmake2 -rootdir .. -mpi \
    -of ../tools/build_options/linux_amd64_ifort+mpi_ice_nas \
    -mo "../../ecco_darwin/regions/nares/v05/code_darwin ../../ecco_darwin/regions/nares/v05/code"

# NOTE THE ORDER: genmake2 resolves -mo left to right and the FIRST match
# wins. code/ is the collapsed physics build (see PROVENANCE.md section 1); it
# already carries packages.conf (darwin, gchem, ptracers, radtrans, oasim,
# shelfice, iceplume, exch2, obcs ...), the size headers and an OLD Darwin
# source set that is inactive in physics runs. code_darwin/ replaces 18 Darwin
# files with the v5 ecosystem set. Reversing the order gives a BGC binary with
# the stale (invalid) 17-argument darwin_plankton.F and no error.
#
# NOTE 2: the BGC combination code_darwin + code was test-built on
# 2026-11-16 (build_v05test: clean; sources = build_bgc + swfix + drag10
# exactly) and again with the 20x20 SIZE.h (OPEN_QUESTIONS.md item 4). It has
# NOT been run. build_bgc (2026-10-29), which ran the 2015 BGC tests,
# predates code_swfix and code_drag10.

# AMD (rom_ait) portability and the >2 GB static arrays (mitgcmuv has a
# 1.93 GB bss): required on Pleiades, harmless otherwise:
  sed -i '/^LINK = /{/-mcmodel=medium/!s/$/ -mcmodel=medium -Wl,--no-relax/}' Makefile
  sed -i '/^FFLAGS = /{/-mcmodel=medium/!s/$/ -mcmodel=medium/}' Makefile
  sed -i '/^F90FLAGS = /{/-mcmodel=medium/!s/$/ -mcmodel=medium/}' Makefile
  sed -i 's/ -axCORE-AVX2 -xSSE4\.2//' Makefile

  make depend
  make -j 16

# Or use the packaged job scripts, which do all of the above plus pre-flight
# checks (33/33 tracers, 678 ranks = 1008 - 330 blank tiles, SIZE.h vs
# data.exch2 vs blankList, every OBCS namelist switch set .TRUE. has its CPP
# flag) and post-build asserts (iceplume, shelfice, JFNK, side drag, radtrans
# objects; which override directory each file came from):
#   qsub -W group_list=<group> -v P=<tree>/regions/nares/v05,B=<dir>,BUILD=build_v05 \
#        scripts/job_build_v05.sh                             (production, 20x20)
# The as-run build scripts of the 2015 tuning runs (job_build_drag10.sh,
# job_build_bgc.sh ...) stop unless SIZE.h gives 192 ranks (40x40 tiles).
# Both wipe their build directory first: genmake2 does not replace existing
# source links, and a build directory that has seen another package set ends
# up as a mixture of both.

==============
# 3. Stage the run directory
  cd ../run
  ln -sf ../build/mitgcmuv .
  cp ../../ecco_darwin/regions/nares/v05/input/* .
  cp ../../ecco_darwin/regions/nares/v05/input_darwin/* .      # BGC runs only
# input_darwin is copied SECOND and overwrites EXACTLY TWO files from input/:
# data.pkg (usePTRACERS/useGCHEM/useRADTRANS/useOASIM = .TRUE.) and
# data.diagnostics (BGC streams added). The other six input_darwin files have
# no physics counterpart. For a physics-only run skip the second cp.

# Grid files (generated, see PIPELINE.md). BOTH are needed:
  ln -sf <config>/grid/Nares.mitgrid tile001.mitgrid        # exch2 facet geometry (16 fields)
  ln -sf <config>/grid/Nares_grid18.face001.bin .           # 16 fields + AngleCS/AngleSN (data: horizGridFile = 'Nares_grid18')
  ln -sf <config>/grid/Nares_bathymetry_walled.bin Nares_bathymetry.bin   # data names Nares_bathymetry.bin
  ln -sf <config>/grid/Nares_shelficeTopo.bin .
  ln -sf <config>/grid/Nares_shelficeLoadAnomaly.bin .

# Atmospheric forcing. data.exf refers to the ERA5 files through the directory
# name, so link the DIRECTORY, not its contents:
  ln -sf /nobackupp17/hzhang1/forcing/era5 ERA5

# Boundary conditions and initial conditions (generated, large):
  ln -sf <config>/forcings/OBCS2_rot_tf_JanAug/*.bin .      # 40 files, 5856 hourly records
  ln -sf <config>/forcings/OBCS_tides_2015Jan_c2/*.bin .    # 16 tidal files (phased for 2015-01-01)
  ln -sf <config>/forcings/pickups/* .                      # 22 files, iteration 604872
  ln -sf <config>/iceplume/Humboldt_plumes_baseline_*.bin . # iceplume mask + length
  ln -sf <config>/forcings/iceplume/runoffQsg_* .           # yearly discharge (read even with iceplume off? not verified)

# BGC runs add:
  ln -sf <config_bgc>/forcings/OBCS_bgc_JanAug/*.bin .      # 152 files (TRACnn, rDOC, CDOM; about 206 GB)
  ln -sf <config_bgc>/forcings/pickups_bgc/* .              # 82 files (pTr01-31, rDOC, CDOM, darwin pH)
  ln -sf /nobackup/dcarrol2/forcing/apCO2/NOAA_MBL/apCO2_* .
  ln -sf /nobackup/dcarrol2/v05_V5r1/ecco_darwin/v05/3deg/data_darwin/3deg_Mahowald_2009_soluble_iron_dust.bin .
  ln -sf /nobackup/ojahn/forcing/oasim oasim
  ln -sf <optics dir> optics        # the AO_ecosystem optics set (see stage_bgc_test.sh)

  mkdir -p diags/3hourly diags/6hourly diags/daily diags/monthly diags/5daily

# Or use the stage scripts (scripts/stage_arch_tc4.sh, stage_bgc_prod.sh ...),
# which do all of the above on Pleiades with absolute paths.

==============
# 4. Run
# Production (20x20 tiles, 678 ranks): 18 sky_ele nodes, 38 ranks per node
# (684 slots), as the 20x20 scaling test ran (run_bgc_sc_20x20, 2026-11-16):
  qsub -q long -l select=18:ncpus=40:mpiprocs=38:model=sky_ele -l walltime=120:00:00 \
       -W group_list=<group> -v RUN=<run dir>,NP=678 scripts/job_run_scaling.sh
# The script checks that every OBCS file holds enough hourly records for
# nIter0 + nTimeSteps.
#
# The 2015 physics tuning runs (40x40, 192 ranks) used 5 sky_ele nodes, 39 ranks
# per node, and scripts/job_run_arch_long.sh. Measured: 0.334 s per step
# (run_arch_Pe16), about 7.5 model days per wall hour; Jan 1 to Sep 1 2015 took
# 27.4 h (run_arch_Y69m_sw) and 28.2 h (run_arch_D10).

# BGC on 40x40 (192 ranks): measured 2.337 s per step = 7.0x the physics, 1.07 model
# days per wall hour (run_bgc_test2d). Jan 1 to Sep 1 2015 does not fit one
# 120 h long-queue window, so scripts/stage_bgc_prod.sh stages 3 segments at
# calendar-month pickups (Jan 1 -> Apr 1, 129600 steps; Apr 1 -> Jul 1,
# 131040; Jul 1 -> Sep 1, 89280) submitted with afterok. Measured: 84 h,
# about 85 h, 58 h. BGC on 20x20 (678 ranks): measured 0.962 s per step
# (2.43x faster). The 2016-2025 production estimate is about 1405 wall hours,
# 12 long-queue segments of 120 h (OPEN_QUESTIONS.md item 3).

==============
# 5. Notes on this configuration

# Grid angles. The model MUST read Nares_grid18.face001.bin (horizGridFile =
# 'Nares_grid18' in data). With only tile001.mitgrid the model computes the
# angles itself (calc_grid_angles.F) and gets the east (i=1120) and north
# (j=360) edge rows wrong by up to 167 degrees, which rotates the EXF wind
# there. Found 2026-10-15 (run_arch_G18 test).

# Tides. useOBCStides with 8 constituents. The am/ph files are phased for
# startDate 2015-01-01 (nodal corrections at 2015-01-16): regenerate them for
# another start year (scripts/pipeline/tides). Tangential tidal amplitudes are
# 0 and the normal amplitude is tapered to 0 over the 3 cells next to each
# coast (variant "C2", 2026-10-15), which removed spurious 2-5 m/s currents in
# 10 m cells at the SE corner. The sponge (10 cells, 6 h -> 10 d) came in with
# the tides.

# Throughflow. The LLC270 parent carries about half the observed southward
# Nares throughflow (0.39 Sv vs about 0.85 Sv). A constant barotropic offset
# of +0.331 Sv southward is added to the OBCS velocities over the two strait
# ends (east and west edges), net zero so useOBCSbalance leaves it alone.
# The OBCS velocities were also rotated from geographic east/north to the
# grid axes (the grid is rotated 43-83 degrees); the first OBCS set
# (2026-10-07) lacked this rotation. See PIPELINE.md.

# Sea-ice tuning (frozen at run_arch_D10, 2026-11-20). The 2015 ice-arch timing was tuned against MODIS,
# OSI SAF drift and CS2SMOS thickness over about 70 short and 20 full runs
# (DECISIONS.md, not in this repository). Changes with a code side:
#   code_capload  2026-10-24  ice load on the sea surface capped at 6 m of water
#                             (default 0.2 * drF(1) = 0.2 m with a 1 m top layer
#                             would remove most of the load; CALC_R_STAR crash)
#   code_swfix    2026-11-08  physics-only runs had NO ocean shortwave heating:
#                             CDOMdecay is set only by darwin_forcing (useGCHEM),
#                             so the SW path in apply_forcing.F and seaice_growth.F
#                             is now guarded by IF (useGCHEM). No effect on BGC runs.
#   code_drag10   2026-11-15  ice-ocean drag uses the 0-10 m mean ocean velocity
#                             (hFac-weighted, levels 1-7) instead of the 1 m top
#                             layer, which flowed about 30 % slower with the parent
#                             GGL90 set and left the strait mouth ice locked.
#                             run_arch_D10 (full year) is the test.

# iceplume is compiled in but OFF (useICEPLUME = .FALSE. in data.pkg) in every
# run so far, including the BGC chain. To switch it on: set useICEPLUME =
# .TRUE. in data.pkg; keep selectAddFluid = 1 in data (ICEPLUME_CHECK needs
# it); link the Humboldt plume mask and length files and the yearly
# runoffQsg_YYYY files (data.iceplume names them; PIPELINE.md STEP 3). The
# plume ptracer concentrations (ICEPLUME_PARM02/03) are unset. Untested in a
# full run (a 1-day test ran on 2026-10-07 with iceplume on).

# Cold start. input/data ships nIter0 = 0 with hydrogThetaFile /
# hydrogSaltFile / uVelInitFile / vVelInitFile / pSurfInitFile pointing at the
# generated pickup fields at parent iteration 604872 (2015-01-01); data.seaice
# reads the split sea-ice IC files (siAREA, siHEFF_ICT, siHSNOW, siUICE,
# siVICE); data.ggl90 reads pickup_ggl90. Pickups are written every 30 days
# (pChkptFreq = chkptFreq = 2592000 s), which the BGC segment plan relies on.

# Diagnostics. The 3 hourly lists were dropped on 2026-11-13 (inode limit on
# /nobackupp19). data.diagnostics writes 3-hourly SIdyn, 6-hourly SIforce,
# daily SIdeform / surf / surfTS / SHI / SImelt and monthly state.
