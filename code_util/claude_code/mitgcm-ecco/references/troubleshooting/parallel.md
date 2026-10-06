# Troubleshooting: MPI / tiling / SIZE.h / exch2 / blank tiles / OpenMP / scaling
Distilled from mitgcm-support (2003-2026) and MITgcm GitHub issues/PRs. Names verified against origin/master (eesupp/inc/CPP_EEOPTIONS.h, eesupp/src, pkg/exch2, doc/tag-index) as of 2026-10. Thread URLs are the month index pages (find the subject inside). GitHub items are cited by number.

## Launch / "No. of procs" / MPI build
### EEBOOT_MINIMAL: No. of procs= 1 not equal to nPx*nPy= N  (also "No. of processes not equal to nPx*nPy 12 2", "EEDIE: earlier error in multi-proc/thread setting")
- Cause: executable sees a different process count than SIZE.h (first number = procs MPI gave it, second = nPx*nPy). Usual reasons: (a) built without `genmake2 -mpi`; (b) run as `./mitgcmuv` instead of `mpirun -np N`; (c) mpirun/mpiexec from a different MPI than the one compiled against (Intel MPI vs OpenMPI); (d) batch script requested fewer cores than nPx*nPy; (e) stale SIZE.h / wrong executable in run dir (use `ln -s ../build/mitgcmuv`, not cp).
- Fix: `genmake2 -mpi -mods ../code ...`, `make CLEAN` after changing SIZE.h, run `mpirun -np $((nPx*nPy)) ./mitgcmuv`; check `which mpirun mpif77` come from the same install; confirm `usingMPI = T` in STDOUT.0000 header. Test MPI with a hello-world first. Note current upstream wording when usingMPI=F: "INI_PROCS: needs MPI for multi-procs (nPx*nPy=..) setup".
- Era: 2009-2024, perennial; message text still in eesupp/src/eeboot_minimal.F / ini_procs.F.
- Src: mitgcm-support 2009-April 'No. of processes not equal to nPx*nPy' http://mailman.mitgcm.org/pipermail/mitgcm-support/2009-April/thread.html ; 2015-April 'error with "No. of procs not equal to nPx*nPy"' http://mailman.mitgcm.org/pipermail/mitgcm-support/2015-April/thread.html ; 2018-August 'An error about procs when I submitted a MITgcm job' http://mailman.mitgcm.org/pipermail/mitgcm-support/2018-August/thread.html ; 2014-March 'Independent Tiling'

### SIZE.h edits ignored with -mpi / "SIZE.h_mpi" silently overrides SIZE.h
- Cause: with `genmake2 -mpi`, files named `*_mpi` in the -mods dir (SIZE.h_mpi, CPP_EEOPTIONS.h_mpi) are linked in and renamed without the suffix, taking priority over SIZE.h.
- Fix: delete or edit `SIZE.h_mpi` in your code dir (verification experiments ship one; copy-over is automatic only with -mpi).
- Era: behaviour added by JMC ~April 2009; still current.
- Src: mitgcm-support 2009-April 'No. of processes not equal to nPx*nPy' (JMC) http://mailman.mitgcm.org/pipermail/mitgcm-support/2009-April/thread.html

### EXCH-1 useCubedSphereExchange unsafe with usingMPI=.TRUE.
- Cause: `useCubedSphereExchange=.TRUE.` in eedata but pkg/exch2 not compiled (often because packages.conf was not picked up: forgot `-mods=../code`). Check `#define ALLOW_EXCH2` in PACKAGES_CONFIG.h.
- Fix: add `exch2` to packages.conf, pass `-mods`, rebuild from clean build dir. useCubedSphereExchange is now set automatically by exch2 for cs/llc topologies.
- Era: 2014; check still valid for old checkpoints.
- Src: mitgcm-support 2014-July 'cubeSphereExchange and MPI error?' http://mailman.mitgcm.org/pipermail/mitgcm-support/2014-July/thread.html

### EESET_PARMS: Error reading parameter file "eedata" / EEDIE
- Cause: namelist syntax. In eedata/data comment lines must have `#` in column 1 (not column 2); unknown parameter name also triggers it; real message is in STDOUT.0000 (last echoed namelist line).
- Fix: fix the namelist; read STDOUT.0000 to see how far parsing got ("Cannot match namelist object name" / "syntax error in NAMELIST input" are the same family).
- Era: 2014-2015.
- Src: mitgcm-support 2014-July 'cubeSphereExchange and MPI error?' ; 2015-February 'restart .meta & .data files from a different processor decomposition?' (Martin Losch) http://mailman.mitgcm.org/pipermail/mitgcm-support/2015-February/thread.html

### Compile/link: "cannot find include file mpif.h", "mpif77: command not found", "/usr/bin/ld: cannot find -lmpi"
- Cause: genmake2/cpp does not use mpif77 for preprocessing, so MPI include dir must be given explicitly; -lmpi in LIBS is redundant when FC=mpif77/mpif90.
- Fix: `export MPI_INC_DIR=/path/to/mpi/include` (honoured by most optfiles; see header of build_options/linux_amd64_gfortran) or edit INCLUDES/INCLUDEDIRS in the optfile; remove `-lmpi` from LIBS; make sure the MPI wrapper is on PATH on compute nodes. `mpif77 --show` prints the include dir. Same family: "mpif.h not found" on Derecho with cray-mpich (use intel-mpi module + MPI_ROOT paths, or `FC=ftn`).
- Era: 2005-2023 (Derecho 2023-Nov).
- Src: mitgcm-support 2009-February 'Problem with MPI (cannot find -lmpi)' http://mailman.mitgcm.org/pipermail/mitgcm-support/2009-February/thread.html ; 2011-November 'Problem with parallel build' ; 2015-October 'build MITgcm with mpi' ; 2023-November 'Derecho environment?' http://mailman.mitgcm.org/pipermail/mitgcm-support/2023-November/thread.html ; GitHub issue #82

### Fatal error in PMPI_Comm_rank: Invalid communicator / symbol lookup error libopen-pal.so ... undefined symbol __intel_sse2_strlen
- Cause: compiled against one MPI (or Intel runtime) and launched with another; with `-shared-intel` the Intel runtime libs must be on LD_LIBRARY_PATH on every node.
- Fix: use mpirun/mpiexec from the same MPI used to build (`which mpirun`, `mpif90 --show`); source ifortvars.sh in the job script or switch `-shared-intel` to `-static-intel` in the optfile. Test a hello-world MPI program first.
- Era: 2011, 2016, 2018.
- Src: mitgcm-support 2018-September 'Fatal error in PMPI_Comm_rank: Invalid communicator, error stack' http://mailman.mitgcm.org/pipermail/mitgcm-support/2018-September/thread.html ; 2011-March 'lib64/libopen-pal.so.0: undefined symbol: __intel_sse2_strlen' ; 2016-January 'Troubleshooting OpenMPI Issues with mpiexec for Jasper'

### MPI_Recv: Invalid tag (coupled runs on Cray XC30 / ARCHER)
- Cause: pkg/compon_communic/generate_tag.F built MPI tags larger than the machine's MPI_TAG_UB (4194303 on Cray XC30).
- Fix: replace generate_tag.F with Jean-Michel's modified version (tags below the limit); Chris Hill suggested plain `iarg1+iarg2` would suffice. Only matters for coupled atm/ocn/seaice setups.
- Era: 2014-October; check whether current generate_tag.F still has the hash (file exists in upstream).
- Src: mitgcm-support 2014-October 'MPI problem on Archer (CRAY XC30)' http://mailman.mitgcm.org/pipermail/mitgcm-support/2014-October/thread.html

## Hangs, aborts, output flushing
### Job hangs until walltime after a model STOP (e.g. "MOM_IMPLICIT_R: error when solving 3-Diag problem", missing pickup)
- Cause: model errors end with Fortran `STOP`, which on some systems (Cray, BlueGene) leaves the other ranks waiting in MPI. `ALL_PROC_DIE` (calls MPI_FINALIZE) is only used when all ranks reach it; ~1400 STOP statements remain, no MPI_ABORT path.
- Fix: no model-side fix. Make the batch job robust (wall-clock watchdog, `mpiexec` kill-on-first-exit option, check STDERR in the job script). Do not patch STOP globally: ALL_PROC_DIE only works if every rank calls it; core devs propose a separate emergency MPI_ABORT routine.
- Era: 2007 and 2021; GitHub issue #439 still open.
- Src: GitHub issue #439 https://github.com/MITgcm/MITgcm/issues/439 ; mitgcm-support 2007-October 'Clean exit from errors during MPI runs' http://mailman.mitgcm.org/pipermail/mitgcm-support/2007-October/thread.html

### STDOUT.* / STDERR.* not updating until job ends (Cray compiler, ARCHER)
- Cause: Fortran buffering of the stdout unit.
- Fix: in eedata set `debugMode=.TRUE.,` (flushes STDOUT after each write) and also set `debugLevel=2,` explicitly because debugMode raises the default debugLevel to 4. Dan Jones' extra tip: reduce per-core memory (smaller tiles) if low memory distorts output behaviour.
- Era: 2019-February; still valid (debugMode exists in eeset_parms.F).
- Src: mitgcm-support 2019-February 'ARCHER, cray, mpi and STDOUT' http://mailman.mitgcm.org/pipermail/mitgcm-support/2019-February/thread.html

### No STDOUT.* files; all ranks print to screen in a jumble
- Cause: file naming is by the model but redirection is by the MPI launcher / batch system; also without MPI stdout goes to the terminal (use `./mitgcmuv > output.txt`).
- Fix: use launcher options for per-rank output (or leave the default STDOUT.NNNN files that eeboot_minimal.F writes with usingMPI=T); `SINGLE_DISK_IO` in CPP_EEOPTIONS.h keeps only rank 0 STDOUT/STDERR (all messages from other ranks are lost, use only on a proven set-up).
- Era: 2004, 2015, 2020.
- Src: mitgcm-support 2004-August 'file size issue with mpi' http://mailman.mitgcm.org/pipermail/mitgcm-support/2004-August/thread.html ; 2015-July 'Poor cpu usage percentage'

## Stack / memory / size limits
### Segmentation fault at start or in first big routine (do_oceanic_phys, mom_calc_visc, dynamics) after raising Nr or tile size
- Cause: large local (automatic) arrays exceed the stack limit, or static executable too big for node memory.
- Fix: `ulimit -s unlimited` (csh: `limit stacksize unlimited`) in the batch script on compute nodes; ifort `-mcmodel=medium` (add `-shared-intel`); `size mitgcmuv` for a memory estimate; reduce tile size/diagnostics/packages; rebuild with `genmake2 -devel` for bounds checks. For OpenMP also `OMP_STACKSIZE`/`MP_STACKSIZE`.
- Era: 2003-2018, perennial.
- Src: mitgcm-support 2018-May 'segmentation fault' http://mailman.mitgcm.org/pipermail/mitgcm-support/2018-May/thread.html ; 2007-August 'Chasing a seg fault' http://mailman.mitgcm.org/pipermail/mitgcm-support/2007-August/thread.html ; 2003-October 'OpenMP Model Code.'

### Link error "relocation truncated to fit: R_X86_64_PC32 / R_X86_64_32S" or "32-bit RIP relative reference out of range" with big sNx*sNy
- Cause: static arrays larger than 2 GB (small code model); more procs (smaller tiles) makes it vanish.
- Fix: `-mcmodel=medium` (or `large`) in FFLAGS and CFLAGS (compiler-specific; gfortran/ifort); otherwise increase nPx*nPy, or split with tiles per process (e.g. sNy=800, nSy=2 instead of sNy=1600, nSy=1). Patrick Heimbach: same message appears when executable size per processor exceeds the node limit.
- Era: 2007, 2013.
- Src: mitgcm-support 2007-April 'Limitation for Tile's dimensions' http://mailman.mitgcm.org/pipermail/mitgcm-support/2007-April/thread.html ; 2013-December 'domain size limit?' http://mailman.mitgcm.org/pipermail/mitgcm-support/2013-December/thread.html

### Reduce memory per process
- Cause: static allocation of everything: diagnostics qdiag, timeave, second-order-moment advection arrays, ptracers etc.
- Fix: more MPI procs (smaller tiles); drop unused packages; make sure GAD second-order-moments advection is not compiled (GAD_OPTIONS.h); reduce numbers of diagnostics; hybrid MPI+OpenMP. Virtual size shown by top/ps is not a good indicator; what matters is not swapping. A steadily growing RSS in an MPI run was traced to an OpenMPI/netcdf install, not MITgcm (MITgcm has no dynamic allocation apart from parts of ptracers).
- Era: 2006-2015.
- Src: mitgcm-support 2014-August 'Memory issues with regional simulation' http://mailman.mitgcm.org/pipermail/mitgcm-support/2014-August/thread.html ; 2015-July 'Memory usage keeps increasing' ; 2006-September 'MITgcm virtual memory usage'

## Decomposition-dependent results (reproducibility)
### Results differ when only nPx/nPy/sNx/sNy (or processor count) change
- Cause: roundoff. Global sums in the pressure solver (cg2d and cg3d) add in a different order for a different tiling; other global sums (OBCS balance, seaice LSR solver, forcing balance) too; aggressive compiler optimisation changes operation order. In a chaotic/eddying flow small differences grow to O(1) locally while statistics stay the same.
- Fix: for identical tiling but different process mapping (nSx,nSy,nPx,nPy) keep `#define GLOBAL_SUM_ORDER_TILES` (default today; added Aug 2015). For tiling-independent sums use `#define CG2D_SINGLECPU_SUM` in eesupp/inc/CPP_EEOPTIONS.h (slow; hydrostatic cg2d only, no OBCS balance, no seaice dynamics). Test with -O0 if differences are suspiciously large. seaice: LSR solver is tiling dependent; EVP assumed independent.
- Era: 2005-2024; GLOBAL_SUM_ORDER_TILES since Aug 2015.
- Src: mitgcm-support 2014-November 'Zonal filter with nPx>1' http://mailman.mitgcm.org/pipermail/mitgcm-support/2014-November/thread.html ; 2016-April 'results quite differents depending on number of procs used' http://mailman.mitgcm.org/pipermail/mitgcm-support/2016-April/thread.html ; 2024-November 'results differences with number of processors' http://mailman.mitgcm.org/pipermail/mitgcm-support/2024-November/thread.html

### Non-hydrostatic run still not decomposition-reproducible with CG2D_SINGLECPU_SUM
- Cause: CG2D_SINGLECPU_SUM only touches cg2d; the cg3d global sums are not covered (and the run is >10x slower anyway).
- Fix: none available in code; treat as roundoff (Martin: "possible, but not yet done").
- Era: 2024-May; CG3D_SINGLECPU_SUM does not exist upstream.
- Src: mitgcm-support 2024-May 'Domain decompositions affecting simulations outcome' http://mailman.mitgcm.org/pipermail/mitgcm-support/2024-May/thread.html

### OBCS sponge layer differs with tile size (strip of non-zero values at inner edge of sponge)
- Cause: OBCS operations (including sponge) are carried out only on tiles that contain an open boundary; if a tile is smaller than the sponge thickness (e.g. sNy=12 vs 15-cell sponge) the next tile is not an OBCS tile and is not touched. Core dev calls it a bug.
- Fix: choose tile size larger than the sponge width (spongeThickness). Not fixed in code as far as these sources say; a warning/stop was suggested.
- Era: 2024-May.
- Src: mitgcm-support 2024-May 'Domain decompositions affecting simulations outcome' (Martin Losch) http://mailman.mitgcm.org/pipermail/mitgcm-support/2024-May/thread.html

### Home-made sponge in external_forcing.F repeats at every tile boundary
- Cause: i,j are local to the tile; hard-coded `i=1,100` applies to every tile.
- Fix: use global indices (myXGlobalLo-1+i) or masks from a file; see how pkg/obcs applies OB_Jn/OB_Ie only at global boundaries; or use pkg/rbcs.
- Era: 2004; still a classic bug.
- Src: mitgcm-support 2004-August 'question on sponge-layer/mpi' http://mailman.mitgcm.org/pipermail/mitgcm-support/2004-August/thread.html

### Spurious standing wave / instability with wavenumber nPx, or blow-up on a tile boundary
- Cause: three separate real cases. (1) Perfectly symmetric flat-bottom set-up where -O2 roundoff triggers/doesn't trigger baroclinic instability differently; adding noise to initial T gave identical MPI/non-MPI results. (2) ifort 13.0 miscompile with odd tile width (sNx=47 blew up at iteration 68; sNx=48 fine; gfortran and older ifort fine). (3) Odd/even tile-size sensitivities with variable dx.
- Fix: perturb initial conditions; compare -O0/-O1/-O2; try another compiler version; use even tile sizes as a workaround; compile with pkg/exch2 to see if it changes anything.
- Era: 2013 (c63o), 2014-2015 (c64u).
- Src: mitgcm-support 2013-June 'problems at tile bdys?' http://mailman.mitgcm.org/pipermail/mitgcm-support/2013-June/thread.html ; 2015-January 'Baroclinic instability with MPI run' http://mailman.mitgcm.org/pipermail/mitgcm-support/2015-January/thread.html ; 2014-March 'Independent Tiling'

### cg2d converges for one processor layout but hits cg2dMaxIters on another (cg2d_res stuck ~1e-6)
- Cause: ill-conditioned problem (point-source forcing in a large domain) where roundoff in the partial sums matters; not an MPI bug.
- Fix: tighten/loosen `cg2dTargetResidual`, raise `cg2dMaxIters`, rule out roundoff via GLOBAL_SUM_ORDER_TILES comparison; check the setup is not nearly unstable.
- Era: 2005.
- Src: mitgcm-support 2005-July 'cg convergence vs processor count' http://mailman.mitgcm.org/pipermail/mitgcm-support/2005-July/thread.html

### cg2dTargetResWunit / cg3dTargetResWunit tolerance scaled wrongly (off by sqrt(Nxy))
- Cause: conversion of the W-unit target to solver tolerance ignored the number of wet surface points.
- Fix: PR #959 (merged, see doc/tag-index "fix scaling from cg2d/cg3d TargetResWunit"): ini_cg2d.F now multiplies by globalArea/sqrt(n2dWetPts). Existing configs using *TargetResWunit need the value corrected to reproduce old results. Non-dimensional `cg2dTargetResidual` (RHS normalisation) is unaffected; recommended for adjoint runs because the adjoint CG problem does not match W scaling.
- Era: fixed Dec 2025 - Jan 2026 (PR #959 / issue #956); older checkpoints (incl. ECCO v4r4 era) keep the old scaling.
- Src: GitHub PR #959 https://github.com/MITgcm/MITgcm/pull/959 ; issue #956

## Large processor counts
### Tile suffix shows ".***.001.data" or "scratch.****", STDOUT.**** (nPx*nSx > 999 or > 9999 procs)
- Cause: fixed-width formats (I3.3 for tile index in mdsio, I4.4 for STDOUT/STDERR names, I3 in ini_procs) in older code.
- Fix: update code (PR #345 fixed ini_procs.F formats, issue #343; scratch files already use FMT_PROC_ID 'I9.9' in eesupp/inc/CPP_EEMACROS.h; ini_procs.F now uses I6). On older checkpoints (c65x) edit eeboot_minimal.F / eeset_parms.F / open_copy_data_file.F I4.4 -> I5.5, or `#define SINGLE_DISK_IO`. For tile indices >999 with exch2, put tiles in both directions (nPx*nSx<1000 and nPy*nSy<1000) or use single-CPU I/O.
- Era: 2014-2020; fixed upstream by PR #345 (2020).
- Src: GitHub issue #343 https://github.com/MITgcm/MITgcm/issues/343 ; mitgcm-support 2020-March 'mpi run with cpu more than 9999' http://mailman.mitgcm.org/pipermail/mitgcm-support/2020-March/thread.html ; 2014-July 'nPx*nSx > 1000' ; 2014-December 'exch2 with lat lon grid' http://mailman.mitgcm.org/pipermail/mitgcm-support/2014-December/thread.html

### EEBOOT_MINIMAL stops when requesting >MAX_NO_PROCS processes (old: 128 -> 1024 -> 2048)
- Cause: hard limit in eesupp/inc/EEPARAMS.h.
- Fix: raise MAX_NO_PROCS (no other consequence).
- Era: obsolete: MAX_NO_PROCS no longer appears in eesupp (only in doc/tag-index history); not an issue in current code.
- Src: mitgcm-support 2010-March 'max number of processors' http://mailman.mitgcm.org/pipermail/mitgcm-support/2010-March/thread.html

### MNC: more than 999 tiles / MNC_MAX_ID too small / gluemncbig "found 999 tiles, need N"
- Cause: pkg/mnc/MNC_SIZE.h index limits (MNC_MAX_ID now 3000); model writes correct t1012.nc style names but the glue utility mis-counted.
- Fix: raise MNC_MAX_ID/MNC_MAX_FID in MNC_SIZE.h if hit; use a recent gluemncbig (utils/python/MITgcmutils/scripts) or read tiles individually; for big runs prefer MDS + exch2 I/O.
- Era: 2004 (MNC_MAX_ID=1000), 2013.
- Src: mitgcm-support 2004-December 'MNC_MAX_ID too small for 24 tiles' ; 2013-January 'netcdf with more than 999 tiles' http://mailman.mitgcm.org/pipermail/mitgcm-support/2013-January/thread.html

### useMNC unsafe with multi-threads (OpenMP+mnc writes mdsio instead)
- Cause: pkg/mnc is not thread safe; mnc_readparms.F switches it off if nThreads>1 ("useMNC unsafe with multi-threads").
- Fix: use MPI instead of OpenMP, or mdsio.
- Era: 2011, still in code.
- Src: mitgcm-support 2011-October 'Run mitgcm with openMP and mnc' http://mailman.mitgcm.org/pipermail/mitgcm-support/2011-October/thread.html

## exch2, cubed-sphere/llc, blank tiles
### How to use blank (land-only) tiles to save cores
- Cause/recipe: pkg/exch2 must be compiled (no run-time switch; compiled = used).
- Fix: (1) add exch2 to packages.conf; (2) SIZE.h: nPx = number of ACTIVE tiles (procs), nPy=1 (tiles may be put in either direction for non-cs grids); (3) data.exch2 (not eedata): `blankList = ...` and `dimsFacets = Nx_full, Ny_full` giving the FULL domain size; for plain lat-lon also `W2_mapIO = 1`; (4) get the list by running one time-step with an un-blanked executable of same sNx,sNy and grep "Empty tile: #" in STDOUT.* (or utils/exch2/matlab-topology-generator/generate_blanklist.m); (5) run with mpirun -np = number of active tiles. Input files keep full-domain size. Examples: verification/global_ocean.90x40x15 (input/data.exch2.mpi, code/SIZE.h_mpi), verification/adjustment.cs-32x32x1.
- Era: 2011-2022 (non-curvilinear + blank tiles supported since Dec 2011).
- Src: mitgcm-support 2022-November '"blanklist" parameter in MITgcm' (Jean-Michel Campin) http://mailman.mitgcm.org/pipermail/mitgcm-support/2022-November/thread.html ; 2019-September 'Excluding Tiles that are All Land' ; 2014-December 'exch2 with lat lon grid' ; GitHub PR #382

### S/R W2_SET_MAP_TILES: Domain Total # of tiles = 468 does not match (SIZE.h+blankList)= 492 (also "ABNORMAL END: S/R W2_SET_MAP_TILES")
- Cause: with exch2, nPx*nPy*nSx*nSy + number of blank tiles must equal the total number of tiles in the topology. Halving sNx/sNy does NOT simply double process count because the blankList changes with tile size.
- Fix: use the blankList that matches your sNx,sNy (comment/uncomment per-size lists in data.exch2), and set nPx so active + blank = total (ECCOv4 LLC90: sNx=sNy=15 -> nPx=360, nPy=1; 360+108 blank = 468). Also tile size must divide the facet dimensions (e.g. sNx=75/sNy=17 does not divide 510).
- Era: 2015-2016.
- Src: mitgcm-support 2016-December 'Changing the number of cores in ECCOv4' http://mailman.mitgcm.org/pipermail/mitgcm-support/2016-December/thread.html ; 2015-March 'crash with a new processor / grid size setup' http://mailman.mitgcm.org/pipermail/mitgcm-support/2015-March/thread.html

### Patchy values / noise at corners of tiles on cubed-sphere with small tiles
- Cause: exch2 bug when GCD(sNx,sNy) < 2*OLx (initially thought to be only for cs-32); not an nPy>1 problem.
- Fix: fixed in pkg/exch2 on 2012-03-26, available from checkpoint63l on. Old checkpoints: keep GCD(sNx,sNy) >= 2*OLx.
- Era: 2012; fixed in checkpoint63l.
- Src: mitgcm-support 2012-April 'arbitrary tiles in cube-sphere grid' http://mailman.mitgcm.org/pipermail/mitgcm-support/2012-April/thread.html

### Blank tiles: output files contain zeros / are not smaller
- Cause: blank tiles are not computed but remain part of the global MDS output (zeros); MITgcm does no compression.
- Fix: compress after the run, or redefine the domain (move grid origin) so the region of interest has fewer blank facets.
- Era: 2025-March.
- Src: mitgcm-support 2025-March 'Zeros in inactive tiles in the regional model' http://mailman.mitgcm.org/pipermail/mitgcm-support/2025-March/thread.html

### Blank tiles on a lat-lon grid give wrong YC/rA/Coriolis (all tiles in lowest row), or global file too short when last tile is blank
- Cause: old exch2 handled grid coordinates only for the cs/curvilinear case.
- Fix: fixed Dec 2011 for non-curvilinear grids (global_ocean.90x40x15 as example). Before that, workaround: selectCoriMap=3 with fCoriC.bin, fCoriG.bin, fCorCs.bin files.
- Era: 2011; fixed Dec 2011.
- Src: mitgcm-support 2011-December 'EXCH2 grid parameter question clarification' http://mailman.mitgcm.org/pipermail/mitgcm-support/2011-December/thread.html

### w2_e2setup.f: This name does not have a type (nTiles = W2_maxNbTiles); old matlab driver.m topology files
- Cause: the generated w2_e2setup.f / matlab-topology-generator technology was retired; exch2 topology is now computed at run time from SIZE.h + data.exch2.
- Fix: delete the old w2_*.f / w2_*.F files from the code dir, then `make makefile && make CLEAN && make depend && make`.
- Era: 2015-March (c65j); manual page at the time was out of date.
- Src: mitgcm-support 2015-March 'w2_e2setup.f: This name does not have a type nTiles = W2_maxNbTiles' http://mailman.mitgcm.org/pipermail/mitgcm-support/2015-March/thread.html ; 2015-October 'exch2 and autogenerated grids (llc, cs, etc.)'

### forrtl: severe (36): attempt to access non-existent record ... tile001.mitgrid
- Cause: grid file shorter than SIZE.h/data.exch2 expects.
- Fix: cs grid file size must be 8*(Nc+1)*(Nc+1)*16 bytes per face (Nc = face size, 16 fields, always real*8 even if readBinaryPrec=32). Check SIZE.h, useCubedSphereExchange, exch2 compiled and data.exch2.
- Era: 2013.
- Src: mitgcm-support 2013-August 'forrtl: severe (36) with tile file' http://mailman.mitgcm.org/pipermail/mitgcm-support/2013-August/thread.html

### Using one global tile001.mitgrid with multiple tiles on a Cartesian/horizGridFile grid
- Cause: without exch2 each tile needs its own tileNNN.mitgrid (EXCH1 only supports one tile/file mapping).
- Fix: compile pkg/exch2 with default parameters (no data.exch2 needed) and supply one full-domain tile001.mitgrid; works for Nx x Ny with nSx/nPx tiling.
- Era: 2016-June.
- Src: mitgcm-support 2016-June 'Internal-Wave: Whether I could enable MPI with using only one tile001.mitgrid (LLC grid)' http://mailman.mitgcm.org/pipermail/mitgcm-support/2016-June/thread.html

### pkg/flt (floats) with exch2/LLC
- Cause: float initialisation (per-tile) and FLT_MAP_XY2IJLOCAL for cs/llc incomplete.
- Fix: no complete solution; see PR #282 and Spencer Jones' llc90 float set-up (cspencerjones/float_stuff) as starting points.
- Era: 2019-2022, issues #283 and #672 open.
- Src: GitHub issue #283 https://github.com/MITgcm/MITgcm/issues/283 ; issue #672

## Changing processor count / restart / I/O layout
### Restart with different nPx/nPy (pickup from other decomposition)
- Cause: per-tile pickup files (pickup.*.001.001.data) are tied to the tiling.
- Fix: works directly if pickups are global (written with `useSingleCpuIO=.TRUE.` or globalFiles). Otherwise join tiles first (rdmds then fwrite 'real*8' to pickup.NNNNNNNNNN.data; keep .meta consistent, per Jean-Michel) or use utils/matlab/interpickups.m (mnc pickups). Nx,Ny must stay the same; recompile (static arrays). The "No. of processes not equal" line seen in this scenario is just the old executable.
- Era: 2006-2015.
- Src: mitgcm-support 2015-February 'restart .meta & .data files from a different processor decomposition?' http://mailman.mitgcm.org/pipermail/mitgcm-support/2015-February/thread.html ; 2006-March 'changing the number of tiles' ; 2009-December 'pickup after changing number of processors'

### Output split in files like U.0000000010.001.002.data (or "output files not separated into 48 subfiles")
- Cause: default is one file per tile (useSingleCpuIO=.FALSE., globalFiles=.FALSE.). If you want per-tile files but get one, useSingleCpuIO=.TRUE. is set.
- Fix: rdmds/read_mds glue tiles using .meta (read the base name, not a tile suffix). `useSingleCpuIO=.TRUE.` gathers to rank 0 and writes one file per field (robust; many times faster than globalFiles on SGI). `globalFiles=.TRUE.` (all ranks write one file) is unsafe on many MPI/filesystem combos (wrong size or contents); avoid it.
- Era: 2004-2024.
- Src: mitgcm-support 2004-August 'file size issue with mpi' ; 2014-July 'Outputting too many files' http://mailman.mitgcm.org/pipermail/mitgcm-support/2014-July/thread.html ; 2024-August 'Raw binary output error when using MPI' http://mailman.mitgcm.org/pipermail/mitgcm-support/2024-August/thread.html

### rdmds('XG.001.001') returns a larger array than the tile (108x30 instead of 54x30)
- Cause: rdmds assembles tiles from .meta info; asking for one tile name pads the global domain with zeros.
- Fix: read the base name (rdmds('XG')); read raw tile with fread of [sNx sNy] if needed.
- Era: 2008.
- Src: mitgcm-support 2008-October 'tile size oddity? (2nd attempt)' http://mailman.mitgcm.org/pipermail/mitgcm-support/2008-October/thread.html

### ZONAL_FILT_INIT: Multi-tiles ( nSx*nPx= N ) in Zonal (X) dir. not implemented in Zonal-Filter code
- Cause: pkg/zonal_filt cannot be decomposed in x.
- Fix: nPx=1 and nSx=1 (put all processes in y); or write a gather/filter/scatter version (Roland Young did, not upstream).
- Era: 2012-2014; message still in zonal_filt_init.F.
- Src: mitgcm-support 2012-July 'nPx vs. nPy' http://mailman.mitgcm.org/pipermail/mitgcm-support/2012-July/thread.html ; 2014-November 'Zonal filter with nPx>1'

## Performance and scaling
### MPI_Allreduce dominates runtime with thousands of tiles (global_sum_tile / cg2d), e.g. llc540 on Pleiades
- Cause: default `#define GLOBAL_SUM_ORDER_TILES` does an Allreduce on a vector with one slot per tile (mostly zeros) to keep sum order independent of processors; gets slow for huge nPx*nPy.
- Fix: in eesupp/inc/CPP_EEOPTIONS.h `#undef GLOBAL_SUM_ORDER_TILES` and `#undef GLOBAL_SUM_SEND_RECV` (plain scalar MPI_Allreduce per process; loses processor-count reproducibility). Dan Kokron measured comm time 342 s -> 99 s (llc_540, 5 days). A SPEED_AT_ALL_COST-style option was agreed in the thread; it is not present in upstream today, so edit the two flags directly. Allgather in place of Allreduce was also tested.
- Era: 2020-November (checkpoint67q); issue closed.
- Src: GitHub issue #385 https://github.com/MITgcm/MITgcm/issues/385

### Same Allreduce problem in a package: thousands of GLOBAL_SUM_RL calls (profiles_make_ncfile.F, obsfit) make MPI adjoint 4-10x slower than serial
- Cause: profiles_make_ncfile.F does NFILESPROFMAX*NVARMAX*NOBSGLOB*NLEVELMAX*2 (35 M in global_oce_biogeo_bling) scalar global sums.
- Fix: replace by one MPI_Allreduce on the whole prof_buff/prof_mask_buff (Martin's hack: 155 s -> 1.5 s); or skip empty profilesfiles, or use a gather to rank 0. Not yet merged upstream (issue open as of 2026-07).
- Era: 2026; issue #1006 open. Lesson: never call GLOBAL_SUM_* inside loops over obs.
- Src: GitHub issue #1006 https://github.com/MITgcm/MITgcm/issues/1006

### Need a global sum over only some tiles, or vector (element-wise) global sum, in a new package
- Cause: GLOBAL_SUM_TILE_* sums scalars over all tiles.
- Fix: multiply by a mask zero outside the tiles of interest and use global_sum_tile (inefficient but simple); for vectors use eesupp/src/global_sum_vector.F (the thread names it global_vec_sum.F, which does not exist upstream); or MPI_Gather to root and sum there (Samar Khatiwala). There is no routine that sums a subset of tiles.
- Era: 2023-June.
- Src: mitgcm-support 2023-June 'Question about MPI sums' http://mailman.mitgcm.org/pipermail/mitgcm-support/2023-June/thread.html

### Adjoint/Pleiades slow, Lustre load: many small file open/close in a tight loop (PLEIADES_LUSTRE_OPT1/2/3)
- Cause: repeated open/read/close of small files and many rank-0 writes on Lustre.
- Fix: no code needed for part: `adTapeDir` (PARM05) to a fast local dir (about 8-15% gain), `mdsioLocalDir` (but it is ignored by mdsio_write_field/tape when useSingleCpuIO=.TRUE.), `profilesDir`, `diagMdsDir`, run executable from node-local /tmp with inputs copied; `SINGLE_DISK_IO` brought ~0.1% (namelist reading is not the bottleneck). Dan Kokron's dkokron/MITgcm checkpoint67l_fileHash branch (PLEIADES_LUSTRE_OPT, OPT2 = MDS_READ_FIELD fast path with open action='read', OPT3) gave ~33% on 3-year adjoint runs (OPT2 ~17%, OPT ~14%, OPT3 ~7%). These macros are not in upstream.
- Era: 2021-2022, issue #535 still open.
- Src: GitHub issue #535 https://github.com/MITgcm/MITgcm/issues/535

### BLOCKING_EXCHANGES takes 80% of time after enabling pkg/ptracers (8 tracers)
- Cause: PTRACERS_FIELDS_BLOCKING_EXCH exchanges each tracer 3-D field separately; on Columbia this was ~100x slower per field than theta/salt.
- Fix: Martin's workaround: copy all ptracers into one 4D/5D array with k index of length nPtracers*Nr and do a single exchange (needs larger exch buffers); commenting the call out only for tests. Status upstream not stated in thread.
- Era: 2010-March.
- Src: mitgcm-support 2010-March 'BLOCKING_EXCHANGES slowdown when using pkg/ptracers on Columbia' http://mailman.mitgcm.org/pipermail/mitgcm-support/2010-March/thread.html

### Implicit vertical solver slow (solve_*diagonal.F)
- Cause: default variant accesses more 3-D arrays; options in model/inc/CPP_OPTIONS.h.
- Fix: `#define SOLVE_DIAGONAL_LOWMEMORY` was ~2x faster for llc_540 in forward runs; not suitable for AD. `SOLVE_DIAGONAL_KINNER` is AD-suitable but slowest (non-vectorising). Martin offered a branch (mjlosch/MITgcm diagonal_lowmemory) hoisting if-statements out of loops.
- Era: 2021-January.
- Src: mitgcm-support 2021-January 'SOLVE_DIAGONAL_LOWMEMORY' http://mailman.mitgcm.org/pipermail/mitgcm-support/2021-January/thread.html

### Poor scaling, wall-clock >> user time (CG2D, BLOCKING_EXCHANGES dominate); rules of thumb
- Cause: tile too small (rule: below ~30x30 gridpoints, with OLx,OLy=3-5 even ~50 for 2D solvers, overlap and the cg2d exchanges dominate); slow interconnect or MPI not using shared-memory comms; memory-bandwidth limited multi-core chips; I/O; shared nodes; swap.
- Fix: read the timing summary at the end of STDOUT.0000 (SOLVE_FOR_PRESSURE, BLOCKING_EXCHANGES, DO_THE_MODEL_IO); keep tiles near-square and >=30x30 (include overlaps and blank tiles in scaling estimates); if system time is within 10x of user time something is wrong with MPI/interconnect; test with I/O off (monitorFreq large, debugLevel=-1); use fewer cores per socket if bandwidth bound. Super-linear speedup is normal (cache effects, better until ~250 procs for a 600x800x50 channel). ECCOv4 LLC90 adjoint on ARCHER scaled linearly 96/192/360 procs (15.6/7.8/3.8 h); nchecklev tuning gave only ~5%.
- Era: 2005-2020.
- Src: mitgcm-support 2012-May 'speedup for cs64 on a linux cluster' http://mailman.mitgcm.org/pipermail/mitgcm-support/2012-May/thread.html ; 2008-October 'scaling on SGI altix, AMD quad core cpus is terrible' ; 2013-June 'MPI speed issues' ; 2020-October 'Better than expected HPC Scaling' ; 2016-December 'Changing the number of cores in ECCOv4'

### Mac laptop/iMac MITgcm run slow, CPU 2-35%
- Cause: not MPI: I/O (forcing reads, MONITOR) and compiler/optfile differences (e.g. -O0 vs -O2).
- Fix: compare timings at end of STDOUT; set `monitorFreq=0`/larger; put forcing on fast drive; compare Makefile flags. (`make -j` only speeds the build.)
- Era: 2015.
- Src: mitgcm-support 2015-July 'Poor cpu usage percentage' http://mailman.mitgcm.org/pipermail/mitgcm-support/2015-July/thread.html

## Domain-size and SIZE.h rules
### Array out of bounds / wrong results with OLx or OLy too small
- Cause: OLx,OLy must cover the advection scheme stencil and at least 2 (1 only offline); old checks were lax (OLx=1 caused out-of-bounds under -devel).
- Fix: PR #473 (merged June 2021) improved the minimum-overlap check in the config check; use OLx,OLy >= 2 and at least what the advection scheme needs (the check reports it). Larger overlaps raise overhead and hurt scaling (Martin saw seaice solver scaling taper near 50-point tiles with OLx=5).
- Era: fixed in PR #473 (2021); older checkpoints silent.
- Src: GitHub issue #125 https://github.com/MITgcm/MITgcm/issues/125

### Nx/Ny, tile terminology (what to set for a given core count)
- Cause: confusion between tiles and processes; SIZE.h: Nx = sNx*nSx*nPx, Ny = sNy*nSy*nPy; nPx*nPy = MPI processes; nSx*nSy = tiles per process (one thread per tile with OpenMP); sNx,sNy = tile size.
- Fix: pick factors such that products equal the full grid; multi-tile per process is mostly for threads, cubed-sphere EXCH1 history, and testing (many verification experiments use more tiles than sensible). Always `make CLEAN`/rebuild when changing SIZE.h (static allocation). With exch2 these identities hold only with blank tiles added (see W2_SET_MAP_TILES entry).
- Era: all; wording cleaned in issue #100.
- Src: GitHub issue #100 https://github.com/MITgcm/MITgcm/issues/100 ; mitgcm-support 2012-July 'nPx vs. nPy' http://mailman.mitgcm.org/pipermail/mitgcm-support/2012-July/thread.html ; 2005-February 'tiles vs processors'

## OpenMP / threading
### Threaded (OpenMP) runs: little gain vs MPI, packages break threading
- Cause: threaded code is less supported; many packages contain constructs that break threading; mnc is disabled; NOTE for threads each of nSx*nSy tiles per process is handled by a thread (nTx*nTy threads = nSx*nSy).
- Fix: prefer MPI even on shared memory (SGI Altix/Origin experience); if using OpenMP set `nTx,nTy` in eedata consistent with nSx,nSy, genmake2 `-omp`; hybrid MPI+OpenMP can lower memory/process. Verification test: eedata.mth in an experiment (PR #834 notes multi-threaded tests must keep nTx,nTy consistent with MPI secondary tests).
- Era: 2003-2024.
- Src: mitgcm-support 2007-October 'Clean exit from errors during MPI runs' (Dimitris Menemenlis) http://mailman.mitgcm.org/pipermail/mitgcm-support/2007-October/thread.html ; 2003-October 'OpenMP Model Code.' ; GitHub PR #834 https://github.com/MITgcm/MITgcm/pull/834

### Array bound error only in multi-threaded configs: "bi loop upper bound = myByLo" in coupling-interface init
- Cause: loop `DO bi = myBxLo, myByLo` typo in atm_compon_interf/ocn_compon_interf init routines.
- Fix: PR #271 (use myBxHi), merged Aug 2019. General lesson: check bi/bj loops use myBxLo..myBxHi, myByLo..myByHi, and that overlap loops use full tile+OL range (PR #256 salt_plume_vol computed only in the overlap).
- Era: fixed upstream in PR #271 / #256 (2019).
- Src: GitHub PR #271 https://github.com/MITgcm/MITgcm/pull/271 ; PR #256

## Legacy / obsolete
### mds_byteswapi4.f: Operands of logical operator '.and.' are INTEGER(4)/BOZ (gfortran >=10, OpenMPI) with -DFAST_BYTESWAP
- Cause: FAST_BYTESWAP code was written for the retired Altix/Columbia; never used with -fconvert=big-endian.
- Fix: do not define FAST_BYTESWAP (leave undefined). Code removed upstream (doc/tag-index: "remove FAST_BYTESWAP code and option").
- Era: obsolete: removed after issue #658 (2022-2023).
- Src: GitHub issue #658 https://github.com/MITgcm/MITgcm/issues/658

### Why MPI_COMM_WORLD still appears in eedie.F
- Cause: deliberate: COMM_WORLD barrier in EEDIE keeps an ocean component from calling MPI_FINALIZE before a coupled atmosphere finishes its sub-steps; model code otherwise uses MPI_COMM_MODEL.
- Fix: in a custom coupler driver add a matching MPI_Barrier(MPI_COMM_WORLD) at the very end. exch_jam/gsum_jam are not used on normal hardware.
- Era: 2006.
- Src: mitgcm-support 2006-January 'MPI_COMM_WORLD question' http://mailman.mitgcm.org/pipermail/mitgcm-support/2006-January/thread.html
