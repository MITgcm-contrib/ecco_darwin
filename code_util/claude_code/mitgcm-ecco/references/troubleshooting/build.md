# Troubleshooting: Build / compile / link / start-up (genmake2, optfiles, MPI, NetCDF, TAF/Tapenade)
Distilled from mitgcm-support 2003-2024 + GitHub issues/PRs (2018-2026). Names verified against origin/master of the local MITgcm clone (Oct 2026); "Era" says when a fix is obsolete or upstream. Mailing-list threads cited by month; thread index URL = http://mailman.mitgcm.org/pipermail/mitgcm-support/YYYY-Month/thread.html.

## genmake2 / Makefile basics

### "These lines are here to deliberately cause a compile-time error" (chksum_tiled.f, "Non-numeric character in statement label")
- Cause: the stub model/inc/SIZE.h was picked up (still starts with those lines) because your code dir with SIZE.h was not in the include path: no `-mods`, wrong dir (building in `input/`), or stale links in build/.
- Fix: build in `build/`, run `../../../tools/genmake2 -mods ../code ...` (code/ must hold SIZE.h; list several dirs as `-mods "dir1 dir2"`), then `make CLEAN; make depend; make`. Check `ls -l SIZE.h` shows a link to ../code/SIZE.h.
- Era: 2003-2024, same message every year (Derecho thread 2024-Aug, 2024-Aug issue #853).
- Src: mitgcm-support 2024-August 'Ask for MITgcm optfile...'; 2014-September 'Compile errors ... ifort on Ubuntu 14.04'; 2007-May 'Porting MITgcm to Solaris x64'

### Edited/new header or source in code/ is ignored; "This name does not have a type" for a new/renamed variable (diagMdsDir, sideDragFactor, viscC4smag, inAdMode, MNC_TAG_ID); "OBCS.h: No such file"
- Cause: stale symlinks in build/ (`make clean` does NOT remove *.h links; `make CLEAN` does), or a local copy in code/ of a header/.F (PARAMS.h, ini_parms.F, the_main_loop.F, gad_calc_rhs.F, DIAGNOSTICS_SIZE.h ...) from an older checkpoint mixed with newer model code. Local copies are never updated by git/cvs.
- Fix: `make CLEAN` (or `make makefile && make CLEAN && make depend && make`) after adding/moving/renaming any file in code/ or changing checkpoint; diff every file in code/ against the current pkg/model version and merge (e.g. inAdMode moved to AUTODIFF_PARAMS.h; OBCS.h was retired). `make depend` re-creates links; also clear stray non-link *.F left in build/.
- Era: all years; 2021-July, 2024-Aug (diagnostics_ini_io: Jean-Michel pointed to issue #853).
- Src: mitgcm-support 2021-July 'Custom CPP_OPTIONS.h ignored'; 2017-March "Can't compile diagnostics_ini_io.F"; 2005-December 'bug with viscC4smag ?'; 2013-February 'make errors for contrib model'; 2015-July 'about compilation'; https://github.com/MITgcm/MITgcm/issues/853

### Fortran errors: "_d" / "WORDLENGTH" not defined, "Missing kind-parameter" `0.5 _d 0`, "Illegal use of symbol d0 - KIND", unexpanded `_RL`/`_EXCH_XYZ_R8`, "missing terminating ' character"
- Cause: the .F -> .f cpp step did not run properly. Makefile rule is `cat x.F | cpp -traditional -P $(DEFINES) $(INCLUDES) | tools/set64bitConst.sh > x.f` (set64bitConst AFTER cpp). Typical triggers: forgot `make depend`; CPP in optfile lacks `-traditional -P` (or is pgcpp/gcc that mangles); compiled .F directly; broken Makefile (CPPCMD order reversed); DEFINES missing -DWORDLENGTH=4; code/ copy of an old file still using a macro that was renamed (`_EXCH_XYZ_R8` -> `_EXCH_XYZ_RL`, Apr 2009).
- Fix: sequence genmake2 -> `make depend` -> `make`; optfile `CPP='cpp -traditional -P'` (macOS: add `-xassembler-with-cpp` for pragma lines); run the cpp command by hand from the make log to see what is left; fresh clone if Makefile looks hand-damaged; stop using `.f` from a manual compile.
- Era: 2003-2020.
- Src: mitgcm-support 2013-April 'Compiling MITgcm'; 2013-May 'strange FORTRAN error'; 2020-July 'Compiling issue'; 2010-July 'a quick compile question'; 2004-August 'no -mpi option for genmake2'; https://github.com/MITgcm/MITgcm/issues/344

### "Error: No Fortran compilers were found in your path" / "can't read OPTFILE=.../linux_aarch64_mpif77|darwin_arm64_mpif77|irix64_ip27_g77" / "no options file found that matches this platform"
- Cause: genmake2 auto-picks `${os}_${arch}_${compiler}`; no file for your combo, or compiler not on PATH (module not loaded; ifort not in cluster PATH).
- Fix: give `-of path/optfile` (or `MITGCM_OF` env), `-fc NAME`, or `FC=` env; check `which gfortran mpif90`; for a new machine pick the closest file by `uname -a` + compiler (e.g. linux_amd64_gfortran) and copy/rename it (`tools/suggest_optfile_names`). Current genmake2 maps x86_64->amd64 and aarch64->arm64 (so linux_arm64_gfortran / darwin_arm64_gfortran are found).
- Era: always; aarch64 mapping added mid-2024 (issue #847).
- Src: mitgcm-support 2008-March 'fortran compiler'; 2020-September 'MITgcm Support Request'; 2021-June 'problem compiling MITgcm with openmpi on mac big sur'; https://github.com/MITgcm/MITgcm/issues/847

### "gfortran: error: unrecognized argument in option '-mcmodel=medium'" on aarch64 Linux (Docker/VM on Apple M-series)
- Cause: genmake2 picked linux_amd64_gfortran for an arm64 CPU; gcc on aarch64 does not support `medium`.
- Fix: use linux_arm64_gfortran (genmake2 now maps aarch64->arm64 automatically); temporary hack `-mcmodel=small` in CFLAGS/FFLAGS.
- Era: 2024 (issue #847, genmake2 commit 105ccd9). Same arch bug seen with MITgcm.jl on M3.
- Src: https://github.com/MITgcm/MITgcm/issues/847

### "Warning: FORTRAN compiler test failed" / "Makefile might be unusable" though everything is fine
- Cause: genmake2 test run used `mpirun`; system only has `mpiexec`.
- Fix: fixed upstream (issue #182 / PR #190, Jan 2019); update tools/genmake2.
- Era: 2018-12 .. 2019-01.
- Src: https://github.com/MITgcm/MITgcm/issues/182

## MPI / NetCDF library detection and linking

### `makedepend: warning ... cannot find include file "mpif.h"` / `EESUPPORT.h: error: mpif.h: No such file` / `MPI_SUCCESS has no type` / `MPI_REAL has no type`
- Cause: cpp (not mpif77) does the include; MPI include dir not passed. `-mpi` takes no argument (`-mpi=PATH` was removed from use).
- Fix: before genmake2: `export MPI_INC_DIR=/path/containing/mpif.h` (or `MPI_HOME`; modern optfiles also try `pkg-config`), run genmake2 with `-mpi`, `make depend`. On modules systems load the MPI module (e.g. `module load impi` on Yellowstone/Derecho). Also set FC=mpif77/mpif90 wrapper, not the bare compiler.
- Era: 2005-2025.
- Src: mitgcm-support 2013-September 'mpi=PATH option for genmake2'; 2014-March 'building on yellowstone'; 2017-March 'help for building the model'; 2023-May 'How to make the executable with Intel compiler?'; https://github.com/MITgcm/MITgcm/issues/635

### "MPI_HOME is not set and pkg-config not available, aborting" (darwin_amd64/arm64_gfortran, linux_amd64_gfortran with -mpi)
- Cause: optfile looks for MPI_INC_DIR, then MPI_HOME, then `pkg-config ompi|mpich`; Homebrew/HPC often provide none.
- Fix: `export MPI_HOME=/path/to/mpi` (or MPI_INC_DIR) before genmake2.
- Era: since PR #158 (Jan 2019); still in optfiles.
- Src: https://github.com/MITgcm/MITgcm/issues/195 ; mitgcm-support 2020-August 'Unable to compile'

### Link step: `undefined reference to mpi_wait_ / mpi_recv_ / mpi_init_` / `_gfortran_st_write` / `multiple definition of main (for_main.o)` / `cannot find -lnetcdf`
- Cause: compiler/MPI/NetCDF mismatch (MPI built with gfortran but linking with ifort; Intel `-lmpi` extras; wrong LIBS), MPI libs not found.
- Fix: first prove a hello-world MPI Fortran program builds with the same wrapper; use matching compiler for MPI and NetCDF (Intel: linux_amd64_ifort+impi); don't add -lmpi by hand when FC is the wrapper; start without NetCDF, then add it.
- Era: all years.
- Src: mitgcm-support 2015-February 'link errors during compiling'; 2021-February "Compiling problem by ifort: undefined reference to '_gfortran_st_write'"; 2013-December 'problems with compile with MPI'

### genmake2: "Can we create NetCDF-enabled binaries... no" / WARNING: the "mnc" package was enabled but tests failed ... now DISABLED
- Cause: check_netcdf_libs() could not compile/link the little test (genmake_tnc.F, `#include "netcdf.inc"`, nf_create). Usual: missing `-lnetcdff -lnetcdf` in LIBS (NetCDF4 splits C and Fortran libs), include/lib dir not set, NetCDF built with a different compiler (name mangling), hdf5 libs missing.
- Fix: read `genmake.log`, section `running: check_netcdf_libs()`; re-run its cpp/gfortran commands by hand and fix flags until link works. Point to the install with `NETCDF_ROOT` (needs $NETCDF_ROOT/include + /lib), or NETCDF_HOME, NETCDF_INC/NETCDF_LIB, else `nf-config --includedir --flibs`; or hard-wire `INCLUDES`/`LIBS` in a copy of the optfile. Model still builds/runs without NetCDF (just no mnc); `SKIP_NETCDF_CHECK=t` in env skips the test.
- Era: always. NetCDF is NOT required to run (mds output).
- Src: mitgcm-support 2023-April 'Loading Netcdf libraries for MITGCM' (Martin's genmake_tnc.F recipe); 2015-April 'Recommended NetCDF libraries'; 2015-June 'netcdf v3.x'; 2022-May 'Errors in the installation of MITgcm'; 2012-July 'mnc package was enabled but tests failed'

### macOS/Homebrew arm64: `ld: library 'netcdf' not found` or NetCDF check silently disabled (nf-config gives `-lnetcdff -lnetcdf -lnetcdf`, two Cellar prefixes)
- Cause: netcdf-c and netcdf-fortran are separate Homebrew prefixes; darwin_arm64_gfortran uses NETCDF_ROOT/include+lib or `nf-config` and then needs the C libs in the same -L; with NETCDF_ROOT pointing at one prefix HAVE_NETCDF ends up empty. Not an optfile bug per Martin/Oliver (conflicting NetCDF installs; check `nf-config --flibs` lib dir really contains libnetcdf).
- Fix: pick one install tree that has BOTH include/netcdf.inc and lib/libnetcdf+libnetcdff (or add the C `-L` to LIBS in a local optfile copy); verify with the manual genmake_tnc steps above. No generic fix was committed.
- Era: 2024 (#846 closed); #982 (open, 2026) proposes skipping all nf-config logic in optfiles when SKIP_NETCDF_CHECK is set.
- Src: https://github.com/MITgcm/MITgcm/issues/846 ; https://github.com/MITgcm/MITgcm/issues/982

### Run time: `error while loading shared libraries: libnetcdff.so.5 | libifport.so.5 | libmpi_f77.so.1` / `symbol lookup error ... undefined symbol: __intel_sse2_strlen` / `Trace/BPT trap` after NetCDF upgrade
- Cause: executable links dynamic libs; compute node/env lacks them on LD_LIBRARY_PATH (DYLD_LIBRARY_PATH on macOS), or `-shared-intel` without the Intel lib dir in path; stale older NetCDF found first.
- Fix: put the lib dirs on LD_LIBRARY_PATH in the run script (all nodes); source the compiler's `*vars.sh`; or link static (`-static-intel`, `-static-libgfortran`); keep run-time env identical to build env (same `which mpif77`/`mpirun`).
- Era: all years.
- Src: mitgcm-support 2012-July 'libnetcdff.so.5 not found!'; 2011-March 'lib64/libopen-pal.so.0: undefined symbol: __intel_sse2_strlen'; 2013-February 'Executable fails, message: Trace/BPT trap: 5'; 2011-November 'Problem with parallel build'

## Large executables

### `relocation truncated to fit: R_X86_64_PC32 / R_X86_64_32S against symbol ..._ defined in COMMON section` ; `failed to convert GOTPCREL relocation; relink with --no-relax` ; macOS `ld: 32-bit RIP relative reference out of range`
- Cause: static arrays (SIZE.h tile * nSx*nSy, plus package sizes) exceed 2 GB with the default x86_64 small code model. Appears when you enlarge sNx/sNy or add packages (e.g. diagnostics) with the same tile.
- Fix: (1) preferred: more/smaller tiles (increase nPx/nPy, decrease sNx/sNy; keep nSx,nSy small) - check with `size mitgcmuv`; (2) add `-mcmodel=medium` to BOTH FFLAGS and CFLAGS (ifort: `-mcmodel=medium -shared-intel`; remove -fPIC which cancels it; link libs built -fPIC/mcmodel too); linux_amd64_gfortran already sets `-mcmodel=medium`; ice_nas has it commented-out for >2 GB. Odd C-file symbols (sigreg.o, tim.o) mean CFLAGS lacks it; a run on more cores often removes the need. Also reduce `numDiags/numlists` in code/DIAGNOSTICS_SIZE.h if diagnostics bloats it.
- Era: 2005-2020 (still hit by large single-tile LLC/high-res tests).
- Src: mitgcm-support 2020-October 'compiling problems'; 2019-July 'Link problem on OSX'; 2013-December 'domain size limit?'; 2008-August 'Compiling error for high resolution configuration'; 2015-December 'Compiler error when doubling horizontal resolution'; 2010-April 'convincing compiler of node memory size'

### Compiles fine, runs out of memory at load ("There is not enough memory", "Out of memory") with diagnostics
- Cause: pkg/diagnostics static arrays from code/DIAGNOSTICS_SIZE.h (numlists, numperlist, numLevels, numDiags) are compile-time; `useDiagnostics=.FALSE.` at run time does not shrink the executable.
- Fix: copy pkg/diagnostics/DIAGNOSTICS_SIZE.h to code/ and reduce `numlists`, `numperlist`, `numLevels`, `numDiags` (e.g. numDiags=1*Nr) to what you need, or drop `diagnostics` from packages.conf, then rebuild.
- Era: all.
- Src: mitgcm-support 2014-September 'Out of memory: optfile for IBM AIX'; 2008-July 'ifort' (Matt: reduce numdiags in DIAGNOSTICS_SIZE.h)

## Compiler-specific

### gfortran >= 10: `Error: Type mismatch between actual argument at (1) and actual argument at (2)` (MPI_RECV, MPI_Allreduce in gather/global_max/cumulsum)
- Cause: gfortran 10 made argument-type mismatch an error (MPI Fortran-77 interface).
- Fix: `FFLAGS="$FFLAGS -fallow-argument-mismatch"` (present in linux_amd64_gfortran, darwin_arm64_gfortran); `#define DISABLE_MPI_READY_TO_RECEIVE` does NOT fix it.
- Era: gcc/gfortran 10+ (2020+); PRs #420, #480; issue #354.
- Src: https://github.com/MITgcm/MITgcm/issues/354 ; https://github.com/MITgcm/MITgcm/pull/480 ; mitgcm-support 2021-June 'problem compiling MITgcm with openmpi on mac big surr'

### clang/gcc >= 10 on macOS: `implicitly declaring library function 'sprintf'` in setdir.c; `clang: error: invalid version number in '-mmacosx-version-min=11.2'`
- Cause: setdir.c missed `#include <stdio.h>`; Big Sur SDK-version string with Homebrew gcc.
- Fix: setdir.c fixed upstream (now includes stdio.h); for the SDK message see the Stack Overflow "Solution no. 1" linked in the issue (not reproduced in the digest). Also add `-w -fallow-argument-mismatch` on old trees (Dustin's 2021 fix).
- Era: 2020-2021; fixed upstream.
- Src: https://github.com/MITgcm/MITgcm/issues/389 ; https://github.com/MITgcm/MITgcm/issues/479 ; mitgcm-support 2021-June (macbook M1/openmpi)

### macOS (case-insensitive FS): "Your file system cannot distinguish between *.F and *.f" ; gfortran "linker input file unused because linking not done" for .fr9 free-form files
- Cause: genmake2 check_for_broken_Ff() switches suffixes to `.for` / `.fr9`; gfortran does not recognise .fr9 as Fortran so no .o/.mod is made; older f90mkdepend/xmakedepend assumed one-letter suffix so module dependencies were lost.
- Fix: current darwin_arm64_gfortran sets `F90FLAGS="$FFLAGS -x f95 -ffree-form"`; use current genmake2/xmakedepend (PR #865, #879). On very old trees put `FS='for'; FS90='fr9'` in the optfile.
- Era: old FS=/FS90= trick 2009; .fr9/-x f95 fix 2023-2024 (issues #754, #850, PR #879, #865).
- Src: https://github.com/MITgcm/MITgcm/issues/754 ; https://github.com/MITgcm/MITgcm/issues/850 ; https://github.com/MITgcm/MITgcm/pull/865 ; mitgcm-support 2009-July 'MITgcm on intel mac 10.5.7'

### Wrong results / NaN / segfault only with optimization or only with one compiler (ifort 8-13, PGI 13.4, ifort vs g77)
- Cause: compiler optimization bugs in specific routines (examples: ifort vectorisation of loops in find_rho (FIND_RHOP0/FIND_BULKMOD/FIND_RHODEN), PGI 13.4.0 `__fsd_exp` floating exception in swfrac.F, segfault in mom_calc_visc.F, obcs_init_fixed.F segfault on Pleiades).
- Fix: rebuild with `-O0` to confirm; bisect with NOOPTFILES/NOOPTFLAGS in the optfile (linux_amd64_ifort+mpi_ice_nas sets NOOPTFILES incl. obcs_init_fixed.F); try other compiler version; ifort `-xN`/`-x` vectorisation was the culprit for find_rho; genmake2 `-devel` for run-time checks.
- Era: 2005-2019; Pleiades obcs_init_fixed.F entry (PR #215/#225, 2019) is already in ice_nas optfile.
- Src: mitgcm-support 2005-March 'ifort optimization puzzle'; 2013-August 'Mysterious initialisation problem'; 2018-June 'segmentation fault'; https://github.com/MITgcm/MITgcm/pull/225 ; https://github.com/MITgcm/MITgcm/pull/215

### Model blows up within a few steps with ifort (or at certain tile sizes) but not gfortran; results differ with tile size
- Cause: overlap too small for the advection scheme (e.g. PTRACERS_advScheme=33 / 7-point schemes need OLx,OLy >= 3, cube-sphere needs more); odd sNx with some ifort/-O versions; optimisation.
- Fix: raise OLx/OLy in SIZE.h (Samar: 3 -> 4 cured), use even sNx (JMC reproduced a blow-up with sNx=47, fine with 48 on ifort 13), lower optimisation. For tile-size dependence: identical results require same tile size (`GLOBAL_SUM_ORDER_TILES` is default define since Aug 2015, so same sNx,sNy with different nPx/nSx matches); small differences otherwise are expected.
- Era: 2005-2016.
- Src: mitgcm-support 2008-September 'problem running tutorial_global_oce_latlon'; 2008-January 'ifort/g77 issue'; 2016-April 'results quite differents depending on number of procs used'

## Namelists, run directory, start-up

### Fortran runtime: "namelist not terminated with / or &end", "Cannot match namelist object name X", "End of file" in ini_parms / packages_boot / ptracers_readparms / mnc_readparms
- Cause: data.* files use `&` as terminator; some compilers/MPI builds want `/` (or `&end`). Also integer set as `0.` (niter0), continuation lines in namelist arrays (PGI), unnamed/typo parameter (e.g. user package flag not added to the PACKAGES namelist), missing terminator.
- Fix: add `-DNML_TERMINATOR` to DEFINES in the optfile (current gfortran optfiles do; code rewrites & -> / on the fly via nml_set_terminator.F; pass a leading blank `-DNML_TERMINATOR=" /"` only if needed), or use `&end`/`/` consistently; strip data.* to bisect; for a new package flag define it in PARAMS.h/packages_boot.F namelist and keep data.pkg consistent.
- Era: all years (gfortran 2008+, PGI, xlf, Mac).
- Src: mitgcm-support 2019-January 'Issue with lagrangian float (FLT) package'; 2014-April 'MITgcm, gfortran and Archer'; 2014-February 'Fortran problem with MITgcm'; 2015-October 'mismatch in namelist'; 2008-June 'Compilation with pgi fortran'; 2008-September 'Problems with tutorial_global_oce_biogeo'

### scratch1.0000#### / scratch2.0000#### files left in run dir; "Permission denied trying to open file /tmp/gfortrantmp..."; "namelist read ... end of file reached without finding group" (many MPI ranks)
- Cause: since checkpoint66j (2017) namelists are read via scratch files named per process (comment lines stripped scratch1->scratch2); left over after a crash in a namelist (look inside to find the bad line); older code: ranks overwrote shared scratch; /tmp not writable for compiler-created temp files.
- Fix: debug namelist by stripping data.* and re-adding lines; `#define USE_FORTRAN_SCRATCH_FILES` in CPP_EEOPTIONS.h for the old STATUS='SCRATCH' behaviour; `#define SINGLE_DISK_IO` (only proc 0 writes STDOUT/STDERR/scratch; hides errors from other ranks); `export TMPDIR=<writable dir>`; TARGET_BGL / TARGET_CRAYXT in DEFINES change scratch handling on old BG/Cray.
- Era: 2009-2023 (USE_FORTRAN_SCRATCH_FILES and SINGLE_DISK_IO verified in current CPP_EEOPTIONS.h).
- Src: mitgcm-support 2018-September 'scratch1.00000#### in run directory'; 2023-March 'MITgcm-support Digest, Vol 237, Issue 1'; 2014-July 'eeset_parms.f and permissions error'; 2009-July 'namelist; end of file reached without finding group'

### `INI_PROCS: needs MPI for multi-procs (nPx*nPy=N) ... usingMPI=False` / `EEBOOT_MINIMAL: No. of procs= 1 not equal to nPx*nPy` / `S/R LOAD_GRID_SPACING ABNORMAL END` / identical output on every tile
- Cause: executable built with a SIZE.h for N procs but run serial, built without `-mpi`, or `mpirun -np` != nPx*nPy; changing nPx/nPy without keeping Nx=sNx*nSx*nPx constant breaks input-file sizes.
- Fix: rebuild with `genmake2 -mpi` (+ matching SIZE.h from the same code dir), run `mpirun -np $((nPx*nPy))` (not nSx*nSy: those are tiles per process), keep global Nx,Ny unchanged. One executable per SIZE.h (static allocation; use separate code dirs).
- Era: all years.
- Src: mitgcm-support 2020-July 'Errors with eedata and eeset_parms'; 2018-August 'An error about procs when I submitted a MITgcm job'; 2011-June 'Question about parallel run'; 2013-April 'Problem for the MITgcm parallelized compilation'

### `forrtl: severe (36): attempt to access non-existent record, unit 9, file bathy.bin|topog.*|OB*.bin` / "Non-existing record number" / "do_ud: end of file" (MDS_READ_FIELD)
- Cause: input file smaller than the model expects: wrong global size vs Nx*Ny*Nr (SIZE.h changed), wrong precision (default readBinaryPrec=32, tutorials with 64; exf_iprec), OBCS files must span the FULL boundary length (Nx or Ny, e.g. 43 not the 13 open points) x Nr x time records; forcing records fewer than externForcingPeriod/externForcingCycle imply; EXF with a different-size field but no nlon/nlat for it; ifort 8+ record length units (needs `-assume byterecl` or -DWORDLENGTH=1).
- Fix: size check `ls -l`; set readBinaryPrec / exf_iprec to match how the file was written; EXF per-field `*_nlon, *_nlat, *_lon0, *_lat0, *_lon_inc, *_lat_inc` in data.exf; keep `externForcingCycle` <= records x period; make sure the file is in the run dir (symlink with `ln -s ../input/* .`); ifort `-assume byterecl` is already in current optfiles.
- Era: all years; also seen as "OBzonalV.bin" with forcing cycle too long.
- Src: mitgcm-support 2013-October 'problems with verification/tutorial_barotropic_gyre'; 2019-August 'Reading errors'; 2016-June '(no subject)' (exf nlon/nlat); 2005-January 'Intel Fortran 8.1: -assume byterecl'; 2004-November 'reading error'

### Binary input/output looks like garbage / NaN / "inf"; python-written files do not match
- Cause: MITgcm mds I/O is always big-endian (unless NetCDF); python writes native little-endian or wrong dtype/dimension order; ifort/pgi flags missing.
- Fix: write with `arr.astype('>f4').tofile(f)` (never mix with byteswap(); read back with dtype='>f4'; `meshgrid(..., indexing='ij')` order is (y,x)); build with `-convert big_endian` (ifort), `-fconvert=big-endian` (gfortran), `-byteswapio` (pgi) OR `-D_BYTESWAPIO`, not both (double swap); match readBinaryPrec.
- Era: all years.
- Src: mitgcm-support 2021-April 'Barotropic ocean gyre tutorial'; 2010-March 'Error in reading the input files'; 2008-October 'Error reading dx file with pathscale'; 2004-November 'More problems'

### STDOUT/STDERR empty or not flushed after crash or queue-kill; "output.txt empty" under MPI
- Cause: output buffering; under MPI stdout goes to STDOUT.0000..., errors to STDERR.0000...
- Fix: `debugMode=.TRUE.` in eedata (extra prints + flush each write, needs HAVE_FLUSH='t' in genmake.log); `debugLevel=4` in data; compile flags `-g -traceback` / gfortran `-fbacktrace`; look at the end of ALL STDOUT.* and STDERR.*; PGI line buffering via SETVBUF3F (old).
- Era: all years (debugMode verified in eeset_parms.F).
- Src: mitgcm-support 2013-August 'Mysterious initialisation problem - debugging options?'; 2014-August 'compiler flags to make output human readable'; 2008-June 'output.txt empty'; 2016-January 'line-buffer the STDOUT'

### Run script needs to detect failed runs; no non-zero exit code on crash
- Cause: MITgcm prints "STOP ABNORMAL END" but historically did not return an error code/MPI_Abort (see #439).
- Fix: grep STDOUT.0000 for `ABNORMAL END`; test last line `Execution ended Normally`; use `writePickupAtEnd=.TRUE.` (or pChkptFreq hitting final step) and check for the pickup.
- Era: 2019-2023.
- Src: https://github.com/MITgcm/MITgcm/issues/197

### Segmentation fault at/near start (before first step, in a package init, with huge arrays)
- Cause: stack limit too small (static arrays on stack, OpenMP/MPI threads), not a model bug. Seen on AIX, Cray, fizhi-ifort, adjoint (mitgcmuv_ad) runs, `CLOSE error unit 11` cascades.
- Fix: `ulimit -s unlimited` (bash) / `limit stacksize unlimited` (csh) in the job script (and on compute nodes); OMP: set `OMP_STACKSIZE`/`GOMP_STACKSIZE`/`KMP_STACKSIZE`; on macOS add `-Wl,-stack_size,...`.
- Era: all years (2005-2023).
- Src: mitgcm-support 2007-August 'Chasing a seg fault'; 2023-March 'MITgcm-support Digest, Vol 237, Issue 4'; 2010-December 'Fizhi-ifort'; https://github.com/MITgcm/MITgcm/issues/427

### "MON_SOLUTION: STOPPED DUE TO EXTREME VALUES OF SOLUTION" / "SOLUTION IS HEADING OUT OF BOUNDS"; "EEDIE: Only 0 threads have completed"
- Cause: model instability (CFL, deltaT too large vs viscAh, diffusion, vertical grid spacing, bad forcing); the EEDIE message is a side effect of the stop, not a threading problem.
- Fix: read monitor output (advcfl_*) first; need CFL<~0.5 and deltaT < dx^2/viscAh; reduce deltaT; then look at compiler/overlap issues above.
- Era: all.
- Src: mitgcm-support 2006-February 'simple question on internal wave simulation'; 2005-October 'Anybody know what's up with this?' (Martin on vertical CFL)

## Performance

### MPI time dominated by GLOBAL_SUM_TILE (MPI_Allreduce of a mostly-zero vector) at many tiles (e.g. llc540, 2819 ranks)
- Cause: GLOBAL_SUM_ORDER_TILES (default #define in eesupp/inc/CPP_EEOPTIONS.h) makes every rank Allreduce a nTile-long vector to keep sums independent of tile count; slows down with very many tiles.
- Fix: `#undef GLOBAL_SUM_ORDER_TILES` (and keep `#undef GLOBAL_SUM_SEND_RECV`) in CPP_EEOPTIONS.h; D. Kokron measured Allreduce time 342 s -> 99 s for a 5-day run; a SPEED_AT_ALL_COST option was agreed but is not in master. Sums then change at round-off level with tile count.
- Era: checkpoint67q, Nov 2020 (issue #385).
- Src: https://github.com/MITgcm/MITgcm/issues/385

### Run much slower than expected: huge Wall vs User time, monitor/file I/O dominates
- Cause: `chkptFreq` tiny (rolling pickup every step), frequent `monitorFreq`, useSingleCpuIO/globalFiles slow on some machines (SGI), diagnostics/I-O on slow disk; `make -j` only speeds compilation.
- Fix: chkptFreq ~ 1/10 of pChkptFreq; monitorFreq 20-50*deltaT or 0; use `useSingleCpuIO=.TRUE.` not `globalFiles` (globalFiles also crashed pickup writing on a Cray: `lib-5058`); look at the timing table at end of STDOUT.0000.
- Era: 2005-2017.
- Src: mitgcm-support 2007-March 'NaNQ' (chkptFreq advice); 2015-July 'Poor cpu usage percentage'; 2017-April 'error while writing pickup files with Cray compilers'; 2005-June 'bluesky build'

## Packages / CPP-option combinations

### "*** ERROR *** from PACKAGES_CHECK: run-time control flag useX is set but pkg/X was not compiled" ; "CONFIG_CHECK: #undef NONLIN_FRSURF and nonlinFreeSurf is non-zero" ; "DIAGSTATS_SET_REGIONS: #define DIAGSTATS_REGION_MASK missing in DIAG_OPTIONS.h"
- Cause: runtime flag (data.pkg/data/data.diagnostics) requires a compile-time package or CPP option that is off; packages.conf/CPP_OPTIONS.h/DIAG_OPTIONS.h in code/ missing it.
- Fix: add the package to code/packages.conf (or `-enable=pkg`); define the option in the local *_OPTIONS.h (e.g. `#define DIAGSTATS_REGION_MASK` in DIAG_OPTIONS.h plus `sizRegMsk` in DIAGNOSTICS_SIZE.h, else "COMMON block data object must not be an automatic object"); `#define NONLIN_FRSURF` for rStar/nonlin FS; then `make CLEAN`. If you do not use the package, remove the run-time flag.
- Era: all years (DIAGSTATS_REGION_MASK, sizRegMsk verified).
- Src: mitgcm-support 2011-January 'oasis ????'; 2021-July 'Custom CPP_OPTIONS.h ignored'; 2020-July 'Errors with eedata and eeset_parms'; 2010-September 'compile problems with regional statistics'

### Compile error when `ALLOW_CTRL_OBCS{N,S,E,W}` defined but `ALLOW_OBCS_{NORTH,SOUTH,EAST,WEST}` undefined (OBNt/OBNu undeclared in ctrl_getobcsn.F); obcs_cost_driver.F missing useObcsCostContribution
- Cause: default CTRL_OPTIONS.h/OBCS_OPTIONS.h allowed the inconsistent combination.
- Fix: keep the pairs consistent; upstream now compiles and STOPs cleanly in CTRL_CHECK ("... but CPP-flag ALLOW_OBCS_NORTH is not defined") (PR #889); obcs_cost_driver fixed in PR #984 (Apr 2026).
- Era: fixed Nov 2024 (#889) / Apr 2026 (#984); older ECCO-style OBCS-control builds hit it.
- Src: https://github.com/MITgcm/MITgcm/pull/889 ; https://github.com/MITgcm/MITgcm/issues/888 ; https://github.com/MITgcm/MITgcm/pull/984

### Type-mismatch compile errors with `genmake2 -devel` / `-ur4` (use _RS=real*4) when writing or updating packages
- Cause: strict S/R-argument checks (ifort -check/-devel) reveal wrong _RS/_RL use: DIAGNOSTICS_FILL_RS called with an _RL array, MDS_READVEC_LOC passed a scalar instead of a length-1 vector, missing argument (ADEXCH_3D_RL needs Nr), `_RL` declared where `_RS` expected (obcs_check summary), locals initialised only `IF (useDiagnostics)` but used always.
- Fix: match types: `DIAGNOSTICS_FILL` (RL) vs `DIAGNOSTICS_FILL_RS`; pass `tmpVal(1)` arrays to MDS_READVEC_LOC; use the right print-summary routine for RS; guard with `#ifdef ALLOW_DIAGNOSTICS` + `useDiagnostics`; test with `-ur4` (`-use_real4`) and `-devel` on a quick testreport. Never use "ALLOW_AUTODIFF_TAMC" for logic that is not TAF-specific (use ALLOW_AUTODIFF).
- Era: 2021-2026 (PRs #428, #838, #708, #1036, #731).
- Src: https://github.com/MITgcm/MITgcm/pull/428 ; https://github.com/MITgcm/MITgcm/pull/838 ; https://github.com/MITgcm/MITgcm/pull/708 ; https://github.com/MITgcm/MITgcm/pull/1036 ; https://github.com/MITgcm/MITgcm/issues/730

### Compiling pkg/seaice without pkg/exf: "seaice_growth.f ... snowPrecip / EVAP / PRECIP has no type"; STOP in seaice_check "need to define pkg/exf ALLOW_ATM_TEMP"
- Cause: since May 2013 pkg_depend no longer enables exf with seaice; seaice needs bulk-formula (ALLOW_ATM_TEMP) fluxes, not pre-computed fluxes alone.
- Fix: add `exf` (and `cal`) to packages.conf; define ALLOW_ATM_TEMP in EXF_OPTIONS.h; start from lab_sea / global_ocean.cs32x15 SEAICE_OPTIONS.h; use SEAICE_CGRID (B-grid code barely tested). Pure-seaice builds without exf were repaired for offline_cheapaml (PR #247).
- Era: 2013-2019.
- Src: mitgcm-support 2013-August 'Error Messages for seaice_growth.f during compile'; https://github.com/MITgcm/MITgcm/pull/247

### F90 module code: "Can't open module file X.mod" / compiled in wrong order with `make -j` / free-format source rejected by a fixed-format compiler
- Cause: genmake2/f90mkdepend only orders modules for files named `*_mod.F` (fixed) or `*_mod.F90` (free) and treats `.F90` as free-format, `.F` as fixed-format; modules in other files are not seen by the dependency scan.
- Fix: rename module files `*_mod.F[90]` (as pkg/ptracers, pkg/atm_phys); for free-format F90 sources use suffix .F90 and, depending on compiler, `ALWAYS_USE_F90=1` in genmake_local; keep header include dependencies via current xmakedepend (F90 includes detected since PR #865); with TAF use -ncad / topological ordering of modules.
- Era: 2021-2025.
- Src: mitgcm-support 2023-July 'configuring compilation dependencies in Makefile'; 2021-August 'using a module instead of a common block?'; https://github.com/MITgcm/MITgcm/issues/475 ; https://github.com/MITgcm/MITgcm/pull/865

### KPP_ESTIMATE_UREF: dbloc / work1 not declared in kpp_forcing_surf.F; LOG of negative zref in vermix
- Cause: option never compiled since 2007; fix passes dbloc from kpp_calc.F and uses abs(rF).
- Fix: current kpp_forcing_surf.F takes `dbloc` as an argument under KPP_ESTIMATE_UREF; physical correctness of zref (sign of rF) was left unresolved in the thread.
- Era: issue #336 (2020); option still `#undef` by default.
- Src: https://github.com/MITgcm/MITgcm/issues/336

### NetCDF-less build: pkg/obsfit / pkg/profiles / mnc fail without NetCDF
- Cause: pkg needs NetCDF but gets enabled by packages.conf.
- Fix: pkg/obsfit is disabled automatically if no NetCDF (PR #957, 2025, like mnc/profiles). Keep `useMNC`/mnc consistent with HAVE_NETCDF.
- Era: fixed Dec 2025.
- Src: https://github.com/MITgcm/MITgcm/pull/957

### mnc output wrong grid (X,Y = 1..N) with pkg/exch2 and a lat-lon grid
- Cause: pkg/mnc assumes exch2 implies curvilinear grid (usingCurvilinearGrid); JMC's diagnosis, not fixed in thread.
- Fix: recompile without pkg/exch2 (use non-mpi SIZE.h, not SIZE.h_mpi with blank tile) or use mds output; glue tiles with utils/python/MITgcmutils/scripts/gluemncbig (needs all tiles; useSingleCpuIO/globalFiles do not apply to mnc).
- Era: 2024-10 (global_ocean.90x40x15).
- Src: mitgcm-support 2024-October 'Issue with MITgcm Global Ocean Example'; 2014-June 'Problems'; 2016-April 'Read NetCDF data (Internal Waves)'

### How to use blank tiles (exch2) - blankList / dimsFacets
- Cause: users unsure how to set data.exch2 (cs and lat-lon).
- Fix: pkg/exch2 must be compiled (no run-time switch); list blank tiles in data.exch2 `blankList`; set `dimsFacets` to the FULL domain incl. blank tiles; find empty tile numbers by running 1 step without blank tiles and grepping `Empty tile:` in STDOUT; examples verification/adjustment.cs-32x32x1/input/data.exch2.mpi and global_ocean.90x40x15/input/data.exch2.mpi; mpirun -np = (non-blank tiles)/nSx*nSy; input files stay full size. PR #382 also relevant.
- Era: 2022.
- Src: mitgcm-support 2022-November '"blanklist" parameter in MITgcm'

## Adjoint / TAF / Tapenade / OpenAD builds

### TAF adjoint link: `undefined reference to mdthe_main_loop_` (or `adthe_main_loop_`) in the_model_main.o
- Cause: DIVA (divided adjoint) code path (`mdthe_main_loop`) only exists when TAF is run with `-pure`; lab_sea (and anything with ALLOW_DIVIDED_ADJOINT / genmake_local pointing at the DIVA adoptfile) needs it. Building in a dir without that genmake_local gives the default AD_TAF_FLAGS (no -pure). Also: building an adjoint-only experiment (obcs_ctrl, global1x1) as a forward model gives `adthe_main_loop`.
- Fix: use the build dir's genmake_local (testreport does) or set `USE_DIVA=1`/DIVA flags: current adjoint_default adds `-pure` only if USE_DIVA=1 (adjoint_diva no longer exists); undef ALLOW_DIVIDED_ADJOINT in code_ad/AUTODIFF_OPTIONS.h if you don't want DIVA; don't compile adjoint-only experiments with plain `make`; check TAF flags in the Makefile `AD_TAF_FLAGS`.
- Era: 2008-2019 reports; adjoint_diva removed, DIVA via USE_DIVA now (verified).
- Src: https://github.com/MITgcm/MITgcm/issues/232 ; https://github.com/MITgcm/MITgcm/issues/804 ; mitgcm-support 2012-August 'OBCS_ctrl package'; 2008-May 'help'; 2015-November 'Issue compiling offline model with TAF on ARCHER'

### TAF: `undefined reference to adexch_3d_rl_ / adexch_uv_3d_rl_ / adexch_uv_xy_rs_ / adexch_xy_rs_` (addummy_in_stepping.o, monitor_ad.o, copy_ad_uv_outp.o)
- Cause: TAF only generates adjoints of EXCH routines for fields that are active; in simple control/cost setups some are never used, but the hand-written AD-variable output (ALLOW_AUTODIFF_MONITOR) calls them. `AUTODIFF_EXCLUDE_ADEXCH_RS` only covers the 2D RS versions, not adexch_3d_rl.
- Fix: `#undef ALLOW_AUTODIFF_MONITOR` in AUTODIFF_OPTIONS.h (or ECCO_CPPOPTIONS.h) - no AD-variable output; or comment the calls in addummy_in_stepping.F/monitor_ad.F/copy_ad_uv_outp.F; make sure exch of a basic state variable is called in do_stagger_fields_exchanges.F so TAF generates it. Martin offers to inspect code_ad if it fails in a verification experiment.
- Era: 2015-2023.
- Src: mitgcm-support 2015-November 'Issue compiling offline model with TAF on ARCHER'; 2023-March 'undefined reference to `adexch_3d_rl_'

### TAF errors: `CADJ STORE ... *ERROR* identifier not defined` (StoreDynVars3D, tapelev4); `keyword REC, KIND, or SHAPE expected` with -f08
- Cause: 4-level checkpointing enabled but AUTODIFF.h not included in (old) the_main_loop.F (version mix); TAF < 6.8.10 mis-parses `kind = isbyte` spacing in STORE directives under F08.
- Fix: start from a working code_ad (lab_sea, natl_box_adjoint) of the same checkpoint and re-copy the_main_loop.F; use TAF >= 6.8.10 (or write `kind=isbyte` without spaces as workaround) when `TAF_FORTRAN_VERS='F08'`/`USE_EXTENDED_SRC`.
- Era: 2008 (c59); 2025 (#940, fixed in TAF 6.8.10).
- Src: mitgcm-support 2008-February 'error in MITgcm compilation'; https://github.com/MITgcm/MITgcm/issues/940

### Parallel make of adjoint: race with `make -j N adtaf` (ad_config.template removed by another make instance)
- Cause: all TAF stages used the same temp file name ad_config.template.
- Fix: fixed in PR #915 (ad_config.template0/1/2, Apr 2025); update tools/genmake2; or run the adjoint build without -j.
- Era: fixed April 2025.
- Src: https://github.com/MITgcm/MITgcm/pull/915

### Tapenade: `undefined reference to the_main_loop_b_`; "docker: command not found" fallback; hangs/ fails at fortranParser
- Cause: Tapenade's precompiled fortranParser (C) incompatible with the system glibc; Tapenade silently falls back to a docker environment which is absent on Linux.
- Fix: compile fortranParser from the tarball provided by the Tapenade developers with the `compile` script (gcc only) and copy to `<tapenade>/bin/linux/fortranParser` (no Tapenade rebuild); test with `tapenade program.f` on a hello-world; for diagnosing add `-tracelevel 15 -traceparser` to the Tapenade command inside tools/genmake2 (not a genmake2 flag). One user's build then stalled at `Created ./autodiff_init_varia_b.msg` (unresolved). Mac: Tapenade needs the docker recipe + darwin_arm64_gfortran (PR #935).
- Era: 2024-2025; see umbrella issue #735.
- Src: mitgcm-support 2024-September 'Compiling Problems' (Shreyas/Laurent); mitgcm-support 2024-October 'Compiling Problems with Tapenade'; https://github.com/MITgcm/MITgcm/pull/935 ; https://github.com/MITgcm/MITgcm/issues/735

### Tapenade makefile: "flow_tap" file not found; TOOLSDIR undefined in adjoint_tap
- Cause: bash variable `${TOOLSDIR}` expands to empty when the ad-optfile is sourced; the Makefile macro is defined later.
- Fix: PR #973: in tools/adjoint_options/adjoint_tap use the escaped Makefile macro `\$(TOOLSDIR)/TAP_support/flow_tap` (verified in master); same for diffsizes.F90 path in genmake2.
- Era: Feb 2026.
- Src: https://github.com/MITgcm/MITgcm/pull/973

### OpenAD (old): `make adAll` redoes everything / dontTransform, dontCompile files ignored
- Cause: *.xsd links break dependency (issue #618, open 2022); dontTransform/dontCompile are only read by genmake2 for AD tools, not documented, and don't override `*_ad_diff.list` entries.
- Fix: read the genmake2 logic (dontTransform ~line 2883); remove a file from `<pkg>_ad_diff.list` to stop it being transformed; keep manual relink workflow. OpenAD is effectively superseded by Tapenade.
- Era: 2014, 2022.
- Src: mitgcm-support 2014-October 'openad/genmake2 configuration files'; https://github.com/MITgcm/MITgcm/issues/618

## Misc build-tool notes

### `makedepend: error: out of space: increase MAXFILES` (adjoint/ECCO builds with thousands of files)
- Cause: system makedepend has a hard file limit.
- Fix: build the shipped one: `cd tools/cyrus-imapd-makedepend; ./configure; make` (do not edit def.h), then `genmake2 -makedepend ../../../tools/cyrus-imapd-makedepend/makedepend ...` (or `-md`).
- Era: 2007-2008; makedepend options still in genmake2.
- Src: mitgcm-support 2007-February 'makedepend: error'

### Harmless genmake2 / compiler messages
- `Do we have etime()/LAPACK ... no`, `Can we register a signal handler ... no`, `cloc() ... no`, `f90mkdepend: no source file found for module this`, ifort remark `LOOP WAS VECTORIZED`, `-O vs -O3` option-override warnings, `Unused dummy argument` warnings, `echo: No match`: ignore if an executable is produced; the only effect is missing timing/signal features. If every probe says "no", the compiler/cpp is broken: look at genmake.log. `undefined _timenow/_system_time/_user_time` at link (old g95/xlf/pgf name mangling): use `genmake2 -ignoretime` (option still in genmake2) or set FC_NAMEMANGLE.
- Era: 2005-2015.
- Src: mitgcm-support 2015-November 'Problems on gfortran in compiling'; 2012-March '(no subject)'; 2010-July 'LOOP WAS VECTORIZED'; 2006-September 'error in exp. s12t_16x32'; 2010-July 'a quick compile question'

### Windows: genmake2 "No Fortran compilers", "C pre-processor failed the test case", sigreg.c ucontext.h missing
- Cause: MITgcm is not supported natively on Windows (Git Bash/MinGW has no matching optfile); Cygwin lacked ucontext.h (fixed by excluding sigreg.c).
- Fix: use WSL, a Linux VM (VirtualBox notes: MITgcm-contrib/ecco_darwin/MITgcm_VirtualBox.txt; coessing-mitgcm-2023 docs Ubuntu_on_Windows.txt), or Docker (MITgcm.jl / ECCO-Docker).
- Era: 2004-2024.
- Src: mitgcm-support 2021-May 'Fortran compilers for Windows'; 2020-September 'MITgcm Support Request'; 2018-August 'Issue setting up GCM'; 2024-August 'Ask for MITgcm optfile and generating make'

### "make: ... multiple definition" / duplicate objects after cvs/git update; or `ln: ./CVS: Operation not permitted` in testreport
- Cause: leftover auto-generated files from templates (exch_*.rl.F) or stale dirs; external filesystem not supporting symlinks (building on an external drive/NTFS).
- Fix: `make CLEAN`, build on a local Unix filesystem (symlinks required), fresh checkout if the tree is mixed.
- Era: 2003-2013.
- Src: mitgcm-support 2004-August 'no -mpi option for genmake2'; 2006-June 'teething.. genmake2 with gfortran'; 2013-November 'Make Depend: does not symbolic link the needed files'
