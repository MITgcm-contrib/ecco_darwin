# Building MITgcm / darwin3

## Standard sequence

```
cd build
../tools/genmake2 -rootdir=<MITgcm or darwin3 root> -mods="<code dir(s)>" -of=<optfile> [-mpi] [-omp] [-ieee|-devel]
make depend
make -j 8
```

- `-mods` takes several dirs, searched left to right before the stock source, e.g.
  `-mods="../code_offline_ggl90 ../code_6+4+0_llc90_ggl90"`. Leaving out a patch dir
  silently builds the stock version of whatever it overrides.
- darwin3 builds need an explicit `-rootdir`, or genmake2 fails with "Cannot determine MITgcm root directory".
- `-ieee` / `-devel` turn off aggressive optimisation; use them for testreport comparisons,
  adjoint gradient checks and anything you'll compare bit-for-bit.
- Rebuild from clean (`make CLEAN`, rerun genmake2) after editing any header with
  common blocks, after adding a file to a `code/` dir, or after changing `packages.conf`.
  genmake2 symlinks sources into the build dir and never replaces an existing link —
  confirm overrides with `ls -la build/<file>.F`.
- `packages.conf` lists packages or groups (`gfd`, `oceanic`, …; see `pkg/pkg_groups`);
  dependencies come from `pkg/pkg_depend`. genmake2 writes `PACKAGES_CONFIG.h` with the
  `ALLOW_<PKG>` macros. A compiled-in package also needs `use<PKG>=.TRUE.` in `data.pkg`
  to run, and often fails if its `data.<pkg>` is missing.

## Local Mac optfiles (arm64 Mac)

| Optfile | Use |
|---|---|
| `ECCO/{wetting_drying,BBL,sea_ice_BCs}/build_options/darwin_x86_gfortran_local` | x86_64 Homebrew gfortran/gcc + open-mpi + netcdf in `/usr/local`, runs under Rosetta. Has `-fconvert=big-endian -fallow-argument-mismatch`. |
| `ECCO/wetting_drying/build_options/darwin_arm64_gfortran_conda` | Native arm64, micromamba env `~/.local/micromamba/envs/mitgcm`. Use `genmake2 ... -make=/usr/bin/make` (the env's GNU make skips cpp). 2–5× faster than Rosetta; results agree to ~1e-14. Preferred for new work. |
| `research/debug/darwin_local_gfortran` | Older local darwin3 build. |

Known Mac failures and fixes:
- `ld: library 'netcdf' not found` → set `NETCDF_INC=/usr/local/include NETCDF_LIB=/usr/local/lib`
  (or add `-L$(nc-config --libdir)`); `nf-config --flibs` omits the C lib dir.
- Undefined `_timenow_` → arm64 C objects mixed with x86 Fortran: `CFLAGS += -arch x86_64`.
- Link failure on huge static arrays (>~2 GB) → fewer tracers/levels per rank, keep `numDiags` ≈ 4·Nr, or use MPI.
- Apparent hang that is really a SIGSEGV in DIAGNOSTICS_WRITE when `numLevels` is large → stack size.
- Segfault at start → `ulimit -s 65520`; with OpenMP also `OMP_STACKSIZE=512M`.
- Check disk first: `df -h /System/Volumes/Data` (often near full; Time Machine local
  snapshots hold deleted files).

## Pleiades (NAS) optfiles

- `linux_amd64_ifort+mpi_ice_nas` — general; and a `…electra_skylake_…` variant with
  `-xCORE-AVX512` that dies with *illegal instruction* on non-Skylake nodes. Match the
  optfile to the `model=` you request.
- On the toss4 OS, older optfiles point at a removed MPT path; use an `optfile_toss4`.
- Build modules (same as run): `module purge; module load comp-intel/2020.4.304 mpi-hpe/mpt
  hdf4/4.2.12 hdf5/1.8.18_mpt netcdf/4.4.1.1_mpt python3/3.9.12`.
- darwin3's `cog` calls `hashlib.md5()`, which FIPS-mode OpenSSL on NAS blocks. Patch
  `tools/darwin/cogapp/cogapp.py` to `hashlib.md5(..., usedforsecurity=False)` or move
  `pkg/darwin/Makefile` aside so `make depend` doesn't run cog.
- AWS optfile example: `ECCO/offline/code_v4r6_forward/linux_ifort_impi_aws_sysmodule`.

## Adjoint / tangent linear

- **TAF** (licensed, approved for use): `genmake2 ... -ieee` then `make adall` →
  `mitgcmuv_ad`; `make ftlall` → `mitgcmuv_ftl`. Gradient check output: `grep grad-res`.
- **Tapenade 3.16** is unpacked in `ECCO/BBL/tools/tapenade_3.16` (also `/nobackup/<nas_user>/bbl_tap`
  on NAS; there `module load jvm/jdk11`), and the build is `genmake2 -tap` then `make tap_adj`, but
  **it does not work yet on either machine**: its Fortran parser falls back to Docker because the
  distribution ships no macOS parser and the Linux one needs glibc 2.34 (pfe has 2.28). Symptom:
  "Parsing error in <every file>" plus "docker: command not found". Fix: build the parser from
  source (Tapenade git clone, `./gradlew frontf`, ~0.5 GB of downloads). Packages need no
  Tapenade-specific code; Tapenade reads the same `<pkg>_ad_diff.list` as TAF.
- **ECCO v4r5–r7 adjoints are built on MITgcm checkpoint68g.** pkg/bbl for them lives on branch
  `bbl-c68g` (worktree `ECCO/BBL/MITgcm_c68g`); ECCO+BBL build in `ECCO/BBL/ecco_ad/`. A
  checkpoint68g TAF build with today's `staf` needs, in the generated `Makefile`: remove
  `-server fastopt.net` (fails with "Host key verification failed") and add `-fixed` to
  `AD_TAF_FLAGS`/`FTL_TAF_FLAGS`. On the Mac, clang `cpp -traditional` does not expand a macro
  called with a space before the bracket (`_GLOBAL_SUM_RL (` / `_EXCH_XY_RS (`): remove the space
  in the code dir. The ECCO LLC90 adjoint has ~3 GB of static data and cannot link on macOS
  (> 2 GB); link it on Pleiades with `-mcmodel=medium` and MPI.
- Lessons from pkg/bbl adjoint work (`ECCO/BBL/plans/BBL_adjoint_scope.md`,
  `ECCO/BBL/tests/adjoint/README.md`):
  - Every new routine must be in `<pkg>_ad_diff.list`. Routines left out are treated as passive and the
    gradient is **silently wrong**; uncompilable code was only a special case (bbl_calc_rho overwriting rhoInSitu).
  - Packages with state need `<pkg>_ad_check_lev{1..4}_dir.h` store directives.
  - Raise `nWh` in `MDSIO_BUFF_WH.h` when the tape outgrows it.
  - A package's init must run after `CTRL_INIT_VARIABLES`.
  - Use FD eps ~1e-4 (not 1e-2) near model switches/thresholds.
- ECCO v4 gradient check from `data.ctrl.iter0` (`doinitxx=.TRUE.`): also set `optimcycle=0` in
  `data.optim`. The ctrl pkg writes zero `xx_*` files only when `optimcycle=0`; with v4r6's 47 it
  stops at once with "MDS_READ_FIELD: xx_etan.0000000047.data: File does not exist".
  Zero every `mult_gencost`/`mult_profiles*`/`mult_genarr*` so a test cost isn't buried in
  ~2e6 of data misfit (FD then sits at round-off). grdchk uses `iGloPos/jGloPos/kGloPos` only
  when `nbeg=0` (then `nend` = extra points after it); with `nbeg>=1` it counts wet points from
  the tile's first surface cell and silently ignores the position. LLC90 v4r6: adjoint of 24 steps
  ≈ 1.7 h on 5 bro_ele nodes, each FD pair ≈ 7 min → ask for ≥ 4 h.
