# Developing MITgcm and Darwin packages

The best templates are the user's own packages, which already follow MITgcm conventions:
`pkg/wad` in `ECCO/MITgcm_wad_checkin` (cleanest: check-in branch), `pkg/bbl` in
`ECCO/BBL/MITgcm` (with adjoint support), and the OBCS sea-ice sponge in
`ECCO/sea_ice_BCs/MITgcm`. Upstream also ships **`pkg/mypackage`**, the official
skeleton for a new package (readparms, check, init, diagnostics, pickups, tendencies; it's
exercised by `verification/hs94.1x64x5`). Start a brand-new package by copying it and
renaming. For a feature similar to an existing one, copy that package's pattern
(`mitgcm-index/pkg/<name>.md` lists its call sites into the core) and read `doc/` alongside it.
Note: upstream `pkg/bbl` already exists (a simple BBL); the user's branch extends it. Read the core routine you are
hooking into before editing it. Most bugs come from misunderstanding where in the step a
field is valid, not from Fortran.

## Anatomy of a package `pkg/<name>/`

| File | Role | Called from |
|---|---|---|
| `<NAME>_OPTIONS.h` | package CPP flags, wrapped in `#ifdef ALLOW_<NAME>`; includes `PACKAGES_CONFIG.h` + `CPP_OPTIONS.h` | every `.F` in the package (first line) |
| `<NAME>.h` | common blocks: parameters (`<NAME>_PARM_*`) and fields; `_RL`/`_RS`/`LOGICAL` blocks kept separate | |
| `<name>_readparms.F` | reads `data.<name>` (`NAMELIST /<NAME>_PARM01/`) via `OPEN_COPY_DATA_FILE`, inside `_BEGIN_MASTER`/`_END_MASTER`, then `_BARRIER` | `model/src/packages_readparms.F` |
| `<name>_check.F` | parameter consistency + incompatible-option stops | `model/src/packages_check.F` |
| `<name>_init_fixed.F` | time-invariant setup; calls `<name>_diagnostics_init` under `ALLOW_DIAGNOSTICS` | `model/src/packages_init_fixed.F` |
| `<name>_init_varia.F` | initial state; calls `<name>_read_pickup(nIter0)` when restarting | `model/src/packages_init_variables.F` |
| `<name>_diagnostics_init.F` | `DIAGNOSTICS_ADDTOLIST` for each field | `<name>_init_fixed` |
| `<name>_read_pickup.F` / `<name>_write_pickup.F` | own `pickup_<name>.*` file | init_varia / `model/src/packages_write_pickup.F` |
| feature routines | the physics, called from the core with guards | e.g. `forward_step.F`, `external_forcing_surf.F`, `do_oceanic_phys.F` |
| `<name>_ad_diff.list`, `<name>_ad_check_lev*_dir.h` | adjoint: routines TAF differentiates; tape store directives | autodiff |

**Registering the package in the core** (this is what `git grep -n useWAD` shows):
1. `model/inc/PARAMS.h`: `LOGICAL use<NAME>` added to the package-flag declaration and common block.
2. `model/src/packages_boot.F`: add to the `PACKAGES` namelist, default `.FALSE.`, and
   `CALL PACKAGES_PRINT_MSG( use<NAME>, '<NAME>', ' ' )` under `#ifdef ALLOW_<NAME>`.
3. `packages_readparms.F`, `packages_check.F` (incl. the `PACKAGES_ERROR_MSG` path when
   `use<NAME>` is set but the package isn't compiled), `packages_init_fixed.F`,
   `packages_init_variables.F`, `packages_write_pickup.F`.
4. Every call site: `#ifdef ALLOW_<NAME>` / `IF ( use<NAME> ) CALL ...` / `#endif`.
   Wrap in `TIMER_START/STOP('<NAME>_X   [CALLER]', myThid)` for timing output.
5. If needed: compile-time dependencies in `pkg/pkg_depend` (`+req`/`-excl`) and groups in `pkg/pkg_groups`.
   Run-time incompatibilities go in `packages_check.F` (e.g. `useDOWN_SLOPE .AND. useBBL` stops the run there).
6. Document: `doc/phys_pkgs/<name>.rst` (+ add to `phys_pkgs.rst`) and a `doc/tag-index` entry.

## Fortran conventions (match them exactly)

- Fixed-form F77 with `.F` (cpp'd): code in columns 7–72, `C` comments in col 1,
  continuation `     &`. Header blocks `CBOP` / `C !ROUTINE:` / `C !INTERFACE:` /
  `C !DESCRIPTION:` / `CEOP` (protex). Argument comments `C     name :: meaning`.
- Types: `_RL` for real state/params, `_RS` for grid/mask fields; `INTEGER myThid` is the
  last argument almost everywhere; `myTime, myIter, myThid` for time-dependent routines.
- Tile loops:
  ```
        DO bj=myByLo(myThid),myByHi(myThid)
         DO bi=myBxLo(myThid),myBxHi(myThid)
          DO k=1,Nr
           DO j=1-OLy,sNy+OLy
            DO i=1-OLx,sNx+OLx
  ```
  Loop over the interior (`1,sNx`) when you will exchange afterwards; loop over overlaps
  only if every input is valid there. After changing a field that neighbours read, call
  `_EXCH_XY_RL(fld, myThid)` / `_EXCH_XYZ_RL` (or the `EXCH_UV_*` vector forms on
  cube-sphere/LLC where u/v rotate across faces).
- Masks: `maskC/W/S`, `hFacC/W/S`, `kLowC`, `kSurfC`. Under z* (`select_rStar`) and with
  pkg/wad, the effective thickness is `hFac*rStarFac`; don't assume static hFac.
- Messages: `WRITE(msgBuf,'(A)') '...'` then
  `CALL PRINT_MESSAGE(msgBuf, standardMessageUnit, SQUEEZE_RIGHT, myThid)`; errors via
  `CALL PRINT_ERROR(msgBuf, myThid)` then `STOP 'ABNORMAL END: S/R <NAME>_CHECK'`.
- I/O: `READ_FLD_XY_RL`, `READ_REC_3D_RL`, `WRITE_FLD_XY_RL`; namelist I/O only on master.
- Global sums: `GLOBAL_SUM_TILE_RL` (tiling-reproducible) instead of hand-rolled MPI.
- Keep new behaviour off by default (a parameter or CPP flag) so stock verification
  results don't change; that's the bar for upstreaming.

## Diagnostics

```
        diagName  = 'WADdryC '        ! exactly 8 chars, padded
        diagTitle = 'WAD dry-column flag (1=at film floor)'
        diagUnits = '1               ' ! 16 chars
        diagCode  = 'SM      L1      ' ! 16-char code: col1 S=scalar/U/V/W, col2 M=cell centre,
                                       ! col9 vertical location, col10 '1'=single level / 'R'=Nr levels
                                       ! (full legend: header of pkg/diagnostics/diagnostics_main_init.F)
        CALL DIAGNOSTICS_ADDTOLIST( diagNum, diagName, diagCode, diagUnits, diagTitle, myThid )
```
Fill during the step with `CALL DIAGNOSTICS_FILL(fld, 'WADdryC ', kLev, nLevs, bibjFlg, bi, bj, myThid)`,
guarded by `IF ( useDiagnostics )`. Check `available_diagnostics.log` after a short run.

## Darwin (darwin3) specifics

- Darwin is called through **pkg/gchem** (`gchem_readparms`, `gchem_init_fixed`,
  `gchem_init_vari`, `gchem_tr_register`, `gchem_forcing_sep` → `DARWIN_FORCING`,
  `gchem_output` → `DARWIN_DIAGS`, `gchem_write_pickup`). darwin3 defines
  `GCHEM_SEPARATE_FORCING`, so tendencies are applied after advection, not inside it.
- Sizes in `DARWIN_SIZE.h`; tracer index layout in `DARWIN_INDICES.h`; parameters in
  `DARWIN_PARAMS.h` (read from `data.darwin`), trait parameters in `DARWIN_TRAITPARAMS.h`
  and per-type traits in `DARWIN_TRAITS.h` (`data.traits`). Allometric traits are generated
  in `darwin_generate_allometric.F`; random traits in `darwin_generate_random.F`.
- `pkg/darwin/darwin_check.F` contains **cog-generated** blocks (`CCOG[[[cog ... ]]]`)
  built from `tools/darwin/checkindices.py` by `pkg/darwin/Makefile` during `make depend`.
  Edit the generator or the cog template, not the generated lines, and rerun cog
  (`tools/darwin/cog -c -r darwin_check.F`). On NAS, cog needs the FIPS md5 patch.
- Tools in `tools/darwin/`: `mkdarwintracers` (data.ptracers block), `mkdiagnosticsdata`,
  `conscheck`/`conscheckall` (conservation check of a run), `checkindices.py`.
- Light: `pkg/radtrans` (spectral, `nlam` bands, OASIM forcing) feeds
  `darwin_light_radtrans.F`; `darwin_light.F` is the non-spectral path.
- Carbonate chemistry: `darwin_carbon_chem.F`; `darwin_solvesaphe.F` when SOLVESAPHE is on (v06).
- Branches: `backport_ckpt68y` (offline LLC90 V4r6), `backport_ckpt68g`; local clones at
  `~/Documents/GitHub/darwin3` (branch `darwin`, the default), `research/debug/darwin3` (`backport_ckpt68y`).
  The OIF kerguelen_eco config pins darwin3 `3f0529872` (Jan 2021); its old cogapp needs `python`/`imp` shims on the Mac.
  Legacy darwin2 macros (e.g. `DARWIN_ALLOW_DIAZ`) are inert in darwin3 group-based setups.
- Adding a tracer or process: update `DARWIN_SIZE.h`/`DARWIN_INDICES.h` logic, register in
  `darwin_tr_register.F`, add parameters + readparms + check, add diagnostics, conservation
  terms in `darwin_cons.F`, then regenerate `data.ptracers` with `mkdarwintracers` and
  run `conscheck`.

## Adjoint-readiness

- Don't use `IF` branches on state that switch discontinuously unless you accept a
  non-differentiable point; provide smooth alternatives behind a flag (cf. `WAD_SMOOTH_MASK`).
- List every routine in `<name>_ad_diff.list`: routines left out are treated as passive and
  the gradient is silently wrong (only some cases fail to compile); add store directives for state that's
  overwritten within a step (`<name>_ad_check_lev{1..4}_dir.h`, included from
  `pkg/autodiff/checkpoint_lev*_directives.h`).
- Recompute rather than store where cheap; check the tape size (`nWh` in `MDSIO_BUFF_WH.h`).
- Test: `make adall`, gradient check vs finite differences (`grdchk`), eps ~1e-4.

## Testing a new feature (minimum before calling it done)

1. Builds with the feature off and on, `-devel` (bounds checking) and optimised.
2. Stock verification experiments that don't use the package are bit-identical (testreport).
3. A new verification experiment `verification/<exp>/{code,input,results/output.txt}`
   that exercises it (see verification-upstream.md).
4. Restart test (N+N = 2N steps), MPI vs serial, and (if threaded) OpenMP bit-identity.
5. Conservation/budget closure (e.g. `BBL_CHECK_BUDGET`-style flag) and a movie of the test.
