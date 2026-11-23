# Provenance of every staged file (2026-11-16, updated 2026-11-22)

Source root on Pleiades: `/nobackup/dcarrol2/downscaling/nares/` (called `nares/`
below). Builds and runs are PBS jobs under `nares/regional_run/`. "Wins" means
genmake2 `-mo` resolution: directories are searched left to right and the
FIRST match is compiled.

Physics frozen 2026-11-20 at `run_arch_D10` (built as `build_drag10`).
`code/` reproduces `build_drag10` file for file, with ONE exception: `SIZE.h`
is the 20x20 production file (sNx = sNy = 20, nPx = 678; from
`nares/regional_run/code_size20x20/`, as built in `build_bgc_s20x20` for the
scaling test). `input/data.exch2` is its partner (from
`nares/prep/scaling/data.exch2_20x20`: 1008 tiles, 330 blank). Every 2015
tuning run used the 40x40 pair (192 ranks).

## 1. `code/` = the physics build `build_drag10`, collapsed

Build script: `scripts/job_build_drag10.sh` (built 2026-11-16, PBS 25394931,
`DRAG10_BUILD_OK`). Its `-mo` line, verified in the script:

```
-mo "$P/code_drag10 $P/code_swfix $P/code_capload $P/code_darwin $P/code"
```

So the collapse rule for `code/` is: `code_drag10` > `code_swfix` >
`code_capload` > `code_darwin` (the 2026-10-07 set) > `code`. The same `code/`
minus the 6 `code_drag10` files = `build_swfix` (used by `run_arch_Y69m_sw`,
`run_arch_Y69c_sw`, `run_arch_Y69m_swD`, `run_arch_Y16sw`). Minus also the 2
`code_swfix` files = `build_capload` (used by `run_arch_Y16` ... `run_arch_Y69m`,
and, with `code_bgc`, by `build_bgc`).

What each source directory is:

| source dir | what it is | used by |
|---|---|---|
| `code` | pan_AO `code/` + Nares changes: `SIZE.h` (40x40, nPx 192), `SEAICE_OPTIONS.h` (JFNK, TEM, side + bottom drag, cap), `CPP_OPTIONS.h` (`ALLOW_ADDFLUID`), `seaice_model.F` (pan_AO 2017 override), 11 core files patched by the iceplume installer (2026-10-07; pkg/iceplume = https://github.com/mhwood/iceplume at e4d4ac9 + local fixes) | every Nares build |
| `code_darwin` | the 2026-10-07 Darwin set (pan_AO code_darwin of that date + 5 radtrans under-ice light overrides + iceplume blocks in `apply_forcing.F`). 17 of its 34 files are stale for BGC (17-argument `darwin_plankton.F`). Compiled into every physics binary but inactive (`useGCHEM = .FALSE.`) | `build`, `build_capload`, `build_swfix`, `build_drag10` |
| `code_capload` | 2026-10-24: `SEAICE_CAP_ICELOAD` + `heffTooHeavy = 6 m` (ice load capped at 6 m of water) | `build_capload` and every later build |
| `code_swfix` | 2026-11-08: `IF (useGCHEM)` guard around the CDOMdecay shortwave path in `apply_forcing.F` and `seaice_growth.F` (physics-only runs had no ocean SW heating) | `build_swfix`, `build_drag10` |
| `code_drag10` | 2026-11-15 (`prep/make_drag10.py`): ice-ocean drag on the 0-10 m mean ocean velocity (`SEAICE_OCNVEL.h`, `seaice_ocnvel_avg.F`, 4 patched seaice routines) | `build_drag10` only (run_arch_D10t passed a 10-day test; run_arch_D10 = the frozen physics) |

File table (61 files; the last two rows are the SEAICE_SWFrac backport, taken from the darwin3 `nares` branch commit ce5f626cf, not from an override directory). "lost" = the same-named files further right in the `-mo` list.

| file | source dir (wins) | lost |
|---|---|---|
| `CPP_EEOPTIONS.h` | `code` | - |
| `CPP_OPTIONS.h` | `code` | - |
| `DARWIN_EXF_FIELDS.h` | `code_darwin` | - |
| `DARWIN_EXF_PARAMS.h` | `code_darwin` | - |
| `DARWIN_FIELDS.h` | `code_darwin` | - |
| `DARWIN_INDICES.h` | `code_darwin` | - |
| `DARWIN_OPTIONS.h` | `code_darwin` | - |
| `DARWIN_PARAMS.h` | `code_darwin` | - |
| `DARWIN_SIZE.h` | `code_darwin` | - |
| `DIAGNOSTICS_SIZE.h` | `code_darwin` | - |
| `EXF_OPTIONS.h` | `code` | - |
| `KPP_OPTIONS.h` | `code` | - |
| `MOM_COMMON_OPTIONS.h` | `code` | - |
| `OBCS_OPTIONS.h` | `code` | - |
| `PARAMS.h` | `code` | - |
| `PTRACERS_SIZE.h` | `code_darwin` | - |
| `RADTRANS_PARAMS.h` | `code_darwin` | - |
| `RADTRANS_SIZE.h` | `code_darwin` | - |
| `SEAICE_OCNVEL.h` | `code_drag10` | - |
| `SEAICE_OPTIONS.h` | `code_capload` | code |
| `SIZE.h` | `code_size20x20` (staged; `build_drag10` compiled the 40x40 `code/SIZE.h`) | - |
| `apply_forcing.F` | `code_swfix` | code_darwin |
| `darwin_add_surfforc.F` | `code_darwin` | - |
| `darwin_check.F` | `code_darwin` | - |
| `darwin_cons.F` | `code_darwin` | - |
| `darwin_diagnostics_init.F` | `code_darwin` | - |
| `darwin_exf_init_fixed.F` | `code_darwin` | - |
| `darwin_exf_init_varia.F` | `code_darwin` | - |
| `darwin_exf_load.F` | `code_darwin` | - |
| `darwin_exf_readparms.F` | `code_darwin` | - |
| `darwin_forcing.F` | `code_darwin` | - |
| `darwin_light.F` | `code_darwin` | - |
| `darwin_monitor.F` | `code_darwin` | - |
| `darwin_plankton.F` | `code_darwin` | - |
| `darwin_read_params.F` | `code_darwin` | - |
| `do_oceanic_phys.F` | `code` | - |
| `exf_getffields.F` | `code` | - |
| `external_fields_load.F` | `code` | - |
| `external_forcing.F` | `code_darwin` | - |
| `gchem_forcing_sep.F` | `code` | - |
| `obcs_balance_flow.F` | `code` | - |
| `packages.conf` | `code_darwin` | - |
| `packages_boot.F` | `code` | - |
| `packages_check.F` | `code` | - |
| `packages_init_fixed.F` | `code` | - |
| `packages_init_variables.F` | `code` | - |
| `packages_readparms.F` | `code` | - |
| `ptracers_apply_forcing.F` | `code` | - |
| `radtrans_calc.F` | `code_darwin` | - |
| `radtrans_check.F` | `code_darwin` | - |
| `radtrans_diagnostics_init.F` | `code_darwin` | - |
| `radtrans_readparms.F` | `code_darwin` | - |
| `seaice_dynsolver.F` | `code_drag10` | - |
| `seaice_growth.F` | `code_swfix` | code_capload code_darwin |
| `seaice_lsr.F` | `code_drag10` | - |
| `seaice_model.F` | `code` | - |
| `seaice_ocean_stress.F` | `code_drag10` | - |
| `seaice_oceandrag_coeffs.F` | `code_drag10` | - |
| `seaice_ocnvel_avg.F` | `code_drag10` | - |
| `SEAICE.h` | darwin3 `nares` branch, `pkg/seaice/` (ce5f626cf) | darwin_ckpt69k `pkg/seaice/SEAICE.h` |
| `seaice_init_fixed.F` | darwin3 `nares` branch, `pkg/seaice/` (ce5f626cf) | darwin_ckpt69k `pkg/seaice/seaice_init_fixed.F` |

## 2. `code_darwin/` = the BGC build `build_bgc`, collapsed, minus what `code/` already has

Build script: `scripts/job_build_bgc.sh` (built 2026-10-29, PBS 25320487). Its
`-mo` line, verified in the script:

```
-mo "$P/code_capload $P/code_bgc $P/code"
```

So the BGC collapse rule is `code_capload` > `code_bgc` > `code`. `code_bgc`
(34 files, its own `PROVENANCE.txt` on pfe) = pan_AO `code_darwin_v5` (= the
west_AO `code_baseline` Darwin files, identical md5) + 2 Nares deltas:
`apply_forcing.F` (v5 + iceplume blocks) and `packages.conf` (+ shelfice,
iceplume). Rule from BGC_PLAN.md: Nares uses the SAME ecosystem code and
namelists as pan_AO run_v5nut77 and west_AO run_prod4_2015H1.

Staging rule: a BGC-chain winner is staged in `code_darwin/` only if its
content differs from the file of the same name in `code/`. Two exceptions
(`apply_forcing.F`, `seaice_growth.F`) are explained in the table. Result:
18 files. The intended BGC build is therefore `-mo "code_darwin code"`, which
= `build_bgc` + `code_swfix` + `code_drag10`. Built once as a test on
2026-11-16 (`regional_run/build_v05test`, PBS 25399785): clean build, and the
compiled sources differ from `build_bgc` only in the swfix and drag10 files
and from `build_drag10` only in the 18 `code_darwin` files. **Not run** (see
OPEN_QUESTIONS.md item 4).

| file | source dir (wins in the BGC chain) | lost | staged? |
|---|---|---|---|
| `CPP_EEOPTIONS.h` | `code` | - | same as code/ (not staged) |
| `CPP_OPTIONS.h` | `code` | - | same as code/ (not staged) |
| `DARWIN_EXF_FIELDS.h` | `code_bgc` | - | staged in code_darwin/ |
| `DARWIN_EXF_PARAMS.h` | `code_bgc` | - | staged in code_darwin/ |
| `DARWIN_FIELDS.h` | `code_bgc` | - | staged in code_darwin/ |
| `DARWIN_INDICES.h` | `code_bgc` | - | same as code/ (not staged) |
| `DARWIN_OPTIONS.h` | `code_bgc` | - | staged in code_darwin/ |
| `DARWIN_PARAMS.h` | `code_bgc` | - | staged in code_darwin/ |
| `DARWIN_SIZE.h` | `code_bgc` | - | same as code/ (not staged) |
| `DIAGNOSTICS_SIZE.h` | `code_bgc` | - | same as code/ (not staged) |
| `EXF_OPTIONS.h` | `code` | - | same as code/ (not staged) |
| `KPP_OPTIONS.h` | `code` | - | same as code/ (not staged) |
| `MOM_COMMON_OPTIONS.h` | `code` | - | same as code/ (not staged) |
| `OBCS_OPTIONS.h` | `code` | - | same as code/ (not staged) |
| `PARAMS.h` | `code` | - | same as code/ (not staged) |
| `PTRACERS_SIZE.h` | `code_bgc` | - | same as code/ (not staged) |
| `RADTRANS_PARAMS.h` | `code_bgc` | - | same as code/ (not staged) |
| `RADTRANS_SIZE.h` | `code_bgc` | - | same as code/ (not staged) |
| `SEAICE_OPTIONS.h` | `code_capload` | code | same as code/ (not staged) |
| `SIZE.h` | `code` | - | same as code/ (not staged) |
| `apply_forcing.F` | `code_bgc` | - | NOT staged: code/apply_forcing.F (code_swfix) = this file + useGCHEM guard |
| `darwin_add_surfforc.F` | `code_bgc` | - | same as code/ (not staged) |
| `darwin_benthic_emulator.F` | `code_bgc` | - | staged in code_darwin/ |
| `darwin_check.F` | `code_bgc` | - | staged in code_darwin/ |
| `darwin_cons.F` | `code_bgc` | - | staged in code_darwin/ |
| `darwin_cons_reset.F` | `code_bgc` | - | staged in code_darwin/ |
| `darwin_diagnostics_init.F` | `code_bgc` | - | staged in code_darwin/ |
| `darwin_exf_init_fixed.F` | `code_bgc` | - | staged in code_darwin/ |
| `darwin_exf_init_varia.F` | `code_bgc` | - | staged in code_darwin/ |
| `darwin_exf_load.F` | `code_bgc` | - | staged in code_darwin/ |
| `darwin_exf_readparms.F` | `code_bgc` | - | staged in code_darwin/ |
| `darwin_forcing.F` | `code_bgc` | - | staged in code_darwin/ |
| `darwin_light.F` | `code_bgc` | - | same as code/ (not staged) |
| `darwin_monitor.F` | `code_bgc` | - | staged in code_darwin/ |
| `darwin_plankton.F` | `code_bgc` | - | staged in code_darwin/ |
| `darwin_read_params.F` | `code_bgc` | - | staged in code_darwin/ |
| `do_oceanic_phys.F` | `code` | - | same as code/ (not staged) |
| `exf_getffields.F` | `code` | - | same as code/ (not staged) |
| `external_fields_load.F` | `code` | - | same as code/ (not staged) |
| `external_forcing.F` | `code_bgc` | - | same as code/ (not staged) |
| `gchem_forcing_sep.F` | `code` | - | same as code/ (not staged) |
| `obcs_balance_flow.F` | `code` | - | same as code/ (not staged) |
| `packages.conf` | `code_bgc` | - | same as code/ (not staged) |
| `packages_boot.F` | `code` | - | same as code/ (not staged) |
| `packages_check.F` | `code` | - | same as code/ (not staged) |
| `packages_init_fixed.F` | `code` | - | same as code/ (not staged) |
| `packages_init_variables.F` | `code` | - | same as code/ (not staged) |
| `packages_readparms.F` | `code` | - | same as code/ (not staged) |
| `ptracers_apply_forcing.F` | `code` | - | same as code/ (not staged) |
| `radtrans_calc.F` | `code_bgc` | - | same as code/ (not staged) |
| `radtrans_check.F` | `code_bgc` | - | same as code/ (not staged) |
| `radtrans_diagnostics_init.F` | `code_bgc` | - | same as code/ (not staged) |
| `radtrans_readparms.F` | `code_bgc` | - | same as code/ (not staged) |
| `seaice_growth.F` | `code_capload` | code_bgc | NOT staged: code/seaice_growth.F (code_swfix) = this file + useGCHEM guard |
| `seaice_model.F` | `code` | - | same as code/ (not staged) |

The 13 `code_bgc` files that are "same as code/" are byte-identical to the
2026-10-07 `code_darwin` files that `code/` carries (the size headers,
`DARWIN_INDICES.h`, `darwin_add_surfforc.F`, `darwin_light.F`,
`external_forcing.F`, the 4 radtrans files, `packages.conf`).

## 3. `input/` = `nares/regional_run/run_arch_Y69m_sw/` namelists

`run_arch_Y69m_sw` (PBS 25386502, 2026-11-13 to 11-14, exit 0, 27.4 h on 5
sky_ele / 192 ranks) is the last scored full Jan 1 to Sep 1 2015 physics run and
the base of `run_arch_D10`. Copied: `data`, `data.cal`, `data.diagnostics`,
`data.exch2`, `data.exf`, `data.ggl90`, `data.iceplume`, `data.obcs`,
`data.pkg`, `data.salt_plume`, `data.seaice`, `data.shelfice`, `eedata`.
`data.exch2` was then replaced by the 20x20 production file (section 1).

Compared with the frozen run `run_arch_D10` on 2026-11-22: every staged
physics namelist is byte-identical to D10's, except `data.exch2` (20x20, see
above) and `data.diagnostics`. D10 adds one evaluation-only list 12 (daily
`uv10`: `UVEL`, `VVEL`, levels 1-7). It is not staged; it does not change the
solution. Not copied: `data.diagnostics.bak_with_hourly` (the 3 hourly lists dropped on
2026-11-13 to save inodes) and `data.obcs_tides.snippet` (a fragment already
inside `data.obcs`).

`run_arch_D10` differs from this set only by one added daily list in
`data.diagnostics` (list 12: `UVEL`, `VVEL`, levels 1-7, `diags/daily/uv10`).

### Differences from `nares/regional_run/input/` (the 2026-10-07 draft set)

`regional_run/input/` is the stale first draft (DECISIONS 2026-10-15: "input/
is NOT this baseline"). The staged files differ from it as follows
(comment-stripped diff, 2026-11-16):

| file | `regional_run/input/` (2026-10-07) | staged (`run_arch_Y69m_sw`) |
|---|---|---|
| `data` | `tRef`/`sRef` placeholder profiles; `deltaT = 120`; `nTimeSteps = 720` | `tRef`/`sRef` from the 2015-01-01 IC (`nares_tref_sref.py`); `deltaT = 60`; `nTimeSteps = 349920`; `horizGridFile = 'Nares_grid18'` |
| `data.pkg` | ptracers, gchem, radtrans, oasim, iceplume ON | all five OFF (physics-only run) |
| `data.obcs` | no sponge, `useOBCStides = .FALSE.` | `useOBCSsponge` (10 cells, 6 h to 10 d), `useOBCStides` with 8 constituents and 16 tidal files, `OB_Ieast` j = 12..30 zeroed |
| `data.seaice` | `presPow0/1 = 1/3`, `area_max 0.95`, `advScheme 30`, `diffKh 400`, `nonLinIterMax 20`, `tensilFac 0.05`, `basalDragK2 15`, `sideDrag 7.2e-3`, default P* | `HO 2.0`, `presPow0/1 = 2/1`, `area_max 0.98`, `SEAICE_strength 6.875e4`, `SEAICE_eccen 1.6`, `advScheme 77`, `diffKh 0`, tensile/basal/side drag 0, `HeffFile = pickup_seaice_siHEFF_ICT` (thinned mouth) |
| `data.ggl90` | TKE file only | + parent GGL90 set: `alpha 30`, `TKEmin 1e-7`, `TKEbottom 1e-6`, `mxlMaxFlag 2`, `mxlSurfFlag` |
| `data.diagnostics` | BGC daily surface lists | sea-ice 3-hourly `SIdyn`, daily `SIdeform`, `surf`, `surfTS`, `SHI`, `SImelt`, 6-hourly `SIforce`, monthly `state` |
| `data.exch2` | 40x40 blankList (60 of 252) | staged: 20x20 blankList (330 of 1008), the production tiling (run_arch_Y69m_sw and D10 ran the 40x40 file) |
| `data.cal`, `data.exf`, `data.iceplume`, `data.salt_plume`, `data.shelfice`, `eedata` | identical | identical |
| `ob_index_namelist.txt` | present (helper) | not copied |

`regional_run/input_baseline/` (2026-10-15 to 10-21 frozen baseline, `run_arch_HO2`
then AM98 + ICT) and `regional_run/input_arch/` are intermediate sets and are
not staged.

## 4. `input_darwin/` = `nares/regional_run/run_bgc15_s1/` BGC namelists

`run_bgc15_s1` (PBS 25321613, 2026-10-30 to 11-02, exit 0, 84 h) is segment 1
of the only full BGC chain (Jan 1 to Apr 1; s2, s3 continue to Sep 1 with
only `nIter0`/`nTimeSteps` changed). Copied: `data.darwin`, `data.gchem`,
`data.oasim`, `data.ptracers`, `data.radtrans`, `data.traits` (all
byte-identical to `nares/regional_run/input_darwin/`), plus the two files that
override `input/`: `data.pkg` (BGC packages on) and `data.diagnostics`.
`data.diagnostics` is DERIVED (2026-11-16, Dustin's decision): the staged
`input/data.diagnostics` (Y69m_sw lists 1-7, 11) plus the BGC streams of
`run_bgc15_s1` (its lists 8-10 kept, its monthly tracer lists 11-43 renumbered
12-44). It has not been run in this form. The as-run `run_bgc15_s1` file had
the older `run_arch_Y16` physics lists (no `SImelt`).

Physics files of `run_bgc15_s1` that are NOT staged, because they are the
older `run_arch_Y16` physics (superseded by Y69m_sw): `data.seaice`
(`SEAICE_strength 8.25e4`, `presPow0/1 = 1/3`), `data.ggl90` (no parent set),
`data` (`nTimeSteps 129600`). `data.obcs`, `data.exf`,
`data.iceplume`, `data.shelfice`, `data.salt_plume`, `data.cal`, `eedata` are
identical to the staged `input/`. Its `data.exch2` is the 40x40 file. (`input_darwin/data.diagnostics` carries the `SImelt`
list; see above.)

## 5. `scripts/`

| file | source | role |
|---|---|---|
| `job_build.sh` | `regional_run/` | first build (2026-10-07), `-mo "code_darwin code"` |
| `job_build_capload.sh` | `regional_run/` | + `code_capload` (2026-10-24) |
| `job_build_swfix.sh` | `regional_run/` | + `code_swfix` (2026-11-08) |
| `job_build_drag10.sh` | `regional_run/` | + `code_drag10` (2026-11-15); = the staged `code/` |
| `job_build_bgc.sh` | `regional_run/` | BGC build (2026-10-29) |
| `job_run_arch.sh`, `job_run_arch_long.sh` | `prep/` | run scripts (192 ranks; the long one checks OBCS record counts) |
| `job_build_v05.sh` | written 2026-11-22 from `regional_run/job_build_v05test.sh` | production build of this tree (`-mo "code_darwin code"`), checks 678 ranks; paths from `P`, `B`, `BUILD`, `DARWIN3` |
| `job_run_scaling.sh` | `prep/` | `job_run_arch_long.sh` with the rank count from `NP`; ran the 20x20 scaling test (678 ranks, 18 sky_ele nodes) |
| `stage_arch_tc4.sh` | `prep/` | stages a cold-start run on the baseline namelists (walled bathymetry, tides, grid18) |
| `stage_arch_restart_long.sh` | `prep/` | continues a run from a pickup with the Jan to Sep OBCS set |
| `stage_y69m_sw.sh`, `stage_drag10_test.sh`, `stage_drag10_full.sh` | `prep/` | how Y69m_sw, D10t and D10 were staged |
| `stage_bgc_test.sh`, `stage_bgc_prod.sh` | `prep/` | the 2-day BGC test and the 3-segment BGC chain |
| `make_swfix.py`, `make_drag10.py` | `prep/` | generators of `code_swfix` and `code_drag10` from the darwin3 sources |
| `check_arch_run.py`, `check_bgc_test.py` | `prep/` | run health checks |

Not staged: `job_build_fix.sh`, `job_build_emul.sh` (dead-end `code_fix` builds),
`job_build_bgc_s20x20.sh`, `job_build_bgc_s40x20.sh`, `make_build_scaling.py`,
`stage_scaling_1d.sh`, `job_run_scaling.sh` (pending scaling test), the 20
other `stage_arch_*.sh` sensitivity scripts, movie/analysis/disk-scan jobs.

## 6. `scripts/pipeline/`

See `scripts/pipeline/README.md`. Sources: `nares/grid/*.py`, `nares/gen_ncgrid/`,
`nares/prep/`, `nares/iceplume/`, `nares/forcings/iceplume/job_runoffQsg.sh`,
the shared parent extraction `/nobackupp19/dcarrol2/downscaling/parent_V5r1_westAO_nares/`
(`code_darwin/`, `input/`, `job_build.sh`, `job_ECCO_darwin`, `chain_next.sh`,
`check_dv.py`), and four scripts that live only on the Mac
(`~/Documents/research/downscaling/python/nares_strait/{gen_exch2_blanklist,gen_ob_index,nares_ob_points,check_tide_files}.py`,
`python/tides/gen_obcs_tides.py`, `model_setup/nares_strait/analysis/make_ic_thin_mouth.py`).

## 7. `pkg/iceplume/`

15 files = `nares/darwin3/pkg/iceplume/` at HEAD 8eb545d38 (byte-identical
check 2026-11-16 against the local fixed fork
`model_setup/nares_strait/iceplume_pkg_fixes/iceplume/pkg/iceplume`, branch
`nares_fixes` c2e6dbd). Upstream: https://github.com/mhwood/iceplume at
e4d4ac9. `LICENSE` = the fork's BSD-3 license file. `README.md` written here.
