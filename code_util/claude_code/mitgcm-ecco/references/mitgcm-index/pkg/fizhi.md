# pkg/fizhi

Fizhi atmospheric physics (NASA GEOS-like) for atmosphere configurations.

**pkg_depend:** +gridalt +diagnostics -aim +atm_common  (`+` requires, `-` excludes)
**runtime switch:** `useFIZHI`-style flag in `data.pkg` (check exact name in packages_boot.F)
**reads:** `data.fizhi`
**manual:** `doc/outp_pkgs/outp_pkgs.rst`

## Namelist parameters
### FIZHI_LIST
- `nymdbegin`
- `nhmsbegin`
- `fizhi_mnc_write_pickup`
- `fizhi_mnc_read_pickup`
- `runlength`
- `climsst`
- `climsice`

## CPP options (defaults as shipped)
- `FIZHI_USE_FIXED_DAY` (undef, FIZHI_OPTIONS.h) — use fixed day in the year:
- `TRY_NEW_GETPWHERE` (define, FIZHI_OPTIONS.h) — use new version of S/R GETPWHERE
- `FIZHI_F77_COMPIL` (undef, FIZHI_OPTIONS.h) — Compiler and Processor specific code
- `FIZHI_CRAY` (undef, FIZHI_OPTIONS.h)
- `FIZHI_SGI` (undef, FIZHI_OPTIONS.h)
- `FIZHI_TRBFLX_OLD_BUG` (undef, FIZHI_OPTIONS.h) — Bring back original/old bug in S/R TRBFLX

## Headers
- `FIZHI_OPTIONS.h` — BOP
- `cah-dat.h` — 
- `cai-dat.h` — 
- `chronos.h` — *****                   Clock Variables                          *****
- `co2-tran3.h` — 
- `fizhi_SHP.h` — The physics state uses the dynamics dimensions in the horizontal and the land dimensions in the horizontal for turbulence variables Secret Hiding Plac
- `fizhi_SIZE.h` — Physics Grid Vertical Dimension
- `fizhi_chemistry_coms.h` — Chemistry Variables Dimensions
- `fizhi_coms.h` — The physics state uses the dynamics dimensions in the horizontal and the land dimensions in the horizontal for turbulence variables Fizhi State Common
- `fizhi_earth_coms.h` — Solid-Earth State Variables
- `fizhi_io_comms.h` — FIZHI I/O flags
- `fizhi_land_SIZE.h` — Land Grid Horizontal Dimension (Number of Tiles)
- `fizhi_land_coms.h` — Land State Common
- `fizhi_ocean_coms.h` — Ocean Parameters
- `h2o-tran3.h` — 
- `o3-tran3.h` — 
- `sibber.h` — **** CHIP HEADER FILE
- `snwmid.h` — Note: SNWMID and SNWALB parameters modified to obtain improved albedo and radswg

## Routines (156)
`AtoC.F`, `CtoA.F`, `do_fizhi.F`, `fizhi_alarms.F`, `fizhi_clockstuff.F`, `fizhi_diagalarms.F`, `fizhi_diagnostics_init.F`, `fizhi_driver.F`, `fizhi_fillnegs.F`, `fizhi_gwdrag.F`, `fizhi_init_chem.F`, `fizhi_init_fixed.F`, `fizhi_init_vars.F`, `fizhi_init_veg.F`, `fizhi_init_vegsurftiles.F`, `fizhi_lsm.F`, `fizhi_lwrad.F`, `fizhi_mnc_init.F`, `fizhi_moist.F`, `fizhi_mpistuff.F`, `fizhi_rayleigh.F`, `fizhi_read_pickup.F`, `fizhi_readparms.F`, `fizhi_readwrite_vegtiles.F`, `fizhi_step_diag.F`, `fizhi_swrad.F`, `fizhi_tendency_apply.F`, `fizhi_turb.F`, `fizhi_update_time.F`, `fizhi_utils.F`, `fizhi_wrapper.F`, `fizhi_write_datetime.F`, `fizhi_write_pickup.F`, `fizhi_write_state.F`, `getcon.F`, `getpwhere.F`, `slprs.F`, `step_fizhi_corr.F`, `step_fizhi_fg.F`, `step_physics.F`, `update_chemistry_exports.F`, `update_earth_exports.F`, `update_ocean_exports.F`

## Called from outside the package
- `FIZHI_TENDENCY_APPLY_S` ← `model/src/apply_forcing.F:868`
- `FIZHI_TENDENCY_APPLY_T` ← `model/src/apply_forcing.F:498`
- `FIZHI_TENDENCY_APPLY_U` ← `model/src/apply_forcing.F:120`
- `FIZHI_TENDENCY_APPLY_V` ← `model/src/apply_forcing.F:310`
- `FIZHI_UPDATE_TIME` ← `model/src/do_atmospheric_phys.F:125`
- `FIZHI_WRAPPER` ← `model/src/do_atmospheric_phys.F:123`
- `STEP_FIZHI_FG` ← `model/src/do_atmospheric_phys.F:124`
- `UPDATE_CHEMISTRY_EXPORTS` ← `model/src/do_atmospheric_phys.F:122`
- `UPDATE_EARTH_EXPORTS` ← `model/src/do_atmospheric_phys.F:121`
- `UPDATE_OCEAN_EXPORTS` ← `model/src/do_atmospheric_phys.F:120`
- `FIZHI_WRITE_STATE` ← `model/src/do_the_model_io.F:122`
- `FIZHI_TENDENCY_APPLY_S` ← `model/src/external_forcing.F:698`
- `FIZHI_TENDENCY_APPLY_T` ← `model/src/external_forcing.F:376`
- `FIZHI_TENDENCY_APPLY_U` ← `model/src/external_forcing.F:83`
- `FIZHI_TENDENCY_APPLY_V` ← `model/src/external_forcing.F:223`
- `STEP_FIZHI_CORR` ← `model/src/forward_step.F:1124`
- `FIZHI_INIT_FIXED` ← `model/src/packages_init_fixed.F:581`
- `FIZHI_INIT_VARS` ← `model/src/packages_init_variables.F:491`
- `FIZHI_READPARMS` ← `model/src/packages_readparms.F:368`
- `FIZHI_WRITE_DATETIME` ← `model/src/packages_write_pickup.F:217`
- `FIZHI_WRITE_PICKUP` ← `model/src/packages_write_pickup.F:215`
- `FIZHI_WRITE_VEGTILES` ← `model/src/packages_write_pickup.F:216`
- `MY_EXIT` ← `pkg/chronos/chronos.F:81`
- `MY_FINALIZE` ← `pkg/chronos/chronos.F:80`
- `QSAT` ← `pkg/diagnostics/diagnostics_fill_state.F:468`
- `FIZHI_DIAGALARMS` ← `pkg/diagnostics/diagnostics_init_fixed.F:51`

## Verification experiments compiling it (3)
`fizhi-cs-32x32x40` `fizhi-cs-aqualev20` `fizhi-gridalt-hs`
