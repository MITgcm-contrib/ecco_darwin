# pkg/land

Simple land model (bucket hydrology, soil temperature) for atmospheric configurations.

**runtime switch:** `useLAND`-style flag in `data.pkg` (check exact name in packages_boot.F)
**reads:** `data.land`
**manual:** `doc/outp_pkgs/outp_pkgs.rst`

## Namelist parameters
### LAND_MODEL_PAR
- `land_calc_grT` — step forward ground Temperature
- `land_calc_grW` — step forward soil moiture
- `land_impl_grT` — solve ground Temperature implicitly
- `land_calc_snow` — step forward snow thickness
- `land_calc_alb` — compute albedo of snow over land
- `land_oldPickup` — restart from an old pickup (= before checkpoint 52l_pre)
- `land_grT_iniFile` — File containing initial ground Temp.
- `land_grW_iniFile` — File containing initial ground Water.
- `land_snow_iniFile` — File containing initial snow thickness.
- `land_deltaT` — land model time-step
- `land_taveFreq`
- `land_diagFreq` — Frequency^-1 for diagnostic output (s)
- `land_monFreq` — Frequency^-1 for monitor    output (s)
- `land_dzF` — layer thickness
- `land_timeave_mnc`
- `land_snapshot_mnc`
- `land_mon_mnc`
- `land_pickup_write_mnc`
- `land_pickup_read_mnc`
### LAND_PHYS_PAR
- `land_grdLambda` — Thermal conductivity of the ground (W/m/K)
- `land_heatCs` — Heat capacity of dry soil (J/m3/K)
- `land_CpWater` — Heat capacity of water    (J/kg/K)
- `land_wTauDiff` — soil moisture diffusion time scale (s)
- `land_waterCap` — field capacity per meter of soil (1)
- `land_fractRunOff` — fraction of water in excess which run-off (1)
- `land_rhoLiqW` — density of liquid water (kg/m3)
- `land_rhoSnow` — density of snow (kg/m3)
- `land_Lfreez` — Latent heat of freezing (J/kg)
- `land_hMaxSnow` — Maximum snow-thickness  (m)
- `diffKsnow` — thermal conductivity of snow (W/m/K)
- `timeSnowAge` — snow aging time scale   (s)
- `hNewSnowAge` — new snow thickness that refreshes snow-age (by 1/e)
- `albColdSnow` — albedo of cold (=dry) new snow (Tsfc < tempSnowAlbL)
- `albWarmSnow` — albedo of warm (=wet) new snow (Tsfc = 0)
- `tempSnowAlbL` — temperature transition from ColdSnow to WarmSnow Alb. (oC)
- `albOldSnow` — albedo of old snow (snowAge > 35.d)
- `hAlbSnow` — snow thickness for albedo transition: snow/ground

## CPP options (defaults as shipped)
- `LAND_DEBUG` (undef, LAND_OPTIONS.h) — to write debugging diagnostics:
- `LAND_OLD_VERSION` (undef, LAND_OPTIONS.h) — to reproduce results from version.1 (not conserving heat)

## Headers
- `LAND_OPTIONS.h` — CPP options file for Land package
- `LAND_PARAMS.h` — Header file for LAND package parameters: - basic parameter ( I/O frequency, etc ...) - physical constants - vertical discretization
- `LAND_SIZE.h` — MITgcm declaration of grid size.
- `LAND_VARS.h` — Land model variables: - prognostic variables - forcing fields - diagnostic variables

## Routines (15)
`land_albedo.F`, `land_check.F`, `land_diagnostics_init.F`, `land_diagnostics_state.F`, `land_do_diags.F`, `land_impl_temp.F`, `land_ini_vars.F`, `land_init_fixed.F`, `land_mnc_init.F`, `land_monitor.F`, `land_output.F`, `land_read_pickup.F`, `land_readparms.F`, `land_stepfwd.F`, `land_write_pickup.F`

## Called from outside the package
- `LAND_DIAGNOSTICS_STATE` ← `model/src/do_statevars_diags.F:114`
- `LAND_OUTPUT` ← `model/src/do_the_model_io.F:128`
- `LAND_CHECK` ← `model/src/packages_check.F:386`
- `LAND_INIT_FIXED` ← `model/src/packages_init_fixed.F:542`
- `LAND_INI_VARS` ← `model/src/packages_init_variables.F:463`
- `LAND_READPARMS` ← `model/src/packages_readparms.F:311`
- `LAND_WRITE_PICKUP` ← `model/src/packages_write_pickup.F:200`
- `LAND_DO_DIAGS` ← `pkg/aim_v23/aim_do_physics.F:154`
- `LAND_STEPFWD` ← `pkg/aim_v23/aim_do_physics.F:150`
- `LAND_ALBEDO` ← `pkg/aim_v23/aim_land2aim.F:178`
- `LAND_IMPL_TEMP` ← `pkg/aim_v23/aim_land_impl.F:111`

## Verification experiments compiling it (2)
`aim.5l_cs` `cpl_aim+ocn`
