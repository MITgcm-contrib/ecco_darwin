# pkg/thsice

Thermodynamic sea ice (Winton 2000 two-layer) usable alone or with seaice dynamics.

**runtime switch:** `useTHSICE`-style flag in `data.pkg` (check exact name in packages_boot.F)
**reads:** `data.ice`
**manual:** `doc/examples/examples.rst`, `doc/outp_pkgs/outp_pkgs.rst`, `doc/phys_pkgs/bulk_force.rst`, `doc/phys_pkgs/seaice.rst`
**adjoint support files:** thsice_ad_check_lev1_dir.h, thsice_ad_check_lev2_dir.h, thsice_ad_check_lev3_dir.h, thsice_ad_check_lev4_dir.h, thsice_ad_diff.list

## Namelist parameters
### THSICE_CONST
- `rhos` — density of snow [kg/m^3]
- `rhoi` — density of ice [kg/m^3]
- `rhosw` — density of seawater [kg/m^3]
- `rhofw` — density of fresh water [kg/m^3]
- `cpIce` — specific heat of fresh ice [J/kg/K]
- `cpWater` — specific heat of water [J/kg/K]
- `kIce` — thermal conductivity of pure ice [W/m/K]
- `kSnow` — thermal conductivity of snow [W/m/K]
- `bMeltCoef` — base-melting heat transfer coefficient (between ice & water) [no unit]
- `Lfresh` — latent heat of melting of pure ice [J/kg]
- `qsnow` — snow enthalpy [J/kg]
- `albColdSnow` — albedo of cold (=dry) new snow (Tsfc < tempSnowAlb)
- `albWarmSnow` — albedo of warm (=wet) new snow (Tsfc = 0)
- `tempSnowAlb` — temperature transition from ColdSnow to WarmSnow Alb. [oC]
- `albOldSnow` — albedo of old snow (snowAge > 35.d)
- `hNewSnowAge` — new snow thickness that refresh the snow-age (by 1/e)
- `snowAgTime` — snow aging time scale [s]
- `albIceMax` — max albedo of bare ice (thick ice)
- `albIceMin` — minimum ice albedo (very thin ice)
- `hAlbIce` — ice thickness for albedo transition: thin/thick ice albedo
- `hAlbSnow` — snow thickness for albedo transition: snow/ice albedo
- `i0swFrac` — fraction of penetrating solar rad
- `ksolar` — bulk solar abs coeff of sea ice [m^-1]
- `dhSnowLin` — half slope of linear distribution of snow thickness within the grid-cell (from hSnow-dhSnow to hSnow+dhSnow, if full ice & snow cover) [m] ; (only used for SW radiation).
- `saltIce` — salinity of ice [g/kg]
- `S_winton` — Winton salinity of ice [g/kg]
- `mu_Tf` — linear dependence of melting temperature on salinity [oC/(g/kg)] Tf(sea-water) = -mu_Tf * S
- `Tf0kel` — Freezing temp of fresh water in Kelvin = 273.15
- `Terrmax` — Temperature convergence criteria [oC]
- `nitMaxTsf` — maximum Nb of iter to find Surface Temp (Trsf)
- `hIceMin` — Minimum ice  thickness [m]
- `hiMax` — Maximum ice  thickness [m]
- `hsMax` — Maximum snow thickness [m]
- `iceMaskMax` — maximum Ice fraction (=1 for no fractional ice)
- `iceMaskMin` — mimimum Ice fraction (=1 for no fractional ice)
- `fracEnMelt` — fraction of energy going to lateral melting (vs height decrease) (=0 for no fract. ice)
- `fracEnFreez` — fraction of energy going to lateral freezing (vs height increase)
- `hThinIce` — ice height above which fracEnMelt/Freez are applied [m] (=hIceMin for no fractional ice)
- `hThickIce` — ice height below which fracEnMelt/Freez are applied [m] (=large for no fractional ice)
- `hNewIceMax` — new ice maximum thickness [m]
### THSICE_PARM01
- `startIceModel` — =1 : start ice model at nIter0 ; =0 : use pickup files
- `stepFwd_oceMxL` — step forward mixed-layer T & S (slab-ocean)
- `thSIce_calc_albNIR` — calculate Near Infra-Red Albedo
- `thSIce_skipThermo` — by-pass seaice thermodynamics
- `thSIce_deltaT` — ice model time-step, seaice thicken/extend [s]
- `thSIce_dtTemp` — ice model time-step, solve4temp [s]
- `ocean_deltaT` — ocean mixed-layer time-step [s]
- `tauRelax_MxL` — Relaxation time scale for MixLayer T [s]
- `tauRelax_MxL_salt` — Relaxation time scale for MixLayer S [s]
- `hMxL_default` — default value for ocean MixLayer thickness [m]
- `sMxL_default` — default value for salinity in MixLayer [g/kg]
- `vMxL_default` — default value for ocean current velocity in MxL [m/s]
- `thSIce_diffK` — thickness (horizontal) diffusivity [m^2/s]
- `thSIceAdvScheme` — thSIce Advection scheme selector
- `stressReduction` — reduction factor for wind-stress under sea-ice [0-1]
- `thSIceBalanceAtmFW` — select balancing Fresh-Water flux from Atm+Land
- `thSIce_taveFreq`
- `thSIce_diagFreq` — Frequency^-1 for diagnostic output [s]
- `thSIce_monFreq` — Frequency^-1 for monitor    output [s]
- `thSIce_tave_mnc`
- `thSIce_snapshot_mnc` — write snap-shot output   using MNC
- `thSIce_mon_mnc` — write monitor to netcdf file
- `thSIce_pickup_read_mnc` — pickup read w/ MNC
- `thSIce_pickup_write_mnc` — pickup write w/ MNC
- `thSIceFract_InitFile` — File name for initial ice fraction
- `thSIceThick_InitFile` — File name for initial ice thickness
- `thSIceSnowH_InitFile` — File name for initial snow thickness
- `thSIceSnowA_InitFile` — File name for initial snow Age
- `thSIceEnthp_InitFile` — File name for initial ice enthalpy
- `thSIceTsurf_InitFile` — File name for initial surf. temp
### THSICE_COST
- `mult_thsice`  _[ifdef ALLOW_COST]_
- `thsice_cost_ice_flag`  _[ifdef ALLOW_COST]_

## CPP options (defaults as shipped)
- `THSICE_FRACEN_POWERLAW` (define, THSICE_OPTIONS.h) — - use continuous power-law function for partition of energy between lateral melting/freezing and thinning/thickening ; otherwise, use step function.
- `ALLOW_DBUG_THSICE` (define, THSICE_OPTIONS.h) — - allow single grid-point debugging write to standard-output
- `CHECK_ENERGY_CONSERV` (undef, THSICE_OPTIONS.h) — - only to check conservation (change content of ICE_qleft,fresh,salFx-T files)
- `THSICE_REGULARIZE_CALC_THICKN` (undef, THSICE_OPTIONS.h) — - replace MIN/MAX by smooth functions, avoid divisions by zero and sqrt of zero mostly to help AD code generation, changes results

## Headers
- `THSICE_COST.h` — /==========================================================\ Sea ice cost terms. \==========================================================/ objf_ths
- `THSICE_DEBUG.h` — BOP
- `THSICE_OPTIONS.h` — Package-specific Options & Macros go here
- `THSICE_PARAMS.h` — Header file for Therm_SeaIce package parameters: - basic parameter ( I/O frequency, etc ...) - physical constants (used in therm_SeaIce pkg)
- `THSICE_SIZE.h` — .. number layers of ice nlyr   ::   maximum number of ice layers
- `THSICE_VARS.h` — variable for thermodynamics - Sea-Ice model
- `thsice_ad_check_lev1_dir.h` — ADJ STORE iceMask = comlev1, key = ikey_dynamics ADJ STORE iceHeight  = comlev1, key = ikey_dynamics ADJ STORE snowHeight = comlev1, key = ikey_dynami
- `thsice_ad_check_lev2_dir.h` — ADJ STORE iceMask    = tapelev2, key = ilev_2 ADJ STORE iceHeight  = tapelev2, key = ilev_2 ADJ STORE snowHeight = tapelev2, key = ilev_2 ADJ STORE sn
- `thsice_ad_check_lev3_dir.h` — ADJ STORE iceMask    = tapelev3, key = ilev_3 ADJ STORE iceHeight  = tapelev3, key = ilev_3 ADJ STORE snowHeight = tapelev3, key = ilev_3 ADJ STORE sn
- `thsice_ad_check_lev4_dir.h` — ADJ STORE iceMask    = tapelev4, key = ilev_4 ADJ STORE iceHeight  = tapelev4, key = ilev_4 ADJ STORE snowHeight = tapelev4, key = ilev_4 ADJ STORE sn
- `thsice_test_addfluid.h` — this is a small piece of code to add to S/R thsice_main.F to check AddFluid implementation :

## Routines (41)
`thsice_advdiff.F`, `thsice_advection.F`, `thsice_albedo.F`, `thsice_ave.F`, `thsice_balance_frw.F`, `thsice_calc_thickn.F`, `thsice_check.F`, `thsice_check_conserv.F`, `thsice_cost_driver.F`, `thsice_cost_final.F`, `thsice_cost_init_varia.F`, `thsice_cost_test.F`, `thsice_diagnostics_init.F`, `thsice_diagnostics_state.F`, `thsice_diffusion.F`, `thsice_do_advect.F`, `thsice_do_exch.F`, `thsice_extend.F`, `thsice_get_bulkf.F`, `thsice_get_exf.F`, `thsice_get_ocean.F`, `thsice_get_precip.F`, `thsice_get_velocity.F`, `thsice_impl_temp.F`, `thsice_ini_vars.F`, `thsice_init_fixed.F`, `thsice_main.F`, `thsice_map_exf.F`, `thsice_mnc_init.F`, `thsice_monitor.F`, `thsice_output.F`, `thsice_read_pickup.F`, `thsice_readparms.F`, `thsice_reshape_layers.F`, `thsice_salt_plume.F`, `thsice_slab_ocean.F`, `thsice_solve4temp.F`, `thsice_step_fwd.F`, `thsice_step_temp.F`, `thsice_turnoff_io.F`, `thsice_write_pickup.F`

## Called from outside the package
- `THSICE_MAIN` ← `model/src/do_oceanic_phys.F:402`
- `THSICE_DIAGNOSTICS_STATE` ← `model/src/do_statevars_diags.F:102`
- `THSICE_OUTPUT` ← `model/src/do_the_model_io.F:207`
- `THSICE_CHECK` ← `model/src/packages_check.F:380`
- `THSICE_INIT_FIXED` ← `model/src/packages_init_fixed.F:532`
- `THSICE_INI_VARS` ← `model/src/packages_init_variables.F:454`
- `THSICE_READPARMS` ← `model/src/packages_readparms.F:306`
- `THSICE_WRITE_PICKUP` ← `model/src/packages_write_pickup.F:192`
- `THSICE_TURNOFF_IO` ← `model/src/turnoff_model_io.F:113`
- `THSICE_DO_ADVECT` ← `pkg/aim_v23/aim_do_physics.F:192`
- `THSICE_DO_EXCH` ← `pkg/aim_v23/aim_do_physics.F:190`
- `THSICE_SLAB_OCEAN` ← `pkg/aim_v23/aim_do_physics.F:197`
- `THSICE_STEP_FWD` ← `pkg/aim_v23/aim_do_physics.F:174`
- `THSICE_ALBEDO` ← `pkg/aim_v23/aim_sice2aim.F:100`
- `THSICE_IMPL_TEMP` ← `pkg/aim_v23/aim_sice_impl.F:89`
- `THSICE_ALBEDO` ← `pkg/atm2d/calc_zonal_means.F:58`
- `THSICE_AVE` ← `pkg/atm2d/forward_step_atm2d.F:212`
- `THSICE_IMPL_TEMP` ← `pkg/atm2d/forward_step_atm2d.F:229`
- `THSICE_STEP_FWD` ← `pkg/atm2d/forward_step_atm2d.F:209`
- `THSICE_ALBEDO` ← `pkg/cheapaml/cheapaml_seaice.F:160`
- `THSICE_GET_OCEAN` ← `pkg/cheapaml/cheapaml_seaice.F:130`
- `THSICE_IMPL_TEMP` ← `pkg/cheapaml/cheapaml_seaice.F:221`
- `THSICE_COST_FINAL` ← `pkg/cost/cost_final.F:110`
- `THSICE_COST_INIT_VARIA` ← `pkg/cost/cost_init_varia.F:85`
- `THSICE_COST_DRIVER` ← `pkg/cost/cost_tile.F:136`
- `THSICE_DO_ADVECT` ← `pkg/seaice/seaice_model.F:213`

## Verification experiments compiling it (4)
`aim.5l_cs` `cpl_aim+ocn` `global_ocean.cs32x15` `offline_exf_seaice`
