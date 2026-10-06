# pkg/ptracers

Passive tracers: arbitrary number of tracers with their own advection schemes, diffusivities, initial files, surface/relaxation; carrier for BGC.

**pkg_depend:** +generic_advdiff  (`+` requires, `-` excludes)
**runtime switch:** `usePTRACERS`-style flag in `data.pkg` (check exact name in packages_boot.F)
**reads:** `data.ptracers`
**manual:** `doc/phys_pkgs/ptracers.rst`, `doc/autodiff/autodiff.rst`, `doc/examples/cfc_offline/cfc_offline.rst`, `doc/examples/examples.rst`, `doc/examples/global_oce_biogeo/global_oce_biogeo.rst`
**adjoint support files:** ptracers_ad_check_lev1_dir.h, ptracers_ad_check_lev2_dir.h, ptracers_ad_check_lev3_dir.h, ptracers_ad_check_lev4_dir.h, ptracers_ad_diff.list

## Namelist parameters
### PTRACERS_PARM01
- `tauTr1ClimRelax`
- `PTRACERS_numInUse` — number of tracers to use
- `PTRACERS_Iter0` — timestep number when tracers are initialized
- `PTRACERS_doAB_onGpTr` — if Adams-Bashforth time stepping is used, apply AB on tracer tendencies (rather than on Tracers)
- `PTRACERS_addSrelax2EmP` — add Salt relaxation to EmP
- `PTRACERS_startStepFwd` — time to start stepping forward this tracer
- `PTRACERS_resetFreq` — Frequency (s) to reset ptracers to original val
- `PTRACERS_resetPhase` — Phase (s) to reset ptracers
- `PTRACERS_advScheme`
- `PTRACERS_ImplVertAdv` — use Implicit Vertical Advection for this tracer
- `PTRACERS_diffKh`
- `PTRACERS_diffK4`
- `PTRACERS_diffKr`
- `PTRACERS_diffKrNr`
- `PTRACERS_ref` — vertical profile for passive tracers, in analogy to tRef and sRef, hence the name
- `PTRACERS_EvPrRn` — tracer concentration in Rain, Evap & RunOff
- `PTRACERS_useGMRedi`
- `PTRACERS_useDWNSLP`
- `PTRACERS_useKPP`
- `PTRACERS_linFSConserve` — apply mean Free-Surf source/sink at surface
- `PTRACERS_stayPositive` — use Smolarkiewicz Hack to ensure Tracer stays >0
- `PTRACERS_initialFile`
- `PTRACERS_names`
- `PTRACERS_long_names`
- `PTRACERS_units`
- `PTRACERS_useRecords` — snap-shot output: put all pTracers in one file
- `PTRACERS_dumpFreq`
- `PTRACERS_taveFreq`
- `PTRACERS_monitorFreq`
- `PTRACERS_timeave_mnc`
- `PTRACERS_snapshot_mnc`
- `PTRACERS_monitor_mnc`
- `PTRACERS_pickup_write_mnc`
- `PTRACERS_pickup_read_mnc`

## CPP options (defaults as shipped)
- `PTRACERS_ALLOW_DYN_STATE` (undef, PTRACERS_OPTIONS.h) — This enables the dynamically allocated internal state data structures for PTracers.  Needed for PTRACERS_SOM_Advection. This requires a Fortran 90 compiler!

## Headers
- `PTRACERS_FIELDS.h` — BOP
- `PTRACERS_MOD.h` — BOP
- `PTRACERS_OPTIONS.h` — CPP options file for PTRACERS package Use this file for selecting options within the PTRACERS package
- `PTRACERS_PARAMS.h` — BOP
- `PTRACERS_SIZE.h` — BOP
- `PTRACERS_START.h` — BOP
- `ptracers_ad_check_lev1_dir.h` — ADJ STORE pTracer   = comlev1, key = ikey_dynamics, kind = isbyte ADJ STORE gpTrNm1   = comlev1, key = ikey_dynamics, kind = isbyte
- `ptracers_ad_check_lev2_dir.h` — ADJ STORE pTracer = tapelev2, key = ilev_2 ADJ STORE gpTrNm1 = tapelev2, key = ilev_2
- `ptracers_ad_check_lev3_dir.h` — ADJ STORE pTracer = tapelev3, key = ilev_3 ADJ STORE gpTrNm1 = tapelev3, key = ilev_3
- `ptracers_ad_check_lev4_dir.h` — ADJ STORE pTracer = tapelev4, key = ilev_4 ADJ STORE gpTrNm1 = tapelev4, key = ilev_4
- `ptracers_adcommon.h` — --   These common blocks are extracted from the --   automatically created tangent linear code. --   You need to make sure that they are up-to-date --

## Routines (32)
`ptracers_ad_dump.F`, `ptracers_apply_forcing.F`, `ptracers_calc_wsurf_tr.F`, `ptracers_check.F`, `ptracers_check_pickup.F`, `ptracers_convect.F`, `ptracers_debug.F`, `ptracers_diagnostics_init.F`, `ptracers_diagnostics_state.F`, `ptracers_dyn_state_data_mod.F`, `ptracers_dyn_state_mod.F`, `ptracers_fields_blocking_exch.F`, `ptracers_forcing_surf.F`, `ptracers_init_fixed.F`, `ptracers_init_varia.F`, `ptracers_integrate.F`, `ptracers_mnc_init.F`, `ptracers_monitor.F`, `ptracers_monitor_ad.F`, `ptracers_output.F`, `ptracers_read_pickup.F`, `ptracers_readparms.F`, `ptracers_reset.F`, `ptracers_set_iolabel.F`, `ptracers_switch_onoff.F`, `ptracers_timeave.F`, `ptracers_turnoff_io.F`, `ptracers_write_pickup.F`, `ptracers_write_state.F`, `ptracers_write_timeave.F`, `ptracers_zonal_filt_apply.F`

## Called from outside the package
- `PTRACERS_CONVECT` ← `model/src/convective_adjustment.F:163`
- `PTRACERS_CONVECT` ← `model/src/convective_adjustment_ini.F:172`
- `PTRACERS_FIELDS_BLOCKING_EXCH` ← `model/src/do_fields_blocking_exchanges.F:97`
- `PTRACERS_DIAGNOSTICS_STATE` ← `model/src/do_statevars_diags.F:78`
- `PTRACERS_OUTPUT` ← `model/src/do_the_model_io.F:213`
- `PTRACERS_FORCING_SURF` ← `model/src/external_forcing_surf.F:191`
- `PTRACERS_RESET` ← `model/src/forward_step.F:1189`
- `PTRACERS_SWITCH_ONOFF` ← `model/src/forward_step.F:501`
- `PTRACERS_CHECK` ← `model/src/packages_check.F:301`
- `PTRACERS_INIT_FIXED` ← `model/src/packages_init_fixed.F:429`
- `PTRACERS_INIT_VARIA` ← `model/src/packages_init_variables.F:334`
- `PTRACERS_READPARMS` ← `model/src/packages_readparms.F:251`
- `PTRACERS_WRITE_PICKUP` ← `model/src/packages_write_pickup.F:156`
- `PTRACERS_CALC_WSURF_TR` ← `model/src/thermodynamics.F:177`
- `PTRACERS_DEBUG` ← `model/src/thermodynamics.F:406`
- `PTRACERS_INTEGRATE` ← `model/src/thermodynamics.F:348`
- `PTRACERS_ZONAL_FILT_APPLY` ← `model/src/tracers_correction_step.F:85`
- `PTRACERS_TURNOFF_IO` ← `model/src/turnoff_model_io.F:119`
- `PTRACERS_AD_DUMP` ← `pkg/autodiff/addummy_in_stepping.F:563`
- `PTRACERS_DEBUG` ← `pkg/longstep/longstep_thermodynamics.F:211`
- `PTRACERS_INTEGRATE` ← `pkg/longstep/longstep_thermodynamics.F:189`
- `ADPTRACERS_MONITOR` ← `pkg/monitor/monitor_ad.F:253`
- `PTRACERS_FIELDS_BLOCKING_EXCH` ← `pkg/obcs/obcs_init_variables.F:461`

## Verification experiments compiling it (12)
`cfc_example` `exp4` `global_oce_biogeo_bling` `global_ocean.90x40x15` `lab_sea` `so_box_biogeo` `tutorial_advection_in_gyre` `tutorial_cfc_offline` `tutorial_dic_adjoffline` `tutorial_global_oce_biogeo` `tutorial_global_oce_latlon` `tutorial_tracer_adjsens`
