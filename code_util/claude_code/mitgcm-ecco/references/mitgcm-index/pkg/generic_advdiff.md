# pkg/generic_advdiff

Generic advection-diffusion operators (2nd/3rd/4th order, DST, flux-limited, OS7MP, Prather SOM, multi-dim) used for T, S and ptracers.

**in groups:** gfd
**always-on utility package** (no data.pkg switch)
**manual:** `doc/phys_pkgs/generic_advdiff.rst`, `doc/algorithm/adv-schemes.rst`, `doc/algorithm/algorithm.rst`, `doc/examples/held_suarez_cs/held_suarez_cs.rst`, `doc/getting_started/getting_started.rst`
**adjoint support files:** gad_ad_check_lev1_dir.h, gad_ad_check_lev2_dir.h, gad_ad_check_lev3_dir.h, gad_ad_check_lev4_dir.h, gad_ad_diff.list, gad_check.F

## CPP options (defaults as shipped)
- `COSINEMETH_III` (define, GAD_OPTIONS.h) — bi-harmonic diffusivity -- only on lat-lon grid. Setting this flag here only affects tracer diffusivity; to use it in the momentum equations it needs to be set in MOM_COMMON_OPTIONS.h
- `ISOTROPIC_COS_SCALING` (undef, GAD_OPTIONS.h) — diffusivity when using the COSINE(lat) scaling -- only on lat-lon grid. Setting this flag here only affects tracer diffusivity; to use it in the momentum equations it needs to be set in MOM_COMMON_OPTIONS.h
- `DISABLE_MULTIDIM_ADVECTION` (undef, GAD_OPTIONS.h) — As of checkpoint41, the inclusion of multi-dimensional advection introduces excessive recomputation/storage for the adjoint. We can disable it here using CPP because run-time flags are insufficient.
- `GAD_MULTIDIM_COMPRESSIBLE` (undef, GAD_OPTIONS.h) — Use compressible flow method for multi-dim advection instead of old, less accurate jmc method. Note: option has no effect on SOM advection which always use compressible flow method.
- `GAD_ALLOW_TS_SOM_ADV` (undef, GAD_OPTIONS.h) — This enable the use of 2nd-Order Moment advection scheme (Prather, 1986) for Temperature and Salinity ; due to large memory space (10 times more / tracer) requirement, by default, this part of the code is not compiled.
- `GAD_SMOLARKIEWICZ_HACK` (undef, GAD_OPTIONS.h) — This hack applies to all tracers except temperature and salinity! Do not use with Adams-Bashforth (for ptracers)! Do not use with OBCS!
- `DISABLE_MULTIDIM_ADVECTION` (define, GAD_OPTIONS.h) — If GAD is disabled then so is multi-dimensional advection

## Headers
- `GAD.h` — BOP
- `GAD_FLUX_LIMITER.h` — BOP
- `GAD_OPTIONS.h` — BOP
- `GAD_SOM_VARS.h` — BOP
- `gad_ad_check_lev1_dir.h` — ADJ STORE som_S      = comlev1, key = ikey_dynamics, kind = isbyte ADJ STORE som_T      = comlev1, key = ikey_dynamics, kind = isbyte
- `gad_ad_check_lev2_dir.h` — ADJ STORE som_S = tapelev2, key = ilev_2 ADJ STORE som_T = tapelev2, key = ilev_2
- `gad_ad_check_lev3_dir.h` — ADJ STORE som_S = tapelev3, key = ilev_3 ADJ STORE som_T = tapelev3, key = ilev_3
- `gad_ad_check_lev4_dir.h` — ADJ STORE som_S = tapelev4, key = ilev_4 ADJ STORE som_T = tapelev4, key = ilev_4

## Routines (101)
`gad_advection.F`, `gad_advscheme.F`, `gad_biharm_r.F`, `gad_biharm_x.F`, `gad_biharm_y.F`, `gad_c2_adv_r.F`, `gad_c2_adv_x.F`, `gad_c2_adv_y.F`, `gad_c2_impl_r.F`, `gad_c4_adv_r.F`, `gad_c4_adv_x.F`, `gad_c4_adv_y.F`, `gad_calc_rhs.F`, `gad_check.F`, `gad_del2.F`, `gad_diagnostics_init.F`, `gad_diagnostics_state.F`, `gad_diff_r.F`, `gad_diff_x.F`, `gad_diff_y.F`, `gad_dst2u1_adv_r.F`, `gad_dst2u1_adv_x.F`, `gad_dst2u1_adv_y.F`, `gad_dst2u1_impl_r.F`, `gad_dst3_adv_r.F`, `gad_dst3_adv_x.F`, `gad_dst3_adv_y.F`, `gad_dst3fl_adv_r.F`, `gad_dst3fl_adv_x.F`, `gad_dst3fl_adv_y.F`, `gad_dst3fl_impl_r.F`, `gad_exch_som.F`, `gad_fluxlimit_adv_r.F`, `gad_fluxlimit_adv_x.F`, `gad_fluxlimit_adv_y.F`, `gad_fluxlimit_impl_r.F`, `gad_grad_x.F`, `gad_grad_y.F`, `gad_implicit_r.F`, `gad_init_fixed.F`, `gad_init_varia.F`, `gad_os7mp_adv_r.F`, `gad_os7mp_adv_x.F`, `gad_os7mp_adv_y.F`, `gad_osc_hat_r.F`, `gad_osc_hat_x.F`, `gad_osc_hat_y.F`, `gad_osc_mul_r.F`, `gad_osc_mul_x.F`, `gad_osc_mul_y.F`, `gad_plm_fun.F`, `gad_ppm_adv_r.F`, `gad_ppm_adv_x.F`, `gad_ppm_adv_y.F`, `gad_ppm_flx_r.F`, `gad_ppm_flx_x.F`, `gad_ppm_flx_y.F`, `gad_ppm_fun.F`, `gad_ppm_hat_r.F`, `gad_ppm_hat_x.F`, `gad_ppm_hat_y.F`, `gad_ppm_p3e_r.F`, `gad_ppm_p3e_x.F`, `gad_ppm_p3e_y.F`, `gad_pqm_adv_r.F`, `gad_pqm_adv_x.F`, `gad_pqm_adv_y.F`, `gad_pqm_flx_r.F`, `gad_pqm_flx_x.F`, `gad_pqm_flx_y.F`, `gad_pqm_fun.F`, `gad_pqm_hat_r.F`, `gad_pqm_hat_x.F`, `gad_pqm_hat_y.F`, `gad_pqm_p5e_r.F`, `gad_pqm_p5e_x.F`, `gad_pqm_p5e_y.F`, `gad_read_pickup.F`, `gad_som_adv_r.F`, `gad_som_adv_x.F`, `gad_som_adv_y.F`, `gad_som_advect.F`, `gad_som_exchanges.F`, `gad_som_fill_cs_corner.F`, `gad_som_lim_r.F`, `gad_som_prep_cs_corner.F`, `gad_u3_adv_r.F`, `gad_u3_adv_x.F`, `gad_u3_adv_y.F`, `gad_u3c4_impl_r.F`, `gad_write_pickup.F`, `salt_fill.F`

## Called from outside the package
- `GAD_SOM_EXCHANGES` ← `model/src/do_fields_blocking_exchanges.F:79`
- `GAD_DIAGNOSTICS_STATE` ← `model/src/do_statevars_diags.F:71`
- `GAD_CHECK` ← `model/src/packages_check.F:517`
- `GAD_INIT_FIXED` ← `model/src/packages_init_fixed.F:195`
- `GAD_INIT_VARIA` ← `model/src/packages_init_variables.F:196`
- `GAD_WRITE_PICKUP` ← `model/src/packages_write_pickup.F:97`
- `GAD_ADVECTION` ← `model/src/salt_integrate.F:275`
- `GAD_CALC_RHS` ← `model/src/salt_integrate.F:334,349`
- `GAD_IMPLICIT_R` ← `model/src/salt_integrate.F:487`
- `GAD_SOM_ADVECT` ← `model/src/salt_integrate.F:261`
- `GAD_ADVECTION` ← `model/src/temp_integrate.F:277`
- `GAD_CALC_RHS` ← `model/src/temp_integrate.F:336,351`
- `GAD_IMPLICIT_R` ← `model/src/temp_integrate.F:489`
- `GAD_SOM_ADVECT` ← `model/src/temp_integrate.F:263`
- `SALT_FILL` ← `model/src/tracers_correction_step.F:94`
- `GAD_C2_ADV_X` ← `pkg/cheapaml/cheapaml_calc_rhs.F:168`
- `GAD_C2_ADV_Y` ← `pkg/cheapaml/cheapaml_calc_rhs.F:245`
- `GAD_DIFF_X` ← `pkg/cheapaml/cheapaml_calc_rhs.F:173`
- `GAD_DIFF_Y` ← `pkg/cheapaml/cheapaml_calc_rhs.F:250`
- `GAD_DST3FL_ADV_X` ← `pkg/cheapaml/cheapaml_calc_rhs.F:161`
- `GAD_DST3FL_ADV_Y` ← `pkg/cheapaml/cheapaml_calc_rhs.F:238`
- `GAD_EXCH_SOM` ← `pkg/ptracers/ptracers_fields_blocking_exch.F:51`
- `GAD_ADVECTION` ← `pkg/ptracers/ptracers_integrate.F:256`
- `GAD_CALC_RHS` ← `pkg/ptracers/ptracers_integrate.F:315`
- `GAD_IMPLICIT_R` ← `pkg/ptracers/ptracers_integrate.F:457`
- `GAD_SOM_ADVECT` ← `pkg/ptracers/ptracers_integrate.F:239`
- `GAD_EXCH_SOM` ← `pkg/ptracers/ptracers_read_pickup.F:319`
- `GAD_DST2U1_ADV_X` ← `pkg/seaice/seaice_advection.F:372`
- `GAD_DST2U1_ADV_Y` ← `pkg/seaice/seaice_advection.F:586`
- `GAD_DST3FL_ADV_X` ← `pkg/seaice/seaice_advection.F:390`
- `GAD_DST3FL_ADV_Y` ← `pkg/seaice/seaice_advection.F:604`
- `GAD_DST3_ADV_X` ← `pkg/seaice/seaice_advection.F:386`
- `GAD_DST3_ADV_Y` ← `pkg/seaice/seaice_advection.F:600`
- `GAD_FLUXLIMIT_ADV_X` ← `pkg/seaice/seaice_advection.F:382`
- `GAD_FLUXLIMIT_ADV_Y` ← `pkg/seaice/seaice_advection.F:596`
- `GAD_OS7MP_ADV_X` ← `pkg/seaice/seaice_advection.F:394`
- `GAD_OS7MP_ADV_Y` ← `pkg/seaice/seaice_advection.F:608`
- `GAD_PPM_ADV_X` ← `pkg/seaice/seaice_advection.F:401`
- `GAD_PPM_ADV_Y` ← `pkg/seaice/seaice_advection.F:615`
- `GAD_PQM_ADV_X` ← `pkg/seaice/seaice_advection.F:407`
- `GAD_PQM_ADV_Y` ← `pkg/seaice/seaice_advection.F:621`
- `GAD_DIFF_X` ← `pkg/seaice/seaice_diffusion.F:79`
- `GAD_DIFF_Y` ← `pkg/seaice/seaice_diffusion.F:81`
- `GAD_DST2U1_ADV_X` ← `pkg/thsice/thsice_advection.F:353`
- `GAD_DST2U1_ADV_Y` ← `pkg/thsice/thsice_advection.F:574`
- `GAD_DST3FL_ADV_X` ← `pkg/thsice/thsice_advection.F:370`
- `GAD_DST3FL_ADV_Y` ← `pkg/thsice/thsice_advection.F:591`
- `GAD_DST3_ADV_X` ← `pkg/thsice/thsice_advection.F:366`
- `GAD_DST3_ADV_Y` ← `pkg/thsice/thsice_advection.F:587`
- `GAD_FLUXLIMIT_ADV_X` ← `pkg/thsice/thsice_advection.F:362`
- `GAD_FLUXLIMIT_ADV_Y` ← `pkg/thsice/thsice_advection.F:583`

## Verification experiments compiling it (56)
`1D_ocean_ice_column` `MLAdjust` `advect_cs` `advect_xz` `aim.5l_Equatorial_Channel` `aim.5l_LatLon` `aim.5l_cs` `atm_gray` `bottom_ctrl_5x5` `cfc_example` `cheapAML_box` `cpl_aim+ocn` `deep_anelastic` `dome` `exp2` `exp4` `fizhi-cs-32x32x40` `fizhi-cs-aqualev20` `fizhi-gridalt-hs` `front_relax` `global_oce_biogeo_bling` `global_oce_latlon` `global_ocean.90x40x15` `global_ocean.cs32x15` `hs94.128x64x5` `hs94.1x64x5` `hs94.cs-32x32x5` `ideal_2D_oce` `internal_wave` `inverted_barometer` `isomip` `lab_sea` `matrix_example` `obcs_ctrl` `offline_exf_seaice` `seaice_itd` `seaice_obcs` `shelfice_2d_remesh` `short_surf_wave` `so_box_biogeo` `solid-body.cs-32x32x1` `tutorial_advection_in_gyre` `tutorial_baroclinic_gyre` `tutorial_cfc_offline` `tutorial_deep_convection` `tutorial_dic_adjoffline` `tutorial_global_oce_biogeo` `tutorial_global_oce_in_p` `tutorial_global_oce_latlon` `tutorial_global_oce_optim` `tutorial_held_suarez_cs` `tutorial_plume_on_slope` `tutorial_reentrant_channel` `tutorial_rotating_tank` `tutorial_tracer_adjsens` `vermix`
