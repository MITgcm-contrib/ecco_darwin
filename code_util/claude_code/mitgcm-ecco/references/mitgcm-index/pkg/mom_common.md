# pkg/mom_common

Momentum code shared by flux-form and vector-invariant: viscosity, bottom/side drag, Smagorinsky/Leith, metric terms.

**in groups:** gfd
**always-on utility package** (no data.pkg switch)
**manual:** `doc/algorithm/algorithm.rst`, `doc/getting_started/getting_started.rst`, `doc/ocean_state_est/ocean_state_est.rst`, `doc/outp_pkgs/outp_pkgs.rst`, `doc/phys_pkgs/mom_packages.rst`
**adjoint support files:** mom_common_ad_diff.list

## CPP options (defaults as shipped)
- `COSINEMETH_III` (define, MOM_COMMON_OPTIONS.h) — bi-harmonic viscosity -- only on lat-lon grid. Setting this flag here only affects momentum viscosity; to use it in the tracer equations it needs to be set in GAD_OPTIONS.h
- `ISOTROPIC_COS_SCALING` (undef, MOM_COMMON_OPTIONS.h) — viscosity when using the COSINE(lat) scaling -- only on lat-lon grid. Setting this flag here only affects momentum viscosity; to use it in the tracer equations it needs to be set in GAD_OPTIONS.h
- `ALLOW_LEITH_QG` (undef, MOM_COMMON_OPTIONS.h) — allow LeithQG coefficient to be calculated
- `ALLOW_SMAG_3D` (undef, MOM_COMMON_OPTIONS.h) — allow isotropic 3-D Smagorinsky viscosity
- `ALLOW_3D_VISCAH` (undef, MOM_COMMON_OPTIONS.h) — allow full 3D specification of horizontal Laplacian Viscosity
- `ALLOW_3D_VISCA4` (undef, MOM_COMMON_OPTIONS.h) — allow full 3D specification of horizontal Biharmonic Viscosity
- `ALLOW_BOTTOMDRAG_ROUGHNESS` (undef, MOM_COMMON_OPTIONS.h) — Compute bottom drag coefficents, following the logarithmic law of the wall, as a function of grid cell thickness and roughness length zRoughBot (order 0.01m), assuming a von Karman constant = 0.4.
- `MOM_USE_OLD_DEEP_VERT_ADV` (undef, MOM_COMMON_OPTIONS.h) — -- Deep-Model: use the original vertical advection (with all NH metric terms) instead of updated version that advects the product deepFac x (u,v) thus removing the need for NH-metric terms: w.(u,v)/r
- `ALLOW_MOM_TEND_EXTRA_DIAGS` (undef, MOM_COMMON_OPTIONS.h) — drag (selectImplicitDrag = 2). For non-r* vertical coordinate, these tendency diagnostics can be derived from existing diagnostics for frictional stress (botTauX etc.) and therefore are not necessarily needed.

## Headers
- `MOM_COMMON_OPTIONS.h` — CPP options file for mom_common package Use this file for selecting CPP options within the mom_common package
- `MOM_VISC.h` — - Common file for length scales

## Routines (33)
`mom_calc_3d_strain.F`, `mom_calc_absvort3.F`, `mom_calc_hdiv.F`, `mom_calc_hfacz.F`, `mom_calc_ke.F`, `mom_calc_relvort3.F`, `mom_calc_smag_3d.F`, `mom_calc_strain.F`, `mom_calc_tension.F`, `mom_calc_visc.F`, `mom_diagnostics_init.F`, `mom_hdissip.F`, `mom_init_fixed.F`, `mom_quasihydrostatic.F`, `mom_u_botdrag_coeff.F`, `mom_u_coriolis_nh.F`, `mom_u_implicit_r.F`, `mom_u_metric_nh.F`, `mom_u_rviscflux.F`, `mom_u_sidedrag.F`, `mom_uv_smag_3d.F`, `mom_v_botdrag_coeff.F`, `mom_v_coriolis_nh.F`, `mom_v_implicit_r.F`, `mom_v_metric_nh.F`, `mom_v_rviscflux.F`, `mom_v_sidedrag.F`, `mom_visc_qgl_limit.F`, `mom_visc_qgl_stretch.F`, `mom_w_coriolis_nh.F`, `mom_w_metric_nh.F`, `mom_w_sidedrag.F`, `mom_w_smag_3d.F`

## Called from outside the package
- `MOM_W_CORIOLIS_NH` ← `model/src/calc_gw.F:629`
- `MOM_W_METRIC_NH` ← `model/src/calc_gw.F:608`
- `MOM_W_SIDEDRAG` ← `model/src/calc_gw.F:459`
- `MOM_W_SMAG_3D` ← `model/src/calc_gw.F:475`
- `MOM_QUASIHYDROSTATIC` ← `model/src/calc_phi_hyd.F:182,339,450`
- `MOM_CALC_3D_STRAIN` ← `model/src/dynamics.F:395`
- `MOM_CALC_SMAG_3D` ← `model/src/dynamics.F:539`
- `MOM_UV_SMAG_3D` ← `model/src/dynamics.F:544`
- `MOM_U_IMPLICIT_R` ← `model/src/dynamics.F:576`
- `MOM_V_IMPLICIT_R` ← `model/src/dynamics.F:578`
- `MOM_INIT_FIXED` ← `model/src/packages_init_fixed.F:204`
- `MOM_CALC_HDIV` ← `pkg/gmredi/gmredi_calc_qgleith.F:111`
- `MOM_CALC_HFACZ` ← `pkg/gmredi/gmredi_calc_qgleith.F:100`
- `MOM_CALC_RELVORT3` ← `pkg/gmredi/gmredi_calc_qgleith.F:110`
- `MOM_VISC_QGL_LIMIT` ← `pkg/gmredi/gmredi_calc_qgleith.F:131`
- `MOM_VISC_QGL_STRETCH` ← `pkg/gmredi/gmredi_calc_qgleith.F:128`
- `MOM_CALC_HDIV` ← `pkg/mom_fluxform/mom_fluxform.F:331`
- `MOM_CALC_HFACZ` ← `pkg/mom_fluxform/mom_fluxform.F:280`
- `MOM_CALC_KE` ← `pkg/mom_fluxform/mom_fluxform.F:329`
- `MOM_CALC_RELVORT3` ← `pkg/mom_fluxform/mom_fluxform.F:332`
- `MOM_CALC_STRAIN` ← `pkg/mom_fluxform/mom_fluxform.F:334`
- `MOM_CALC_TENSION` ← `pkg/mom_fluxform/mom_fluxform.F:333`
- `MOM_CALC_VISC` ← `pkg/mom_fluxform/mom_fluxform.F:454`
- `MOM_U_BOTDRAG_COEFF` ← `pkg/mom_fluxform/mom_fluxform.F:675`
- `MOM_U_CORIOLIS_NH` ← `pkg/mom_fluxform/mom_fluxform.F:1113`
- `MOM_U_METRIC_NH` ← `pkg/mom_fluxform/mom_fluxform.F:738`
- `MOM_U_RVISCFLUX` ← `pkg/mom_fluxform/mom_fluxform.F:622,623`
- `MOM_U_SIDEDRAG` ← `pkg/mom_fluxform/mom_fluxform.F:658`
- `MOM_VISC_QGL_LIMIT` ← `pkg/mom_fluxform/mom_fluxform.F:340`
- `MOM_VISC_QGL_STRETCH` ← `pkg/mom_fluxform/mom_fluxform.F:337`
- `MOM_V_BOTDRAG_COEFF` ← `pkg/mom_fluxform/mom_fluxform.F:970`
- `MOM_V_CORIOLIS_NH` ← `pkg/mom_fluxform/mom_fluxform.F:1121`
- `MOM_V_METRIC_NH` ← `pkg/mom_fluxform/mom_fluxform.F:1033`
- `MOM_V_RVISCFLUX` ← `pkg/mom_fluxform/mom_fluxform.F:917,918`
- `MOM_V_SIDEDRAG` ← `pkg/mom_fluxform/mom_fluxform.F:953`
- `MOM_CALC_ABSVORT3` ← `pkg/mom_vecinv/mom_vecinv.F:670`
- `MOM_CALC_HDIV` ← `pkg/mom_vecinv/mom_vecinv.F:329,404`
- `MOM_CALC_HFACZ` ← `pkg/mom_vecinv/mom_vecinv.F:263`
- `MOM_CALC_KE` ← `pkg/mom_vecinv/mom_vecinv.F:284`
- `MOM_CALC_RELVORT3` ← `pkg/mom_vecinv/mom_vecinv.F:286,405`
- `MOM_CALC_STRAIN` ← `pkg/mom_vecinv/mom_vecinv.F:333`
- `MOM_CALC_TENSION` ← `pkg/mom_vecinv/mom_vecinv.F:332`
- `MOM_CALC_VISC` ← `pkg/mom_vecinv/mom_vecinv.F:378`
- `MOM_HDISSIP` ← `pkg/mom_vecinv/mom_vecinv.F:421`
- `MOM_U_BOTDRAG_COEFF` ← `pkg/mom_vecinv/mom_vecinv.F:488`
- `MOM_U_CORIOLIS_NH` ← `pkg/mom_vecinv/mom_vecinv.F:879`
- `MOM_U_METRIC_NH` ← `pkg/mom_vecinv/mom_vecinv.F:902`
- `MOM_U_RVISCFLUX` ← `pkg/mom_vecinv/mom_vecinv.F:442`
- `MOM_U_SIDEDRAG` ← `pkg/mom_vecinv/mom_vecinv.F:470`
- `MOM_VISC_QGL_LIMIT` ← `pkg/mom_vecinv/mom_vecinv.F:349`
- `MOM_VISC_QGL_STRETCH` ← `pkg/mom_vecinv/mom_vecinv.F:346`
- `MOM_V_BOTDRAG_COEFF` ← `pkg/mom_vecinv/mom_vecinv.F:598`
- `MOM_V_CORIOLIS_NH` ← `pkg/mom_vecinv/mom_vecinv.F:887`
- `MOM_V_METRIC_NH` ← `pkg/mom_vecinv/mom_vecinv.F:903`
- `MOM_V_RVISCFLUX` ← `pkg/mom_vecinv/mom_vecinv.F:554`
- `MOM_V_SIDEDRAG` ← `pkg/mom_vecinv/mom_vecinv.F:580`
- `MOM_CALC_HFACZ` ← `pkg/seaice/seaice_mom_advection.F:101`
- `MOM_CALC_KE` ← `pkg/seaice/seaice_mom_advection.F:111`
- `MOM_CALC_RELVORT3` ← `pkg/seaice/seaice_mom_advection.F:113`
- `MOM_CALC_HDIV` ← `pkg/shap_filt/shap_filt_uv_s2.F:133`
- `MOM_CALC_HFACZ` ← `pkg/shap_filt/shap_filt_uv_s2.F:129`
- `MOM_CALC_RELVORT3` ← `pkg/shap_filt/shap_filt_uv_s2.F:143`

## Verification experiments compiling it (54)
`1D_ocean_ice_column` `MLAdjust` `adjustment.cs-32x32x1` `advect_cs` `aim.5l_Equatorial_Channel` `aim.5l_LatLon` `aim.5l_cs` `atm_gray` `bottom_ctrl_5x5` `cfc_example` `cheapAML_box` `cpl_aim+ocn` `deep_anelastic` `dome` `exp2` `exp4` `fizhi-cs-32x32x40` `fizhi-cs-aqualev20` `fizhi-gridalt-hs` `front_relax` `global_oce_biogeo_bling` `global_oce_latlon` `global_ocean.90x40x15` `global_ocean.cs32x15` `hs94.128x64x5` `hs94.1x64x5` `hs94.cs-32x32x5` `ideal_2D_oce` `internal_wave` `inverted_barometer` `isomip` `lab_sea` `matrix_example` `obcs_ctrl` `offline_exf_seaice` `seaice_itd` `seaice_obcs` `shelfice_2d_remesh` `short_surf_wave` `so_box_biogeo` `solid-body.cs-32x32x1` `tutorial_advection_in_gyre` `tutorial_baroclinic_gyre` `tutorial_deep_convection` `tutorial_global_oce_biogeo` `tutorial_global_oce_in_p` `tutorial_global_oce_latlon` `tutorial_global_oce_optim` `tutorial_held_suarez_cs` `tutorial_plume_on_slope` `tutorial_reentrant_channel` `tutorial_rotating_tank` `tutorial_tracer_adjsens` `vermix`
