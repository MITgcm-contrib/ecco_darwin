# pkg/mom_vecinv

Vector-invariant momentum equations (used by ECCO/cubed-sphere/LLC configs; enstrophy/energy-conserving Coriolis).

**pkg_depend:** +mom_common  (`+` requires, `-` excludes)
**in groups:** gfd
**runtime switch:** `useMOM_VECINV`-style flag in `data.pkg` (check exact name in packages_boot.F)
**manual:** `doc/algorithm/algorithm.rst`, `doc/examples/held_suarez_cs/held_suarez_cs.rst`, `doc/getting_started/getting_started.rst`, `doc/outp_pkgs/outp_pkgs.rst`
**adjoint support files:** mom_vecinv_ad_diff.list

## CPP options (defaults as shipped)
- `MOM_VI_ORIGINAL_VISCA4` (undef, MOM_VECINV_OPTIONS.h) — use the original discretization (not recommended) for biharmonic viscosity that was in mom_vi_hdissip.F, version 1.1.2.1

## Headers
- `MOM_VECINV_OPTIONS.h` — CPP options file for mom_vecinv package Use this file for selecting CPP options within the mom_vecinv package

## Routines (12)
`mom_vecinv.F`, `mom_vi_coriolis.F`, `mom_vi_del2uv.F`, `mom_vi_hdissip.F`, `mom_vi_u_coriolis.F`, `mom_vi_u_coriolis_c4.F`, `mom_vi_u_grad_ke.F`, `mom_vi_u_vertshear.F`, `mom_vi_v_coriolis.F`, `mom_vi_v_coriolis_c4.F`, `mom_vi_v_grad_ke.F`, `mom_vi_v_vertshear.F`

## Called from outside the package
- `MOM_VECINV` ← `model/src/dynamics.F:527`
- `MOM_VI_U_CORIOLIS` ← `pkg/seaice/seaice_mom_advection.F:135`
- `MOM_VI_U_CORIOLIS_C4` ← `pkg/seaice/seaice_mom_advection.F:129`
- `MOM_VI_U_GRAD_KE` ← `pkg/seaice/seaice_mom_advection.F:171`
- `MOM_VI_V_CORIOLIS` ← `pkg/seaice/seaice_mom_advection.F:152`
- `MOM_VI_V_CORIOLIS_C4` ← `pkg/seaice/seaice_mom_advection.F:146`
- `MOM_VI_V_GRAD_KE` ← `pkg/seaice/seaice_mom_advection.F:177`
- `MOM_VI_DEL2UV` ← `pkg/shap_filt/shap_filt_uv_s2.F:178`

## Verification experiments compiling it (53)
`1D_ocean_ice_column` `MLAdjust` `adjustment.cs-32x32x1` `advect_cs` `aim.5l_Equatorial_Channel` `aim.5l_LatLon` `aim.5l_cs` `atm_gray` `bottom_ctrl_5x5` `cfc_example` `cheapAML_box` `cpl_aim+ocn` `deep_anelastic` `dome` `exp2` `exp4` `fizhi-cs-32x32x40` `fizhi-cs-aqualev20` `fizhi-gridalt-hs` `front_relax` `global_oce_biogeo_bling` `global_oce_latlon` `global_ocean.90x40x15` `global_ocean.cs32x15` `hs94.128x64x5` `hs94.1x64x5` `hs94.cs-32x32x5` `ideal_2D_oce` `internal_wave` `inverted_barometer` `isomip` `lab_sea` `matrix_example` `obcs_ctrl` `offline_exf_seaice` `seaice_itd` `seaice_obcs` `shelfice_2d_remesh` `short_surf_wave` `so_box_biogeo` `solid-body.cs-32x32x1` `tutorial_advection_in_gyre` `tutorial_baroclinic_gyre` `tutorial_deep_convection` `tutorial_global_oce_biogeo` `tutorial_global_oce_latlon` `tutorial_global_oce_optim` `tutorial_held_suarez_cs` `tutorial_plume_on_slope` `tutorial_reentrant_channel` `tutorial_rotating_tank` `tutorial_tracer_adjsens` `vermix`
