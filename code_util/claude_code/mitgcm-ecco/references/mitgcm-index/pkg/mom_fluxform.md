# pkg/mom_fluxform

Flux-form momentum equations (default unless vectorInvariantMomentum=.TRUE.).

**pkg_depend:** +mom_common  (`+` requires, `-` excludes)
**in groups:** gfd
**runtime switch:** `useMOM_FLUXFORM`-style flag in `data.pkg` (check exact name in packages_boot.F)
**manual:** `doc/algorithm/algorithm.rst`, `doc/getting_started/getting_started.rst`, `doc/phys_pkgs/mom_packages.rst`
**adjoint support files:** mom_fluxform_ad_diff.list

## CPP options (defaults as shipped)
- `MOM_BOUNDARY_CONSERVE` (undef, MOM_FLUXFORM_OPTIONS.h) — A trick to conserve U,V momemtum next to a step (vertical plane) or a coastline edge (horizontal plane).

## Headers
- `MOM_FLUXFORM.h` — BOP
- `MOM_FLUXFORM_OPTIONS.h` — CPP options file for mom_fluxform package Use this file for selecting CPP options within the mom_fluxform package

## Routines (21)
`mom_calc_rtrans.F`, `mom_fluxform.F`, `mom_u_adv_uu.F`, `mom_u_adv_vu.F`, `mom_u_adv_wu.F`, `mom_u_coriolis.F`, `mom_u_del2u.F`, `mom_u_metric_cylinder.F`, `mom_u_metric_sphere.F`, `mom_u_xviscflux.F`, `mom_u_yviscflux.F`, `mom_uv_boundary.F`, `mom_v_adv_uv.F`, `mom_v_adv_vv.F`, `mom_v_adv_wv.F`, `mom_v_coriolis.F`, `mom_v_del2v.F`, `mom_v_metric_cylinder.F`, `mom_v_metric_sphere.F`, `mom_v_xviscflux.F`, `mom_v_yviscflux.F`

## Called from outside the package
- `MOM_FLUXFORM` ← `model/src/dynamics.F:517`

## Verification experiments compiling it (48)
`1D_ocean_ice_column` `MLAdjust` `adjustment.cs-32x32x1` `advect_cs` `aim.5l_Equatorial_Channel` `aim.5l_LatLon` `atm_gray` `bottom_ctrl_5x5` `cfc_example` `cheapAML_box` `cpl_aim+ocn` `deep_anelastic` `dome` `exp2` `exp4` `front_relax` `global_oce_biogeo_bling` `global_oce_latlon` `global_ocean.90x40x15` `global_ocean.cs32x15` `hs94.128x64x5` `hs94.1x64x5` `hs94.cs-32x32x5` `ideal_2D_oce` `internal_wave` `inverted_barometer` `isomip` `lab_sea` `matrix_example` `offline_exf_seaice` `seaice_itd` `seaice_obcs` `shelfice_2d_remesh` `short_surf_wave` `so_box_biogeo` `tutorial_advection_in_gyre` `tutorial_baroclinic_gyre` `tutorial_deep_convection` `tutorial_global_oce_biogeo` `tutorial_global_oce_in_p` `tutorial_global_oce_latlon` `tutorial_global_oce_optim` `tutorial_held_suarez_cs` `tutorial_plume_on_slope` `tutorial_reentrant_channel` `tutorial_rotating_tank` `tutorial_tracer_adjsens` `vermix`
