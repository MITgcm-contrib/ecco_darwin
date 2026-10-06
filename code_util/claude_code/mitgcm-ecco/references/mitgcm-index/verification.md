# Verification experiments

Source: `verification/` in the indexed tree. Grid = global nx×ny×Nr from `code/SIZE.h` (sNx·nSx·nPx etc.). "tests" counts reference outputs in `results/` (incl. adjoint/TLM/secondary).
Per-experiment details: `verification/<exp>.md`. Packages listed as written in packages.conf (groups like `gfd`/`oceanic` unexpanded; `-pkg` = excluded).

| experiment | grid | packages (code/packages.conf) | tests | description |
|---|---|---|---|---|
| `1D_ocean_ice_column` | 1x1x23 | gfd kpp exf seaice | 3 |  |
| `MLAdjust` | 50x26x40 | exch2 gfd gmredi flt diagnostics mnc | 7 | Simple set-up to test flow-dependent horizontal viscosity implementation. |
| `adjustment.128x64x1` | 128x64x1 |  | 1 |  |
| `adjustment.cs-32x32x1` | 768x8x1 | exch2 gfd -generic_advdiff diagnostics | 2 | Simple 1 layer, Barotropic adjustment on the Sphere, using the |
| `advect_cs` | 64x96x1 | exch2 gfd diagnostics | 1 |  |
| `advect_xy` | 20x20x1 |  | 2 |  |
| `advect_xz` | 20x1x20 | gfd -mom_common -mom_fluxform -mom_vecinv diagnostics | 3 |  |
| `aim.5l_Equatorial_Channel` | 128x23x5 | atmospheric aim_v23 | 1 | Five level intermediate atmospheric physics example. |
| `aim.5l_LatLon` | 128x64x5 | atmospheric aim_v23 zonal_filt | 1 | Intermediate Atmospheric physics, 5 layers Molteni Physics package. |
| `aim.5l_cs` | 192x32x5 | exch2 atmospheric -mom_fluxform aim_v23 land thsice diagnostics mnc | 2 |  |
| `atm_gray` | 192x32x26 | exch2 gfd shap_filt atm_phys diagnostics | 2 | Gray atmosphere physics example on Cubed-Sphere grid |
| `bottom_ctrl_5x5` | 5x5x4 | gfd cd_code adjoint | 4 |  |
| `cfc_example` | 128x64x15 | gfd cd_code gmredi ptracers gchem cfc layers diagnostics | 1 |  |
| `cheapAML_box` | 100x100x1 | gfd cheapaml diagnostics mnc | 2 |  |
| `cpl_aim+ocn` | 192x32x5 | exch2 atm_compon_interf compon_communic gfd -mom_fluxform shap_filt aim_v23 land thsice diagnostics | 0 | Atmosphere-Ocean coupled set-up example "cpl_aim+ocn" |
| `deep_anelastic` | 1x160x120 | gfd diagnostics | 2 |  |
| `dome` | 200x45x25 | gfd obcs diagnostics | 1 |  |
| `exp2` | 90x40x20 | gfd cd_code | 2 | Example: "4x4 Steady Global Simulation" |
| `exp4` | 80x42x8 | gfd obcs rbcs ptracers layers flt | 4 | Example: "Flow over a bump with Open Boundaries and passive tracers" |
| `fizhi-cs-32x32x40` | 192x32x40 | exch2 gfd -mom_fluxform shap_filt fizhi gridalt diagnostics | 1 |  |
| `fizhi-cs-aqualev20` | 192x32x20 | exch2 gfd -mom_fluxform shap_filt fizhi gridalt diagnostics mnc | 1 |  |
| `fizhi-gridalt-hs` | 192x32x10 | exch2 gfd -mom_fluxform shap_filt fizhi gridalt diagnostics | 1 |  |
| `front_relax` | 1x32x25 | gfd gmredi diagnostics | 5 | # Relaxation of a front in a channel : simplest example that uses GM-Redi parameterization |
| `global_oce_biogeo_bling` | 128x64x15 | obsfit cal gfd cd_code gmredi ptracers gchem bling diagnostics mnc | 6 |  |
| `global_oce_latlon` | 90x40x15 | gfd cd_code gmredi bbl ebm exf frazil profiles diagnostics | 13 | # Global Ocean Simulation at 4 degree Resolution, including Adjoint Set-Up |
| `global_ocean.90x40x15` | 90x40x15 | exch2 oceanic -kpp cd_code down_slope ggl90 ptracers sbo diagnostics | 11 | Example: "4x4 Global Simulation with Seasonal Forcing" |
| `global_ocean.cs32x15` | 384x16x15 | exch2 gfd gmredi ggl90 bulk_force exf -cal seaice thsice diagnostics mnc | 16 | global ocean using the cubed-sphere grid 32x32x32 with 15 levels |
| `halfpipe_streamice` | 40x20x1 | gfd -mom_common -mom_fluxform -mom_vecinv -generic_advdiff streamice diagnostics | 4 |  |
| `hs94.128x64x5` | 128x64x5 | gfd shap_filt zonal_filt | 1 |  |
| `hs94.1x64x5` | 1x64x5 | atmospheric diagnostics mnc mypackage | 3 | Held-Suarez zonal average config. - no eddy param |
| `hs94.cs-32x32x5` | 192x32x5 | exch2 atmospheric diagnostics | 2 |  |
| `ideal_2D_oce` | 1x56x15 | gfd cd_code gmredi diagnostics | 2 | test - experiment : ideal_2D_oce |
| `internal_wave` | 60x1x20 | gfd obcs kl10 mnc | 2 | Example: "Internal Wave Forced by Open Boundary" |
| `inverted_barometer` | 60x60x4 | gfd diagnostics mnc | 1 | Example: "4 Layer Double Gyre with Pressure Loading" |
| `isomip` | 50x100x30 | gfd ggl90 cd_code obcs shelfice steep_icecavity icefront diagnostics mnc | 10 |  |
| `lab_sea` | 20x16x23 | oceanic cd_code exf seaice salt_plume ptracers longstep diagnostics mnc | 12 | Labrador Sea Region with Sea-Ice |
| `matrix_example` | 32x32x1 | gfd matrix | 1 | Example to test pkg/matrix |
| `obcs_ctrl` | 64x64x8 | gfd -mom_fluxform obcs exf cal diagnostics ecco autodiff cost ctrl grdchk | 1 |  |
| `offline_exf_seaice` | 80x42x1 | gfd exf -cal seaice thsice diagnostics | 15 | Seaice-only verification experiment in idealized periodic channel |
| `seaice_itd` | 80x42x1 | gfd exf seaice diagnostics | 3 | Seaice-only verification experiment in idealized periodic channel with |
| `seaice_obcs` | 10x8x23 | oceanic obcs exf seaice salt_plume | 4 | Test set-up for seaice pkg with Open-Boundary Conditions |
| `shelfice_2d_remesh` | 1x200x90 | gfd obcs shelfice diagnostics | 1 | Simplified experiment to test pkg/shelfice vertical remeshing code |
| `short_surf_wave` | 52x1x50 | gfd diagnostics | 1 |  |
| `so_box_biogeo` | 42x20x15 | gfd cd_code obcs gmredi ptracers gchem dic diagnostics | 4 | Southern-Ocean box with Biochemistry, using Open-Boundary Conditions |
| `solid-body.cs-32x32x1` | 192x32x1 | exch2 gfd -mom_fluxform diagnostics | 1 | Simple solid-body rotation test on cubed-sphere grid |
| `tutorial_advection_in_gyre` | 60x60x1 | gfd ptracers diagnostics mnc | 1 |  |
| `tutorial_baroclinic_gyre` | 62x62x15 | gfd diagnostics mnc | 1 | Tutorial Example: "Baroclinic gyre" |
| `tutorial_barotropic_gyre` | 62x62x1 |  | 1 | Tutorial Example: "Barotropic gyre" |
| `tutorial_cfc_offline` | 128x64x15 | gfd -mom_common -mom_fluxform -mom_vecinv gmredi offline ptracers gchem cfc | 1 | Tutorial Example: "Offline CFC Experiments" |
| `tutorial_deep_convection` | 100x100x50 | gfd diagnostics | 2 | Tutorial Example: "Surface Driven (Deep) Convection" |
| `tutorial_dic_adjoffline` | 128x64x15 | gfd -mom_common -mom_fluxform -mom_vecinv gmredi offline ptracers gchem dic mnc autodiff cost ctrl grdchk | 2 |  |
| `tutorial_global_oce_biogeo` | 128x64x15 | gfd cd_code gmredi ptracers gchem dic diagnostics mnc | 5 | Tutorial Example: "Biochemistry Tutorial" |
| `tutorial_global_oce_in_p` | 90x40x15 | gfd -mom_vecinv | 1 | Tutorial Example: "P coordinate Global Ocean" |
| `tutorial_global_oce_latlon` | 90x40x15 | gfd cd_code gmredi ptracers mnc | 1 | Tutorial Example: "Global ocean" |
| `tutorial_global_oce_optim` | 90x40x15 | gfd cd_code gmredi autodiff cost ctrl grdchk | 1 | Tutorial Example: "Global Ocean State Estimation" |
| `tutorial_held_suarez_cs` | 192x32x20 | exch2 gfd shap_filt diagnostics mnc | 1 |  |
| `tutorial_plume_on_slope` | 320x1x60 | gfd obcs | 2 | Tutorial Example: "Gravity plume on a continental slope" |
| `tutorial_reentrant_channel` | 20x40x49 | gfd gmredi rbcs layers diagnostics | 1 | Tutorial Example: "Reentrant channel" |
| `tutorial_rotating_tank` | 120x23x29 | gfd diagnostics mnc | 1 |  |
| `tutorial_tracer_adjsens` | 90x40x20 | gfd cd_code gmredi kpp ptracers autodiff cost ctrl grdchk | 6 | Tutorial Example: "Centennial Time Scale Tracer Injection" |
| `vermix` | 1x1x26 | gfd kpp pp81 my82 ggl90 opps diagnostics mnc | 7 |  |

## Experiments in the user's branches (new or modified vs upstream)

`name@label` — label = branch clone (see README.md).

| experiment | grid | packages | tests | description |
|---|---|---|---|---|
| `global_oce_latlon@bbl` | 90x40x15 | gfd cd_code gmredi bbl ebm exf frazil profiles diagnostics | 14 | # Global Ocean Simulation at 4 degree Resolution, including Adjoint Set-Up |
| `wad_balzano@wad` | 140x1x1 | gfd obcs wad diagnostics | 4 |  |
| `wad_estuary_3d@wad` | 40x30x5 | gfd obcs kpp gmredi ptracers diagnostics exf cal seaice wad | 4 |  |
| `wad_flat_xz@wad` | 100x1x10 | gfd obcs kpp wad | 2 |  |
| `wad_iceground@wad` | 60x10x5 | gfd diagnostics exf cal seaice wad | 2 |  |
| `wad_mangrove@wad` | 100x50x5 | gfd obcs kpp ptracers diagnostics wad sediment mangrove | 6 |  |
| `wad_mudflat@wad` | 100x60x5 | gfd obcs kpp ptracers rbcs diagnostics exf cal seaice wad sediment | 6 |  |
| `wad_overflood@wad` | 80x40x5 | gfd diagnostics exf cal seaice wad overflood | 3 |  |
| `wad_thacker_1d@wad` | 250x1x1 | gfd wad diagnostics | 3 |  |
| `wad_balzano@wadcheckin` | 140x1x1 | gfd obcs wad diagnostics | 4 | # wad_balzano: Balzano (1998) tidal-flat tests: slope, step, pool |
| `wad_estuary_3d@wadcheckin` | 40x30x5 | gfd obcs kpp gmredi ptracers diagnostics wad | 4 | # wad_estuary_3d: 3-D tidal estuary with drying flats |
| `wad_flat_xz@wadcheckin` | 100x1x10 | gfd obcs kpp wad | 2 | # wad_flat_xz: Stratified beach at rest, and under a tide |
| `wad_mudflat@wadcheckin` | 100x60x5 | gfd obcs kpp ptracers rbcs diagnostics wad | 3 | # wad_mudflat: Macrotidal mudflat with tidal creeks and a river |
| `wad_thacker_1d@wadcheckin` | 250x1x1 | gfd wad diagnostics | 3 | # wad_thacker_1d: Thacker (1981) oscillating basin |
