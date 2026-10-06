# pkg/monitor

Monitor statistics (%MON lines in STDOUT: CFL, KE, field min/max/mean/sd) used by testreport.

**in groups:** gfd
**always-on utility package** (no data.pkg switch)
**manual:** `doc/contributing/contributing.rst`, `doc/examples/baroclinic_gyre/baroclinic_gyre.rst`, `doc/getting_started/getting_started.rst`, `doc/outp_pkgs/outp_pkgs.rst`

## CPP options (defaults as shipped)
- `MONITOR_TEST_HFACZ` (undef, MONITOR_OPTIONS.h) — Disable use of hFacZ

## Headers
- `MONITOR.h` — BOP
- `MONITOR_OPTIONS.h` — CPP options file for monitor package Use this file for selecting options within the monitor package

## Routines (29)
`mon_advcfl.F`, `mon_advcflw.F`, `mon_advcflw2.F`, `mon_calc_advcfl.F`, `mon_calc_stats_rl.F`, `mon_calc_stats_rs.F`, `mon_init.F`, `mon_ke.F`, `mon_out.F`, `mon_printstats_rl.F`, `mon_printstats_rs.F`, `mon_set_iounit.F`, `mon_set_pref.F`, `mon_solution.F`, `mon_stats_latbnd_rl.F`, `mon_stats_rl.F`, `mon_stats_rs.F`, `mon_surfcor.F`, `mon_vort3.F`, `mon_writestats_rl.F`, `mon_writestats_rs.F`, `monitor.F`, `monitor_ad.F`, `monitor_g.F`

## Called from outside the package
- `MONITOR` ← `model/src/forward_step.F:1154`
- `MON_PRINTSTATS_RS` ← `model/src/ini_cori.F:208,209,210`
- `MON_SET_PREF` ← `model/src/ini_cori.F:207`
- `MON_PRINTSTATS_RS` ← `model/src/ini_forcing.F:208`
- `MON_PRINTSTATS_RS` ← `model/src/ini_grid.F:209,210,211,212`
- `MON_INIT` ← `model/src/ini_model_io.F:249`
- `MONITOR` ← `model/src/initialise_varia.F:383`
- `MON_CALC_ADVCFL_GLOB` ← `model/src/thermodynamics.F:388`
- `MON_CALC_ADVCFL_TILE` ← `model/src/thermodynamics.F:279`
- `MON_OUT_I` ← `pkg/exf/exf_monitor.F:108`
- `MON_OUT_RL` ← `pkg/exf/exf_monitor.F:109`
- `MON_SET_PREF` ← `pkg/exf/exf_monitor.F:107`
- `MON_WRITESTATS_RL` ← `pkg/exf/exf_monitor.F:113,115,118,120`
- `MON_OUT_I` ← `pkg/exf/exf_monitor_ad.F:114`
- `MON_OUT_RL` ← `pkg/exf/exf_monitor_ad.F:115`
- `MON_SET_PREF` ← `pkg/exf/exf_monitor_ad.F:113`
- `MON_WRITESTATS_RL` ← `pkg/exf/exf_monitor_ad.F:121,123,126,128`
- `MON_WRITESTATS_RS` ← `pkg/exf/exf_monitor_ad.F:203,205,207,209`
- `MON_OUT_RL` ← `pkg/land/land_monitor.F:122,144,145,146`
- `MON_SET_PREF` ← `pkg/land/land_monitor.F:121`
- `MON_STATS_LATBND_RL` ← `pkg/land/land_monitor.F:129,165,193,221`
- `MON_OUT_I` ← `pkg/obcs/obcs_monitor.F:95`
- `MON_OUT_RL` ← `pkg/obcs/obcs_monitor.F:96,258,259,260`
- `MON_SET_PREF` ← `pkg/obcs/obcs_monitor.F:94,99`
- `MON_OUT_I` ← `pkg/ptracers/ptracers_monitor.F:103`
- `MON_OUT_RL` ← `pkg/ptracers/ptracers_monitor.F:104`
- `MON_SET_PREF` ← `pkg/ptracers/ptracers_monitor.F:102,107`
- `MON_WRITESTATS_RL` ← `pkg/ptracers/ptracers_monitor.F:111`
- `MON_SET_PREF` ← `pkg/ptracers/ptracers_monitor_ad.F:105`
- `MON_WRITESTATS_RL` ← `pkg/ptracers/ptracers_monitor_ad.F:108`
- `MON_OUT_I` ← `pkg/seaice/seaice_monitor.F:99`
- `MON_OUT_RL` ← `pkg/seaice/seaice_monitor.F:100`
- `MON_SET_PREF` ← `pkg/seaice/seaice_monitor.F:98`
- `MON_WRITESTATS_RL` ← `pkg/seaice/seaice_monitor.F:104,106,110,112`
- `MON_OUT_I` ← `pkg/seaice/seaice_monitor_ad.F:105`
- `MON_OUT_RL` ← `pkg/seaice/seaice_monitor_ad.F:106`
- `MON_SET_PREF` ← `pkg/seaice/seaice_monitor_ad.F:104`
- `MON_WRITESTATS_RL` ← `pkg/seaice/seaice_monitor_ad.F:110,112,116,118`
- `MON_CALC_STATS_RL` ← `pkg/thsice/thsice_monitor.F:245,249`
- `MON_OUT_RL` ← `pkg/thsice/thsice_monitor.F:124,137,138,139`
- `MON_SET_PREF` ← `pkg/thsice/thsice_monitor.F:123`
- `MON_STATS_LATBND_RL` ← `pkg/thsice/thsice_monitor.F:127,149,168,200`

## Verification experiments compiling it (58)
`1D_ocean_ice_column` `MLAdjust` `adjustment.cs-32x32x1` `advect_cs` `advect_xz` `aim.5l_Equatorial_Channel` `aim.5l_LatLon` `aim.5l_cs` `atm_gray` `bottom_ctrl_5x5` `cfc_example` `cheapAML_box` `cpl_aim+ocn` `deep_anelastic` `dome` `exp2` `exp4` `fizhi-cs-32x32x40` `fizhi-cs-aqualev20` `fizhi-gridalt-hs` `front_relax` `global_oce_biogeo_bling` `global_oce_latlon` `global_ocean.90x40x15` `global_ocean.cs32x15` `halfpipe_streamice` `hs94.128x64x5` `hs94.1x64x5` `hs94.cs-32x32x5` `ideal_2D_oce` `internal_wave` `inverted_barometer` `isomip` `lab_sea` `matrix_example` `obcs_ctrl` `offline_exf_seaice` `seaice_itd` `seaice_obcs` `shelfice_2d_remesh` `short_surf_wave` `so_box_biogeo` `solid-body.cs-32x32x1` `tutorial_advection_in_gyre` `tutorial_baroclinic_gyre` `tutorial_cfc_offline` `tutorial_deep_convection` `tutorial_dic_adjoffline` `tutorial_global_oce_biogeo` `tutorial_global_oce_in_p` `tutorial_global_oce_latlon` `tutorial_global_oce_optim` `tutorial_held_suarez_cs` `tutorial_plume_on_slope` `tutorial_reentrant_channel` `tutorial_rotating_tank` `tutorial_tracer_adjsens` `vermix`
