# pkg/openad

OpenAD automatic-differentiation support (legacy).

**runtime switch:** `useOPENAD`-style flag in `data.pkg` (check exact name in packages_boot.F)
**adjoint support files:** openad_ad_diff.list, openad_checkpointInit.F

## CPP options (defaults as shipped)
- `ALLOW_OPENAD_ACTIVE_READ_XY` (undef, OPENAD_OPTIONS.h)
- `ALLOW_OPENAD_ACTIVE_READ_XYZ` (undef, OPENAD_OPTIONS.h)
- `ALLOW_OPENAD_ACTIVE_WRITE` (undef, OPENAD_OPTIONS.h)
- `ALLOW_OPENAD_DIVA` (undef, OPENAD_OPTIONS.h)
- `OAD_DEBUG` (undef, OPENAD_OPTIONS.h)

## Headers
- `OPENAD_OPTIONS.h` — BOP

## Routines (11)
`ad_s_different_multiple.F`, `ad_s_ifnblnk.F`, `ad_s_ilnblnk.F`, `ad_s_master_cpu_thread.F`, `ad_s_mds_reclen.F`, `externalDummies.F`, `inner_do_loop.F`, `openad_checkpointInit.F`, `openad_dumpAdjoint.F`

## Called from outside the package
- `EXCH1_RL` ← `eesupp/src/exch_tap_b.F:328`
- `EXCH1_RL` ← `eesupp/src/exch_tap_d.F:197,199`
- `GLOBAL_SUM_TILE_RL` ← `model/src/calc_wsurf_tr.F:91,92`
- `GLOBAL_SUM_TILE_RL` ← `model/src/cg2d.F:185,187,243,295`
- `GLOBAL_SUM_TILE_RL` ← `model/src/cg2d_nsa.F:216,217,277,330`
- `GLOBAL_SUM_TILE_RL` ← `model/src/cg2d_sr.F:198,199,246,273`
- `GLOBAL_SUM_TILE_RL` ← `model/src/cg3d.F:239,240,330,470`
- `GLOBAL_SUM_TILE_RL` ← `model/src/forcing_surf_relax.F:180,216`
- `GLOBAL_SUM_TILE_RL` ← `model/src/ini_global_domain.F:88,89,110,111`
- `INNER_DO_LOOP` ← `model/src/main_do_loop.F:228`
- `GLOBAL_SUM_TILE_RL` ← `model/src/remove_mean.F:101,102,258,259`
- `GLOBAL_SUM_TILE_RL` ← `model/src/solve_for_pressure.F:164`
- `GLOBAL_SUM_TILE_RL` ← `pkg/aim_v23/aim_do_co2.F:91`
- `GLOBAL_SUM_TILE_RL` ← `pkg/cost/cost_final.F:215`
- `GLOBAL_SUM_TILE_RL` ← `pkg/ctrl/ctrl_map_gentim2d.F:117`
- `GLOBAL_SUM_TILE_RL` ← `pkg/debug/debug_fld_stats_rl.F:90,91,123`
- `GLOBAL_SUM_TILE_RL` ← `pkg/debug/debug_fld_stats_rs.F:90,91,123`
- `GLOBAL_SUM_TILE_RL` ← `pkg/diagnostics/diag_cg2d.F:197,199,250,287`
- `GLOBAL_SUM_TILE_RL` ← `pkg/dic/dic_atmos.F:124,125`
- `GLOBAL_SUM_TILE_RL` ← `pkg/dic/dic_cost.F:54`
- `GLOBAL_SUM_TILE_RL` ← `pkg/dic/dic_ini_atmos.F:121`
- `GLOBAL_SUM_TILE_RL` ← `pkg/ebm/ebm_area_t.F:106,107,108,109`
- `GLOBAL_SUM_TILE_RL` ← `pkg/ebm/ebm_zonalmean.F:71,72`
- `GLOBAL_SUM_TILE_RL` ← `pkg/ecco/cost_gencost_boxmean.F:152,164`
- `GLOBAL_SUM_TILE_RL` ← `pkg/ecco/cost_gencost_moc.F:222`
- `GLOBAL_SUM_TILE_RL` ← `pkg/ecco/cost_gencost_transp.F:244,245`
- `GLOBAL_SUM_TILE_RL` ← `pkg/ecco/ecco_cost_final.F:150`
- `GLOBAL_SUM_TILE_RL` ← `pkg/ecco/ecco_phys.F:155,156,198,199`
- `GLOBAL_SUM_TILE_RL` ← `pkg/ecco/ecco_toolbox.F:772,773`
- `GLOBAL_SUM_TILE_RL` ← `pkg/gchem/gchem_surfmean.F:52`
- `GLOBAL_SUM_TILE_RL` ← `pkg/monitor/mon_calc_stats_rl.F:111,112,113,114`
- `GLOBAL_SUM_TILE_RL` ← `pkg/monitor/mon_calc_stats_rs.F:111,112,113,114`
- `GLOBAL_SUM_TILE_RL` ← `pkg/monitor/mon_ke.F:155,156,157,158`
- `GLOBAL_SUM_TILE_RL` ← `pkg/monitor/mon_stats_latbnd_rl.F:122,123,124`
- `GLOBAL_SUM_TILE_RL` ← `pkg/monitor/mon_stats_rl.F:104,105,106,107`
- `GLOBAL_SUM_TILE_RL` ← `pkg/monitor/mon_stats_rs.F:77,78,113`
- `GLOBAL_SUM_TILE_RL` ← `pkg/monitor/mon_surfcor.F:174,175,176,178`
- `GLOBAL_SUM_TILE_RL` ← `pkg/monitor/mon_vort3.F:327,328,329,330`
- `GLOBAL_SUM_TILE_RL` ← `pkg/obcs/obcs_balance_flow.F:138,140,187,189`
- `GLOBAL_SUM_TILE_RL` ← `pkg/obcs/obcs_diag_balance.F:191,193`
- `GLOBAL_SUM_TILE_RL` ← `pkg/obcs/obcs_mon_stats.F:140,143,144,213`
- `GLOBAL_SUM_TILE_RL` ← `pkg/profiles/profiles_cost.F:609,611,639,641`
- `GLOBAL_SUM_TILE_RL` ← `pkg/ptracers/ptracers_calc_wsurf_tr.F:77`
- `GLOBAL_SUM_TILE_RL` ← `pkg/sbo/sbo_calc.F:241,242,243,424`
- `GLOBAL_SUM_TILE_RL` ← `pkg/seaice/seaice_cost_final.F:76,77`
- `GLOBAL_SUM_TILE_RL` ← `pkg/seaice/seaice_evp.F:732,1006`
- `GLOBAL_SUM_TILE_RL` ← `pkg/seaice/seaice_fgmres.F:669`
- `GLOBAL_SUM_TILE_RL` ← `pkg/seaice/seaice_growth.F:2583,2585,2592`
- `GLOBAL_SUM_TILE_RL` ← `pkg/shelfice/shelfice_cost_final.F:92`
- `GLOBAL_SUM_TILE_RL` ← `pkg/shelfice/shelfice_step_icemass.F:228`
- `STREAMICE_SMOOTH_ADJOINT_FIELD` ← `pkg/streamice/streamice_advect_thickness.F:73`
- `GLOBAL_SUM_TILE_RL` ← `pkg/streamice/streamice_cg_solve.F:199,355,356,422`
- `GLOBAL_SUM_TILE_RL` ← `pkg/streamice/streamice_cg_solve_matfree.F:141,281,282,348`
- `GLOBAL_SUM_TILE_RL` ← `pkg/streamice/streamice_cost_final.F:314,316,318,320`
- `GLOBAL_SUM_TILE_RL` ← `pkg/streamice/streamice_get_vel_fp_err.F:120`
- `GLOBAL_SUM_TILE_RL` ← `pkg/streamice/streamice_get_vel_resid_err.F:131`
- `GLOBAL_SUM_TILE_RL` ← `pkg/streamice/streamice_tap_fixedpoint_notreduced.F:42`
- `GLOBAL_SUM_TILE_RL` ← `pkg/streamice/streamice_timestep.F:164`
- `GLOBAL_SUM_TILE_RL` ← `pkg/tapenade/stubs_tap_adj.F:75`
- `GLOBAL_SUM_TILE_RL` ← `pkg/tapenade/stubs_tap_tlm.F:15,16`
- `GLOBAL_SUM_TILE_RL` ← `pkg/thsice/thsice_balance_frw.F:90,91`
