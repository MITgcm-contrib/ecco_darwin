# pkg/obcs

Open boundary conditions: prescribed/relaxed boundary values (T,S,U,V,eta,seaice, ptracers), sponges, Orlanski radiation, Stevens and Flather schemes, tidal forcing, balancing.

**runtime switch:** `useOBCS`-style flag in `data.pkg` (check exact name in packages_boot.F)
**reads:** `data.obcs`
**manual:** `doc/phys_pkgs/obcs.rst`, `doc/examples/plume_on_slope/plume_on_slope.rst`, `doc/ocean_state_est/ocean_state_est.rst`, `doc/outp_pkgs/outp_pkgs.rst`
**adjoint support files:** obcs_ad_check_lev1_dir.h, obcs_ad_check_lev2_dir.h, obcs_ad_check_lev3_dir.h, obcs_ad_check_lev4_dir.h, obcs_ad_diff.list

## Namelist parameters
### OBCS_PARM01
- `insideOBmaskFile` — File to specify Inside OB region mask (zero beyond OB)
- `OBNconnectFile`
- `OBSconnectFile`
- `OBEconnectFile`
- `OBWconnectFile`
- `OB_Jnorth`
- `OB_Jsouth`
- `OB_Ieast`
- `OB_Iwest`
- `OB_singleJnorth`
- `OB_singleJsouth`
- `OB_singleIeast`
- `OB_singleIwest`
- `useOrlanskiNorth`
- `useOrlanskiSouth`
- `useOrlanskiEast`
- `useOrlanskiWest`
- `useStevensNorth`
- `useStevensSouth`
- `useStevensEast`
- `useStevensWest`
- `useOBCSprescribe` — read boundary conditions from a file (overrides Orlanski and other boundary values)
- `useOBCStides` — modify OB velocity by adding barotropic tidal component
- `OBCS_tidalPeriod` — tidal period (s) for each tidal component
- `OBCS_u1_adv_T` — >0: use 1rst O. upwind adv-scheme @ OB (=1: only if outflow)
- `OBCS_u1_adv_S` — >0: use 1rst O. upwind adv-scheme @ OB (=1: only if outflow)
- `OBNuFile`
- `OBNvFile`
- `OBNtFile`
- `OBNsFile`
- `OBNetaFile`
- `OBSuFile`
- `OBSvFile`
- `OBStFile`
- `OBSsFile`
- `OBSetaFile`
- `OBEuFile`
- `OBEvFile`
- `OBEtFile`
- `OBEsFile`
- `OBEetaFile`
- `OBWuFile`
- `OBWvFile`
- `OBWtFile`
- `OBWsFile`
- `OBWetaFile`
- `OBNwFile`
- `OBSwFile`
- `OBEwFile`
- `OBWwFile`
- `OBN_uTidAmFile`
- `OBS_uTidAmFile`
- `OBE_uTidAmFile`
- `OBW_uTidAmFile`
- `OBN_vTidAmFile`
- `OBS_vTidAmFile`
- `OBE_vTidAmFile`
- `OBW_vTidAmFile`
- `OBN_uTidPhFile`
- `OBS_uTidPhFile`
- `OBE_uTidPhFile`
- `OBW_uTidPhFile`
- `OBN_vTidPhFile`
- `OBS_vTidPhFile`
- `OBE_vTidPhFile`
- `OBW_vTidPhFile`
- `OBNaFile`
- `OBSaFile`
- `OBEaFile`
- `OBWaFile`
- `OBNhFile`
- `OBShFile`
- `OBEhFile`
- `OBWhFile`
- `OBNslFile`
- `OBSslFile`
- `OBEslFile`
- `OBWslFile`
- `OBNsnFile`
- `OBSsnFile`
- `OBEsnFile`
- `OBWsnFile`
- `OBNuiceFile`
- `OBSuiceFile`
- `OBEuiceFile`
- `OBWuiceFile`
- `OBNviceFile`
- `OBSviceFile`
- `OBEviceFile`
- `OBWviceFile`
- `OBCS_u1_adv_Tr` — >0: use 1rst O. upwind adv-scheme @ OB (=1: only if outflow)  _[ifdef ALLOW_PTRACERS]_
- `OBNptrFile`  _[ifdef ALLOW_PTRACERS]_
- `OBSptrFile`  _[ifdef ALLOW_PTRACERS]_
- `OBEptrFile`  _[ifdef ALLOW_PTRACERS]_
- `OBWptrFile`  _[ifdef ALLOW_PTRACERS]_
- `useOBCSsponge` — turns on sponge layer along boundaries (def=false)
- `useSeaiceSponge` — turns on seaice sponge layer along boundary (def=false)
- `useSeaiceNeumann` — use Neumann conditions for sea ice variables (def=false)
- `OBCSsponge_N` — turns on sponge layer along North boundary (def=true)
- `OBCSsponge_S` — turns on sponge layer along South boundary (def=true)
- `OBCSsponge_E` — turns on sponge layer along East boundary (def=true)
- `OBCSsponge_W` — turns on sponge layer along West boundary (def=true)
- `OBCSsponge_UatNS` — turns on uVel sponge at North/South boundaries (def=true)
- `OBCSsponge_UatEW` — turns on uVel sponge at East/West boundaries (def=true)
- `OBCSsponge_VatNS` — turns on vVel sponge at North/South boundaries (def=true)
- `OBCSsponge_VatEW` — turns on vVel sponge at East/West boundaries (def=true)
- `OBCSsponge_Theta` — turns on Theta sponge along boundaries (def=true)
- `OBCSsponge_Salt` — turns on Salt sponge along boundaries (def=true)
- `useLinearSponge` — use linear instead of exponential sponge (def=false)
- `useOBCSbalance` — balance the volume flux through boundary at every time step
- `OBCSbalanceSurf` — also include surface flux of mass into balance
- `OBCS_balanceFacN`
- `OBCS_balanceFacS`
- `OBCS_balanceFacE`
- `OBCS_balanceFacW`
- `OBCSfixTopo` — check and adjust topography for problematic gradients across boundaries (def=true)
- `OBCS_uvApplyFac` — multiplying factor to U,V normal comp. when applying OBC to 2nd column/row (for backward compatibility).
- `OBCS_monitorFreq` — monitor output frequency (s) for OB statistics
- `OBCS_monSelect` — select group of variables to monitor
- `OBCSprintDiags` — print boundary values to STDOUT (def=true)
- `useOBCSYearlyFields`
- `OBNAmFile`
- `OBSAmFile`
- `OBEAmFile`
- `OBWAmFile`
- `OBNPhFile`
- `OBSPhFile`
- `OBEPhFile`
- `OBWPhFile`
- `tidalPeriod`
### OBCS_PARM02
- `CMAX`  _[ifdef ALLOW_ORLANSKI]_
- `cvelTimeScale`  _[ifdef ALLOW_ORLANSKI]_
- `CFIX`  _[ifdef ALLOW_ORLANSKI]_
- `useFixedCEast`  _[ifdef ALLOW_ORLANSKI]_
- `useFixedCWest`  _[ifdef ALLOW_ORLANSKI]_
### OBCS_PARM03
- `Urelaxobcsinner`  _[ifdef ALLOW_OBCS_SPONGE]_
- `Urelaxobcsbound`  _[ifdef ALLOW_OBCS_SPONGE]_
- `Vrelaxobcsinner`  _[ifdef ALLOW_OBCS_SPONGE]_
- `Vrelaxobcsbound`  _[ifdef ALLOW_OBCS_SPONGE]_
- `spongeThickness` — number grid points that make up the sponge layer (def=0)  _[ifdef ALLOW_OBCS_SPONGE]_
### OBCS_PARM04
- `TrelaxStevens`  _[ifdef ALLOW_OBCS_STEVENS]_
- `SrelaxStevens`  _[ifdef ALLOW_OBCS_STEVENS]_
- `useStevensPhaseVel`  _[ifdef ALLOW_OBCS_STEVENS]_
- `useStevensAdvection`  _[ifdef ALLOW_OBCS_STEVENS]_
### OBCS_PARM05
- `Arelaxobcsinner`  _[ifdef ALLOW_OBCS_SEAICE_SPONGE]_
- `Arelaxobcsbound`  _[ifdef ALLOW_OBCS_SEAICE_SPONGE]_
- `Hrelaxobcsinner`  _[ifdef ALLOW_OBCS_SEAICE_SPONGE]_
- `Hrelaxobcsbound`  _[ifdef ALLOW_OBCS_SEAICE_SPONGE]_
- `SLrelaxobcsinner`  _[ifdef ALLOW_OBCS_SEAICE_SPONGE]_
- `SLrelaxobcsbound`  _[ifdef ALLOW_OBCS_SEAICE_SPONGE]_
- `SNrelaxobcsinner`  _[ifdef ALLOW_OBCS_SEAICE_SPONGE]_
- `SNrelaxobcsbound`  _[ifdef ALLOW_OBCS_SEAICE_SPONGE]_
- `seaiceSpongeThickness` — number grid points that make up the sponge layer (def=0)  _[ifdef ALLOW_OBCS_SEAICE_SPONGE]_

## CPP options (defaults as shipped)
- `ALLOW_OBCS_NORTH` (define, OBCS_OPTIONS.h) — Enable individual open boundaries
- `ALLOW_OBCS_SOUTH` (define, OBCS_OPTIONS.h)
- `ALLOW_OBCS_EAST` (define, OBCS_OPTIONS.h)
- `ALLOW_OBCS_WEST` (define, OBCS_OPTIONS.h)
- `ALLOW_ORLANSKI` (define, OBCS_OPTIONS.h) — This include hooks to the Orlanski Open Boundary Radiation code
- `ALLOW_OBCS_PRESCRIBE` (define, OBCS_OPTIONS.h) — Enable OB values to be prescribed via external fields that are read from a file
- `ALLOW_OBCS_STEVENS` (undef, OBCS_OPTIONS.h) — Enable OB conditions following Stevens (1990)
- `ALLOW_OBCS_SPONGE` (undef, OBCS_OPTIONS.h) — Allow sponge layer treatment of open boundary conditions
- `ALLOW_OBCS_SEAICE_SPONGE` (undef, OBCS_OPTIONS.h) — Include hooks to sponge layer treatment of pkg/seaice variables
- `ALLOW_OBCS_BALANCE` (define, OBCS_OPTIONS.h) — balance barotropic velocity
- `ALLOW_OBCS_TIDES` (undef, OBCS_OPTIONS.h) — Allow to add barotropic tidal contributions to OB velocity
- `OBCS_UVICE_OLD` (undef, OBCS_OPTIONS.h) — Use older implementation of obcs in seaice-dynamics note: most of the "experimental" options listed below have not yet been implementated in new version.
- `OBCS_SEAICE_AVOID_CONVERGENCE` (undef, OBCS_OPTIONS.h) — Ice convergence at edges can cause model to blow up.  The following CPP option fixes this problem at the expense of less accurate boundary conditions.
- `OBCS_SEAICE_SMOOTH_UVICE_PERP` (undef, OBCS_OPTIONS.h) — Smooth the component of sea-ice velocity perpendicular to the edge.
- `OBCS_SEAICE_SMOOTH_UVICE_PAR` (undef, OBCS_OPTIONS.h) — Smooth the component of sea ice velocity parallel to the edge.
- `OBCS_SEAICE_COMPUTE_UVICE` (undef, OBCS_OPTIONS.h) — Compute rather than specify seaice velocities at the edges.
- `OBCS_SEAICE_SMOOTH_EDGE` (undef, OBCS_OPTIONS.h) — Smooth the tracer sea-ice variables near the edges.
- `OBCS_AGEOS_COST_CONTRIBUTION` (undef, OBCS_OPTIONS.h) — o Flags related to Open-Boundary cost contributions o these flags refer to untested and potentially broken code
- `OBCS_VOLFLUX_COST_CONTRIBUTION` (undef, OBCS_OPTIONS.h)

## Headers
- `OBCS_FIELDS.h` — BOP
- `OBCS_GRID.h` — BOP
- `OBCS_OPTIONS.h` — CPP options file for OBCS package Use this file for selecting options within the OBCS package
- `OBCS_PARAMS.h` — BOP
- `OBCS_PTRACERS.h` — -- Fields and files for OBCS-support for passive tracers package PTRACERS
- `OBCS_SEAICE.h` — BOP
- `ORLANSKI.h` — SPK 6/2/00: Added storage arrays for salinity. Removed some unneeded arrays. SPK 7/18/00: Added dimensional phase speed arrays CVEL_**, where **=Varia
- `obcs_ad_check_lev1_dir.h` — ADJ STORE OBNu    = comlev1, key = ikey_dynamics ADJ STORE OBNv    = comlev1, key = ikey_dynamics ADJ STORE OBNt    = comlev1, key = ikey_dynamics ADJ
- `obcs_ad_check_lev2_dir.h` — ADJ STORE StoreOBCSN     = tapelev2, key = ilev_2 ADJ STORE StoreOBCSS     = tapelev2, key = ilev_2 ADJ STORE StoreOBCSE     = tapelev2, key = ilev_2 
- `obcs_ad_check_lev3_dir.h` — ADJ STORE StoreOBCSN     = tapelev3, key = ilev_3 ADJ STORE StoreOBCSS     = tapelev3, key = ilev_3 ADJ STORE StoreOBCSE     = tapelev3, key = ilev_3 
- `obcs_ad_check_lev4_dir.h` — ADJ STORE StoreOBCSN     = tapelev4, key = ilev_4 ADJ STORE StoreOBCSS     = tapelev4, key = ilev_4 ADJ STORE StoreOBCSE     = tapelev4, key = ilev_4 

## Routines (70)
`obcs_add_tides.F`, `obcs_adjust.F`, `obcs_adjust_uvice.F`, `obcs_apply_eta.F`, `obcs_apply_ptracer.F`, `obcs_apply_r_star.F`, `obcs_apply_seaice.F`, `obcs_apply_surf_dr.F`, `obcs_apply_ts.F`, `obcs_apply_uv.F`, `obcs_apply_uvice.F`, `obcs_apply_w.F`, `obcs_balance_flow.F`, `obcs_calc.F`, `obcs_calc_stevens.F`, `obcs_check.F`, `obcs_check_depths.F`, `obcs_copy_tracer.F`, `obcs_copy_uv_n.F`, `obcs_cost_ageos.F`, `obcs_cost_driver.F`, `obcs_cost_final.F`, `obcs_cost_ob_e.F`, `obcs_cost_ob_n.F`, `obcs_cost_ob_s.F`, `obcs_cost_ob_w.F`, `obcs_cost_vol.F`, `obcs_cost_weights.F`, `obcs_diag_balance.F`, `obcs_exchanges.F`, `obcs_exf_load.F`, `obcs_fields_load.F`, `obcs_init_fixed.F`, `obcs_init_variables.F`, `obcs_mon_stats.F`, `obcs_monitor.F`, `obcs_output.F`, `obcs_prescribe_read.F`, `obcs_read_pickup.F`, `obcs_readparms.F`, `obcs_save_uv_n.F`, `obcs_seaice_buffer_init.F`, `obcs_seaice_sponge.F`, `obcs_seaice_sponge_uvice.F`, `obcs_set_connect.F`, `obcs_sponge.F`, `obcs_u1_adv_tracer.F`, `obcs_write_pickup.F`, `orlanski_east.F`, `orlanski_init.F`, `orlanski_north.F`, `orlanski_south.F`, `orlanski_west.F`

## Called from outside the package
- `OBCS_SPONGE_S` ← `model/src/apply_forcing.F:969`
- `OBCS_SPONGE_T` ← `model/src/apply_forcing.F:737`
- `OBCS_SPONGE_U` ← `model/src/apply_forcing.F:180`
- `OBCS_SPONGE_V` ← `model/src/apply_forcing.F:369`
- `OBCS_APPLY_R_STAR` ← `model/src/calc_r_star.F:170`
- `OBCS_APPLY_SURF_DR` ← `model/src/calc_surf_dr.F:188`
- `OBCS_ADJUST` ← `model/src/do_oceanic_phys.F:595`
- `OBCS_CALC` ← `model/src/do_oceanic_phys.F:322`
- `OBCS_OUTPUT` ← `model/src/do_the_model_io.F:133`
- `OBCS_APPLY_UV` ← `model/src/dynamics.F:610`
- `OBCS_COPY_UV_N` ← `model/src/dynamics.F:410`
- `OBCS_EXCHANGES` ← `model/src/dynamics.F:665`
- `OBCS_SAVE_UV_N` ← `model/src/dynamics.F:607`
- `OBCS_SPONGE_S` ← `model/src/external_forcing.F:799`
- `OBCS_SPONGE_T` ← `model/src/external_forcing.F:595`
- `OBCS_SPONGE_U` ← `model/src/external_forcing.F:131`
- `OBCS_SPONGE_V` ← `model/src/external_forcing.F:270`
- `OBCS_CHECK_DEPTHS` ← `model/src/ini_depths.F:286`
- `OBCS_APPLY_W` ← `model/src/integr_continuity.F:299`
- `OBCS_APPLY_UV` ← `model/src/momentum_correction_step.F:95`
- `OBCS_CHECK` ← `model/src/packages_check.F:200`
- `OBCS_INIT_FIXED` ← `model/src/packages_init_fixed.F:223`
- `OBCS_INIT_VARIABLES` ← `model/src/packages_init_variables.F:631`
- `OBCS_READPARMS` ← `model/src/packages_readparms.F:166`
- `OBCS_WRITE_PICKUP` ← `model/src/packages_write_pickup.F:111`
- `OBCS_APPLY_TS` ← `model/src/thermodynamics.F:359`
- `OBCS_APPLY_ETA` ← `model/src/update_etah.F:75`
- `OBCS_COST_DRIVER` ← `pkg/cost/cost_driver.F:38`
- `OBCS_COST_FINAL` ← `pkg/cost/cost_final.F:105`
- `OBCS_DIAG_BALANCE` ← `pkg/diagnostics/diagnostics_calc_phivel.F:164`
- `OBCS_APPLY_PTRACER` ← `pkg/gchem/gchem_forcing_sep.F:275`
- `OBCS_U1_ADV_TRACER` ← `pkg/generic_advdiff/gad_advection.F:452,673`
- `OBCS_U1_ADV_TRACER` ← `pkg/generic_advdiff/gad_calc_rhs.F:302,431`
- `OBCS_APPLY_PTRACER` ← `pkg/ptracers/ptracers_integrate.F:514`
- `OBCS_APPLY_UVICE` ← `pkg/seaice/seaice_dynsolver.F:327`
- `OBCS_ADJUST_UVICE` ← `pkg/seaice/seaice_model.F:193`
- `OBCS_APPLY_SEAICE` ← `pkg/seaice/seaice_model.F:312`

## Verification experiments compiling it (10)
`dome` `exp4` `internal_wave` `isomip` `obcs_ctrl` `offline_exf_seaice` `seaice_obcs` `shelfice_2d_remesh` `so_box_biogeo` `tutorial_plume_on_slope`
