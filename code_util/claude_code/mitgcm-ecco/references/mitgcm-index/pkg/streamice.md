# pkg/streamice

Ice-stream/shelf dynamics model (shallow-shelf / hybrid), optionally coupled to shelfice.

**runtime switch:** `useSTREAMICE`-style flag in `data.pkg` (check exact name in packages_boot.F)
**reads:** `data.streamice`, `data.strmctrlflux`
**manual:** `doc/phys_pkgs/streamice.rst`, `doc/examples/examples.rst`, `doc/phys_pkgs/remesh.rst`
**adjoint support files:** streamice_ad_check_lev1_dir.h, streamice_ad_check_lev2_dir.h, streamice_ad_check_lev3_dir.h, streamice_ad_check_lev4_dir.h, streamice_ad_diff.list

## Namelist parameters
### STREAMICE_PARM01
- `streamice_density` — average ice density
- `streamice_density_ocean_avg` — average ocean density determining ice floatation
- `streamice_density_firn` — firn density in column
- `B_glen_isothermal` — (sqrt of) uniform ice stiffness coefficient (Pa 1/2 yr 1/6)
- `n_glen` — Glen s law exponent
- `eps_glen_min` — min strain rate in ice viscosity
- `eps_u_min` — min velocity in nonlinear sliding law
- `C_basal_fric_const` — (sqrt of) coefficient in sliding law (Pa 1/2 (m/yr) m/2)
- `n_basal_friction` — exponent in basal sliding law (tau = C u^n)
- `streamice_vel_update` — frequency of velocity solve (s) -- coupled ice-ocean only
- `streamice_cg_tol` — conj gradient tolerance
- `streamice_nonlin_tol` — nonlinear solver tolerance (relative residual, unitless)
- `streamice_nonlin_tol_fp` — fixed point nonlinear solver tolerance(absolute change, m/a)
- `streamice_err_norm`
- `streamice_max_cg_iter` — max conj gradient iterations
- `streamice_max_nl_iter` — max nonlin iterations in vel solve
- `streamice_maxcgiter_cpl` — max CG iters, coupled mode
- `streamice_maxnliter_cpl` — max NL iters, coupled mode
- `STREAMICEthickInit` — mode of thickness initialisation FILE - via STREAMICEthickFile PARAM - from STREAMICE_H_INIT_R common block
- `STREAMICEsigcoordInit`
- `STREAMICEsigcoordFile`
- `STREAMICEthickFile`
- `STREAMICEcalveMaskFile` — calving mask file
- `STREAMICEcostMaskFile` — mask to be used in "custom" cost function
- `STREAMICE_dump_mdsio`
- `STREAMICE_tave_mdsio`
- `STREAMICE_dump_mnc`
- `STREAMICE_tave_mnc`
- `STREAMICE_move_front` — advance ice-shelf front
- `STREAMICE_calve_to_mask` — do not advance front past streamice_calve_mask
- `STREAMICE_diagnostic_only` — do not update thickness
- `STREAMICE_lower_cg_tol` — lower CG tolerance when NL error is lowered by factor of .5e2
- `streamice_CFL_factor` — time step limiting factor
- `streamice_adjDump` — write frequency (s) of adjoint sensitivity fields
- `streamice_bg_surf_slope_x` — uniform surface slope, x-dir
- `streamice_bg_surf_slope_y` — uniform surface slope, y-dir
- `streamice_kx_b_init` — x-wave number for periodically initialised basal friction coeff
- `streamice_ky_b_init` — y-wave number for periodically initialised basal friction coeff
- `STREAMICEbasalTracConfig` — mode of sliding factor init FILE - via STREAMICEbasalTracFile UNIFORM - C_basal_fric_const 1DPERIODIC - varies in x-dir via streamice_kx_b_init and C_basal_fric_const 2DPERIODIC - varies in x- and y-dirs via streamice_kx_b_init
- `STREAMICEBdotConfig` — mode of ice-shelf melt rate init FILE - via STREAMICEBdotFile overridden in coupled mode
- `STREAMICEAdotConfig` — mode of SMB init FILE - via STREAMICEAdotFile o/w streamice_adot_uniform
- `STREAMICEbasalTracFile`
- `STREAMICEBdotFile`
- `STREAMICEAdotFile`
- `STREAMICEBdotTimeDepFile`
- `streamice_bdot_depth_nomelt` — pw linear depth-dependent melt param: depth above which no melt
- `streamice_bdot_depth_maxmelt` — pw linear depth-dependent melt param: depth below which const melt
- `streamice_bdot_maxmelt` — pw linear depth-dependent melt param: maximum melt rate
- `streamice_bdot_exp` — pw linear depth-dependent melt param: melt exponent
- `STREAMICEtopogFile` — bed topography (separate from ocean bathy)
- `STREAMICEhmaskFile` — ice mask file see EXPLANATION OF MASKS below
- `STREAMICEHBCyFile` — upstream thickness at y-boundaries
- `STREAMICEHBCxFile` — upstream thickness at x-boundaries -- to be used only with inhomogen. velocity condition
- `STREAMICEuFaceBdryFile` — streamice_ufacemask_bdry values see EXPLANATION OF MASKS below
- `STREAMICEvFaceBdryFile` — streamice_vfacemask_bdry values see EXPLANATION OF MASKS below
- `STREAMICEuDirichValsFile` — inhomogeneous x-vel dirich values to be set only where bound mask=3 see EXPLANATION OF MASKS below
- `STREAMICEvDirichValsFile` — inhomogeneous y-vel dirich values to be set only where bound mask=3 see EXPLANATION OF MASKS below
- `STREAMICEuMassFluxFile` — file to set u_flux_bdry_SI see EXPLANATION OF MASKS below
- `STREAMICEvMassFluxFile` — file to set v_flux_bdry_SI see EXPLANATION OF MASKS below
- `STREAMICEuNormalStressFile`
- `STREAMICEvNormalStressFile`
- `STREAMICEuShearStressFile`
- `STREAMICEvShearStressFile`
- `STREAMICEuNormalTimeDepFile`
- `STREAMICEvNormalTimeDepFile`
- `STREAMICEuShearTimeDepFile`
- `STREAMICEvShearTimeDepFile`
- `STREAMICEuFluxTimeDepFile`
- `STREAMICEvFluxTimeDepFile`
- `bdotMaxmeltTimeDepFile` — file giving a time and spatially dependent max melt at depth when Bdot_config='PARAM' overrides streamice_bdot_maxmelt and STREAMICEBdotMaxMeltFile
- `bglenTimeDepFile` — file giving a time and spatially dependent rheology param (Bbar) when STREAMICEGlenconstConfig='FILE' overrides STREAMICEGlenConstFile
- `cfricTimeDepFile` — file giving a time and spatially dependent sliding param when STREAMICEbasalTracConfig is 'FILE' overrides STREAMICEbasaltracFile
- `STREAMICEGlenConstFile`
- `STREAMICEGlenConstConfig` — mode of Glen s const init FILE - via STREAMICEGlenConstFile UNIFORM - B_glen_isothermal
- `STREAMICE_ppm_driving_stress` — use partial parabolic method to reconstruct surf slope
- `STREAMICE_h_ctrl_const_surf`
- `streamice_addl_backstress`
- `streamice_smooth_gl_width` — grounding line regularisation width (m)
- `streamice_adot_uniform` — uniform surface mass balance (m/yr)
- `streamice_firn_correction` — air thickness in column (m)
- `STREAMICE_apply_firn_correction`
- `STREAMICE_ADV_SCHEME` — DST3 -- 3rd order direct ST o/w 2nd order flux limited
- `streamice_forcing_period` — forcing freq (s)
- `STREAMICE_chkfixedptconvergence` — terminate velocity iteration based on fp_error
- `STREAMICE_chkresidconvergence` — terminate velocity iteration based on residual error
- `STREAMICE_alt_driving_stress` — use finite difference based driving stress (overrides above option)
- `STREAMICE_allow_reg_coulomb` — rather than using power-law sliding, implements "regularised coulomb" sliding law Asay-Davis et al (2016), Geosci. Model Dev.,
- `STREAMICE_use_log_ctrl` — fields C_basal_friction and Bglen (and initialisation values) given as the *logarithm* of physical values (if false, sqrt is used)
- `STREAMICE_vel_ext` — impose velocity with external files
- `STREAMICE_vel_ext_cgrid` — impose velocity with external files on C grid
- `STREAMICE_uvel_ext_file` — x-velocity file to replace velocity calc
- `STREAMICE_vvel_ext_file` — y-velocity file to replace velocity calc
- `STREAMICEBdotDepthFile` — file giving a spatially dependent depth below which const melt when Bdot_config='PARAM'. overrides streamice_bdot_depth_maxmelt
- `STREAMICEBdotMaxMeltFile` — file giving a spatially dependent max melt at depth when Bdot_config='PARAM' overrides streamice_bdot_maxmelt
- `STREAMICE_shelf_dhdt_ctrl` — option to apply surface elevation constraint to floating ice in cost function
- `streamice_buttr_width` — effective width for parameterisation of buttressing -- flowline mode only  _[ifdef STREAMICE_FLOWLINE_BUTTRESS]_
- `useStreamiceFlowlineButtr`  _[ifdef STREAMICE_FLOWLINE_BUTTRESS]_
- `STREAMICE_allow_cpl` — enable streamice-ocean coupling
- `streamice_smooth_thick_adjoint` — facility to smooth adjoint thickness sensitivity after advect_thickness 0 -> no smoothing  _[ifdef ALLOW_OPENAD]_
### STREAMICE_PARMTRACER
- `STREAMICETrac2DBCxFile`  _[ifdef ALLOW_STREAMICE_2DTRACER]_
- `STREAMICETrac2DBCyFile`  _[ifdef ALLOW_STREAMICE_2DTRACER]_
- `STREAMICETrac2DINITFile`  _[ifdef ALLOW_STREAMICE_2DTRACER]_
### STREAMICE_PARMPETSC
- `PETSC_PRECOND_TYPE` — JACOBI -- a jacobi precond (equiv to no petsc) BLOCKJACOBI -- block incomplete cholesky GAMG -- Algebraic multigrid MUMPS -- Direct ILU -- incomplete ILU (will not work in parallel)  _[ifdef ALLOW_PETSC]_
- `PETSC_SOLVER_TYPE` — CG, BICG, GMRES  _[ifdef ALLOW_PETSC]_
- `streamice_use_petsc`  _[ifdef ALLOW_PETSC]_
- `streamice_maxnliter_petsc` — max NL iters with PETSC unavailable with OpenAD  _[ifdef ALLOW_PETSC]_
- `streamice_petsc_pcfactorlevels` — fill level of incomplete cholesky preconditioner for use with PETSC and BLOCKJACOBI precond ONLY  _[ifdef ALLOW_PETSC]_
### STREAMICE_PARMOAD
- `streamice_nonlin_tol_adjoint` — fixed-point error of adjoint iterative solve (absolute)  _[ifdef ALLOW_STREAMICE_FP_ADJ]_
- `streamice_nonlin_tol_adjoint_rl`  _[ifdef ALLOW_STREAMICE_FP_ADJ]_
- `STREAMICE_OAD_petsc_reuse`  _[ifdef ALLOW_STREAMICE_FP_ADJ & ifdef ALLOW_PETSC]_
- `PETSC_PRECOND_OAD`  _[ifdef ALLOW_STREAMICE_FP_ADJ & ifdef ALLOW_PETSC]_
### STREAMICE_COST
- `STREAMICEvelOptimSnapBasename`  _[ifdef ALLOW_COST_STREAMICE]_
- `STREAMICEvelOptimTCBasename` — file prefix for obs velocities in transient inversion, expects files  _[ifdef ALLOW_COST_STREAMICE]_
- `STREAMICEsurfOptimTCBasename` — file prefix for obs surf elev in transient inversion, expects files  _[ifdef ALLOW_COST_STREAMICE]_
- `STREAMICEBglenCostMaskFile` — prior values for Bglen in transient or snapshot inversion  _[ifdef ALLOW_COST_STREAMICE]_
- `streamice_wgt_drift` — cost function coefficient of drift term  _[ifdef ALLOW_COST_STREAMICE]_
- `streamice_wgt_vel` — cost function coefficient of vel misfit term  _[ifdef ALLOW_COST_STREAMICE]_
- `streamice_wgt_vel_norm` — cost function coefficient of vel misfit normalised by obs. vel  _[ifdef ALLOW_COST_STREAMICE]_
- `streamice_wgt_surf` — cost function coefficient of surface misfit term  _[ifdef ALLOW_COST_STREAMICE]_
- `streamice_wgt_tikh_beta` — cost function coefficient of sq gradient penalty for sliding param  _[ifdef ALLOW_COST_STREAMICE]_
- `streamice_wgt_tikh_bglen` — cost function coefficient of sq gradient penalty for stiffness  _[ifdef ALLOW_COST_STREAMICE]_
- `streamice_wgt_tikh_gen` — cost function coefficient of sq gradient penalty for generic ctrl  _[ifdef ALLOW_COST_STREAMICE]_
- `streamice_wgt_prior_bglen` — cost function coefficient of sq deviation from prior stiffness  _[ifdef ALLOW_COST_STREAMICE]_
- `streamice_wgt_prior_gen` — cost function coefficient of sq deviation from prior, generic ctrl  _[ifdef ALLOW_COST_STREAMICE]_
- `STREAMICE_do_snapshot_cost` — accumulate snapshot cost function at final time step  _[ifdef ALLOW_COST_STREAMICE]_
- `STREAMICE_do_timedep_cost` — accumulate cost at specified time steps  _[ifdef ALLOW_COST_STREAMICE]_
- `STREAMICE_do_verification_cost`  _[ifdef ALLOW_COST_STREAMICE]_
- `STREAMICE_do_vaf_cost` — do cost for volume above floatation  _[ifdef ALLOW_COST_STREAMICE]_
- `streamice_vel_cost_timesteps` — array of time steps where velocity misfit is accumulated, expects files  _[ifdef ALLOW_COST_STREAMICE & ifdef ALLOW_STREAMICE_TC_COST]_
- `streamice_surf_cost_timesteps` — array of time steps where surface misfit is accumulated, expects files  _[ifdef ALLOW_COST_STREAMICE & ifdef ALLOW_STREAMICE_TC_COST]_
### STREAMICE_PARM02
- `shelf_max_draft`
- `shelf_min_draft`
- `shelf_edge_pos`
- `shelf_slope_scale`
- `shelf_flat_width`
- `flow_dir`
### STREAMICE_PARM03
- `min_x_noflow_NORTH`
- `max_x_noflow_NORTH`
- `min_x_noflow_SOUTH`
- `max_x_noflow_SOUTH`
- `min_y_noflow_WEST`
- `max_y_noflow_WEST`
- `min_y_noflow_EAST`
- `max_y_noflow_EAST`
- `min_x_noStress_NORTH`
- `max_x_noStress_NORTH`
- `min_x_noStress_SOUTH`
- `max_x_noStress_SOUTH`
- `min_y_noStress_WEST`
- `max_y_noStress_WEST`
- `min_y_noStress_EAST`
- `max_y_noStress_EAST`
- `min_x_FluxBdry_NORTH`
- `max_x_FluxBdry_NORTH`
- `min_x_FluxBdry_SOUTH`
- `max_x_FluxBdry_SOUTH`
- `min_y_FluxBdry_WEST`
- `max_y_FluxBdry_WEST`
- `min_y_FluxBdry_EAST`
- `max_y_FluxBdry_EAST`
- `min_x_Dirich_NORTH`
- `max_x_Dirich_NORTH`
- `min_x_Dirich_SOUTH`
- `max_x_Dirich_SOUTH`
- `min_y_Dirich_WEST`
- `max_y_Dirich_WEST`
- `min_y_Dirich_EAST`
- `max_y_Dirich_EAST`
- `min_x_CFBC_NORTH`
- `max_x_CFBC_NORTH`
- `min_x_CFBC_SOUTH`
- `max_x_CFBC_SOUTH`
- `min_y_CFBC_WEST`
- `max_y_CFBC_WEST`
- `min_y_CFBC_EAST`
- `max_y_CFBC_EAST`
- `flux_bdry_val_SOUTH`
- `flux_bdry_val_NORTH`
- `flux_bdry_val_WEST`
- `flux_bdry_val_EAST`
- `STREAMICE_NS_periodic`
- `STREAMICE_EW_periodic`

## CPP options (defaults as shipped)
- `STREAMICE_CONSTRUCT_MATRIX` (define, STREAMICE_OPTIONS.h)
- `STREAMICE_HYBRID_STRESS` (undef, STREAMICE_OPTIONS.h)
- `STREAMICE_FLOWLINE_BUTTRESS` (undef, STREAMICE_OPTIONS.h)
- `USE_ALT_RLOW` (undef, STREAMICE_OPTIONS.h)
- `STREAMICE_GEOM_FILE_SETUP` (undef, STREAMICE_OPTIONS.h)
- `STREAMICE_SMOOTH_FLOATATION` (undef, STREAMICE_OPTIONS.h) — on height above floatation, and option (2) will also smooth surface elevation across grounding line; only one should be defined
- `STREAMICE_SMOOTH_FLOATATION2` (undef, STREAMICE_OPTIONS.h)
- `ALLOW_PETSC` (undef, STREAMICE_OPTIONS.h)
- `ALLOW_STREAMICE_2DTRACER` (undef, STREAMICE_OPTIONS.h)
- `STREAMICE_TRACER_AB` (undef, STREAMICE_OPTIONS.h)
- `STREAMICE_SERIAL_TRISOLVE` (undef, STREAMICE_OPTIONS.h)
- `STREAMICE_3D_GLEN_CONST` (undef, STREAMICE_OPTIONS.h) — -  Undocumented Options:
- `STREAMICE_COULOMB_SLIDING` (undef, STREAMICE_OPTIONS.h)
- `STREAMICE_ECSECRYO_DOSUM` (undef, STREAMICE_OPTIONS.h)
- `STREAMICE_FALSE` (undef, STREAMICE_OPTIONS.h)
- `STREAMICE_FIRN_CORRECTION` (undef, STREAMICE_OPTIONS.h)
- `STREAMICE_PETSC_3_8` (undef, STREAMICE_OPTIONS.h)
- `STREAMICE_STRESS_BOUNDARY_CONTROL` (undef, STREAMICE_OPTIONS.h)
- `ALLOW_STREAMICE_TIMEDEP_FORCING` (undef, STREAMICE_OPTIONS.h)
- `ALLOW_STREAMICE_FLUX_CONTROL` (undef, STREAMICE_OPTIONS.h)
- `ALLOW_STREAMICE_TC_COST` (undef, STREAMICE_OPTIONS.h)
- `ALLOW_STREAMICE_FP_ADJ` (undef, STREAMICE_OPTIONS.h) — Christianson et al 1994, Opt. Meth. & Software ; this reduce size of adjoint tape in memory as well as decouple forward and reverse convergence criteria

## Headers
- `STREAMICE.h` — ---+----1--+-+----2----+----3----+----4----+----5----+----6----+----7-|--+----
- `STREAMICE_ADV.h` — ---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----
- `STREAMICE_BDRY.h` — ---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----
- `STREAMICE_CG.h` — ---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----
- `STREAMICE_COST_SIZE.h` — Maximum number of cost levels for STREAMICE transient calibration cost streamiceMaxCostLevel :: max list of timestep levels where cost is applied
- `STREAMICE_CTRL_FLUX.h` — ---+----1--+-+----2----+----3----+----4----+----5----+----6----+----7-|--+----
- `STREAMICE_FP.h` — ---+----1--+-+----2----+----3----+----4----+----5----+----6----+----7-|--+----
- `STREAMICE_OPTIONS.h` — BOP
- `STREAMICE_PETSC_MOD.h` — ---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----
- `streamice_ad_check_lev1_dir.h` — ADJ STORE area_shelf_streamice ADJ &     = comlev1, key=ikey_dynamics, kind=isbyte ADJ STORE streamice_hmask ADJ &     = comlev1, key=ikey_dynamics, k
- `streamice_ad_check_lev2_dir.h` — ADJ STORE area_shelf_streamice = tapelev2, key = ilev_2 ADJ STORE streamice_hmask = tapelev2, key = ilev_2 ADJ STORE u_streamice = tapelev2, key = ile
- `streamice_ad_check_lev3_dir.h` — ADJ STORE area_shelf_streamice = tapelev3, key = ilev_3 ADJ STORE streamice_hmask = tapelev3, key = ilev_3 ADJ STORE u_streamice = tapelev3, key = ile
- `streamice_ad_check_lev4_dir.h` — ADJ STORE area_shelf_streamice = tapelev4, key = ilev_4 ADJ STORE streamice_hmask = tapelev4, key = ilev_4 ADJ STORE u_streamice = tapelev4, key = ile

## Routines (66)
`adstreamice_cg_solve.F`, `adstreamice_invert_surf_forthick.F`, `eta_gl_prime_streamice.F`, `eta_gl_streamice.F`, `eta_gl_streamice_prime.F`, `phi_gl_streamice.F`, `slope_limiter.F`, `streamice_adv_flux_fl_x.F`, `streamice_adv_flux_fl_y.F`, `streamice_adv_front.F`, `streamice_advect_2dtracer.F`, `streamice_advect_thickness.F`, `streamice_apply_flux_ctrl.F`, `streamice_bstress_exponent.F`, `streamice_cg_functions.F`, `streamice_cg_solve.F`, `streamice_cg_solve_matfree.F`, `streamice_cg_solve_petsc.F`, `streamice_cg_wrapper.F`, `streamice_check.F`, `streamice_cost_accum.F`, `streamice_cost_final.F`, `streamice_cost_reg_accum.F`, `streamice_cost_surf_accum.F`, `streamice_cost_vel_accum.F`, `streamice_diagnostics_state.F`, `streamice_driving_stress.F`, `streamice_driving_stress_fd.F`, `streamice_driving_stress_ppm.F`, `streamice_dump.F`, `streamice_dump_ad.F`, `streamice_fields_load.F`, `streamice_finalize_petsc.F`, `streamice_forced_buttress.F`, `streamice_get_fp_err_oad.F`, `streamice_get_vel_fp_err.F`, `streamice_get_vel_resid_err.F`, `streamice_get_vel_resid_err_oad.F`, `streamice_init_diagnostics.F`, `streamice_init_fixed.F`, `streamice_init_phi.F`, `streamice_init_varia.F`, `streamice_initialise_petsc.F`, `streamice_invert_surf_forthick.F`, `streamice_petsc_numerate.F`, `streamice_petscmatdestroy.F`, `streamice_read_pickup.F`, `streamice_readparms.F`, `streamice_tap_fixedpoint_notreduced.F`, `streamice_taub.F`, `streamice_timestep.F`, `streamice_tridiag_solve.F`, `streamice_upd_ffrac_uncoupled.F`, `streamice_vel_phi.F`, `streamice_vel_phistage.F`, `streamice_vel_solve.F`, `streamice_vel_solve_fp.F`, `streamice_vel_solve_openad.F`, `streamice_velmask_upd.F`, `streamice_visc_beta.F`, `streamice_visc_beta_hybrid.F`, `streamice_write_pickup.F`

## Called from outside the package
- `STREAMICE_DIAGNOSTICS_STATE` ← `model/src/do_statevars_diags.F:108`
- `STREAMICE_TIMESTEP` ← `model/src/forward_step.F:662`
- `STREAMICE_CHECK` ← `model/src/packages_check.F:366`
- `STREAMICE_INIT_FIXED` ← `model/src/packages_init_fixed.F:457`
- `STREAMICE_INIT_VARIA` ← `model/src/packages_init_variables.F:372`
- `STREAMICE_READPARMS` ← `model/src/packages_readparms.F:276`
- `STREAMICE_WRITE_PICKUP` ← `model/src/packages_write_pickup.F:171`
- `STREAMICE_FINALIZE_PETSC` ← `model/src/the_model_main.F:754`
- `STREAMICE_COST_FINAL` ← `pkg/cost/cost_final.F:121`
- `STREAMICE_COST_ACCUM` ← `pkg/cost/cost_tile.F:132`
- `STREAMICE_FINALIZE_PETSC` ← `pkg/openad/the_model_main.F:327`
- `STREAMICE_INITIALIZE_PETSC` ← `pkg/openad/the_model_main.F:131`
- `ADSTREAMICE_CG_SOLVE` ← `pkg/tapenade/stubs_tap_adj.F:348,482`
- `STREAMICE_CG_SOLVE` ← `pkg/tapenade/stubs_tap_adj.F:244,295,413`

## Verification experiments compiling it (1)
`halfpipe_streamice`
