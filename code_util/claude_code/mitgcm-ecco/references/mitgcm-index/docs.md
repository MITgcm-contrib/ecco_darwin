# MITgcm manual section map (doc/**/*.rst)

The online manual (https://mitgcm.readthedocs.io) is built from these files. Every heading is listed with
its line number (`L123`) and the parameters / CPP flags / source files the section cites, so
`grep -n <name> docs.md` finds the section that explains a parameter. Read that line range from the
source tree (`git -C <clone> show origin/master:doc/... | sed -n <a>,<b>p`) for equations and intent.


## `doc/algorithm/adv-schemes.rst` — Linear advection schemes
- L1 Linear advection schemes
  - L13 Centered second order advection-diffusion — `uT`, `uTrans`, `tracer`, `vT`, `vTrans`, `wT`, `rTrans`, `gad_c2_adv_x.F`, `gad_c2_adv_y.F`, `gad_c2_adv_r.F`
  - L75 Third order upwind bias advection — `uT`, `uTrans`, `tracer`, `vT`, `vTrans`, `wT`, `rTrans`, `gad_u3_adv_x.F`, `gad_u3_adv_y.F`, `gad_u3_adv_r.F`
  - L120 Centered fourth order advection — `uT`, `uTrans`, `tracer`, `vT`, `vTrans`, `wT`, `rTrans`, `gad_c4_adv_x.F`, `gad_c4_adv_y.F`, `gad_c4_adv_r.F`
  - L165 First order upwind advection
- L190 Non-linear advection schemes
  - L205 Second order flux limiters — `uT`, `uTrans`, `tracer`, `vT`, `vTrans`, `wT`, `rTrans`, `gad_fluxlimit_adv_x.F`, `gad_fluxlimit_adv_y.F`, `gad_fluxlimit_adv_r.F`
  - L262 Third order direct space-time — `uT`, `uTrans`, `tracer`, `vT`, `vTrans`, `wT`, `rTrans`, `gad_dst3_adv_x.F`, `gad_dst3_adv_y.F`, `gad_dst3_adv_r.F`
  - L327 Third order direct space-time with flux limiting — `uT`, `uTrans`, `tracer`, `vT`, `vTrans`, `wT`, `rTrans`, `gad_dst3fl_adv_x.F`, `gad_dst3fl_adv_y.F`, `gad_dst3fl_adv_r.F`
  - L371 Multi-dimensional advection — `tracer`, `gTracer`, `aF`, `uTrans`, `vTrans`, `rTrans`, `gad_advection.F`
- L428 Comparison of advection schemes — `multiDimAdvection`, `OLx`

## `doc/algorithm/algorithm.rst` — Discretization and Algorithm
- L3 Discretization and Algorithm
  - L15 Notation
  - L66 Time-stepping
  - L120 Pressure method with rigid-lid — `implicitViscosity`, `timestep.F`, `solve_for_pressure.F`, `correction_step.F`, `forward_step.F`, `dynamics.F`, `calc_div_ghat.F`, `cg2d.F`, `momentum_correction_step.F`, `calc_grad_phi_surf.F`
  - L268 Pressure method with implicit linear free-surface — `freesurfFac`
  - L357 Explicit time-stepping: Adams-Bashforth
    - L379 Adams-Bashforth II — `forcing_In_AB`, `tracForcingOutAB`, `momForcingOutAB`, `momDissip_In_AB`
    - L424 Adams-Bashforth III — `momDissip_In_AB`, `ALLOW_ADAMSBASHFORTH_3`, `CPP_OPTIONS.h`
  - L504 Implicit time-stepping: backward method — `forward_step.F`, `thermodynamics.F`, `temp_integrate.F`, `gad_calc_rhs.F`, `external_forcing.F`, `adams_bashforth2.F`, `timestep_tracer.F`, `impldiff.F`
  - L568 Synchronous time-stepping: variables co-located in time — `hFac`, `thermodynamics.F`, `dynamics.F`, `solve_for_pressure.F`, `momentum_correction_step.F`, `forward_step.F`, `external_fields_load.F`, `do_atmospheric_phys.F`, `do_oceanic_phys.F`, `calc_gt.F`, `gad_calc_rhs.F`, `external_forcing.F`, `adams_bashforth2.F`, `timestep_tracer.F`, `impldiff.F`, `calc_phi_hyd.F`, `mom_fluxform.F`, `mom_vecinv.F`, `timestep.F`, `update_r_star.F`, `update_surf_dr.F`, `calc_div_ghat.F`, `cg2d.F`, `calc_grad_phi_surf.F`, `correction_step.F`, `tracers_correction_step.F`, `cycle_tracer.F`, `shap_filt_apply_ts.F`, `zonal_filt_apply_ts.F`, `convective_adjustment.F`
  - L692 Staggered baroclinic time-stepping — `staggerTimeStep`, `forward_step.F`, `external_fields_load.F`, `do_atmospheric_phys.F`, `do_oceanic_phys.F`, `dynamics.F`, `calc_phi_hyd.F`, `mom_fluxform.F`, `mom_vecinv.F`, `timestep.F`, `impldiff.F`, `update_r_star.F`, `update_surf_dr.F`, `solve_for_pressure.F`, `calc_div_ghat.F`, `cg2d.F`, `momentum_correction_step.F`, `calc_grad_phi_surf.F`, `correction_step.F`, `thermodynamics.F`, `calc_gt.F`, `gad_calc_rhs.F`, `external_forcing.F`, `adams_bashforth2.F`, `timestep_tracer.F`, `tracers_correction_step.F`, `cycle_tracer.F`, `shap_filt_apply_ts.F`, `zonal_filt_apply_ts.F`, `convective_adjustment.F`
  - L821 Non-hydrostatic formulation
  - L964 Variants on the Free Surface — `gU`, `gV`, `cg2d_b`, `etaN`, `phi_nh`, `uVel`, `vVel`, `solve_for_pressure.F`, `DYNVARS.h`, `SOLVE_FOR_PRESSURE.h`, `correction_step.F`, `NH_VARS.h`, `cg2d.F`, `ini_cg2d.F`, `calc_div_ghat.F`, `cg3d.F`, `ini_cg3d.F`
  - L1078 Spatial discretization of the dynamical equations
  - L1094 Continuity and horizontal pressure gradient term
  - L1145 Hydrostatic balance — `buoyancyRelation`, `calc_phi_hyd.F`
  - L1190 Flux-form momentum equations — `gU`, `gV`, `gW`, `mom_fluxform.F`, `DYNVARS.h`, `NH_VARS.h`
    - L1232 Advection of momentum — `fZon`, `fMer`, `fVerUkp`, `fVerVkp`, `mom_u_adv_uu.F`, `mom_u_adv_vu.F`, `mom_u_adv_wu.F`, `mom_fluxform.F`, `mom_v_adv_uv.F`, `mom_v_adv_vv.F`, `mom_v_adv_wv.F`
    - L1292 Coriolis terms — `selectCoriScheme`, `cF`, `cd_code_scheme.F`, `mom_u_coriolis.F`, `mom_v_coriolis.F`, `mom_fluxform.F`
    - L1363 Curvature metric terms — `mT`, `mom_u_metric_sphere.F`, `mom_v_metric_sphere.F`, `mom_fluxform.F`
    - L1412 Non-hydrostatic metric terms — `mT`, `mom_u_metric_nh.F`, `mom_v_metric_nh.F`, `mom_fluxform.F`
    - L1454 Lateral dissipation — `viscAh`, `viscA4`, `vF`, `v4F`, `mom_u_xviscflux.F`, `mom_u_yviscflux.F`, `mom_fluxform.F`, `mom_v_xviscflux.F`, `mom_v_yviscflux.F`, `mom_u_sidedrag.F`, `mom_v_sidedrag.F`
    - L1561 Vertical dissipation — `fVrUp`, `fVrDw`, `bottomDragLinear`, `bottomDragQuadratic`, `ALLOW_BOTTOMDRAG_ROUGHNESS`, `zRoughBot`, `vF`, `mom_u_rviscflux.F`, `mom_v_rviscflux.F`, `mom_fluxform.F`, `MOM_COMMON_OPTIONS.h`, `mom_u_bottomdrag.F`, `mom_v_bottomdrag.F`
    - L1678 Derivation of discrete energy conservation
    - L1691 Mom Diagnostics
  - L1762 Vector invariant momentum equations — `gU`, `gV`, `gW`, `mom_vecinv.F`, `DYNVARS.h`, `NH_VARS.h`
    - L1820 Relative vorticity — `vort3`, `mom_calc_relvort3.F`, `mom_vecinv.F`
    - L1844 Kinetic energy — `KE`, `mom_calc_KE.F`, `mom_vecinv.F`
    - L1859 Coriolis terms — `uCf`, `vCf`, `mom_vi_coriolis.F`, `mom_vi_u_coriolis.F`, `mom_vi_v_coriolis.F`, `mom_vecinv.F`
    - L1919 Shear terms — `uCf`, `vCf`, `mom_vi_u_vertshear.F`, `mom_vi_v_vertshear.F`, `mom_vecinv.F`
    - L1944 Gradient of Bernoulli function — `uCf`, `vCf`, `mom_vi_u_grad_ke.F`, `mom_vi_v_grad_ke.F`, `mom_vecinv.F`
    - L1964 Horizontal divergence — `hDiv`, `mom_calc_ke.F`, `mom_vecinv.F`
    - L1982 Horizontal dissipation — `uDissip`, `vDissip`, `mom_vi_hdissip.F`
    - L2020 Vertical dissipation — `vrf`, `mom_u_rviscflux.F`, `mom_vecinv.F`
  - L2052 Tracer equations
    - L2067 Time-stepping of tracers: ABII — `tau`, `gTracer`, `fVerT`, `gTrNm1`, `ABeps`, `tracer`, `deltaTtracer`, `calc_gt.F`, `calc_gs.F`, `gad_calc_rhs.F`, `PARAMS.h`, `timestep_tracer.F`
  - L2165 Advection schemes
  - L2173 Shapiro Filter
    - L2224 SHAP Diagnostics
  - L2239 Nonlinear Viscosities for Large Eddy Simulation
    - L2303 Eddy Viscosity
      - L2335 Reynolds-Number Limited Eddy Viscosity — `viscAh`, `viscAhReMax`
      - L2370 Vertical Eddy Viscosities — `viscAr`
      - L2381 Smagorinsky Viscosity — `viscAh`, `viscC2Smag`
      - L2467 Leith Viscosity — `useFullLeith`
      - L2505 Modified Leith Viscosity — `viscC2LeithD`, `viscC2Leith`
      - L2555 Quasi-Geostrophic Leith Viscosity — `ALLOW_LEITH_QG`, `viscC2LeithQG`, `useFullLeith`, `ALLOW_GM_LEITH_QG`, `GM_useLeithQG`, `MOM_COMMON_OPTIONS.h`, `GMREDI_OPTIONS.h`
      - L2624 Courant–Freidrichs–Lewy Constraint on Viscosity — `viscAhGridMax`, `viscA4GridMax`, `viscAhGridMin`, `viscA4GridMin`
      - L2649 Biharmonic Viscosity
      - L2749 Selection of Length Scale — `useAreaViscLength`
    - L2765 Mercator, Nondimensional Equations

## `doc/algorithm/c-grid.rst` — C grid staggering of variables
- L1 C grid staggering of variables
- L24 Grid initialization and data — `ini_grid.F`, `ini_vertical_grid.F`, `ini_cartesian_grid.F`, `ini_spherical_polar_grid.F`, `ini_curvilinear_grid.F`, `ini_masks_etc.F`, `initialise_fixed.F`, `GRID.h`

## `doc/algorithm/crank-nicol.rst` — Crank-Nicolson barotropic time stepping
- L3 Crank-Nicolson barotropic time stepping — `implicSurfPress`, `implicDiv2DFlow`, `useRealFreshWaterFlux`, `implicitNHPress`, `correction_step.F`

## `doc/algorithm/finitevol-meth.rst` — The finite volume method: finite volumes versus finite difference
- L1 The finite volume method: finite volumes versus finite difference

## `doc/algorithm/horiz-grid.rst` — Horizontal grid
- L3 Horizontal grid — `dxG`, `dyG`, `rA`, `dxC`, `dyC`, `rAz`, `dxV`, `dyF`, `rAw`, `dxF`, `dyU`, `rAs`, `ini_cartesian_grid.F`, `ini_spherical_polar_grid.F`, `ini_curvilinear_grid.F`, `GRID.h`
  - L96 Reciprocals of horizontal grid descriptors — `recip_rA`, `recip_rAz`, `recip_rAw`, `recip_rAs`, `recip_dxG`, `recip_dyG`, `recip_dxC`, `recip_dyC`, `recip_dxF`, `recip_dyF`, `recip_dxV`, `recip_dyU`, `recip_`, `GRID.h`, `ini_masks_etc.F`
  - L117 Cartesian coordinates — `usingCartesianGrid`, `dXspacing`, `dYspacing`, `DELX`, `DELY`
  - L128 Spherical-polar coordinates — `usingSphericalPolarGrid`, `dXspacing`, `dYspacing`, `DELX`, `DELY`
  - L139 Curvilinear coordinates — `usingCurvilinearGrid`

## `doc/algorithm/nonlinear-freesurf.rst` — Non-linear free-surface
- L3 Non-linear free-surface
  - L8 Pressure/geo-potential and free surface — `nonlinFreeSurf`, `uniformLin_PhiSurf`, `ini_linear_phisurf.F`
  - L92 Free surface effect on column total thickness (Non-linear free-surface) — `nonlinFreeSurf`, `EXACT_CONSERV`, `select_rStar`
  - L208 Tracer conservation with non-linear free-surface
  - L275 Time stepping implementation of the non-linear free-surface — `exactConserv`, `cg2dTargetResidual`, `etaH`, `etaHnm1`, `dEtaHdt`, `solve_for_pressure.F`, `integr_continuity.F`, `timestep.F`, `calc_gt.F`, `calc_gs.F`, `forward_step.F`, `DYNVARS.h`, `SURFACE.h`
  - L408 Non-linear free-surface and vertical resolution — `hFacInf`

## `doc/algorithm/vert-grid.rst` — Vertical grid
- L2 Vertical grid — `delR`, `delZ`, `delP`, `delRc`, `drF`, `drC`, `recip_drF`, `recip_drC`, `ini_vertical_grid.F`, `GRID.h`
- L59 Topography: partially filled cells — `hFacC`, `hFacW`, `hFacS`, `recip_hFacC`, `recip_hFacW`, `recip_hFacS`, `hFacMin`, `hFacMinDr`, `ini_masks_etc.F`, `GRID.h`

## `doc/autodiff/autodiff.rst` — Automatic Differentiation
- L3 Automatic Differentiation
  - L50 Some basic algebra
    - L81 Forward or direct sensitivity
    - L115 Reverse or adjoint sensitivity
      - L482 Example 1: :math:`{\cal J} = v_{j} (T)`
      - L497 Example 2: :math:`{\cal J} = \langle \, {\cal H}(\vec{v}) - \vec{d} \, , \, {\cal H}(\vec{v}) - \vec{d} \, \rangle`
    - L526 Storing vs. recomputation in reverse mode
  - L614 TLM and ADM generation in general — `ALLOW_AUTODIFF_TAMC`, `ALLOW_ADJOINT_RUN`, `ALLOW_TANGENTLINEAR_RUN`, `ALLOW_GRDCHK`, `the_model_main.F`, `the_main_loop.F`, `ctrl_unpack.F`, `ctrl_pack.F`, `grdchk_main.F`
    - L697 General setup — `AUTODIFF_OPTIONS.h`, `COST_OPTIONS.h`, `CTRL_OPTIONS.h`, `tamc.h`
    - L731 Building the AD code using TAF
    - L794 The AD build process in detail — `ALLOW_ADJOINT_RUN`, `ALLOW_TANGENTLINEAR_RUN`, `myThid`, `AD_FILES`, `AD_FLOW_FILES`
      - L829 The list ``AD_FILES`` and ``*_ad_diff.list`` files — `AD_FILES`
      - L853 The list ``AD_FLOW_FILES`` and ``.flow`` files — `AD_FILES`, `AD_FLOW_FILES`
      - L887 Store directives for 3-level checkpointing
      - L912 Changing the default AD tool flags: ad_options files
      - L915 Hand-written adjoint code
    - L920 The cost function (dependent variable)
      - L964 Enabling the package — `ALLOW_COST`, `ALLOW_COST_TRACER`, `ALLOW_AUTODIFF_TAMC`, `COST_OPTIONS.h`, `AUTODIFF_OPTIONS.h`
      - L981 Initialization — `ALLOW_COST`, `mult_tracer`, `objf_tracer`, `cost_readparms.F`, `cost_init_varia.F`
      - L1001 Accumulation — `ALLOW_COST_TRACER`, `objf_tracer`, `cost_tile.F`, `cost_tracer.F`
      - L1013 Finalize all contributions — `fc`, `cost_final.F`
    - L1084 The control variables (independent variables) — `CTRL_OPTIONS.h`
      - L1136 Initialization — `nvarlength`, `nWetCtile`, `nWetWtile`, `nWetStile`, `ALLOW_NONDIMENSIONAL_CONTROL_IO`, `ctrl_readparms.F`, `ctrl_unpack.F`
      - L1174 Perturbation of the independent variables — `pTracer`, `xx_tr1`, `diffkr`, `kapgm`, `xx_tr1_dummy`, `ptracers_init_varia.F`, `ctrl_map_ini.F`, `ctrl_map_forcing.F`
      - L1266 Output of adjoint variables and gradient — `ALLOW_AUTODIFF_MONITOR`, `addynvars_r`, `addynvars_cd`, `addynvars_diffkr`, `addynvars_kapgm`, `adtr1_r`, `adffields`, `dynvars_r`, `dynvars_cd`, `vector_ctrl`, `vector_grad`, `ctrl_map_ini.F`, `ctrl_map_forcing.F`, `ctrl_unpack.F`, `ctrl_pack.F`, `addummy_in_stepping.F`, `AUTODIFF_OPTIONS.h`, `dummy_in_stepping.F`, `the_main_loop.F`, `adcommon.h`, `dynamics.F`, `addummy_in_dynamics.F`
      - L1336 Control variable handling for optimization applications — `ctrl_unpack.F`, `ctrl_pack.F`
  - L1395 The gradient check package — `grdchk_main.F`, `the_model_main.F`
    - L1424 Code description
    - L1427 Code configuration — `ALLOW_ADJOINT_RUN`, `ALLOW_GRADIENT_CHECK`, `useGrdchk`, `grdchk_eps`, `nbeg`, `nstep`, `nend`, `grdchkvarindex`, `CPP_OPTIONS.h`
  - L1481 Adjoint dump & restart – divided adjoint (DIVA)
    - L1492 Introduction
    - L1583 Recipe for divided adjoint code generation — `ALLOW_DIVIDED_ADJOINT`, `genmake_local`, `adthe_main_loop_ad`, `adthe_main_loop`, `AUTODIFF_OPTIONS.h`
    - L1633 Special considerations for multi processor (MPI) runs — `mpi_comm_world`, `mpi_integer`
  - L1647 Adjoint code generation using OpenAD
    - L1662 Introduction
    - L1689 Downloading and installing OpenAD
    - L1698 Building MITgcm adjoint with OpenAD
    - L1721 Building the MITgcm adjoint using an OpenAD Singularity container
  - L1752 Adjoint code generation using Tapenade
    - L1762 Introduction
    - L1775 Downloading and installing Tapenade
    - L1786 Prerequisites for Linux or Mac OS
    - L1793 Steps for Mac OS — `MITGCM_ROOTDIR`, `genmake_local`
    - L1832 Steps for Linux — `JAVA_HOME`
    - L1854 Prerequisites for Windows
    - L1862 Steps for Windows — `JAVA_HOME`
    - L1907 Prerequisites for Tapenade setup — `code_tap`, `ALLOW_AUTODIFF_MONITOR`
    - L1923 Building MITgcm TLM with Tapenade
    - L1951 Building MITgcm adjoint with Tapenade

## `doc/contributing/contributing.rst` — Contributing to the MITgcm
- L3 Contributing to the MITgcm
  - L9 Bugs and feature requests
  - L30 Using Git and Github
    - L41 Quickstart Guide
    - L129 Detailed guide for those less familiar with Git and GitHub
  - L423 Coding style guide
  - L428 Creating MITgcm packages
    - L448 Package structure — `CPP_OPTIONS.h`, `GMREDI_OPTIONS.h`
    - L538 Package boot sequence — `packages_boot.F`, `packages_readparms.F`, `initialise_fixed.F`, `packages_init_fixed.F`, `packages_check.F`, `packages_init_variables.F`, `initialise_varia.F`
    - L606 Package S/R calls — `gmredi_calc_tensor.F`, `do_oceanic_phys.F`, `do_the_model_io.F`, `packages_write_pickup.F`, `the_model_main.F`, `gmredi_diagnostics_init.F`
    - L659 Package “mypackage”
  - L669 MITgcm code testing protocols
    - L678 Test-experiment directory content — `nPx`, `nPy`, `_mpi`, `mitgcmuv_ad`, `genmake_local`, `prepare_run`, `SIZE.h`
    - L790 The testreport utility — `cg2d_init_res`, `tr_checklist`
      - L945 Reference Output
    - L957 The do_tst_2+2 utility
    - L1020 Daily Testing of MITgcm
    - L1032 Required Testing for MITgcm Code Contributors
      - L1035 Using testreport to check your new code — `global_sum.F`
      - L1084 Using do_tst_2+2 to check your new code
      - L1093 Automatic testing with GitHub Actions
  - L1113 Contributing to the manual
    - L1143 Section headings
    - L1156 Internal document references
    - L1181 Citations
    - L1211 Other embedded links — `dynamics.F`
    - L1245 Symbolic Notation
    - L1284 Figures
    - L1315 Tables
    - L1378 Other text blocks — `subroutine_name.F`, `where_var1_defined.h`, `where_var2_defined.h`, `where_var3_defined.h`
    - L1437 Other style conventions
    - L1461 Building the manual
  - L1513 Reviewing pull requests

## `doc/examples/advection_in_gyre/advection_in_gyre.rst` — Ocean Gyre Advection Schemes
- L3 Ocean Gyre Advection Schemes
  - L40 Advection and tracer transport
  - L63 Introducing a tracer into the flow
  - L84 Selecting an advection scheme — `PTRACERS_ALLOW_DYN_STATE`, `PTRACERS_OPTIONS.h`
  - L95 Comparison of different advection schemes

## `doc/examples/baroclinic_gyre/baroclinic_gyre.rst` — Baroclinic Ocean Gyre
- L4 Baroclinic Ocean Gyre
  - L81 Equations solved — `gBaro`, `gravity`
  - L164 Discrete Numerical Configuration
    - L204 Numerical Stability Criteria
  - L277 Configuration — `SIZE.h`, `DIAGNOSTICS_SIZE.h`
    - L301 Compile-time Configuration
      - L304 File :filelink:`code/packages.conf <verification/tutorial_baroclinic_gyre/code/packages.conf>`
      - L325 File :filelink:`code/SIZE.h <verification/tutorial_baroclinic_gyre/code/SIZE.h>` — `nSx`, `nSy`, `nPx`, `nPy`, `SIZE.h`
      - L398 File :filelink:`code/DIAGNOSTICS_SIZE.h <verification/tutorial_baroclinic_gyre/code/DIAGNOSTICS_SIZE.h>` — `numDiags`, `DIAGNOSTICS_SIZE.h`
    - L417 Run-time Configuration
      - L422 File :filelink:`input/data <verification/tutorial_baroclinic_gyre/input/data>`
      - L650 File :filelink:`input/data.pkg <verification/tutorial_baroclinic_gyre/input/data.pkg>` — `useMNC`, `useDiagnostics`
      - L666 File :filelink:`input/data.mnc <verification/tutorial_baroclinic_gyre/input/data.mnc>` — `monitor_mnc`, `mnc_use_outdir`, `mnc_outdir_str`, `mnc_test_`
      - L695 File :filelink:`input/data.diagnostics <verification/tutorial_baroclinic_gyre/input/data.diagnostics>`
      - L777 File :filelink:`input/eedata <verification/tutorial_baroclinic_gyre/input/eedata>`
      - L788 Files ``input/bathy.bin``, ``input/windx_cosy.bin``
      - L796 File ``input/SST_relax.bin``
  - L805 Building and running the model — `endTime`, `monitorFreq`
    - L817 Output Files — `monitor_mnc`, `XC`, `YC`, `dxC`, `dyC`, `XG`, `YG`, `dxG`, `dyG`, `dyU`, `dxV`, `RC`, `RF`, `drC`, `drF`, `rA`, `rAw`, `rAs`, `rAz`, `Depth`, `HFacC`, `HFacW`, `HFacS`, `nIter0`, `dumpFreq`, `deltaT`, `Nr`, `mnc_test_0001`, `mnc_test_0002`, `mnc_test_`, `DIAG_STATIS_PARMS`, `PHrefC`, `PHrefF`, `RhoRef`, `phiHyd`
  - L967 Running with MPI — `nSx`, `nSy`, `nPx`, `nPy`, `mnc_use_outdir`, `build_mpi`, `mnc_test_00001`, `mnc_test_00004`, `SIZE.h`
  - L1037 Running with OpenMP — `nTy`, `nTx`, `nSx`, `nSy`, `useMNC`, `globalFiles`, `build_openmp`, `OMP_NUM_THREADS`, `surfDiag`, `oceStDiag`, `SIZE.h`
  - L1100 Model solution — `drF`, `dyG`, `TRELAX_ave`, `THETA_lv_avg`, `THETA_lv_std`

## `doc/examples/barotropic_gyre/barotropic_gyre.rst` — Barotropic Ocean Gyre
- L3 Barotropic Ocean Gyre
  - L44 Equations Solved
  - L76 Discrete Numerical Configuration
    - L89 Numerical Stability Criteria
  - L155 Configuration — `SIZE.h`
    - L176 Compile-time Configuration
      - L179 File :filelink:`code/SIZE.h <verification/tutorial_barotropic_gyre/code/SIZE.h>` — `sNx`, `sNy`, `OLx`, `OLy`, `nSx`, `nSy`, `nPx`, `nPy`, `Nr`, `SIZE.h`
    - L242 Run-time Configuration
      - L247 File :filelink:`input/data <verification/tutorial_barotropic_gyre/input/data>`
      - L487 File :filelink:`input/data.pkg <verification/tutorial_barotropic_gyre/input/data.pkg>`
      - L499 File :filelink:`input/eedata <verification/tutorial_barotropic_gyre/input/eedata>`
      - L511 File ``input/bathy.bin``
      - L527 File ``input/windx_cosy.bin``
  - L544 Building and running the model
    - L568 Standard output — `monitorFreq`, `debugLevel`, `SIZE.h`
    - L607 Other output files — `etaN`, `dumpFreq`, `uVel`, `theta`, `salt`, `wVel`, `pChkptFreq`, `ChkptFreq`, `hFacC`, `hFacS`, `hFacW`, `PHrefC`, `PHrefF`
  - L684 Model Solution

## `doc/examples/cfc_offline/cfc_offline.rst` — Offline Experiments
- L3 Offline Experiments
  - L11 Overview
  - L31 Time-stepping of tracers
  - L37 Code Configuration — `PTRACERS_SIZE.h`, `GMREDI_OPTIONS.h`, `SIZE.h`
    - L69 File :filelink:`input_tutorial/data <verification/tutorial_cfc_offline/input_tutorial/data>` — `nIter0`, `nTimesteps`, `deltaTtracer`, `deltatTracer`, `deltaToffline`, `deltaTClock`, `deltatTtracer`, `periodicExternalForcing`, `externForcingPeriod`, `externForcingCycle`
    - L215 File :filelink:`input_tutorial/data.off <verification/tutorial_cfc_offline/input_tutorial/data.off>` — `UvelFile`, `VvelFile`, `WvelFile`, `ConvFile`, `deltaToffline`, `offlineForcingPeriod`, `offlineForcingCycle`, `offlineIter0`, `uVeltave`, `vVeltave`, `wVeltave`, `nIter0`, `deltatToffline`
    - L297 File :filelink:`input_tutorial/data.pkg <verification/tutorial_cfc_offline/input_tutorial/data.pkg>` — `usePTRACERS`
    - L316 File :filelink:`input_tutorial/data.ptracers <verification/tutorial_cfc_offline/input_tutorial/data.ptracers>` — `PTRACERS_numInUse`, `PTRACERS_Iter0`, `PTRACERS_advScheme`, `PTRACERS_diffKh`, `PTRACERS_diffKr`, `PTRACERS_initialFile`, `PTRACERS_SIZE.h`
    - L392 File :filelink:`input_tutorial/eedata <verification/tutorial_cfc_offline/input_tutorial/eedata>`
    - L398 File :filelink:`code/packages.conf <verification/tutorial_cfc_offline/code/packages.conf>`
    - L410 File :filelink:`code/PTRACERS_SIZE.h <verification/tutorial_cfc_offline/code/PTRACERS_SIZE.h>` — `PTRACERS_num`
    - L426 File :filelink:`code/SIZE.h <verification/tutorial_cfc_offline/code/SIZE.h>`
  - L461 Running the Experiment
  - L469 A more complicated example
    - L503 File :filelink:`input/data <verification/tutorial_cfc_offline/input/data>` — `implicitDiffusion`, `nTimeSteps`, `pChkptFreq`, `chkptFreq`, `dumpFreq`, `taveFreq`
    - L532 File :filelink:`input/data.off <verification/tutorial_cfc_offline/input/data.off>`
    - L548 File :filelink:`input/data.pkg <verification/tutorial_cfc_offline/input/data.pkg>`
    - L559 File :filelink:`input/data.ptracers <verification/tutorial_cfc_offline/input/data.ptracers>` — `PTRACERS_useGMRedi`, `PTRACERS_useKPP`, `PTRACERS_initialFile`
    - L609 File :filelink:`input/data.gchem <verification/tutorial_cfc_offline/input/data.gchem>`
    - L619 File :filelink:`input/data.gmredi <verification/tutorial_cfc_offline/input/data.gmredi>`
    - L628 File :filelink:`input/cfc1112.atm <verification/tutorial_cfc_offline/input/cfc1112.atm>`
    - L635 Running the Experiment

## `doc/examples/deep_convection/deep_convection.rst` — Deep Convection
- L3 Deep Convection
  - L33 Overview — `theta`, `SIZE.h`
  - L93 Equations solved
  - L164 Discrete numerical configuration
  - L172 Numerical stability criteria and other considerations
  - L199 Experiment configuration — `CPP_OPTIONS.h`, `SIZE.h`
    - L215 File :filelink:`code/CPP_OPTIONS.h <verification/tutorial_deep_convection/code/CPP_OPTIONS.h>`
    - L221 File :filelink:`code/SIZE.h <verification/tutorial_deep_convection/code/SIZE.h>`
    - L257 File :filelink:`input/data <verification/tutorial_deep_convection/input/data>` — `SIZE.h`
    - L525 File :filelink:`input/data.pkg <verification/tutorial_deep_convection/input/data.pkg>`
    - L531 File :filelink:`input/eedata <verification/tutorial_deep_convection/input/eedata>`
    - L537 File ``input/Qsurf.bin``

## `doc/examples/examples.rst` — MITgcm Tutorial Example Experiments
- L3 MITgcm Tutorial Example Experiments — `SIZE.h`, `CPP_OPTIONS.h`
  - L249 Additional Example Experiments: Forward Model Setups — `useLANGMUIR`
  - L604 Additional Example Experiments: Adjoint Model Setups — `code_ad`, `input_ad`, `code_tap`, `input_tap`, `cg2d_nsa.F`, `cg2d.F`, `cg2d_mad.F`

## `doc/examples/global_oce_biogeo/global_oce_biogeo.rst` — Biogeochemistry Simulation
- L3 Biogeochemistry Simulation
  - L8 Overview
  - L56 Equations Solved
  - L129 Code configuration — `nTimeSteps`, `taveFreq`, `PTRACERS_numInUse`, `PTRACERS_Iter0`, `SIZE.h`, `PTRACERS_SIZE.h`, `DIAGNOSTICS_SIZE.h`
  - L206 Running the example — `taveFreq`, `DIC_Biotave`, `DIC_Cartave`, `DIC_fluxCO2ave`, `DIC_pCO2tave`, `DIC_pHtave`, `DIC_SurOtave`, `DIC_Surtave`

## `doc/examples/global_oce_in_p/global_oce_in_p.rst` — Global Ocean Simulation in Pressure Coordinates
- L3 Global Ocean Simulation in Pressure Coordinates
  - L23 Overview
  - L80 Discrete Numerical Configuration
  - L203 Experiment Configuration — `CPP_OPTIONS.h`, `SIZE.h`
    - L240 Driving Datasets
    - L309 File :filelink:`input/data <verification/tutorial_global_oce_in_p/input/data>` — `startTime`
    - L633 File :filelink:`input/data.pkg <verification/tutorial_global_oce_in_p/input/data.pkg>`
    - L639 File :filelink:`input/eedata <verification/tutorial_global_oce_in_p/input/eedata>`
    - L645 File ``input/topog.bin``
    - L657 File ``input/deltageopotjmd95.box``
    - L667 Files ``input/lev_t.bin`` and ``input/lev_s.bin``
    - L677 Files ``input/trenberth_taux.bin`` and ``input/trenberth_tauy.bin``
    - L686 File ``input/lev_sst.bin``
    - L693 Files ``input/shi_qnet.bin`` and ``input/shi_empmr.bin``
    - L702 File :filelink:`code/SIZE.h <verification/tutorial_global_oce_in_p/code/SIZE.h>` — `SIZE.h`
    - L708 File :filelink:`code/CPP_OPTIONS.h <verification/tutorial_global_oce_in_p/code/CPP_OPTIONS.h>` — `ATMOSPHERIC_LOADING`, `EXACT_CONSERV`, `NONLIN_FRSURF`

## `doc/examples/global_oce_latlon/global_oce_latlon.rst` — Global Ocean Simulation
- L3 Global Ocean Simulation
  - L19 Overview
  - L77 Discrete Numerical Configuration
    - L192 Numerical Stability Criteria
  - L273 Experiment Configuration — `SIZE.h`
    - L304 Driving Datasets
    - L365 File :filelink:`input/data <verification/tutorial_global_oce_latlon/input/data>` — `nIter0`, `startTime`, `endTime`
    - L569 File :filelink:`input/data.pkg <verification/tutorial_global_oce_latlon/input/data.pkg>`
    - L575 File :filelink:`input/eedata <verification/tutorial_global_oce_latlon/input/eedata>`
    - L581 Files ``input/trenberth_taux.bin`` and ``input/trenberth_tauy.bin``
    - L589 File ``input/bathymetry.bin``
    - L600 File :filelink:`code/SIZE.h <verification/tutorial_global_oce_latlon/code/SIZE.h>`

## `doc/examples/global_oce_optim/global_oce_optim.rst` — Global Ocean State Estimation
- L3 Global Ocean State Estimation
  - L8 Overview
  - L116 Implementation of the control variable and the cost function — `ALLOW_HFLUXM_CONTROL`, `cHFLUXM_CONTROL`
    - L128 The control variable — `ALLOW_HFLUXM_CONTROL`, `CTRL_OPTIONS.h`, `FFIELDS.h`, `ini_forcing.F`, `external_forcing_surf.F`, `ctrl_init.F`, `ctrl_pack.F`, `ctrl_unpack.F`, `ctrl_map_forcing.F`
    - L153 Cost functions — `ALLOW_COST_TEMP`, `objf_temp`, `mult_temp`, `ALLOW_COST_HFLUXM`, `objf_hfluxm`, `mult_hflux`, `cost_temp.F`, `COST_OPTIONS.h`, `cost_tile.F`, `cost_hflux.F`, `cost_final.F`, `cost_weights.F`, `cost_readparms.F`
  - L182 Code Configuration
    - L191 Compilation-time customizations in :filelink:`code_ad <verification/tutorial_global_oce_optim/code_ad/>` — `ALLOW_ECCO_OPTIMIZATION`, `CTRL_OPTIONS.h`
    - L200 Running-time customizations in :filelink:`input_ad <verification/tutorial_global_oce_optim/input_ad/>` — `cg2dTargetResidual`, `lastinterval`, `useGrdchk`
  - L222 Compiling — `mitgcmuv_ad`
    - L229 Compilation of MITgcm and its adjoint: ``mitcgmuv_ad`` — `mitgcmuv_ad`, `cost_temp.F`, `cost_hflux.F`
    - L252 Compilation of the line-search algorithm: ``optim.x`` — `INCLUDEDIRS`, `xx_hfluxm_file`, `mitgcm_ad`, `optim_numbmod.F`
  - L282 Running the estimation — `fmin`, `nTimeSteps`, `lastinterval`, `mitgcmuv_ad`

## `doc/examples/held_suarez_cs/held_suarez_cs.rst` — Held-Suarez Atmosphere
- L3 Held-Suarez Atmosphere
  - L15 Overview
  - L51 Forcing
  - L116 Set-up description
    - L186 Numerical Stability Criteria
  - L218 Experiment Configuration — `CPP_OPTIONS.h`, `SIZE.h`, `DIAGNOSTICS_SIZE.h`, `apply_forcing.F`
    - L246 File :filelink:`input/data <verification/tutorial_held_suarez_cs/input/data>` — `nIter0`, `startTime`, `grid_cs32`
    - L540 File :filelink:`input/data.pkg <verification/tutorial_held_suarez_cs/input/data.pkg>`
    - L580 File :filelink:`input/data.shap <verification/tutorial_held_suarez_cs/input/data.shap>` — `nShapUV`, `nShapUVPhys`
    - L639 File :filelink:`input/eedata <verification/tutorial_held_suarez_cs/input/eedata>`
    - L655 File :filelink:`code/SIZE.h <verification/tutorial_held_suarez_cs/code/SIZE.h>`
    - L700 File :filelink:`code/packages.conf <verification/tutorial_held_suarez_cs/code/packages.conf>` — `useEXCH2`
    - L758 File :filelink:`code/CPP_OPTIONS.h <verification/tutorial_held_suarez_cs/code/CPP_OPTIONS.h>`
    - L770 Other Files — `EXTERNAL_FORCING_U`, `EXTERNAL_FORCING_V`, `EXTERNAL_FORCING_T`, `EXTERNAL_FORCING_S`, `apply_forcing.F`, `GRID.h`

## `doc/examples/plume_on_slope/plume_on_slope.rst` — Gravity Plume On a Continental Slope
- L3 Gravity Plume On a Continental Slope
  - L92 Configuration
  - L103 Binary input data
  - L205 Code configuration — `ALLOW_NONHYDROSTATIC`, `nonHydrostatic`, `useOBCS`, `SIZE.h`, `CPP_OPTIONS.h`
  - L230 Model parameters — `eosType`, `sBeta`

## `doc/examples/reentrant_channel/reentrant_channel.rst` — Southern Ocean Reentrant Channel Example
- L3 Southern Ocean Reentrant Channel Example
  - L96 Equations Solved
  - L138 Discrete Numerical Configuration — `tempAdvScheme`
    - L178 Numerical Stability Criteria — `implicitViscosity`
  - L253 Configuration — `SIZE.h`, `LAYERS_SIZE.h`, `DIAGNOSTICS_SIZE.h`
    - L279 Compile-time Configuration
      - L282 File :filelink:`code/packages.conf <verification/tutorial_reentrant_channel/code/packages.conf>` — `diffKh`
      - L316 File :filelink:`code/SIZE.h <verification/tutorial_reentrant_channel/code/SIZE.h>` — `nSy`, `OLx`, `OLy`
      - L333 File :filelink:`code/LAYERS_SIZE.h <verification/tutorial_reentrant_channel/code/LAYERS_SIZE.h>` — `Nlayers`
      - L346 File :filelink:`code/DIAGNOSTICS_SIZE.h <verification/tutorial_reentrant_channel/code/DIAGNOSTICS_SIZE.h>` — `numDiags`
    - L355 Run-time Configuration
      - L360 File :filelink:`input/data <verification/tutorial_reentrant_channel/input/data>`
      - L552 File :filelink:`input/data.pkg <verification/tutorial_reentrant_channel/input/data.pkg>`
      - L576 File :filelink:`input/data.gmredi <verification/tutorial_reentrant_channel/input/data.gmredi>` — `background_K`, `GM_isopycK`
      - L627 File :filelink:`input/data.rbcs <verification/tutorial_reentrant_channel/input/data.rbcs>` — `useRBCtemp`, `tauRelaxT`, `relaxMaskFile`, `relaxTFile`
      - L642 File :filelink:`input/data.layers <verification/tutorial_reentrant_channel/input/data.layers>` — `pkg`, `layers_maxNum`, `layers_name`, `layers_bounds`, `Nlayers`, `LAYERS_SIZE.h`
      - L683 File :filelink:`input/data.diagnostics <verification/tutorial_reentrant_channel/input/data.diagnostics>`
      - L750 File :filelink:`input/eedata <verification/tutorial_reentrant_channel/input/eedata>`
      - L758 File ``input/bathy.50km.bin`` — `rF`, `hFacMin`, `hFacMinDr`, `hFacC`
      - L790 File ``input/zonal_wind.50km.bin``, ``input/SST_relax.50km.bin``
      - L801 File ``input/temperature.50km.bin``
      - L816 File ``input/T_relax_mask.50km.bin``
  - L830 Building and running the model — `nTimeSteps`, `monitorFreq`, `nPy`, `nSy`, `nPx`, `useGMRedi`, `diffKhT`, `SIZE.h`
  - L860 Model Solution
    - L876 Coarse Resolution Solution — `rC`, `stdiags_bylev`, `stdiags_2D`, `LAYERS_SIZE.h`
    - L1106 Eddy Permitting Solution — `useGMRedi`, `DeltaT`, `nTimeSteps`, `viscAh`, `viscC2Leith`, `useFullLeith`, `viscAhGridMax`, `useSingleCpuIO`, `GM_background_K`, `HeatCapacity_Cp`, `rA`, `drF`, `hFacC`, `tauRelaxT`

## `doc/examples/rotating_tank/rotating_tank.rst` — Rotating Tank
- L3 Rotating Tank
  - L18 Equations Solved
  - L21 Discrete Numerical Configuration
  - L37 Code Configuration — `CPP_OPTIONS.h`, `SIZE.h`
    - L55 File :filelink:`input/data <verification/tutorial_rotating_tank/input/data>`
    - L216 File - :filelink:`input/data.pkg <verification/tutorial_rotating_tank/input/data.pkg>`
    - L222 File - :filelink:`input/eedata <verification/tutorial_rotating_tank/input/eedata>`
    - L228 File ``input/thetaPolR.bin``
    - L236 File ``input/bathyPolR.bin``
    - L245 File :filelink:`code/SIZE.h <verification/tutorial_rotating_tank/code/SIZE.h>`
    - L274 File :filelink:`code/CPP_OPTIONS.h <verification/tutorial_rotating_tank/code/CPP_OPTIONS.h>`

## `doc/examples/tracer_adjsens/tracer_adjsens.rst` — Adjoint Sensitivity Analysis for Tracer Injection
- L3 Adjoint Sensitivity Analysis for Tracer Injection
  - L16 Overview of the experiment
    - L23 Passive tracer equation
    - L54 Model configuration
    - L65 Out-gassing cost function
  - L88 Code configuration — `COST_OPTIONS.h`, `CTRL_OPTIONS.h`, `CPP_OPTIONS.h`, `AUTODIFF_OPTIONS.h`, `CTRL_SIZE.h`, `GAD_OPTIONS.h`, `GMREDI_OPTIONS.h`, `SIZE.h`, `tamc.h`, `ctrl_map_ini_genarr.F`, `ptracers_forcing_surf.F`
    - L151 File :filelink:`code_ad/COST_OPTIONS.h /<verification/tutorial_tracer_adjsens/code_ad/COST_OPTIONS.h>`
    - L156 File :filelink:`code_ad/CTRL_OPTIONS.h /<verification/tutorial_tracer_adjsens/code_ad/CTRL_OPTIONS.h>`
    - L161 File :filelink:`code_ad/CPP_OPTIONS.h /<verification/tutorial_tracer_adjsens/code_ad/CPP_OPTIONS.h>` — `the_main_loop.F`
    - L189 File ``ECCO_OPTIONS.h`` — `ALLOW_AUTODIFF_TAMC`, `ALLOW_TAMC_CHECKPOINTING`, `ALLOW_AUTODIFF_MONITOR`, `ALLOW_DIVIDED_ADJOINT`, `ALLOW_COST`, `ALLOW_COST_TRACER`, `ALLOW_THETA0_CONTROL`, `ALLOW_SALT0_CONTROL`, `ALLOW_TR10_CONTROL`, `ALLOW_TAUU0_CONTROL`, `ALLOW_TAUV0_CONTROL`, `ALLOW_SFLUX0_CONTROL`, `ALLOW_HFLUX0_CONTROL`, `ALLOW_DIFFKR_CONTROL`, `ALLOW_KAPGM_CONTROL`, `tamc.h`, `checkpoint_lev3_directives.h`, `checkpoint_lev2_directives.h`
    - L250 File :filelink:`SIZE.h <verification/tutorial_tracer_adjsens/code_ad/SIZE.h>`
    - L260 File :filelink:`/pkg/autodiff/adcommon.h` — `addynvars_r`, `addynvars_cd`, `addynvars_diffkr`, `addynvars_kapgm`, `adtr1_r`, `adffields`, `ALLOW_AUTODIFF_MONITOR`, `addummy_in_stepping.F`, `DYNVARS.h`, `FFIELDS.h`, `adcommon.h`
    - L291 File :filelink:`code_ad/tamc.h <verification/tutorial_tracer_adjsens/code_ad/tamc.h>` — `ALLOW_TAMC_CHECKPOINTING`, `nchklev_3`, `nchklev_2`, `nchklev_1`, `nTimeSteps`, `endTime`, `startTime`, `deltaTClock`, `nchklev_0`, `isbyte`, `maxpass`, `the_main_loop.F`
    - L333 File ``makefile``
      - L410 File ``input/topog.bin``
      - L415 Files ``input/windx.bin``, ``input/windy.bin``, ``input/salt.bin``, ``input/theta.bin``, ``input/SSS.bin``, ``input/SST.bin``
  - L421 Compiling the model and its adjoint
    - L439 Adjoint code generation and compilation – step by step — `myThId`
      - L510 Adjoint code generation and compilation – summary

## `doc/getting_started/getting_started.rst` — Getting Started with MITgcm
- L3 Getting Started with MITgcm
  - L32 Where to find information
  - L44 Obtaining the code
    - L66 Method 1
    - L100 Method 2
  - L117 Updating the code
  - L190 Model and directory structure — `main.F`, `the_model_main.F`
  - L257 Building the model
    - L262 Quickstart Guide — `SIZE.h`
    - L345 Generating a ``Makefile`` using genmake2 — `genmake_local`
      - L529 Command-line options:
      - L672 Optfiles in tools/build_options directory:
    - L858 ``make`` commands
    - L921 Building with MPI — `nPx`, `nPy`, `SIZE.h`
    - L1024 Building  with OpenMP — `EEPARAMS.h`
  - L1053 Running the model — `prepare_run`
    - L1099 Running with MPI
    - L1127 Running with OpenMP — `nTx`, `nTy`, `nSx`, `nSy`, `OMP_NUM_THREADS`, `OMP_STACKSIZE`, `SIZE.h`
    - L1169 Output files
      - L1181 Raw binary output files — `ckptA`, `ckptB`
      - L1243 NetCDF output files — `mnc_output_`
    - L1279 Looking at the output
      - L1282 MATLAB
      - L1327 Python
      - L1372 Bash scripts
  - L1397 Customizing the Model Configuration - Code Parameters and Compilation Options
    - L1400 Model Array Dimensions — `Nx`, `Ny`, `sNx`, `sNy`, `Nr`, `OLx`, `OLy`, `nSx`, `nSy`, `nPx`, `nPy`, `SIZE.h`
    - L1445 C Preprocessor Options — `SHORTWAVE_HEATING`, `ALLOW_GEOTHERMAL_FLUX`, `ALLOW_FRICTION_HEATING`, `ALLOW_ADDFLUID`, `ATMOSPHERIC_LOADING`, `ALLOW_BALANCE_FLUXES`, `ALLOW_BALANCE_RELAX`, `CHECK_SALINITY_FOR_NEGATIVE_VALUES`, `EXCLUDE_FFIELDS_LOAD`, `INCLUDE_PHIHYD_CALCULATION_CODE`, `INCLUDE_CONVECT_CALL`, `INCLUDE_CALC_DIFFUSIVITY_CALL`, `ALLOW_3D_DIFFKR`, `ALLOW_BL79_LAT_VARY`, `EXCLUDE_PCELL_MIX_CODE`, `ALLOW_SMAG_3D_DIFFUSIVITY`, `ALLOW_SOLVE4_PS_AND_DRAG`, `INCLUDE_IMPLVERTADV_CODE`, `ALLOW_ADAMSBASHFORTH_3`, `ALLOW_QHYD_STAGGER_TS`, `EXACT_CONSERV`, `NONLIN_FRSURF`, `ALLOW_NONHYDROSTATIC`, `ALLOW_EDDYPSI`, `ALLOW_CG2D_NSA`, `ALLOW_SRCG`, `SOLVE_DIAGONAL_LOWMEMORY`, `SOLVE_DIAGONAL_KINNER`, `CPP_OPTIONS.h`
    - L1551 Preprocessor Execution Environment Options — `GLOBAL_SUM_ORDER_TILES`, `CG2D_SINGLECPU_SUM`, `SINGLE_DISK_IO`, `USE_FORTRAN_SCRATCH_FILES`, `COMPONENT_MODULE`, `DISCONNECTED_TILES`, `REAL4_IS_SLOW`, `CPP_EEOPTIONS.h`, `global_sum_singlecpu.F`, `cg2d.F`
  - L1678 Customizing the Model Configuration - Runtime Parameters — `PARAMS.h`, `set_defaults.F`, `ini_parms.F`
    - L1711 Parameters: Configuration, Computational Domain, Geometry, and Time-Discretization
      - L1716 Model Configuration — `buoyancyRelation`, `nonHydrostatic`, `quasiHydrostatic`, `rhoRefFile`, `ALLOW_NONHYDROSTATIC`
      - L1741 Grid — `usingCartesianGrid`, `usingSphericalPolarGrid`, `usingCylindricalGrid`, `usingCurvilinearGrid`, `xgOrigin`, `ygOrigin`, `delX`, `delY`, `cosPower`, `delR`, `horizGridFile`, `seaLev_Z`, `top_Pres`, `dxSpacing`, `delXFile`, `dySpacing`, `delYFile`, `delRc`, `delRFile`, `delRcFile`, `rSphere`, `selectFindRoSurf`, `radius_fromHorizGrid`, `phiEuler`, `thetaEuler`, `psiEuler`
      - L1855 Topography - Full and Partial Cells — `bathyFile`, `hFacMin`, `hFacC`, `delR`, `topoFile`, `addWwallFile`, `addSwallFile`, `hFacMinDr`, `hFacInf`, `hFacSup`, `useMin4hFacEdges`, `hFacW`, `hFacS`, `pCellMix_select`, `pCellMix_maxFac`, `pCellMix_delR`
      - L1924 Physical Constants — `rhoConst`, `rhoNil`, `gravity`, `gravityFile`, `gBaro`
      - L1941 Rotation — `f0`, `beta`, `rotationPeriod`, `omega`, `selectCoriMap`, `fPrime`
      - L1977 Free Surface — `rigidLid`, `implicitFreeSurface`, `useRealFreshWaterFlux`, `implicSurfPress`, `implicDiv2Dflow`, `implicitNHPress`, `nonlinFreeSurf`, `NONLIN_FRSURF`, `select_rStar`, `selectNHfreeSurf`, `exactConserv`
      - L2027 Time-Discretization — `deltaTMom`, `deltaTtracer`, `deltaT`, `deltaTClock`, `baseTime`, `deltaTmom`, `dTtracerLev`, `deltaTfreesurf`
    - L2062 Parameters: Main Algorithmic Parameters
      - L2065 Pressure Solver — `cg2dMaxIters`, `cg2dTargetResidual`, `cg3dMaxIters`, `cg3dTargetResidual`, `cg2dTargetResWunit`, `cg2dPreCondFreq`, `cg2dUseMinResSol`, `ALLOW_NONHYDROSTATIC`, `cg3dTargetResWunit`, `useSRCGSolver`, `printResidualFreq`, `debugLevel`, `integr_GeoPot`, `uniformLin_PhiSurf`, `deepAtmosphere`, `nh_Am2`
      - L2122 Time-Stepping Algorithm — `abEps`, `staggerTimeStep`, `alph_AB`, `ALLOW_ADAMSBASHFORTH_3`, `beta_AB`, `multiDimAdvection`, `implicitIntGravWave`, `ALLOW_NONHYDROSTATIC`
    - L2150 Parameters: Equation of State — `eosType`, `tAlpha`, `sBeta`, `tRef`, `sRef`, `implicitIntGravWave`, `selectP_inEOS_Zc`, `salt`, `tRefFile`, `thetaConst`, `sRefFile`, `rhonil`
      - L2282 Thermodynamic Constants — `HeatCapacity_Cp`, `celsius2K`, `atm_Cp`, `atm_Rd`, `atm_Rq`, `atm_Po`
    - L2305 Parameters: Momentum Equations
      - L2308 Configuration — `momViscosity`, `momAdvection`, `useCoriolis`, `momStepping`, `metricTerms`, `momPressureForcing`, `implicitViscosity`, `selectmetricTerms`, `uVel`, `dxC`, `useNHMTerms`, `momImplVertAdv`, `INCLUDE_IMPLVERTADV_CODE`, `interViscAr_pCell`, `momDissip_In_AB`, `selectCoriScheme`, `select3dCoriScheme`, `vectorInvariantMomentum`, `useJamartMomAdv`, `selectVortScheme`, `upwindVorticity`, `useAbsVorticity`, `highOrderVorticity`, `upwindShear`, `selectKEscheme`, `mom_calc_ke.F`
      - L2405 Initialization — `uVelInitFile`, `vVelInitFile`, `pSurfInitFile`
      - L2432 General Dissipation Scheme — `viscAh`, `viscAr`, `viscA4`, `viscAhD`, `viscAhZ`, `viscAhW`, `viscAhDfile`, `ALLOW_3D_VISCAH`, `viscAhZfile`, `viscAhGrid`, `viscAhMax`, `viscAhGridMax`, `viscAhGridMin`, `viscAhReMax`, `viscC2leith`, `viscC2leithD`, `viscC2LeithQG`, `viscC2smag`, `viscA4D`, `viscA4Z`, `viscA4W`, `viscA4Dfile`, `ALLOW_3D_VISCA4`, `viscA4Zfile`, `viscA4Grid`, `viscA4Max`, `viscA4GridMax`, `viscA4GridMin`, `viscA4ReMax`, `viscC4leith`, `viscC4leithD`, `viscC4smag`, `useFullLeith`, `useSmag3D`, `ALLOW_SMAG_3D`, `smag3D_coeff`, `useStrainTensionVisc`, `useAreaViscLength`, `viscArNr`, `pCellMix_viscAr` (+1)
      - L2535 Sidewall/Bottom Dissipation — `no_slip_sides`, `no_slip_bottom`, `bottomDragLinear`, `bottomDragQuadratic`, `sideDragFactor`, `zRoughBot`, `selectBotDragQuadr`, `selectImplicitDrag`, `ALLOW_SOLVE4_PS_AND_DRAG`, `bottomVisc_pCell`
    - L2585 Parameters: Tracer Equations
      - L2592 Configuration — `tempAdvection`, `tempStepping`, `saltAdvection`, `implicitDiffusion`, `tempAdvScheme`, `tempVertAdvScheme`, `tempImplVertAdv`, `addFrictionHeating`, `ALLOW_FRICTION_HEATING`, `temp_stayPositive`, `GAD_SMOLARKIEWICZ_HACK`, `saltStepping`, `saltAdvScheme`, `saltVertAdvScheme`, `saltImplVertAdv`, `salt_stayPositive`, `interDiffKr_pCell`, `linFSConserveTr`, `doAB_onGtGs`, `GAD_OPTIONS.h`
      - L2649 Initialization — `hydrogThetaFile`, `hydrogSaltFile`, `tRef`, `sRef`, `maskIniTemp`, `maskIniSalt`, `checkIniTemp`, `checkIniSalt`
      - L2678 Tracer Diffusivities — `diffKhT`, `diffKhS`, `diffKrT`, `diffKrS`, `diffK4T`, `diffK4S`, `diffKr4T`, `diffKrNrT`, `pCellMix_diffKr`, `diffKrNr`, `diffKr4S`, `diffKrNrS`, `diffKrFile`, `ALLOW_3D_DIFFKR`, `diffKrBL79surf`, `diffKrBL79deep`, `diffKrBL79scl`, `diffKrBL79Ho`, `diffKrBLEQsurf`, `ALLOW_BL79_LAT_VARY`, `diffKrBLEQdeep`, `diffKrBLEQscl`, `diffKrBLEQHo`, `BL79LatVary`
      - L2745 Ocean Convection — `cadjFreq`, `implicitDiffusion`, `ivdc_kappa`, `cAdjFreq`, `deltaTclock`, `hMixCriteria`, `hMixSmooth`
    - L2778 Parameters: Model Forcing
      - L2787 Momentum Forcing — `zonalWindFile`, `meridWindFile`, `momForcing`, `momForcingOutAB`, `momTidalForcing`, `ploadFile`
      - L2823 Tracer Forcing — `surfQnetfile`, `thetaClimFile`, `tauThetaClimRelax`, `EmPmRfile`, `saltClimFile`, `tauSaltClimRelax`, `tempForcing`, `surfQnetFile`, `surfQswFile`, `SHORTWAVE_HEATING`, `lambdaThetaFile`, `ThetaClimFile`, `balanceThetaClimRelax`, `ALLOW_BALANCE_RELAX`, `balanceQnet`, `ALLOW_BALANCE_FLUXES`, `geothermalFile`, `ALLOW_GEOTHERMAL_FLUX`, `temp_EvPrRn`, `allowFreezing`, `saltForcing`, `convertFW2Salt`, `useRealFreshWaterFlux`, `rhoConstFresh`, `rhoConst`, `EmPmRFile`, `saltFluxFile`, `lambdaSaltFile`, `balanceSaltClimRelax`, `selectBalanceEmPmR`, `wghtBalancedFile`, `wghtBalanceFile`, `salt_EvPrRn`, `selectAddFluid`, `ALLOW_ADDFLUID`, `temp_addMass`, `salt_addMass`, `addMassFile`, `balancePrintMean`, `latBandClimRelax` (+1)
      - L2923 Periodic Forcing — `periodicExternalForcing`, `externForcingPeriod`, `externForcingCycle`
    - L2952 Parameters: Simulation Controls
      - L2955 Run Start and Duration — `startTime`, `nIter0`, `endTime`, `nTimeSteps`, `deltaTClock`, `nEndIter`, `baseTime`
      - L2985 Input/Output Files — `readBinaryPrec`, `writeBinaryPrec`, `globalFiles`, `useSingleCpuIO`, `the_run_name`, `outputTypesInclusive`, `rwSuffixType`, `myTime`, `mdsioLocalDir`
      - L3039 Frequency/Amount of Output — `dumpFreq`, `monitorFreq`, `deltaTClock`, `dumpInitAndLast`, `diagFreq`, `monitorSelect`, `debugLevel`, `debugMode`, `plotLevel`
      - L3073 Restart/Pickup Files — `chkPtFreq`, `pchkPtFreq`, `deltaTClock`, `pChkPtFreq`, `pickupSuff`, `nIter0`, `pickupStrictlyMatch`, `writePickupAtEnd`, `usePickupBeforeC54`, `startFromPickupAB2`, `ALLOW_ADAMSBASHFORTH_3`
    - L3101 Parameters Used In Optional Packages
      - L3110 C-D Scheme — `tauCD`, `useCDscheme`, `deltaTMom`, `rCD`, `epsAB_CD`, `abEps`
      - L3134 Automatic Differentiation — `nTimeSteps_l2`, `adjdumpFreq`, `adjMonitorFreq`, `adTapeDir`
    - L3155 Execution Environment Parameters — `nTx`, `nTy`, `useCubedSphereExchange`, `debugMode`, `debugLevel`, `useCoupler`, `useSETRLSTK`, `useSIGREG`, `printMapIncludesZeros`, `maxLengthPrt1D`
  - L3201 MITgcm Input Data File Format — `delXFile`, `tauThetaClimRelax`, `viscAhZfile`, `readBinaryPrec`

## `doc/index.rst` — Welcome to MITgcm's user manual
- L6 Welcome to MITgcm's user manual

## `doc/ocean_state_est/ocean_state_est.rst` — Packages III - Ocean State Estimation
- L3 Packages III - Ocean State Estimation
  - L13 ECCO: model-data comparisons using gridded data sets
    - L105 Generic Cost Function — `gencost_datafile`, `gencost_barfile`, `gencost_outputlevel`, `gencost_name`, `gencost_errfile`, `gencost_avgperiod`, `gencost_preproc`, `gencost_preproc_i`, `gencost_posproc`, `gencost_posproc_c`, `gencost_posproc_i`, `gencost_kLev_select`, `gencost_is3d`, `gencost_mask`, `mult_gencost`, `gencost_preproc_c`, `gencost_preproc_r`, `gencost_posproc_r`, `gencost_spmin`, `gencost_spmax`, `gencost_spzero`, `gencost_startdate1`, `gencost_startdate2`, `gencost_enddate1`, `gencost_enddate2`, `gencost_useDensityMask`, `gencost_sigmaLow`, `gencost_sigmaHigh`, `gencost_refPressure`, `gencost_tanhScale`, `ATMOSPHERIC_LOADING`, `ALLOW_PSBAR_STERIC`, `m_eta`, `m_sst`, `m_sss`, `m_bp`, `m_siarea`, `m_siheff`, `m_sihsnow`, `m_theta` (+19)
    - L435 Generic Integral Function — `gencost_barfile`, `gencost_mask`, `gencost_avgperiod`, `gencost_useDensityMask`, `gencost_sigmaLow`, `gencost_sigmaHigh`, `gencost_refPressure`, `foo_maskC`, `foo_maskK`, `foo_maskT`, `foo_maskW`, `foo_maskS`, `m_boxmean`, `eccoVol_0`, `ECCO_VARIABLE_AREAVOLGLOB`, `m_boxmean_theta`, `m_boxmean_salt`, `m_boxmean_eta`, `m_boxmean_shifwf`, `m_boxmean_shihf`, `m_boxmean_vol`, `m_horflux_vol`, `cost_gencost_boxmean.F`
    - L521 Custom Cost Functions — `gencost_barfile`, `gencost_name`, `m_trVol`, `m_trHeat`, `m_trSalt`, `m_horflux_vol`, `ENUM_CENTERED_2ND`, `transp_trVol`, `transp_trHeat`, `transp_trSalt`, `moc_trVol`, `cost_gencost_bpv4.F`, `cost_gencost_seaicev4.F`, `cost_gencost_sshv4.F`, `cost_gencost_sstv4.F`, `cost_gencost_transp.F`, `cost_gencost_moc.F`
    - L622 Key Routines — `ecco_readparms.F`, `ecco_check.F`, `ecco_summary.F`, `cost_generic.F`, `cost_gencost_boxmean.F`, `ecco_toolbox.F`, `ecco_phys.F`, `cost_gencost_customize.F`, `cost_averagesfields.F`
    - L636 Compile Options — `ALLOW_GENCOST_CONTRIBUTION`, `ALLOW_GENCOST3D`, `ALLOW_PSBAR_STERIC`, `ALLOW_SHALLOW_ALTIMETRY`, `ALLOW_HIGHLAT_ALTIMETRY`, `ALLOW_PROFILES_CONTRIBUTION`, `ALLOW_ECCO_OLD_FC_PRINT`
  - L656 PROFILES: model-data comparisons at observed locations
  - L758 OBSFIT: grid-independent model-data comparisons
    - L763 Introduction
    - L779 Description
      - L793 Observations vs. Samples
      - L812 Sample types — `etaN`
      - L829 Observation duration
      - L846 Interpolation
      - L854 Cost Functions
    - L858 OBSFIT configuration and compiling — `OBSFIT_SIZE.h`
    - L875 Run-time requirements
      - L878 Pre-processing: How to make OBSFIT input files
      - L959 Enabling the package — `useOBSFIT`
    - L965 General flags and parameters — `obsfitDir`, `obsfitFiles`, `mult_obsfit`, `obsfit_facmod`, `obsfitDoNcOutput`
      - L1046 Post-processing
    - L1064 Experiments and tutorials that use OBSFIT
  - L1072 CTRL: Model Parameter Adjustment Capability
    - L1079 Generic Control Parameters — `ALLOW_GENTIM2D_CONTROL`, `ALLOW_GENARR2D_CONTROL`, `ALLOW_GENARR3D_CONTROL`, `EXCLUDE_CTRL_PACK`, `CTRL_OPTIONS.h`
      - L1101 Run-time Parameters — `ctrl_nml_genarr`, `maxCtrlArr2D`, `maxCtrlArr3D`, `maxCtrlTim2D`, `xx_gentim2d_period`, `xx_gentim2d_startdate1`, `startdate_1`, `xx_gentim2d_startdate2`, `startdate_2`, `xx_gentim2d_cumsum`, `xx_gentim2d_glosum`
      - L1150 Generic Control Fields — `ALLOW_DEPTH_CONTROL`, `xx_etan`, `xx_bottomdrag`, `xx_geothermal`, `xx_shicoefft`, `xx_shicoeffs`, `xx_shicdrag`, `xx_depth`, `xx_siheff`, `xx_siarea`, `xx_theta`, `xx_salt`, `xx_uvel`, `xx_vvel`, `xx_kapgm`, `xx_kapredi`, `xx_diffkr`, `xx_atemp`, `xx_aqh`, `xx_swdown`, `xx_lwdown`, `xx_precip`, `xx_runoff`, `xx_uwind`, `xx_vwind`, `xx_tauu`, `xx_tauv`, `xx_gen_precip`, `xx_hflux`, `xx_sflux`, `xx_shifwflx`
      - L1256 Generic Control Processing Options — `xx_gentim2d_period`, `xx_gentim2d_weight`, `xx_gentim2d_cumsum`, `xx_gentim2d_glosum`, `xx_gentim2d`, `adctrl_bound.F`
      - L1370 Generic Control Record Access — `docycle`, `rmcycle`, `noscaling`, `ALLOW_GENTIM2D_CONTROL`, `startrec`, `endrec`, `diffrec`, `startdate_1`, `startdate_2`, `xx_gentim2d_startdate1`, `xx_gentim2d_startdate2`, `nIter0`, `yadprefix`, `xx_gentim2d_startdate`, `xx_gentim2d`, `packages_init_fixed.F`, `ctrl_init.F`, `ctrl_init_rec.F`, `ctrl_init_ctrlvar.F`, `the_model_main.F`, `the_main_loop.F`, `initialise_fixed.F`, `ini_parms.F`, `packages_boot.F`, `packages_readparms.F`, `set_parms.F`, `ini_model_io.F`, `ini_grid.F`, `load_ref_files.F`, `ini_eos.F`, `set_ref_state.F`, `set_grid_factors.F`, `ini_depths.F`, `ini_masks_etc.F`, `cal_init_fixed.F`, `diagnostics_init_early.F`, `diagnostics_main_init.F`, `gad_init_fixed.F`, `mom_init_fixed.F`, `obcs_init_fixed.F` (+24)
    - L1534 Shelfice Control Parameters — `SHI_ALLOW_GAMMAFRICT`, `SHELFICEuseGammaFrict`, `SHELFICEsaltToHeatRatio`, `xx_shicoefft`, `xx_shicoeffs`, `xx_shicdrag`, `SHELFICE_OPTIONS.h`
    - L1565 Logarithmic Control Parameters — `log10InitVal`
  - L1599 SMOOTH: Smoothing And Covariance Model
  - L1608 The line search optimisation algorithm
    - L1613 General features
    - L1621 The online vs. offline version
    - L1650 Number of iterations vs. number of simulations
      - L1656 Summary
      - L1665 Description
      - L1713 The parameter file lsopt.par
      - L1744 OPWARMI, OPWARMD files
      - L1805 Error handling
    - L1978 Alternative code to :filelink:`optim` and :filelink:`lsopt` — `simul_rc`, `m1qn3_offline`, `optim_sub`
  - L2042 Test Cases For Estimation Package Capabilities

## `doc/outp_pkgs/flt.rst` — Introduction
- L2 Introduction
- L17 Compile-time options in `FLT_OPTIONS.h` — `ALLOW_3D_FLT`, `USE_FLT_ALT_NOISE`, `ALLOW_FLT_3D_NOISE`, `FLT_SECOND_ORDER_RUNGE_KUTTA`, `FLT_WITHOUT_X_PERIODICITY`, `FLT_WITHOUT_Y_PERIODICITY`, `DEVEL_FLT_EXCH2`
- L41 Compile-time parameters in `FLT_SIZE.h` include:
- L56 Run-time options in `data.flt` include:
- L95 Input Files
- L135 Output Files
- L144 Verification Experiment
- L159 Algorithm details

## `doc/outp_pkgs/outp_pkgs.rst` — Packages II - Diagnostics and I/O
- L4 Packages II - Diagnostics and I/O
  - L13 pkg/diagnostics – A Flexible Infrastructure
    - L16 Introduction
    - L48 Equations
    - L53 Key Subroutines and Parameters — `gdiag`, `diagnostics_addtolist.F`, `diagnostics_fill.F`, `diagnostics_scale_fill.F`, `diagnostics_fract_fill.F`, `diagnostics_frac_fill.F`, `diagnostics_is_on.F`, `diagnostics_count.F`, `DIAGNOSTICS.h`
    - L283 Usage Notes
      - L286 Using available diagnostics — `useDiagnostics`, `fileFlags`, `writeBinaryPrec`, `numDiags`, `numLists`, `numperList`, `diagSt_size`, `gdiag`, `DIAGNOSTICS_SIZE.h`
      - L434 Adjoint variables
      - L566 Adding new diagnostics to the code — `diagMate`, `diagNum`, `diagCode`, `diagnostics_fill.F`, `diagnostics_addtolist.F`, `diagnostics_main_init.F`
      - L638 MITgcm kernel available diagnostics list:
      - L880 MITgcm packages: available diagnostics lists
  - L920 Fortran Native I/O: pkg/mdsio and pkg/rw
    - L923 pkg/mdsio
      - L926 Introduction
      - L933 Using pkg/mdsio — `precFloat64`, `precFloat32`, `globalFile`
      - L1030 Important considerations — `SAFE_IO`, `ALLOW_WHIO`, `_BYTESWAPIO`
    - L1141 pkg/rw basic binary I/O utilities
      - L1148 Introduction — `RW_SAFE_MFLDS`, `RW_DISABLE_SMALL_OVERLAP`
  - L1187 NetCDF I/O: pkg/mnc
    - L1216 Using pkg/mnc
      - L1219 pkg/mnc configuration: — `MNC_COMMON.h`
      - L1256 pkg/mnc Inputs: — `useMNC`, `outputTypesInclusive`, `mnc_use_outdir`, `mnc_outdir_str`, `mnc_outdir_date`, `mnc_outdir_num`, `pickup_write_mnc`, `pickup_read_mnc`, `mnc_use_indir`, `mnc_indir_str`, `snapshot_mnc`, `monitor_mnc`, `timeave_mnc`, `autodiff_mnc`, `mnc_max_fsize`, `mnc_filefreq`, `readgrid_mnc`, `mnc_echo_gvtypes`
      - L1375 pkg/mnc output: — `dumpFreq`
    - L1465 pkg/mnc Troubleshooting
      - L1468 Build troubleshooting:
      - L1487 Runtime troubleshooting:
    - L1526 pkg/mnc Internals
      - L1572 pkg/mnc grid–tTypes and variable–types: — `DYNVARS.h`
      - L1639 Using pkg/mnc: examples — `ini_model_io.F`, `write_state.F`, `mom_vecinv.F`
  - L1721 Monitor: Simulation State Monitoring Toolkit
    - L1724 Introduction
    - L1763 Using pkg/monitor — `MONITOR_TEST_HFACZ`, `outputTypesInclusive`
  - L1790 Grid Generation
    - L1824 Using SPGrid
      - L1855 SPGrid requirements
      - L1875 Obtaining SPGrid
      - L1881 Building SPGrid
      - L1901 Running SPGrid — `SpF_test_cube_cap`
    - L1933 Example Grids
  - L1938 Pre– and Post–Processing Scripts and Utilities
    - L1949 Utilities Supplied With the Model
      - L1955 utils/scripts
      - L1965 utils/matlab
      - L1988 pkg/mnc utils
    - L2016 Pre-Processing Software
  - L2024 Potential Vorticity Matlab Toolbox
    - L2029 Introduction
    - L2045 Equations
      - L2048 Potential vorticity
      - L2096 Surface vertical potential vorticity fluxes
    - L2198 Key routines — `splPV`
    - L2265 Technical details
      - L2268 File name — `netcdf_UVEL`, `netcdf_domain`, `netcdf_suff`
      - L2305 Path to file
      - L2325 Grids
    - L2336 Notes on the flux form of the PV equation and vertical PV fluxes
      - L2339 Flux form of the PV equation
      - L2397 Determining the PV flux at the ocean’s surface
  - L2561 pkg/flt – Simulation of float / parcel displacements

## `doc/overview/adjoint.rst` — Adjoint
- L1 Adjoint

## `doc/overview/atmosphere.rst` — Atmosphere
- L1 Atmosphere

## `doc/overview/bound_forc_inter_waves.rst` — Boundary forced internal waves
- L2 Boundary forced internal waves

## `doc/overview/coordinate_sys.rst` — Coordinate systems
- L3 Coordinate systems
  - L6 Spherical coordinates

## `doc/overview/cvct_mixing_topo.rst` — Convection and mixing over topography
- L1 Convection and mixing over topography

## `doc/overview/eqn_motion_ocn.rst` — Equations of Motion for the Ocean
- L3 Equations of Motion for the Ocean
  - L91 Compressible z-coordinate equations
  - L135 ‘Anelastic’ z-coordinate equations
  - L221 Incompressible z-coordinate equations
  - L257 Compressible non-divergent equations

## `doc/overview/finding_pressure.rst` — Finding the pressure field
- L3 Finding the pressure field
  - L14 Hydrostatic pressure
  - L36 Surface pressure
  - L77 Non-hydrostatic pressure
    - L96 Boundary Conditions

## `doc/overview/forcing_dissip.rst` — Forcing/dissipation
- L1 Forcing/dissipation
  - L4 Forcing
  - L11 Dissipation
    - L14 Momentum
    - L30 Tracers

## `doc/overview/global_atmos_hs.rst` — Global atmosphere: ‘Held-Suarez’ benchmark
- L1 Global atmosphere: ‘Held-Suarez’ benchmark

## `doc/overview/global_ocean_circ.rst` — Global ocean circulation
- L2 Global ocean circulation

## `doc/overview/global_state_est.rst` — Global state estimation of the ocean
- L1 Global state estimation of the ocean

## `doc/overview/hydro_prim_eqn.rst` — Hydrostatic Primitive Equations for the Atmosphere in Pressure Coordinates
- L3 Hydrostatic Primitive Equations for the Atmosphere in Pressure Coordinates
  - L92 Boundary conditions
  - L111 Splitting the geopotential

## `doc/overview/hydrostatic.rst` — Hydrostatic, Quasi-hydrostatic, Quasi-nonhydrostatic and Non-hydrostatic forms
- L3 Hydrostatic, Quasi-hydrostatic, Quasi-nonhydrostatic and Non-hydrostatic forms
  - L77 Shallow atmosphere approximation
  - L92 Hydrostatic and quasi-hydrostatic forms
  - L125 Non-hydrostatic and quasi-nonhydrostatic forms
    - L131 Non-hydrostatic Ocean
    - L145 Quasi-nonhydrostatic Atmosphere
  - L156 Summary of equation sets supported by model
    - L159 Atmosphere
      - L166 Hydrostatic and quasi-hydrostatic
      - L174 Quasi-nonhydrostatic
    - L179 Ocean
      - L182 Hydrostatic and quasi-hydrostatic
      - L188 Non-hydrostatic

## `doc/overview/kinematic_bound.rst` — Kinematic Boundary conditions
- L1 Kinematic Boundary conditions
  - L4 Vertical
  - L26 Horizontal

## `doc/overview/ocean.rst` — Ocean
- L1 Ocean

## `doc/overview/ocean_biogeo_cyc.rst` — Ocean biogeochemical cycles
- L1 Ocean biogeochemical cycles

## `doc/overview/ocean_gyres.rst` — Ocean gyres
- L2 Ocean gyres

## `doc/overview/overview.rst` — Overview
- L1 Overview
  - L14 Introduction
  - L117 Illustrations of the model in action
  - L145 Continuous equations in ‘r’ coordinates
  - L276 Appendix ATMOSPHERE
  - L285 Appendix OCEAN
  - L294 Appendix OPERATORS

## `doc/overview/parm_sens.rst` — Parameter sensitivity using the adjoint of MITgcm
- L1 Parameter sensitivity using the adjoint of MITgcm

## `doc/overview/sim_lab_exp.rst` — Simulations of laboratory experiments
- L1 Simulations of laboratory experiments

## `doc/overview/soln_strategy.rst` — Solution strategy
- L1 Solution strategy

## `doc/overview/vector_invar.rst` — Vector invariant form
- L1 Vector invariant form

## `doc/phys_pkgs/aim.rst` — Atmospheric Intermediate Physics: AIM
- L3 Atmospheric Intermediate Physics: AIM — `aim_v23`
  - L10 Key subroutines, parameters and files
  - L15 AIM Diagnostics
  - L61 Experiments and tutorials that use aim

## `doc/phys_pkgs/bulk_force.rst` — BULK_FORCE: Bulk Formula Package
- L3 BULK_FORCE: Bulk Formula Package
  - L31 subroutine BULKF_FIELDS_LOAD
  - L53 subroutine BULKF_FORCING
  - L88 subroutine BULKF_FORMULA_LANL
  - L213 Initializing subroutines
  - L221 Diagnostic subroutines
  - L229 Common Blocks
  - L241 Input file DATA.ICE
  - L250 Important Notes
  - L258 References
  - L270 Experiments and tutorials that use bulk\_force

## `doc/phys_pkgs/cal.rst` — CAL: The calendar package
- L3 CAL: The calendar package
  - L19 Basic assumptions for the calendar tool
  - L31 Format of calendar dates
  - L64 Calendar dates and time intervals
  - L82 Using the calendar together with MITgcm
  - L125 The individual calendars
  - L151 Short routine description
  - L259 Experiments and tutorials that use cal

## `doc/phys_pkgs/dic.rst` — DIC Package
- L3 DIC Package
  - L6 Introduction
  - L25 Key subroutines and parameters
  - L123 Do’s and Don’ts
  - L130 Reference Material
  - L146 Experiments and tutorials that use dic

## `doc/phys_pkgs/exch2.rst` — exch2: Extended Cubed Sphere Topology
- L3 exch2: Extended Cubed Sphere Topology
  - L8 Introduction — `W2_EXCH2_TOPOLOGY.h`, `w2_e2setup.F`
  - L33 Invoking exch2
  - L60 Generating Topology Files for exch2
  - L129 exch2, SIZE.h, and Multiprocessing — `sNx`, `sNy`, `OLx`, `OLy`, `Nr`, `nSx`, `nSy`, `nPx`, `nPy`
  - L216 Key Variables
    - L227 Scalars: — `exch2_nTiles`, `W2_maxNeighbours`
    - L250 Arrays indexed to tile number: — `exch2_tnx`, `exch2_tny`, `exch2_tbasex`, `exch2_tbasey`, `exch2_txglobalo`, `exch2_myFace`, `exch2_nNeighbours`, `exch2_tProc`, `exch2_isWedge`, `exch2_isEedge`, `exch2_isSedge`, `exch2_isNedge`
    - L302 Arrays Indexed to Tile Number and Neighbor: — `W2_maxNeighbours`, `exch2_nTiles`, `exch2_pi`, `exch2_pj`, `exch2_oi`, `exch2_oj`, `exch2_oi_f`, `exch2_oj_f`, `exch2_itlo_c`, `exch2_ithi_c`, `exch2_jtlo_c`, `exch2_jthi_c`
  - L453 Key Routines
  - L485 Experiments and tutorials that use exch2

## `doc/phys_pkgs/exf.rst` — EXF: The external forcing package
- L3 EXF: The external forcing package
  - L10 Introduction
  - L33 EXF configuration, compiling & running
    - L36 Compile-time options
  - L103 Run-time parameters
    - L112 Enabling the package
    - L118 General flags and parameters — `_YYYY`
    - L225 Field attributes — `exf_inscal_`, `exf_outscal_`, `EXF_USE_INTERPOLATION`, `_lon0`, `_lon_inc`, `_lat0`, `_lat_inc`, `_nlon`, `_nlat`
    - L289 Example configuration — `global_oce_latlon`
  - L325 EXF bulk formulae
  - L332 EXF input fields and units
  - L470 Key subroutines
  - L553 EXF diagnostics
  - L585 References
  - L588 Experiments and tutorials that use exf

## `doc/phys_pkgs/fizhi.rst` — Fizhi: High-end Atmospheric Physics
- L5 Fizhi: High-end Atmospheric Physics
  - L9 Introduction
  - L20 Equations
    - L28 Sub-grid and Large-scale Convection
    - L149 Cloud Formation
    - L246 Shortwave Radiation
    - L339 Longwave Radiation
    - L423 Cloud-Radiation Interaction
    - L471 Turbulence
    - L677 Atmospheric Boundary Layer
    - L686 Surface Energy Budget
    - L748 Surface Type
    - L809 Surface Roughness
    - L818 Albedo
    - L831 Gravity Wave Drag
    - L869 Boundary Conditions and other Input Data
    - L906 Topography and Topography Variance
    - L921 Upper Level Moisture
  - L937 Fizhi Diagnostics
  - L1192 Fizhi Diagnostic Description
    - L1207 Surface Zonal Wind Stress on the Atmosphere (:math:`Newton/m^2`)
    - L1221 Surface Meridional Wind Stress on the Atmosphere (:math:`Newton/m^2`)
    - L1235 Surface Flux of Sensible Heat (W m\ :sup:`--2`)
    - L1256 Surface Flux of Latent Heat (:math:`Watts/m^2`)
    - L1279 Heat Conduction Through Sea Ice (:math:`Watts/m^2`)
    - L1298 Net upward Longwave Flux at the surface (:math:`Watts/m^2`)
    - L1312 Net downard shortwave Flux at the surface (:math:`Watts/m^2`)
    - L1325 Richardson number (:math:`dimensionless`)
    - L1344 CT - Surface Exchange Coefficient for Temperature and Moisture (dimensionless)
    - L1378 CU - Surface Exchange Coefficient for Momentum (dimensionless)
    - L1398 ET - Diffusivity Coefficient for Temperature and Moisture (m^2/sec)
    - L1440 EU - Diffusivity Coefficient for Momentum (m^2/sec)
    - L1483 TURBU - Zonal U-Momentum changes due to Turbulence (m/sec/day)
    - L1497 TURBV - Meridional V-Momentum changes due to Turbulence (m/sec/day)
    - L1512 TURBT - Temperature changes due to Turbulence (deg/day)
    - L1528 TURBQ - Specific Humidity changes due to Turbulence (g/kg/day)
    - L1543 MOISTT - Temperature Changes Due to Moist Processes (deg/day)
    - L1571 MOISTQ - Specific Humidity Changes Due to Moist Processes (g/kg/day)
    - L1600 RADLW - Heating Rate due to Longwave Radiation (deg/day)
    - L1631 RADSW - Heating Rate due to Shortwave Radiation (deg/day)
    - L1662 PREACC - Total (Large-scale + Convective) Accumulated Precipition (mm/day)
    - L1678 PRECON - Convective Precipition (mm/day)
    - L1693 TUFLUX - Turbulent Flux of U-Momentum (Newton/m^2)
    - L1708 TVFLUX - Turbulent Flux of V-Momentum (Newton/m^2)
    - L1724 TTFLUX - Turbulent Flux of Sensible Heat (Watts/m^2)
    - L1741 TQFLUX - Turbulent Flux of Latent Heat (Watts/m^2)
    - L1757 CN - Neutral Drag Coefficient (dimensionless)
    - L1768 WINDS - Surface Wind Speed (meter/sec)
    - L1797 TG - Ground Temperature (deg K)
    - L1830 TS - Surface Temperature (deg K)
    - L1840 DTG - Surface Temperature Adjustment (deg K)
    - L1854 QG - Ground Specific Humidity (g/kg)
    - L1870 QS - Saturation Surface Specific Humidity (g/kg)
    - L1878 TGRLW - Instantaneous ground temperature used as input to the Longwave radiation subroutine (deg)
    - L1886 ST4 - Upward Longwave flux at the surface (Watts/m^2)
    - L1895 OLR - Net upward Longwave flux at :math:`p=p_{top}` (Watts/m^2)
    - L1904 OLRCLR - Net upward clearsky Longwave flux at :math:`p=p_{top}` (Watts/m^2)
    - L1913 LWGCLR - Net upward clearsky Longwave flux at the surface (Watts/m^2)
    - L1928 LWCLR - Heating Rate due to Clearsky Longwave Radiation (deg/day)
    - L1959 TLW - Instantaneous temperature used as input to the Longwave radiation subroutine (deg)
    - L1968 SHLW - Instantaneous specific humidity used as input to the Longwave radiation subroutine (kg/kg)
    - L1977 OZLW - Instantaneous ozone used as input to the Longwave radiation subroutine (kg/kg)
    - L1986 CLMOLW - Maximum Overlap cloud fraction used in LW Radiation (0-1)
    - L1999 CLDTOT - Total cloud fraction used in LW and SW Radiation (0-1)
    - L2016 CLMOSW - Maximum Overlap cloud fraction used in SW Radiation (0-1)
    - L2029 CLROSW - Random Overlap cloud fraction used in SW Radiation (0-1)
    - L2042 RADSWT - Incident Shortwave radiation at the top of the atmosphere (Watts/m^2)
    - L2055 EVAP - Surface Evaporation (mm/day)
    - L2072 DUDT - Total Zonal U-Wind Tendency  (m/sec/day)
    - L2081 DVDT - Total Zonal V-Wind Tendency  (m/sec/day)
    - L2090 DTDT - Total Temperature Tendency  (deg/day)
    - L2103 DQDT - Total Specific Humidity Tendency  (g/kg/day)
    - L2115 USTAR -  Surface-Stress Velocity (m/sec)
    - L2131 Z0 - Surface Roughness Length (m)
    - L2146 FRQTRB - Frequency of Turbulence (0-1)
    - L2156 PBL - Planetary Boundary Layer Depth (mb)
    - L2170 SWCLR - Clear sky Heating Rate due to Shortwave Radiation (deg/day)
    - L2201 OSR - Net upward Shortwave flux at the top of the model (Watts/m^2)
    - L2210 OSRCLR - Net upward clearsky Shortwave flux at the top of the model (Watts/m^2)
    - L2219 CLDMAS - Convective Cloud Mass Flux (kg/m^2)
    - L2233 UAVE - Time-Averaged Zonal U-Wind (m/sec)
    - L2245 VAVE - Time-Averaged Meridional V-Wind (m/sec)
    - L2258 TAVE - Time-Averaged Temperature (Kelvin)
    - L2268 QAVE - Time-Averaged Specific Humidity (g/kg)
    - L2279 PAVE - Time-Averaged Surface Pressure - PTOP (mb)
    - L2293 QQAVE - Time-Averaged Turbulent Kinetic Energy (m/sec)^2
    - L2308 SWGCLR - Net downward clearsky Shortwave flux at the surface (Watts/m^2)
    - L2324 DIABU - Total Diabatic Zonal U-Wind Tendency  (m/sec/day)
    - L2334 DIABV - Total Diabatic Meridional V-Wind Tendency  (m/sec/day)
    - L2343 DIABT Total Diabatic Temperature Tendency (deg/day)
    - L2374 DIABQ - Total Diabatic Specific Humidity Tendency (g/kg/day)
    - L2397 VINTUQ - Vertically Integrated Moisture Flux (m/sec  g/kg)
    - L2413 VINTVQ - Vertically Integrated Moisture Flux (m/sec g/kg)
    - L2429 VINTUT - Vertically Integrated Heat Flux (m/sec deg)
    - L2443 VINTVT - Vertically Integrated Heat Flux (m/sec deg)
    - L2457 CLDFRC - Total 2-Dimensional Cloud Fracton (0-1)
    - L2509 QINT - Total Precipitable Water (gm/cm^2)
    - L2526 U2M  Zonal U-Wind at 2 Meter Depth (m/sec)
    - L2543 V2M - Meridional V-Wind at 2 Meter Depth (m/sec)
    - L2560 T2M - Temperature at 2 Meter Depth (deg K)
    - L2583 Q2M - Specific Humidity at 2 Meter Depth (g/kg)
    - L2606 U10M - Zonal U-Wind at 10 Meter Depth (m/sec)
    - L2623 V10M - Meridional V-Wind at 10 Meter Depth (m/sec)
    - L2640 T10M - Temperature at 10 Meter Depth (deg K)
    - L2664 Q10M - Specific Humidity at 10 Meter Depth (g/kg)
    - L2688 DTRAIN - Cloud Detrainment Mass Flux (kg/m^2)
    - L2700 QFILL - Filling of negative Specific Humidity (g/kg/day)
  - L2717 Key subroutines, parameters and files
  - L2721 Dos and don'ts
  - L2725 Fizhi Reference
  - L2729 Experiments and tutorials that use fizhi

## `doc/phys_pkgs/gchem.rst` — GCHEM Package
- L3 GCHEM Package
  - L6 Introduction
  - L25 Key subroutines and parameters
  - L114 GCHEM Diagnostics
  - L133 Do’s and Don’ts
  - L142 Reference Material
  - L145 Experiments and tutorials that use gchem

## `doc/phys_pkgs/generic_advdiff.rst` — Generic Advection/Diffusion
- L3 Generic Advection/Diffusion
  - L12 Introduction — `gad_advection.F`
  - L26 Key subroutines, parameters and files — `COSINEMETH_III`, `ISOTROPIC_COS_SCALING`, `DISABLE_MULTIDIM_ADVECTION`, `GAD_MULTIDIM_COMPRESSIBLE`, `GAD_ALLOW_TS_SOM_ADV`, `GAD_SMOLARKIEWICZ_HACK`, `gad_calc_rhs.F`, `gad_advection.F`
  - L63 GAD Diagnostics
  - L86 Experiments and tutorials that use GAD

## `doc/phys_pkgs/ggl90.rst` — GGL90: a TKE vertical mixing scheme
- L3 GGL90: a TKE vertical mixing scheme
  - L11 Key subroutines, parameters and files
  - L18 Experiments and tutorials that use GGL90

## `doc/phys_pkgs/gmredi.rst` — GMREDI: Gent-McWilliams/Redi Eddy Parameterization
- L3 GMREDI: Gent-McWilliams/Redi Eddy Parameterization
  - L6 Introduction
  - L31 Description
    - L59 Redi scheme: Isopycnal diffusion
    - L102 GM parameterization — `GM_AdvForm`, `GM_BOLUS_ADVEC`, `GM_EXTRA_DIAGONAL`, `uVel`, `vVel`, `wVel`, `GM_AdvSeparate`, `GM_PsiX`, `GM_PsiY`
    - L186 Griffies Skew Flux — `gmredi_calc_tensor.F`
    - L275 Redi and GM schemes in pressure coordinate
    - L347 Visbeck et al. 1997 GM diffusivity :math:`\kappa_{GM}(x,y)`
    - L376 Marshall et al. 2012 GM diffusivity :math:`\kappa_{GM}(x,y)` — `GEOM_alpha`, `GEOM_vert_struc`, `GEOM_lmbda`, `GEOM_diffKh_EKE`, `pickup_gmredi`
    - L426 Tapering and stability
      - L435 Slope clipping — `GM_taper_scheme`, `gmredi_slope_limit.F`
      - L501 Tapering: Gerdes, Koberle and Willebrand, 1991 (GKW91) — `GM_taper_scheme`
      - L542 Tapering: Danabasoglu and McWilliams, 1995 (DM95) — `GM_taper_scheme`
      - L559 Tapering: Large, Danabasoglu and Doney, 1997 (LDD97) — `GM_taper_scheme`
  - L580 GMREDI configuration and compiling
    - L583 Compile-time options — `GM_NON_UNITY_DIAGONAL`, `GM_EXTRA_DIAGONAL`, `GM_BOLUS_ADVEC`, `GM_BOLUS_BVP`, `ALLOW_GM_LEITH_QG`, `GM_VISBECK_VARIABLE_K`, `GM_GEOM_VARIABLE_K`, `GMREDI_OPTIONS.h`
  - L620 Run-time parameters — `gmredi_readparms.F`
    - L626 Enabling the package — `useGMREDI`
    - L632 General flags and parameters — `GM_AdvForm`, `GM_AdvSeparate`, `GM_background_K`, `GM_isopycK`, `GM_maxSlope`, `GM_Kmin_horiz`, `GM_Small_Number`, `GM_slopeSqCutoff`, `GM_taper_scheme`, `GM_maxTransLay`, `GM_facTrL2ML`, `GM_facTrL2dz`, `GM_Scrit`, `GM_Sd`, `GM_UseBVP`, `GM_BVP_ModeNumber`, `GM_BVP_cMin`, `GM_UseSubMeso`, `subMeso_Ceff`, `subMeso_invTau`, `subMeso_LfMin`, `subMeso_Lmax`, `GM_Visbeck_alpha`, `GM_Visbeck_length`, `GM_Visbeck_depth`, `GM_Visbeck_maxSlope`, `GM_Visbeck_minVal_K`, `GM_Visbeck_maxVal_K`, `GM_useGEOM`, `GM_GEOM_VARIABLE_K`, `GEOM_alpha`, `GEOM_lmbda`, `GEOM_diffKh_EKE`, `GEOM_ini_EKE`, `GEOM_vert_struc`, `GEOM_vert_struc_min`, `GEOM_vert_struc_max`, `GEOM_minVal_K`, `GEOM_maxVal_K`, `GM_useLeithQG` (+8)
  - L746 GMREDI Diagnostics
  - L788 Experiments and tutorials that use GMREDI

## `doc/phys_pkgs/gridalt.rst` — Gridalt - Alternate Grid Package
- L3 Gridalt - Alternate Grid Package
  - L7 Introduction
  - L72 Equations on Both Grids
  - L108 Time stepping Sequence
  - L125 Interpolation
  - L143 Key subroutines, parameters and files — `gridalt_initialise.F`, `make_phys_grid.F`, `gridalt_mapping.F`, `gridalt_update.F`
  - L256 Gridalt Diagnostics
  - L267 Dos and donts
  - L270 Gridalt Reference
  - L273 Experiments and tutorials that use gridalt

## `doc/phys_pkgs/kl10.rst` — KL10: Vertical Mixing Due to Breaking Internal Waves
- L3 KL10: Vertical Mixing Due to Breaking Internal Waves
  - L13 Introduction
  - L104 KL10 configuration and compiling
  - L125 Run-time parameters
    - L133 Enabling the package
    - L139 Required MITgcm flags
    - L151 Package flags and parameters — `dumpFreq`, `taveFreq`
  - L176 Equations and key routines
    - L179 KL10_CALC: — `viscAz`, `diffKzT`
    - L187 KL10_CALC_VISC:
    - L192 KL10_CALC_DIFF:
  - L199 KL10 diagnostics
  - L218 References
  - L224 Experiments and tutorials that use KL10

## `doc/phys_pkgs/kpp.rst` — KPP: Nonlocal K-Profile Parameterization for Vertical Mixing
- L3 KPP: Nonlocal K-Profile Parameterization for Vertical Mixing
  - L12 Introduction
  - L70 KPP configuration and compiling — `_KPP_RL`, `FRUGAL_KPP`, `KPP_SMOOTH_SHSQ`, `KPP_SMOOTH_DVSQ`, `KPP_SMOOTH_DENS`, `KPP_SMOOTH_VISC`, `KPP_SMOOTH_DIFF`, `KPP_ESTIMATE_UREF`, `INCLUDE_DIAGNOSTICS_INTERFACE_CODE`, `KPP_GHAT`, `EXCLUDE_KPP_SHEAR_MIX`
  - L126 Run-time parameters
    - L134 Enabling the package
    - L140 Required MITgcm flags
    - L152 Package flags and parameters — `deltaTClock`, `dumpFreq`, `taveFreq`
  - L252 Equations and key routines
    - L259 KPP_CALC:
    - L265 KPP_MIX:
    - L271 BLMIX: Mixing in the boundary layer
    - L340 RI\_IWMIX: Mixing in the interior
    - L353 BLDEPTH: Boundary layer depth calculation:
    - L377 KPP\_CALC\_DIFF\_T/\_S, KPP\_CALC\_VISC:
    - L385 KPP\_TRANSPORT\_T/\_S/\_PTR:
    - L394 Implicit time integration
    - L399 Penetration of shortwave radiation
  - L406 Flow chart
  - L436 KPP diagnostics
  - L456 Reference experiments
  - L463 References
  - L466 Experiments and tutorials that use kpp

## `doc/phys_pkgs/land.rst` — Land package
- L3 Land package
  - L7 Introduction
  - L26 Equations and Key Parameters
  - L83 Land diagnostics
  - L104 References
  - L111 Experiments and tutorials that use land

## `doc/phys_pkgs/mom_packages.rst` — Momentum Packages
- L1 Momentum Packages — `COSINEMETH_III`, `ISOTROPIC_COS_SCALING`, `ALLOW_SMAG_3D`, `ALLOW_3D_VISCAH`, `ALLOW_3D_VISCA4`, `ALLOW_BOTTOMDRAG_ROUGHNESS`, `MOM_BOUNDARY_CONSERVE`, `MOM_COMMON_OPTIONS.h`, `MOM_FLUXFORM_OPTIONS.h`

## `doc/phys_pkgs/obcs.rst` — OBCS: Open boundary conditions for regional modeling
- L3 OBCS: Open boundary conditions for regional modeling
  - L11 Introduction — `obcs_calc.F`
  - L29 OBCS configuration and compiling — `obcs_fields_load.F`, `obcs_prescribe_read.F`, `OBCS_OPTIONS.h`
  - L85 Run-time parameters — `packages_readparms.F`, `obcs_readparms.F`, `exf_readparms.F`
    - L102 Enabling the package — `useOBCS`
    - L108 Package flags and parameters — `OB_Jnorth`, `OB_Jsouth`, `OB_Ieast`, `OB_Iwest`, `useOBCSprescribe`, `useOBCSsponge`, `useOBCSbalance`, `OBCS_balanceFacN`, `OBCS_balanceFacS`, `OBCS_balanceFacE`, `OBCS_balanceFacW`, `OBCSbalanceSurf`, `useOrlanskiNorth`, `useOrlanskiSouth`, `useOrlanskiEast`, `useOrlanskiWest`, `useStevensNorth`, `useStevensSouth`, `useStevensEast`, `useStevensWest`, `cvelTimeScale`, `CMAX`, `CFIX`, `useFixedCEast`, `useFixedCWest`, `spongeThickness`, `Urelaxobcsinner`, `Vrelaxobcsinner`, `Urelaxobcsbound`, `Vrelaxobcsbound`, `TrelaxStevens`, `SrelaxStevens`, `useStevensPhaseVel`, `useStevensAdvection`
  - L200 Defining open boundary positions — `OB_Jnorth`, `OB_Jsouth`, `OB_Ieast`, `OB_Iwest`
    - L244 Simple examples — `OB_singleJnorth`, `OB_singleJsouth`, `OB_singleIeast`, `OB_singleIwest`
    - L276 A more complex example — `OB_Ieast`, `OB_Iwest`, `insideOBmaskFile`
  - L313 Equations and key routines — `OBSs`, `tRef`, `sRef`, `useOBCSprescribe`, `useSEAICE`, `HEFF`, `OBNu0`, `OBNu1`, `OBNu`, `OBCSWstartdate1`, `OBCSWstartdate2`, `OBCSWperiod`, `siobWstartdate1`, `siobWstartdate2`, `siobWperiod`, `externForcingPeriod`, `externForcingCycle`, `TrelaxStevens`, `SrelaxStevens`, `useStevensPhaseVel`, `useStevensAdvection`, `nonlinFreeSurf`, `ALLOW_OBCS_BALANCE`, `useOBCSbalance`, `OBCS_balanceFacN`, `OBCS_balanceFacS`, `OBCS_balanceFacE`, `OBCS_balanceFacW`, `OBWu`, `OBCSbalanceSurf`, `EmPmR`, `EXF_NML_OBCS`, `obcs_readparms.F`, `obcs_calc.F`, `ORLANSKI.h`, `external_fields_load.F`, `obcs_calc_stevens.F`, `dynamics.F`, `obcs_save_uv_n.F`, `obcs_balance_flow.F` (+1)
    - L544 OBCS\_APPLY\_*: — `ALLOW_OBCS_SPONGE`, `useOBCSsponge`, `spongeThickness`, `Urelaxobcsbound`, `Vrelaxobcsbound`, `Urelaxobcsinner`, `Vrelaxobcsinner`, `obcs_sponge.F`
    - L579 OB's with nonlinear free surface
    - L582 OB's with sea ice — `ALLOW_OBCS_SEAICE_SPONGE`, `useSeaiceSponge`, `seaiceSpongeThickness`, `SEAICEuseNeumannBC`, `OBCS_SEAICE_SMOOTH_EDGE`
  - L605 Flow chart
  - L660 OBCS diagnostics
  - L670 Experiments and tutorials that use obcs — `obcs_calc.F`

## `doc/phys_pkgs/opps.rst` — OPPS: Ocean Penetrative Plume Scheme
- L4 OPPS: Ocean Penetrative Plume Scheme
  - L12 Key subroutines, parameters and files
  - L20 Experiments and tutorials that use OPPS

## `doc/phys_pkgs/packages_overview.rst` — Using MITgcm Packages
- L3 Using MITgcm Packages — `gmRedi`
  - L21 Package Inclusion/Exclusion — `default_pkg_list`
  - L67 Package Activation — `usePackageName`
  - L85 Package Coding Standards
    - L91 Packages are Not Libraries
    - L99 File Inclusion Rules — `PackB`, `PackA`
    - L134 Conditional Compilation and ``PACKAGES_CONFIG.h`` — `ALLOW_GMREDI`, `useGMRedi`
    - L185 Package Startup or Boot Sequence
    - L244 Adding a package to PARAMS.h and packages\_boot() — `PARM_PACKAGES`

## `doc/phys_pkgs/phys_pkgs.rst` — Packages I - Physical Parameterizations
- L3 Packages I - Physical Parameterizations
  - L56 Overview
  - L65 Packages Related to Hydrodynamical Kernel
  - L78 General purpose numerical infrastructure packages
  - L88 Ocean Packages
  - L103 Atmosphere Packages
  - L113 Ice and Sea Ice Packages
  - L125 Biogeochemistry Packages

## `doc/phys_pkgs/ptracers.rst` — PTRACERS Package
- L3 PTRACERS Package
  - L7 Introduction
  - L23 Equations
  - L26 Key subroutines and parameters — `PTRACERS_num`, `PTRACERS_Iter0`, `PTRACERS_numInUse`, `PTRACERS_dumpFreq`, `dumpFreq`, `PTRACERS_taveFreq`, `taveFreq`, `PTRACERS_monitorFreq`, `monitorFreq`, `PTRACERS_timeave_mnc`, `useMNC`, `timeave_mnc`, `PTRACERS_snapshot_mnc`, `snapshot_mnc`, `PTRACERS_monitor_mnc`, `monitor_mnc`, `PTRACERS_pickup_write_mnc`, `pickup_write_mnc`, `PTRACERS_pickup_read_mnc`, `pickup_read_mnc`, `PTRACERS_useRecords`, `PTRACERS_advScheme`, `saltAdvScheme`, `PTRACERS_ImplVertAdv`, `PTRACERS_diffKh`, `diffKhS`, `PTRACERS_diffK4`, `diffK4S`, `PTRACERS_diffKr`, `PTRACERS_diffKrNr`, `diffKrNrS`, `PTRACERS_ref`, `PTRACERS_EvPrRn`, `convertFW2Salt`, `PTRACERS_useGMRedi`, `useGMREdi`, `PTRACERS_useKPP`, `useKPP`, `PTRACERS_initialFile`, `PTRACERS_names` (+5)
  - L66 PTRACERS Diagnostics
  - L148 Do’s and Don’ts
  - L151 Reference Material

## `doc/phys_pkgs/rbcs.rst` — RBCS Package
- L3 RBCS Package
  - L8 Introduction
  - L32 Key subroutines and parameters — `maskLEN`, `rbcsForcingPeriod`, `rbcsForcingCycle`, `rbcsForcingOffset`, `rbcsSingleTimeFiles`, `deltaTrbcs`, `deltaTclock`, `rbcsVanishingTime`, `myTime`, `rbcsIter0`, `useRBCtemp`, `useRBCsalt`, `useRBCuVel`, `useRBCvVel`, `tauRelaxT`, `tauRelaxS`, `tauRelaxU`, `tauRelaxV`, `relaxMaskFile`, `relaxMaskUFile`, `relaxMaskVFile`, `relaxTFile`, `relaxSFile`, `relaxUFile`, `relaxVFile`, `useRBCpTrNum`, `tauRelaxPTR`, `relaxPtracerFile`, `useRBCxxx`, `RBCS_SIZE.h`
  - L93 Timing of relaxation forcing fields — `rbcsForcingPeriod`, `rbcsSingleTimeFiles`, `rbcsForcingCycle`, `rbcsForcingOffset`, `rbcsIter0`, `deltaTrbcs`
  - L145 Example 1: forcing with time averages starting at :math:`t=0`
    - L148 Cyclic data in a single file — `rbcsSingleTimeFiles`, `rbcsForcingOffset`
    - L157 Non-cyclic data, multiple files — `rbcsForcingCycle`, `rbcsSingleTimeFiles`, `rbcsForcingOffset`, `rbcsIter0`, `deltaTrbcs`, `rbcsForcingPeriod`
  - L167 Example 2: forcing with snapshots starting at :math:`t=0`
    - L170 Cyclic data in a single file — `rbcsSingleTimeFiles`, `rbcsForcingOffset`
    - L177 Non-cyclic data, multiple files — `rbcsForcingCycle`, `rbcsSingleTimeFiles`, `rbcsForcingOffset`, `rbcsIter0`, `deltaTrbcs`, `rbcsForcingPeriod`
  - L188 Do’s and Don’ts
  - L191 Reference Material
  - L194 Experiments and tutorials that use rbcs

## `doc/phys_pkgs/remesh.rst` — SHELFICE Remeshing
- L3 SHELFICE Remeshing
  - L10 Introduction
  - L23 REMESHING configuration and compiling
    - L26 Compile-time options — `NONLIN_FRSURF`, `ALLOW_SHELFICE_REMESHING`, `SHI_ALLOW_GAMMAFRICT`, `SHI_withBL_uStarTopDz`, `CPP_OPTIONS.h`, `SHELFICE_OPTIONS.h`
  - L41 Run-time parameters — `nonlinFreeSurf`, `select_rstar`, `SHI_withBL_realFWflux`, `SHELFICEboundaryLayer`, `SHI_withBL_uStarTopDz`, `SHELFICEmassFile`, `SHELFICEMassStepping`, `SHELFICEMassDynTendFile`, `SHELFICEDynMassOnly`, `shelficeMass`, `SHELFICERemeshFrequency`, `SHELFICESplitThreshold`, `hFacC`, `SHELFICEMergeThreshold`
  - L75 Description — `implicitFreeSurface`, `nonlinFreeSurf`, `SHELFICERemeshFrequency`, `SHELFICESplitThreshold`, `SHELFICEMergeThreshold`, `useRealFreshWaterFlux`, `SHELFICEboundaryLayer`, `SHI_withBL_realFWflux`
  - L126 Alternate boundary layer formulation — `SHELFICEboundaryLayer`, `SHI_ALLOW_GAMMAFRICT`, `SHELFICEuseGammaFrict`, `SHI_withBL_uStarTopDz`
  - L139 Coupling with :filelink:`pkg/streamice` — `SHELFICEMassDynTendFile`, `pkg`, `conserve_ssh`, `OBCS_BALANCE_FLOW`, `useOBCSbalance`
  - L158 Diagnostics
  - L164 Experiments that use Remeshing

## `doc/phys_pkgs/seaice.rst` — SEAICE Package
- L3 SEAICE Package
  - L11 Introduction
  - L26 SEAICE configuration and compiling
    - L29 Compile-time options — `SEAICE_BGRID_DYNAMICS`, `SEAICE_DEBUG`, `SEAICE_CGRID`, `SEAICE_ALLOW_EVP`, `SEAICE_ALLOW_JFNK`, `SEAICE_ALLOW_KRYLOV`, `SEAICE_ALLOW_TEM`, `SEAICE_ALLOW_MCS`, `SEAICE_ALLOW_MCE`, `SEAICE_ALLOW_TD`, `SEAICE_LSR_ZEBRA`, `SEAICE_ALLOW_FREEDRIFT`, `SEAICE_EXTERNAL_FLUXES`, `SEAICE_ZETA_SMOOTHREG`, `SEAICE_DELTA_SMOOTHREG`, `SEAICE_ALLOW_BOTTOMDRAG`, `SEAICE_ALLOW_SIDEDRAG`, `SEAICE_BICE_STRESS`, `EXPLICIT_SSH_SLOPE`, `SEAICE_LSRBNEW`, `SEAICE_ITD`, `SEAICE_VARIABLE_SALINITY`, `SEAICE_CAP_ICELOAD`, `siceLoad`, `ALLOW_SITRACER`, `SEAICE_USE_GROWTH_ADX`, `SEAICE_OPTIONS.h`, `seaice_growth_adx.F`, `seaice_growth.F`
  - L89 Run-time parameters — `seaice_readparms.F`
    - L95 Enabling the package — `useSEAICE`
    - L101 General flags and parameters — `SEAICEwriteState`, `SEAICEuseDYNAMICS`, `SEAICEuseJFNK`, `SEAICEuseTEM`, `SEAICEuseMCS`, `SEAICEuseMCE`, `SEAICEuseTD`, `SEAICEusePL`, `SEAICEuseStrImpCpl`, `SEAICEselectMetricTerms`, `SEAICEuseEVPpickup`, `SEAICEuseFREEDRIFT`, `SEAICEuseFluxForm`, `SEAICErestoreUnderIce`, `SEAICEupdateOceanStress`, `SEAICEscaleSurfStress`, `SEAICEaddSnowMass`, `useHB87stressCoupling`, `usePW79thermodynamics`, `SEAICEadvHeff`, `HEFF`, `SEAICEadvArea`, `AREA`, `SEAICEadvSnow`, `HSNOW`, `SEAICEadvSalt`, `HSALT`, `SEAICEadvScheme`, `SEAICEuseFlooding`, `SINegFac`, `SEAICE_no_slip`, `SEAICE_deltaTtherm`, `dTtracerLev`, `SEAICE_deltaTdyn`, `SEAICE_deltaTevp`, `SEAICEuseEVPstar`, `SEAICEuseEVPrev`, `SEAICEnEVPstarSteps`, `SEAICE_evpAlpha`, `SEAICE_evpTauRelax` (+97)
  - L409 Description — `SEAICE_USE_GROWTH_ADX`, `SINegFac`, `seaice_growth.F`, `seaice_growth_adx.F`, `SEAICE_OPTIONS.h`
    - L474 Compatibility with ice-thermodynamics package :filelink:`pkg/thsice`
    - L499 Surface forcing
  - L514 Dynamics
    - L566 Viscous-Plastic (VP) Rheology — `seaice_sigma1`, `seaice_sigma2`, `PRESS0`, `SEAICE_strength`, `SEAICE_cStar`, `PRESS`, `SEAICEpressReplFac`, `SEAICE_DELTA_SMOOTHREG`, `SEAICE_deltaMin`, `SEAICE_tensilFac`, `SEAICE_eccen`, `SEAICE_eccfr`, `SEAICE_ALLOW_TEM`, `SEAICEuseTEM`, `SEAICEmcMU`, `SEAICE_ALLOW_MCE`, `SEAICEuseMCE`, `SEAICE_ALLOW_MCS`, `SEAICEuseMCS`, `SEAICE_ALLOW_TD`, `SEAICEuseTD`, `SEAICEusePL`
      - L707 Elliptical yield curve with normal flow rule — `SEAICE_eccen`, `SEAICE_deltaMin`, `SEAICE_EPS`, `SEAICE_zetaMaxFac`, `SEAICE_zetaMin`, `SEAICE_ZETA_SMOOTHREG`, `SEAICE_tensilFac`, `SEAICE_OPTIONS.h`
      - L782 Elliptical yield curve with non-normal flow rule — `SEAICE_eccfr`, `SEAICE_eccen`
      - L809 Truncated ellipse method (TEM) for elliptical yield curve — `SEAICE_ALLOW_TEM`, `SEAICEuseTEM`, `SEAICEmcMU`, `SEAICE_tensilFac`, `SEAICE_OPTIONS.h`
      - L841 Mohr-Coulomb yield curve with elliptical plastic potential — `SEAICE_ALLOW_MCE`, `SEAICEuseMCE`, `SEAICEmcMU`, `SEAICE_eccfr`, `SEAICE_tensilFac`, `SEAICE_OPTIONS.h`
      - L859 Mohr-Coulomb yield curve with shear flow rule — `SEAICE_ALLOW_MCS`, `SEAICEuseMCS`, `SEAICEmcMU`, `SEAICE_tensilFac`, `SEAICE_OPTIONS.h`
      - L878 Teardrop yield curve with normal flow rule — `SEAICE_ALLOW_TEARDROP`, `SEAICEuseTD`, `SEAICE_tensFac`, `SEAICE_tensilFac`, `SEAICE_OPTIONS.h`
      - L897 Parabolic lens yield curve with normal flow rule — `SEAICE_ALLOW_TEARDROP`, `SEAICEusePL`, `SEAICE_tensFac`, `SEAICE_tensilFac`, `SEAICE_OPTIONS.h`
    - L916 LSR and JFNK solver — `SEAICEnonLinIterMax`, `SEAICE_JFNK_lsIter`, `SEAICE_JFNK_lsGamma`, `SEAICE_JFNK_lsLmax`, `SEAICE_JFNKepsilon`, `SEAICEuseJFNK`, `SEAICE_ALLOW_JFNK`, `SEAICE_ZETA_SMOOTHREG`, `SEAICEnonLinTol`, `JFNKgamma_lin_max`, `JFNKgamma_lin_min`, `JFNKres_tFac`, `SEAICEnewtonIterMax`, `SEAICEkrylovIterMax`, `SEAICEuseStrImpCpl`, `SEAICE_OPTIONS.h`
    - L1082 Elastic-Viscous-Plastic (EVP) Dynamics — `SEAICE_deltaTevp`, `SEAICE_deltaTdyn`, `SEAICE_CGRID`, `SEAICE_ALLOW_EVP`, `SEAICEuseEVPstar`, `SEAICEuseEVPrev`, `deltaTmom`, `SEAICE_elasticParm`, `SEAICE_evpTauRelax`, `SEAICE_OPTIONS.h`
    - L1161 More stable variants of Elastic-Viscous-Plastic Dynamics: EVP\*, mEVP, and aEVP — `SEAICEuseEVPstar`, `SEAICE_evpAlpha`, `SEAICE_evpBeta`, `SEAICE_deltaTevp`, `SEAICE_elasticParm`, `SEAICE_evpTauRelax`, `SEAICEnEVPstarSteps`, `SEAICE_CGRID`, `SEAICE_ALLOW_EVP`, `SEAICEuseEVPrev`, `SEAICEaEVPcoeff`, `SEAICEaEVPcStar`, `SEAICEaEVPalphaMin`, `SEAICE_OPTIONS.h`
    - L1263 Ice-Ocean stress — `useHB87StressCoupling`
    - L1286 Finite-volume discretization of the stress tensor divergence
  - L1569 Thermodynamics
    - L1576 Zero-layer thermodynamics — `SEAICE_rhoAir`, `SEAICE_dalton`, `SEAICE_lhEvap`, `SEAICE_lhFusion`, `SEAICE_cpAir`, `useMaykutSatVapPoly`, `SEAICE_multDim`, `SEAICE_PDF`, `nITD`, `SEAICE_useMultDimSnow`, `SEAICE_gamma_t`, `HO`, `SEAICEuseFlooding`, `SEAICE_SIZE.h`
    - L1699 Advection of thermodynamic variables — `HEFF`, `AREA`, `HSNOW`, `SEAICEadvScheme`, `DIFF1`, `thSIceAdvScheme`, `thSIce_diffK`
    - L1749 Dynamical Ice Thickness Distribution (ITD)
      - L1763 Distribution, participation and redistribution functions in ridging — `SEAICE_ITD`, `SEAICEshearParm`, `SEAICEsimpleRidging`, `SEAICEpartFunc`, `SEAICEredistFunc`, `SEAICEaStar`, `SEAICEmaxRaft`, `SEAICEmuRidging`, `nITD`, `HO`, `Hlimit`, `Hlimit_c1`, `Hlimit_c2`, `Hlimit_c3`, `SEAICE_OPTIONS.h`, `SEAICE_SIZE.h`
      - L1874 Ice strength parameterization — `SEAICE_strength`, `SEAICE_cStar`, `SEAICE_cf`, `useHibler79IceStrength`
  - L1908 Known issues and work-arounds — `useRealFreshWaterFlux`, `SEAICEpressReplFac`, `SEAICE_CAP_ICELOAD`, `sIceLoad`, `heffTooHeavy`, `seaice_growth.F`
  - L1936 Key subroutines — `seaice_model.F`
  - L1987 SEAICE diagnostics
  - L2086 Experiments and tutorials that use seaice

## `doc/phys_pkgs/shap_filt.rst` — Shapiro Filter
- L1 Shapiro Filter
  - L7 Key subroutines, parameters and files
  - L12 Experiments and tutorials that use shap filter

## `doc/phys_pkgs/shelfice.rst` — SHELFICE Package
- L3 SHELFICE Package
  - L8 Introduction
  - L23 SHELFICE configuration — `ALLOW_SHELFICE_DEBUG`, `ALLOW_ISOMIP_TD`, `SHI_ALLOW_GAMMAFRICT`, `SHI_SALTBAL_FWFLX`, `SHELFICE_OPTIONS.h`
  - L61 SHELFICE run-time parameters — `useSHELFICE`, `SHELFICEtopoFile`, `SHELFICEloadAnomalyFile`, `rhoConst`, `useISOMIPTD`, `SHELFICEconserve`, `SHELFICEboundaryLayer`, `SHI_withBL_realFWflux`, `SHI_withBL_uStarTopDz`, `SHELFICEmassFile`, `SHELFICEMassDynTendFile`, `SHELFICETransCoeffTFile`, `SHELFICElatentHeat`, `SHELFICEHeatCapacity_Cp`, `rhoShelfIce`, `SHELFICEsalinity`, `SHELFICEheatTransCoeff`, `SHELFICEsaltTransCoeff`, `SHELFICEsaltToHeatRatio`, `SHELFICEkappa`, `SHELFICEthetaSurface`, `no_slip_shelfice`, `no_slip_bottom`, `SHELFICEDragLinear`, `bottomDragLinear`, `SHELFICEDragQuadratic`, `bottomDragQuadratic`, `SHELFICEselectDragQuadr`, `SHELFICEMassStepping`, `SHELFICEDynMassOnly`, `SHELFICEmassStepping`, `SHELFICEadvDiffHeatFlux`, `SHELFICEuseGammaFrict`, `SHELFICE_oldCalcUStar`, `SHELFICEwriteState`, `SHELFICE_dumpFreq`, `dumpFreq`, `SHELFICE_dump_mnc`, `snapshot_mnc`, `SHELFICE_PARM01` (+1)
  - L153 SHELFICE description — `SHELFICEloadAnomalyFile`, `SHELFICEboundaryLayer`
    - L292 Three-equations thermodynamics — `rhoConst`, `SHELFICEheatTransCoeff`, `SHELFICEsaltTransCoeff`, `rhoShelfIce`, `SHELFICEkappa`, `SHELFICElatentHeat`, `HeatCapacity_Cp`, `SHELFICEHeatCapacity_Cp`, `SHELFICEthetaSurface`, `SHELFICEadvDiffHeatFlux`, `SHELFICEsaltToHeatRatio`, `SHELFICEconserve`
    - L452 Solving the three-equations system — `cFac`, `SHELFICEconserve`, `rFac`, `useRealFreshWaterFlux`, `dFac`, `SHELFICEadvDiffHeatFlux`, `fwFlxFac`, `rFWinBL`, `SHI_withBL_realFWflux`, `SHELFICEboundaryLayer`, `SHI_SALTBAL_FWFLX`, `shelfice_thermodynamics.F`, `SHELFICE_OPTIONS.h`
    - L596 ISOMIP thermodynamics — `useISOMIPTD`
    - L619 Exchange coefficients — `SHELFICEheatTransCoeff`, `SHELFICEsaltTransCoeff`, `SHELFICEuseGammaFrict`, `shelfice_readparms.F`
    - L633 Remark
  - L643 Key subroutines — `shelfice_thermodynamics.F`
  - L686 SHELFICE diagnostics
  - L705 Experiments and tutorials that use shelfice

## `doc/phys_pkgs/streamice.rst` — STREAMICE Package
- L3 STREAMICE Package
  - L11 Introduction
  - L24 STREAMICE configuration
    - L27 Compile-time options — `STREAMICE_CONSTRUCT_MATRIX`, `STREAMICE_HYBRID_STRESS`, `USE_ALT_RLOW`, `STREAMICE_GEOM_FILE_SETUP`, `STREAMICE_PARM03`, `ALLOW_PETSC`, `STREAMICE_COULOMB_SLIDING`, `STREAMICE_SMOOTH_FLOATATION`, `STREAMICE_OPTIONS.h`
    - L63 Enabling the package — `useSTREAMICE`
    - L68 Runtime parmeters: general flags and parameters — `STREAMICE_PARM01`, `streamice_density`, `streamice_density_ocean_avg`, `n_glen`, `eps_glen_min`, `eps_u_min`, `n_basal_friction`, `streamice_cg_tol`, `streamice_lower_cg_tol`, `streamice_max_cg_iter`, `streamice_maxcgiter_cpl`, `streamice_nonlin_tol`, `streamice_max_nl_iter`, `streamice_maxnliter_cpl`, `streamice_nonlin_tol_fp`, `streamice_err_norm`, `streamice_chkfixedptconvergence`, `streamice_chkresidconvergence`, `streamicethickInit`, `streamicethickFile`, `STREAMICE_PARM03`, `streamice_move_front`, `streamice_calve_to_mask`, `streamice_calve_mask`, `STREAMICE_use_log_ctrl`, `streamicecalveMaskFile`, `streamice_diagnostic_only`, `streamice_CFL_factor`, `streamice_adjDump`, `streamicebasalTracConfig`, `streamicebasalTracFile`, `C_basal_fric_const`, `streamiceGlenConstConfig`, `streamiceGlenConstFile`, `B_glen_isothermal`, `streamiceBdotFile`, `streamiceBdotTimeDepFile`, `streamice_forcing_period`, `streamiceTopogFile`, `USE_ALT_RLOW` (+27)
    - L203 Configuring domain through files — `STREAMICE_GEOM_FILE_SETUP`, `streamice_hmask`, `streamiceHmaskFile`, `streamice_ufacemask_bdry`, `streamice_vfacemask_bdry`, `streamiceuFaceBdryFile`, `streamicevFaceBdryFile`, `u_flux_bdry_SI`, `v_flux_bdry_SI`, `streamiceuMassFluxFile`, `streamicevMassFluxFile`, `streamicethickFile`
    - L232 Configuring domain through parameters — `STREAMICE_GEOM_FILE_SETUP`, `flux_bdry_val_NORTH`, `flux_bdry_val_SOUTH`, `flux_bdry_val_EAST`, `flux_bdry_val_WEST`, `STREAMICE_PARM03`, `min_x_noflow_NORTH`, `max_x_noflow_NORTH`, `min_x_noflow_SOUTH`, `max_x_noflow_SOUTH`, `min_y_noflow_EAST`, `max_y_noflow_EAST`, `min_y_noflow_WEST`, `max_y_noflow_WEST`, `min_x_nostress_NORTH`, `max_x_nostress_NORTH`, `min_x_nostress_SOUTH`, `max_x_nostress_SOUTH`, `min_y_nostress_EAST`, `max_y_nostress_EAST`, `min_y_nostress_WEST`, `max_y_nostress_WEST`, `min_x_fluxbdry_NORTH`, `max_x_fluxbdry_NORTH`, `min_x_fluxbdry_SOUTH`, `max_x_fluxbdry_SOUTH`, `min_y_fluxbdry_EAST`, `max_y_fluxbdry_EAST`, `min_y_fluxbdry_WEST`, `max_y_fluxbdry_WEST`, `min_x_CFBC_NORTH`, `max_x_CFBC_NORTH`, `min_x_CFBC_SOUTH`, `max_x_CFBC_SOUTH`, `min_y_CFBC_EAST`, `max_y_CFBC_EAST`, `min_y_CFBC_WEST`, `max_y_CFBC_WEST`
  - L340 Description
    - L345 Equations Solved — `float_frac_streamice`, `STREAMICE_COULOMB_SLIDING`, `streamice_allow_reg_coulomb`
    - L486 Hybrid SIA-SSA stress balance
    - L519 Ice front advance — `streamice_move_front`, `streamice_calve_to_mask`, `streamice_calve_mask`
    - L543 Units of input files — `streamicebasalTracFile`, `C_basal_fric_const`, `streamiceGlenConstFile`, `B_glen_isothermal`, `n_glen`, `n_basal_friction`
  - L559 Numerical Details — `float_frac_streamice`, `streamice_hmask`, `streamice_umask`, `streamice_vmask`, `streamice_ufacemask_bdry`, `u_flux_bdry_SI`, `streamice_vfacemask_bdry`, `STREAMICE_GEOM_FILE_SETUP`, `SIZE.h`, `GRID.h`, `STREAMICE_OPTIONS.h`, `streamice_vel_solve.F`
  - L631 Additional Features — `float_frac_streamice`, `STREAMICE_SMOOTH_FLOATATION2`, `streamice_smooth_gl_width`
    - L645 PETSc — `ALLOW_PETSC`
    - L657 Boundary Stresses — `streamiceuNormalStressFile`, `streamicevNormalStressFile`, `streamiceuNormalTimeDepFile`, `streamicevNormalTimeDepFile`, `streamiceuShearStressFile`, `streamicevShearStressFile`, `streamiceuShearTimeDepFile`, `streamicevShearTimeDepFile`
  - L694 Adjoint
  - L703 Key Subroutines — `streamice_timestep.F`, `do_oceanic_phys.F`
  - L748 STREAMICE diagnostics
  - L773 Experiments and tutorials that use streamice

## `doc/phys_pkgs/thsice.rst` — THSICE: The Thermodynamic Sea Ice Package
- L3 THSICE: The Thermodynamic Sea Ice Package
  - L23 Key parameters and Routines
    - L56 subroutine ICE_FREEZE
    - L73 subroutine ICE_START
    - L117 subroutine ICE_THERM
    - L283 subroutine SFC_ALBEDO
    - L318 subroutine NEW_LAYERS_WINTON
    - L339 Initializing subroutines
    - L348 Diagnostic subroutines
    - L356 Common Blocks
    - L366 Input file DATA.ICE
  - L381 Important Notes
  - L391 THSICE Diagnostics
  - L420 References
  - L432 Experiments and tutorials that use thsice

## `doc/phys_pkgs/zonal_filt.rst` — FFT Filtering Code
- L1 FFT Filtering Code
  - L7 Key subroutines, parameters and files
  - L10 Experiments and tutorials that use zonal filter

## `doc/references.rst` — 

## `doc/related_projects/related_projects.rst` — Related Projects and Highlighted Papers
- L1 Related Projects and Highlighted Papers
  - L5 Projects Related to MITgcm
    - L8 Estimating the Circulation and Climate of the Ocean (ECCO)
    - L20 Southern Ocean State Estimation (SOSE)
    - L28 MITgcmIS: global ice sheet model for the coupled atm-ocn-sea ice MITgcm
    - L42 MITprof: In-Situ Ocean Data In Matlab And Octave
    - L51 OceanParcels - Lagrangian Particle Tracker
    - L60 Xgcm: General Circulation Model Postprocessing with xarray
    - L73 Xmitgcm
    - L80 Gcmfaces: Gridded Earth Variables In Matlab And Octave
  - L89 Highlighted Papers

## `doc/software_arch/software_arch.rst` — Software Architecture
- L3 Software Architecture
  - L21 Overall architectural goals
  - L78 WRAPPER
    - L106 Target hardware
    - L127 Supporting hardware neutrality
    - L142 WRAPPER machine model
    - L154 Machine model parallelism
      - L178 Tiles
      - L204 Tile layout
    - L229 Communication mechanisms
      - L241 Shared memory communication
      - L348 Distributed memory communication
    - L388 Communication primitives
    - L450 Memory architecture
    - L472 Summary
  - L522 Using the WRAPPER
    - L547 Specifying a domain decomposition — `sNx`, `OLx`, `OLy`, `nSx`, `nSy`, `nPx`, `nPy`, `myThid`, `cg2d_r`, `SIZE.h`, `dynamics.F`
      - L710 Examples of :filelink:`SIZE.h <model/inc/SIZE.h>` specifications — `SIZE.h`
    - L803 Starting the code — `main.F`, `the_model_main.F`
      - L845 Multi-threaded execution — `myThid`, `the_model_main.F`
      - L925 Multi-process execution
      - L980 Environment variables
      - L990 Runtime input parameters — `myProcId`, `MPI_COMM_MODEL`, `myXGlobalLo`, `myYGlobalLo`, `pidW`, `pidE`, `pidS`, `pidN`, `eeboot_minimal.F`, `ini_procs.F`, `EESUPPORT.h`
    - L1045 Controlling communication — `tileNo`, `tileNoN`, `tileNoS`, `tileNoE`, `tileNoW`, `tileCommModeN`, `tileCommModeS`, `tileCommModeE`, `tileCommModeW`, `nThreads`, `nTx`, `nTy`, `exchNeedsMemSync`, `cacheLineSize`, `lShare1`, `lShare4`, `lShare8`, `theSimulationMode`, `MAX_NO_THREADS`, `nSx`, `nSy`, `COMM_NONE`, `COMM_MSG`, `COMM_PUT`, `COMM_GET`, `_BARRIER`, `_GSUM`, `_EXCH`, `_EXCH_XY`, `_EXCH_XYZ`, `_EXCH_XY_R4`, `_EXCH_XYZ_R4`, `_EXCH_XY_R8`, `_EXCH_XYZ_R8`, `REVERSE_SIMULATION`, `ini_communication_patterns.F`, `the_model_main.F`, `main.F`, `MAIN_PDIRECTIVES1.h`, `ini_threading_environment.F` (+8)
      - L1303 Specializing the Communication Code
      - L1313 JAM example — `_EXCH`, `_GSUM`, `LETS_MAKE_JAM`, `eeboot.F`, `CPP_EEMACROS.h`, `cg2d.F`
      - L1341 Cube sphere communication — `useCubedSphereExchange`, `_EXCH`, `exch2_rx`, `exch2_uv_rx`, `exch_uv_rx`, `exch_z_rx`
  - L1369 MITgcm execution under WRAPPER
    - L1375 Annotated call tree for MITgcm and WRAPPER
    - L1416 Measuring and Characterizing Performance
    - L1421 Estimating Resource Requirements
      - L1426 Atlantic 1/6 degree example
      - L1429 Dry Run testing
      - L1432 Adjoint Resource Requirements
      - L1435 State Estimation Environment Resources

## `doc/utilities/utilities.rst` — Utilities
- L1 Utilities
  - L6 MITgcmutils
    - L33 mds
    - L39 mnc
    - L45 diagnostics
    - L51 ptracers
    - L57 density
    - L63 miscellaneous utilities
    - L69 conversion
    - L81 llc
    - L87 examples
    - L95 gluemncbig
  - L107 Bash scripts
    - L121 gluemnc
      - L130 Usage
      - L164 Dependencies
      - L175 Notes
