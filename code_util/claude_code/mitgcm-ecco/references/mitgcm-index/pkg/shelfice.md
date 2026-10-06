# pkg/shelfice

Ice-shelf cavities: pressure loading, three-equation melt, boundary-layer options, conservative remeshing.

**runtime switch:** `useSHELFICE`-style flag in `data.pkg` (check exact name in packages_boot.F)
**reads:** `data.shelfice`
**manual:** `doc/phys_pkgs/shelfice.rst`, `doc/examples/examples.rst`, `doc/ocean_state_est/ocean_state_est.rst`, `doc/phys_pkgs/remesh.rst`, `doc/phys_pkgs/streamice.rst`
**adjoint support files:** shelfice_ad_check_lev1_dir.h, shelfice_ad_check_lev2_dir.h, shelfice_ad_check_lev3_dir.h, shelfice_ad_check_lev4_dir.h, shelfice_ad_diff.list

## Namelist parameters
### SHELFICE_PARM01
- `SHELFICEsaltToHeatRatio` — constant ratio giving SHELFICEsaltTransCoeff/SHELFICEheatTransCoeff (def: 5.05e-3)
- `SHELFICEheatTransCoeff` — constant heat transfer coefficient that determines heat flux into shelfice (def: 1e-4 m/s)
- `SHELFICEsaltTransCoeff` — constant salinity transfer coefficient that determines salt flux into shelfice (def: SHELFICEsaltToHeatRatio * SHELFICEheatTransCoeff)
- `SHELFICEMassStepping` — flag to step forward ice shelf mass/thickness accounts for melting/freezing & dynamics (from file or from coupling), def: F
- `SHI_update_kTopC` — update lateral extension (kTopC) if ice-shelf retreats from or expands to model top level (requires to define ALLOW_SHELFICE_REMESHING)
- `rhoShelfice` — density of ice shelf (def: 917.0 kg/m^3)
- `SHELFICEkappa`
- `SHELFICElatentHeat` — latent heat of fusion (def: 334000 J/kg)
- `SHELFICEHeatCapacity_Cp` — heat capacity of ice shelf (def: 2000 J/K/kg)
- `no_slip_shelfice` — set slip conditions for shelfice separately, (by default the same as no_slip_bottom, but really should be false when there is linear or quadratic drag)
- `SHELFICEDragLinear` — linear drag at bottom shelfice (1/s)
- `SHELFICEDragQuadratic` — quadratic drag at bottom shelfice (default = shiCdrag or bottomDragQuadratic)
- `SHELFICEselectDragQuadr`
- `SHELFICEthetaSurface`
- `SHELFICEsalinity`
- `useISOMIPTD` — use simple ISOMIP thermodynamics, def: F
- `SHELFICEconserve` — use conservative form of H&O-thermodynamics following Jenkins et al. (2001, JPO), def: F
- `SHELFICEboundaryLayer` — turn on vertical merging of cells to for a boundary layer of drF thickness, def: F
- `SHI_withBL_realFWflux` — with above BL, allow to use real-FW flux (and adjust advective flux at boundary accordingly) def: F
- `SHI_withBL_uStarTopDz` — with SHELFICEboundaryLayer, compute uStar from uVel,vVel avergaged over top Dz thickness; def: F
- `SHELFICEwriteState` — enable output
- `SHELFICE_dumpFreq` — analoguous to dumpFreq (= default)
- `SHELFICE_taveFreq`
- `SHELFICE_tave_mnc`
- `SHELFICE_dump_mnc` — use netcdf for snapshot output
- `SHELFICEtopoFile` — File containing the topography of the shelfice draught (unit=m)
- `SHELFICEmassFile` — name of shelfice Mass file
- `SHELFICEloadAnomalyFile` — name of shelfice load anomaly file
- `SHELFICEMassDynTendFile` — file name for other mass tendency (e.g. dynamics)
- `SHELFICETransCoeffTFile`
- `SHELFICEDynMassOnly` — step ice mass ONLY with Shelficemassdyntendency (not melting/freezing) def: F
- `SHELFICEadvDiffHeatFlux` — use advective-diffusive heat flux into the shelf instead of default diffusive heat flux, see Holland and Jenkins (1999), eq.21,22,26,31; def: F
- `SHELFICEuseGammaFrict` — use velocity dependent exchange coefficients, see Holland and Jenkins (1999), eq.11-18, with the following parameters (def: F):
- `SHELFICE_oldCalcUStar` — use old uStar averaging expression
- `shiCdrag` — quadratic drag coefficient to compute uStar (def: 0.0015)
- `shiZetaN` — ??? (def: 0.052)
- `shiRc` — ??? (not used, def: 0.2)
- `shiPrandtl` — constant Prandtl (13.8) and Schmidt (2432.0) numbers used to compute gammaTurb
- `shiSchmidt` — constant Prandtl (13.8) and Schmidt (2432.0) numbers used to compute gammaTurb
- `shiKinVisc` — constant kinetic viscosity used to compute gammaTurb (def: 1.95e-5)
- `mult_shelfice`  _[ifdef ALLOW_COST]_
- `mult_shifwflx`  _[ifdef ALLOW_COST]_
- `wshifwflx0`  _[ifdef ALLOW_COST]_
- `shifwflx_errfile`  _[ifdef ALLOW_COST]_
- `SHELFICEremeshFrequency` — Frequency (in seconds) of call to SHELFICE_REMESHING (def: 0. --> no remeshing)
- `SHELFICEsplitThreshold` — Thickness fraction remeshing threshold above which top-cell splits (no unit)
- `SHELFICEmergeThreshold` — Thickness fraction remeshing threshold below which top-cell merges with below (no unit)

## CPP options (defaults as shipped)
- `ALLOW_ISOMIP_TD` (define, SHELFICE_OPTIONS.h) — allow code for simple ISOMIP thermodynamics
- `SHI_ALLOW_GAMMAFRICT` (define, SHELFICE_OPTIONS.h) — allow friction velocity-dependent transfer coefficient following Holland and Jenkins, JPO, 1999
- `SHI_SALTBAL_FWFLX` (undef, SHELFICE_OPTIONS.h) — balance equation instead of the heat balance equation. Ill-defined when all salinities (ice, ocean, boundary layer) are zero, therefore deprecated.
- `ALLOW_SHELFICE_REMESHING` (undef, SHELFICE_OPTIONS.h) — allow (vertical) remeshing whenever ocean top thickness factor exceeds thresholds
- `SHELFICE_REMESH_PRINT` (define, SHELFICE_OPTIONS.h) — and allow to print message to STDOUT when this happens

## Headers
- `SHELFICE.h` — BOP
- `SHELFICE_COST.h` — BOP
- `SHELFICE_OPTIONS.h` — CPP options file for SHELFICE package. Use this file for selecting options within the SHELFICE package.
- `shelfice_ad_check_lev1_dir.h` — ADJ STORE kTopC            = comlev1, key=ikey_dynamics ADJ STORE shelficeMass     = comlev1, key=ikey_dynamics, kind=isbyte ADJ STORE shelficeForcing
- `shelfice_ad_check_lev2_dir.h` — ADJ STORE kTopC            = tapelvi2, key = ilev_2 ADJ STORE phi0surf         = tapelev2, key = ilev_2 ADJ STORE shelficeMass     = tapelev2, key = i
- `shelfice_ad_check_lev3_dir.h` — ADJ STORE kTopC            = tapelvi3, key = ilev_3 ADJ STORE phi0surf         = tapelev3, key = ilev_3 ADJ STORE shelficeMass     = tapelev3, key = i
- `shelfice_ad_check_lev4_dir.h` — ADJ STORE kTopC            = tapelvi4, key = ilev_4 ADJ STORE phi0surf         = tapelev4, key = ilev_4 ADJ STORE shelficeMass     = tapelev4, key = i

## Routines (27)
`shelfice_check.F`, `shelfice_cost_accumulate.F`, `shelfice_cost_final.F`, `shelfice_cost_shifwflx.F`, `shelfice_diagnostics_drag.F`, `shelfice_forcing.F`, `shelfice_forcing_surf.F`, `shelfice_init_depths.F`, `shelfice_init_fixed.F`, `shelfice_init_varia.F`, `shelfice_mask_seaice.F`, `shelfice_mnc_init.F`, `shelfice_output.F`, `shelfice_read_pickup.F`, `shelfice_readparms.F`, `shelfice_remesh_c_mask.F`, `shelfice_remesh_calc_w.F`, `shelfice_remesh_state.F`, `shelfice_remesh_uv_mask.F`, `shelfice_remeshing.F`, `shelfice_step_icemass.F`, `shelfice_thermodynamics.F`, `shelfice_u_drag_coeff.F`, `shelfice_v_drag_coeff.F`, `shelfice_write_pickup.F`

## Called from outside the package
- `SHELFICE_FORCING_S` ← `model/src/apply_forcing.F:937`
- `SHELFICE_FORCING_T` ← `model/src/apply_forcing.F:705`
- `SHELFICE_DIAGNOSTICS_DRAG` ← `model/src/correction_step.F:291`
- `SHELFICE_THERMODYNAMICS` ← `model/src/do_oceanic_phys.F:523`
- `SHELFICE_OUTPUT` ← `model/src/do_the_model_io.F:200`
- `SHELFICE_FORCING_S` ← `model/src/external_forcing.F:767`
- `SHELFICE_FORCING_T` ← `model/src/external_forcing.F:563`
- `SHELFICE_FORCING_SURF` ← `model/src/external_forcing_surf.F:396`
- `SHELFICE_REMESHING` ← `model/src/forward_step.F:456`
- `SHELFICE_INIT_DEPTHS` ← `model/src/ini_masks_etc.F:56`
- `SHELFICE_CHECK` ← `model/src/packages_check.F:346`
- `SHELFICE_INIT_FIXED` ← `model/src/packages_init_fixed.F:469`
- `SHELFICE_INIT_VARIA` ← `model/src/packages_init_variables.F:386`
- `SHELFICE_READPARMS` ← `model/src/packages_readparms.F:281`
- `SHELFICE_WRITE_PICKUP` ← `model/src/packages_write_pickup.F:178`
- `SHELFICE_COST_FINAL` ← `pkg/cost/cost_final.F:117`
- `SHELFICE_COST_ACCUMULATE` ← `pkg/cost/cost_tile.F:128`
- `SHELFICE_U_DRAG_COEFF` ← `pkg/ggl90/ggl90_calc.F:834`
- `SHELFICE_V_DRAG_COEFF` ← `pkg/ggl90/ggl90_calc.F:838`
- `SHELFICE_U_DRAG_COEFF` ← `pkg/mom_common/mom_u_implicit_r.F:188`
- `SHELFICE_V_DRAG_COEFF` ← `pkg/mom_common/mom_v_implicit_r.F:188`
- `SHELFICE_U_DRAG_COEFF` ← `pkg/mom_fluxform/mom_fluxform.F:704`
- `SHELFICE_V_DRAG_COEFF` ← `pkg/mom_fluxform/mom_fluxform.F:999`
- `SHELFICE_U_DRAG_COEFF` ← `pkg/mom_vecinv/mom_vecinv.F:521`
- `SHELFICE_V_DRAG_COEFF` ← `pkg/mom_vecinv/mom_vecinv.F:631`
- `SHELFICE_NETMASSFLUX_SURF` ← `pkg/obcs/obcs_balance_flow.F:329`
- `SHELFICE_MASK_SEAICE` ← `pkg/seaice/seaice_init_fixed.F:470`
- `SHELFICE_MASK_SEAICE` ← `pkg/thsice/thsice_get_ocean.F:94`

## Verification experiments compiling it (2)
`isomip` `shelfice_2d_remesh`
