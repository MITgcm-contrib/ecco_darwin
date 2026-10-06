# MITgcm packages

One row per `pkg/` directory. Details (namelists with descriptions, CPP options, call sites, experiments): `pkg/<name>.md`. "#nml" = namelist parameters, "#cpp" = CPP options in its *_OPTIONS.h, AD = ships adjoint diff lists. Experiments = verification experiments whose packages.conf (groups expanded) includes it. Rows named `pkg@label` come from the user's branch clones and are new or differ from upstream (see README.md for label → path).

| pkg | data file | #nml | #cpp | AD | experiments | description |
|---|---|---|---|---|---|---|
| `admtlm` |  | 0 | 0 |  | 0 | Adjoint/tangent-linear combined (ADM-TLM) driver for singular-vector / Hessian-type calculations (legacy). |
| `aim_v23` | data.aimphys | 92 | 6 |  | 4: aim.5l_Equatorial_Channel, aim.5l_LatLon, aim.5l_cs, cpl_aim+ocn | AIM intermediate-complexity atmospheric physics (SPEEDY v23: convection, clouds, radiation, surface fluxes) for atmosphere set-ups |
| `atm2d` | data.atm2d | 28 | 17 |  | 0 | 2-D (zonally averaged) atmosphere coupled to an MITgcm ocean (IGSM-style climate coupling). |
| `atm_common` |  | 0 | 0 |  | 0 | Shared diagnostics/variables for atmospheric physics packages. |
| `atm_compon_interf` | data.cpl | 10 | 0 |  | 1: cpl_aim+ocn | Atmosphere-side interface of the coupler (cpl_aim+ocn style atmosphere-ocean coupling). |
| `atm_ocn_coupler` |  | 8 | 0 |  | 1: cpl_aim+ocn | Stand-alone coupler component exchanging fields between atmosphere and ocean MITgcm components. |
| `atm_phys` | data.atm_gray data.atm_phys | 23 | 0 |  | 1: atm_gray | Grey-radiation idealized atmospheric physics (Frierson/O'Gorman-type moist physics). |
| `autodiff` | data.autodiff | 19 | 16 | y | 16: 1D_ocean_ice_column, bottom_ctrl_5x5, global_oce_biogeo_bling, global_oce_latlon, … | Automatic differentiation support (TAF/Tapenade/OpenAD): checkpointing, tape storage directives, adjoint I/O and dumps. |
| `bbl` | data.bbl | 6 | 0 | y | 1: global_oce_latlon | Bottom boundary layer: a thin bottom layer with its own T/S exchanging with the interior and downslope (dense-overflow representat |
| `bling` | data.bling | 134 | 17 | y | 1: global_oce_biogeo_bling | BLING biogeochemistry (Biogeochemistry with Light, Iron, Nutrients and Gases; Galbraith et al.) via gchem/ptracers. |
| `bulk_force` | data.blk | 53 | 2 |  | 1: global_ocean.cs32x15 | Simple bulk-formula surface forcing (older alternative to exf). |
| `cal` | data.cal | 4 | 0 | y | 2: global_oce_biogeo_bling, obcs_ctrl | Calendar tool (Gregorian/360-day/model calendar) used by exf, ecco, ctrl, diagnostics with calendarDumps. |
| `cd_code` |  | 0 | 2 | y | 14: bottom_ctrl_5x5, cfc_example, exp2, global_oce_biogeo_bling, … | C-D scheme: carries D-grid velocities to stabilise Coriolis on the C-grid at coarse resolution (useCDscheme). |
| `cfc` | data.cfc | 10 | 0 | y | 2: cfc_example, tutorial_cfc_offline | CFC-11/CFC-12 ocean tracers with air-sea gas exchange (OCMIP protocol) via gchem/ptracers. |
| `cheapaml` | data.cheapaml | 61 | 1 |  | 1: cheapAML_box | Cheap atmospheric mixed layer: simple prognostic atmospheric boundary layer above the ocean for surface fluxes. |
| `chronos` |  | 0 | 0 |  | 0 | Clock/date utilities (legacy). |
| `compon_communic` |  | 0 | 0 |  | 1: cpl_aim+ocn | MPI component communication layer for coupled (multi-executable) runs. |
| `cost` | data.cost | 13 | 13 | y | 16: 1D_ocean_ice_column, bottom_ctrl_5x5, global_oce_biogeo_bling, global_oce_latlon, … | Cost-function framework for adjoint/state estimation (accumulates and writes the cost). |
| `ctrl` | data.ctrl data.optim | 239 | 33 | y | 16: 1D_ocean_ice_column, bottom_ctrl_5x5, global_oce_biogeo_bling, global_oce_latlon, … | Control-vector handling for adjoint/optimisation: generic 2-D/3-D/time-varying controls, packing/unpacking, preconditioning, smoot |
| `darwin (darwin3)` | data.darwin data.traits | 894 | 47 |  | 0 | darwin3 Darwin ecosystem model: many plankton types/groups (allometric or random traits), nutrients, carbon chemistry, iron, radtr |
| `debug` |  | 0 | 0 |  | 58: 1D_ocean_ice_column, MLAdjust, adjustment.cs-32x32x1, advect_cs, … | Debug utilities (debugMode, debugLevel, field stats printing, DEBUG_ENTER/LEAVE). |
| `diagnostics` | data.diagnostics | 37 | 2 | y | 41: MLAdjust, adjustment.cs-32x32x1, advect_cs, advect_xz, … | Diagnostics framework: data.diagnostics output streams (time averages/snapshots, levels, frequencies) and statistics-diagnostics;  |
| `dic` | data.dic | 55 | 18 | y | 3: so_box_biogeo, tutorial_dic_adjoffline, tutorial_global_oce_biogeo | Simple dissolved inorganic carbon biogeochemistry (DIC, ALK, PO4, DOP, O2, Fe) with carbonate chemistry and air-sea CO2 flux (OCMI |
| `down_slope` | data.down_slope | 5 | 0 | y | 2: global_ocean.90x40x15, lab_sea | Down-slope flow parameterisation (Campin & Goosse 1999) moving dense bottom water downslope (excludes bbl). |
| `ebm` | data.ebm | 3 | 3 | y | 1: global_oce_latlon | Energy-balance atmosphere model coupled to the ocean. |
| `ecco` | data.ecco | 280 | 15 | y | 4: 1D_ocean_ice_column, global_oce_biogeo_bling, lab_sea, obcs_ctrl | ECCO state-estimation package: model-data misfit cost terms, generic cost (gencost), averaging, ECCO-specific I/O. |
| `embed_files` |  | 0 | 0 |  | 0 | Embeds input files into the executable (e.g. for testing/portability). |
| `exch2` | data.exch2 | 7 | 5 | y | 15: MLAdjust, adjustment.cs-32x32x1, advect_cs, aim.5l_cs, … | Generalised tile exchanges for cubed-sphere and LLC grids (facets, blank tiles via data.exch2, W2_mapIO). |
| `exf` | data.exf | 585 | 25 | y | 8: 1D_ocean_ice_column, global_oce_latlon, global_ocean.cs32x15, lab_sea, … | External forcing: reads/interpolates/time-interpolates atmospheric state or fluxes, bulk formulae (Large & Yeager), runoff, open-b |
| `fizhi` | data.fizhi | 7 | 6 |  | 3: fizhi-cs-32x32x40, fizhi-cs-aqualev20, fizhi-gridalt-hs | Fizhi atmospheric physics (NASA GEOS-like) for atmosphere configurations. |
| `flt` | data.flt | 9 | 7 |  | 2: MLAdjust, exp4 | Lagrangian floats/particles advected online (2-D/3-D, profiling floats). |
| `frazil` |  | 0 | 0 | y | 1: global_oce_latlon | Frazil ice formation: removes supercooling and adjusts heat/salt (used with shelfice/seaice set-ups). |
| `gchem` | data.gchem | 25 | 4 | y | 6: cfc_example, global_oce_biogeo_bling, so_box_biogeo, tutorial_cfc_offline, … | Geochemistry driver: interface between ptracers and BGC packages (dic, bling, cfc, darwin); separate forcing step, surface forcing |
| `generic_advdiff` |  | 0 | 7 | y | 56: 1D_ocean_ice_column, MLAdjust, advect_cs, advect_xz, … | Generic advection-diffusion operators (2nd/3rd/4th order, DST, flux-limited, OS7MP, Prather SOM, multi-dim) used for T, S and ptra |
| `ggl90` | data.ggl90 | 40 | 7 | y | 5: global_oce_latlon, global_ocean.90x40x15, global_ocean.cs32x15, isomip, … | Gaspar-Gregoris-Lefevre (1990) TKE vertical mixing scheme (ECCO v4 default), with IDEMIX and Langmuir options. |
| `gmredi` | data.gmredi | 75 | 13 | y | 20: MLAdjust, cfc_example, cpl_aim+ocn, front_relax, … | Gent-McWilliams / Redi isopycnal eddy parameterisation (GM, Redi, tapering, Visbeck, GEOMETRIC, advective/skew forms). |
| `grdchk` | data.grdchk | 17 | 0 |  | 16: 1D_ocean_ice_column, bottom_ctrl_5x5, global_oce_biogeo_bling, global_oce_latlon, … | Gradient check: compares adjoint gradient against finite differences (data.grdchk). |
| `gridalt` |  | 0 | 0 |  | 3: fizhi-cs-32x32x40, fizhi-cs-aqualev20, fizhi-gridalt-hs | Alternative vertical grid support for atmospheric physics (fizhi). |
| `icefront` | data.icefront | 18 | 0 |  | 1: isomip | Vertical ice-front (tidewater glacier face) melt parameterisation at side walls. |
| `kl10` | data.kl10 | 4 | 0 |  | 1: internal_wave | Klymak & Legg (2010) internal-wave-breaking vertical mixing scheme. |
| `kpp` | data.kpp | 48 | 17 | y | 5: 1D_ocean_ice_column, lab_sea, seaice_obcs, tutorial_tracer_adjsens, … | K-Profile Parameterization (Large et al. 1994) vertical mixing with nonlocal transport. |
| `land` | data.land | 37 | 2 |  | 2: aim.5l_cs, cpl_aim+ocn | Simple land model (bucket hydrology, soil temperature) for atmospheric configurations. |
| `layers` | data.layers | 10 | 8 | y | 3: cfc_example, exp4, tutorial_reentrant_channel | Diagnostics of transport in isopycnal/temperature/salinity layers (residual overturning). |
| `longstep` | data.longstep | 2 | 0 |  | 1: lab_sea | Longer time step for passive tracers than for dynamics (ptracers stepped every N steps with averaged velocities). |
| `matrix` | data.matrix | 2 | 0 |  | 1: matrix_example | Transport-matrix method: extracts explicit/implicit tracer transport matrices. |
| `mdsio` |  | 0 | 9 |  | 58: 1D_ocean_ice_column, MLAdjust, adjustment.cs-32x32x1, advect_cs, … | MDS (.data/.meta) binary I/O routines; always compiled. |
| `mnc` |  | 21 | 2 |  | 21: MLAdjust, aim.5l_cs, cheapAML_box, fizhi-cs-aqualev20, … | NetCDF I/O (per-tile files) for state, pickups, diagnostics when diag_mnc is on. |
| `mom_common` |  | 0 | 9 | y | 54: 1D_ocean_ice_column, MLAdjust, adjustment.cs-32x32x1, advect_cs, … | Momentum code shared by flux-form and vector-invariant: viscosity, bottom/side drag, Smagorinsky/Leith, metric terms. |
| `mom_fluxform` |  | 0 | 1 | y | 48: 1D_ocean_ice_column, MLAdjust, adjustment.cs-32x32x1, advect_cs, … | Flux-form momentum equations (default unless vectorInvariantMomentum=.TRUE.). |
| `mom_vecinv` |  | 0 | 1 | y | 53: 1D_ocean_ice_column, MLAdjust, adjustment.cs-32x32x1, advect_cs, … | Vector-invariant momentum equations (used by ECCO/cubed-sphere/LLC configs; enstrophy/energy-conserving Coriolis). |
| `monitor` |  | 0 | 1 |  | 58: 1D_ocean_ice_column, MLAdjust, adjustment.cs-32x32x1, advect_cs, … | Monitor statistics (%MON lines in STDOUT: CFL, KE, field min/max/mean/sd) used by testreport. |
| `my82` | data.my82 | 6 | 1 |  | 1: vermix | Mellor-Yamada (1982) level-2 vertical mixing. |
| `mypackage` | data.mypackage | 21 | 5 | y | 1: hs94.1x64x5 | TEMPLATE package: skeleton showing how to write a new package (readparms, check, init, diagnostics, pickup, tendencies). Start new |
| `obcs` | data.obcs | 153 | 19 | y | 10: dome, exp4, internal_wave, isomip, … | Open boundary conditions: prescribed/relaxed boundary values (T,S,U,V,eta,seaice, ptracers), sponges, Orlanski radiation, Stevens  |
| `obsfit` | data.obsfit | 6 | 2 | y | 1: global_oce_biogeo_bling | Model-data comparison for generic (unstructured) observations in the cost function. |
| `ocn_compon_interf` | data.cpl | 15 | 0 |  | 1: cpl_aim+ocn | Ocean-side interface of the coupler for coupled atmosphere-ocean runs. |
| `offline` | data.off | 22 | 1 | y | 2: tutorial_cfc_offline, tutorial_dic_adjoffline | Offline mode: reads precomputed velocities/diffusivities/forcing (e.g. ECCO archive) to drive passive tracers/BGC without dynamics |
| `openad` |  | 0 | 5 | y | 0 | OpenAD automatic-differentiation support (legacy). |
| `opps` | data.opps | 9 | 2 |  | 1: vermix | OPPS convection (Paluszkiewicz & Romea penetrative plume scheme). |
| `pp81` | data.pp81 | 9 | 2 |  | 1: vermix | Pacanowski & Philander (1981) Richardson-number vertical mixing. |
| `profiles` | data.profiles | 14 | 3 | y | 2: global_oce_biogeo_bling, global_oce_latlon | Model-data comparison for in situ profiles (Argo, CTD, XBT...) for ECCO cost. |
| `ptracers` | data.ptracers | 34 | 1 | y | 12: cfc_example, exp4, global_oce_biogeo_bling, global_ocean.90x40x15, … | Passive tracers: arbitrary number of tracers with their own advection schemes, diffusivities, initial files, surface/relaxation; c |
| `radtrans (darwin3)` | data.radtrans | 52 | 1 |  | 0 | darwin3 spectral radiative transfer (direct/diffuse irradiance in nlam wavebands, OASIM forcing) used by Darwin. |
| `rbcs` | data.rbcs | 27 | 1 | y | 2: exp4, tutorial_reentrant_channel | Relaxation boundary conditions: 3-D relaxation of T, S, ptracers (and U/V) toward prescribed fields with masks/timescales. |
| `regrid` | data.regrid | 5 | 1 |  | 0 | Regridding utilities for diagnostics output. |
| `runclock` | data.runclock | 3 | 1 |  | 0 | Wall-clock based run termination (stop gracefully before walltime). |
| `rw` |  | 0 | 2 | y | 58: 1D_ocean_ice_column, MLAdjust, adjustment.cs-32x32x1, advect_cs, … | Read/write utilities on top of mdsio (READ_FLD_XY_RL etc.); always compiled. |
| `salt_plume` | data.salt_plume | 11 | 3 | y | 3: cpl_aim+ocn, lab_sea, seaice_obcs | Brine rejection from sea-ice growth distributed vertically (salt plume parameterisation, ECCO v4). |
| `sbo` | data.sbo | 2 | 0 |  | 1: global_ocean.90x40x15 | Solid-body ocean diagnostics: angular momentum, centre of mass (Earth rotation studies). |
| `seaice` | data.seaice | 282 | 53 | y | 7: 1D_ocean_ice_column, cpl_aim+ocn, global_ocean.cs32x15, lab_sea, … | Dynamic-thermodynamic sea ice (viscous-plastic/EVP/JFNK/Krylov dynamics, zero-layer/ITD thermodynamics, snow, advection). |
| `shap_filt` | data.shap | 16 | 4 | y | 12: aim.5l_Equatorial_Channel, aim.5l_LatLon, aim.5l_cs, atm_gray, … | Shapiro filter for grid-scale noise. |
| `shelfice` | data.shelfice | 47 | 5 | y | 2: isomip, shelfice_2d_remesh | Ice-shelf cavities: pressure loading, three-equation melt, boundary-layer options, conservative remeshing. |
| `showflops` |  | 0 | 1 |  | 0 | FLOP counting/timing (PAPI). |
| `smooth` | data.smooth | 18 | 0 | y | 0 | Diffusion-operator smoothing for control/covariance (2-D/3-D correlation operators). |
| `sphere` |  | 0 | 0 | y | 0 | Spherical-harmonic utilities. |
| `steep_icecavity` | data.stic | 2 | 1 | y | 1: isomip | Ice cavities with steep ice draft (alternative shelfice treatment). |
| `streamice` | data.streamice data.strmctrlflux | 182 | 22 | y | 1: halfpipe_streamice | Ice-stream/shelf dynamics model (shallow-shelf / hybrid), optionally coupled to shelfice. |
| `tapenade` |  | 0 | 3 |  | 8: global_oce_biogeo_bling, global_oce_latlon, global_ocean.cs32x15, halfpipe_streamice, … | Tapenade automatic-differentiation support. |
| `thsice` | data.ice | 72 | 4 | y | 4: aim.5l_cs, cpl_aim+ocn, global_ocean.cs32x15, offline_exf_seaice | Thermodynamic sea ice (Winton 2000 two-layer) usable alone or with seaice dynamics. |
| `zonal_filt` | data.zonfilt | 6 | 0 |  | 2: aim.5l_LatLon, hs94.128x64x5 | Zonal (polar) Fourier filter for lat-lon grids. |
| `autodiff@bbl` | data.autodiff | 19 | 16 | y | 1: global_oce_latlon@bbl | ad dump record number (used only if dumpAdByRec is true) |
| `autodiff@bbl_c68g` | data.autodiff | 19 | 12 | y | 0 | ad dump record number (used only if dumpAdByRec is true) |
| `bbl@bbl` | data.bbl | 16 | 1 | y | 1: global_oce_latlon@bbl | bbl_wvel    :: default vertical entrainment velocity (m/s) bbl_hvvel   :: default horizontal velocity of BBL (m/s); upper bound of |
| `bbl@bbl_c68g` | data.bbl | 16 | 1 | y | 0 | bbl_wvel    :: default vertical entrainment velocity (m/s) bbl_hvvel   :: default horizontal velocity of BBL (m/s); upper bound of |
| `darwin@d3backport` | data.darwin data.traits | 926 | 48 |  | 0 | --  File darwin_carbon_chem.F: --   Contents --   o DARWIN_CALC_PCO2 --   o DARWIN_CALC_PCO2_APPROX --   o DARWIN_CARBON_COEFFS |
| `exch2@bbl` | data.exch2 | 7 | 5 | y | 0 | ---+----1----+----2----+----3----+----4----+----5----+----6----+----7-/--+---- BOP 0 |
| `gmredi@wad` | data.gmredi | 75 | 13 | y | 1: wad_estuary_3d@wad | Calculates the 3D diffusivity as per Bates et al. (2014) \ev |
| `gmredi@wadcheckin` | data.gmredi | 75 | 13 | y | 1: wad_estuary_3d@wadcheckin | Calculates the 3D diffusivity as per Bates et al. (2014) \ev |
| `kpp@wad` | data.kpp | 48 | 17 | y | 4: wad_estuary_3d@wad, wad_flat_xz@wad, wad_mangrove@wad, wad_mudflat@wad | ---+----1----+----2----+----3----+----4----+----5----+----6----+----7-/--+---- BOP 0 |
| `kpp@wadcheckin` | data.kpp | 48 | 17 | y | 3: wad_estuary_3d@wadcheckin, wad_flat_xz@wadcheckin, wad_mudflat@wadcheckin | ---+----1----+----2----+----3----+----4----+----5----+----6----+----7-/--+---- BOP 0 |
| `mangrove@wad` | data.mangrove | 14 | 0 |  | 1: wad_mangrove@wad | ---+----1----+----2----+----3----+----4----+----5----+----6----+----7-/--+---- BOP |
| `mnc@bbl` |  | 21 | 2 |  | 1: global_oce_latlon@bbl | EH3 package-specific options go here |
| `mom_vecinv@wad` |  | 0 | 1 | y | 8: wad_balzano@wad, wad_estuary_3d@wad, wad_flat_xz@wad, wad_iceground@wad, … | S/R MOM_VECINV Form the right hand-side of the momentum equation. Terms are evaluated one layer at a time working from the bottom  |
| `mom_vecinv@wadcheckin` |  | 0 | 1 | y | 5: wad_balzano@wadcheckin, wad_estuary_3d@wadcheckin, wad_flat_xz@wadcheckin, wad_mudflat@wadcheckin, … | S/R MOM_VECINV Form the right hand-side of the momentum equation. Terms are evaluated one layer at a time working from the bottom  |
| `obcs@seaicebc` | data.obcs | 165 | 19 | y | 0 | -- Fields and files for OBCS-support for passive tracers package PTRACERS |
| `overflood@wad` | data.overflood | 16 | 0 |  | 1: wad_overflood@wad | pkg/overflood: a thin sheet of water on top of floating sea ice (river overflood at breakup), coupled to the ocean beneath and to  |
| `ptracers@bbl` | data.ptracers | 34 | 1 | y | 0 | ---+----1----+----2----+----3----+----4----+----5----+----6----+----7-/--+---- BOP |
| `ptracers@bbl_c68g` | data.ptracers | 34 | 1 | y | 0 | ---+----1----+----2----+----3----+----4----+----5----+----6----+----7-/--+---- BOP |
| `ptracers@wad` | data.ptracers | 34 | 1 | y | 3: wad_estuary_3d@wad, wad_mangrove@wad, wad_mudflat@wad | ---+----1----+----2----+----3----+----4----+----5----+----6----+----7-/--+---- BOP |
| `radtrans@d3backport` | data.radtrans | 52 | 1 |  | 0 | ---+----1----+----2----+----3----+----4----+----5----+----6----+----7-/--+---- BOP |
| `regrid@bbl` | data.regrid | 5 | 1 |  | 0 | ---+----1----+----2----+----3----+----4----+----5----+----6----+----7-/--+---- BOP 0 |
| `seaice@seaicebc` | data.seaice | 282 | 53 | y | 0 | Basic parameter header for sea ice model. |
| `seaice@wad` | data.seaice | 283 | 53 | y | 4: wad_estuary_3d@wad, wad_iceground@wad, wad_mudflat@wad, wad_overflood@wad | Basic parameter header for sea ice model. |
| `sediment@wad` | data.sediment | 45 | 0 |  | 2: wad_mangrove@wad, wad_mudflat@wad | ---+----1----+----2----+----3----+----4----+----5----+----6----+----7-/--+---- BOP |
| `wad@wad` | data.wad | 28 | 2 | y | 8: wad_balzano@wad, wad_estuary_3d@wad, wad_flat_xz@wad, wad_iceground@wad, … | ---+----1----+----2----+----3----+----4----+----5----+----6----+----7-/--+---- BOP |
| `wad@wadcheckin` | data.wad | 20 | 2 | y | 5: wad_balzano@wadcheckin, wad_estuary_3d@wadcheckin, wad_flat_xz@wadcheckin, wad_mudflat@wadcheckin, … | ---+----1----+----2----+----3----+----4----+----5----+----6----+----7-/--+---- BOP |

## Core (model/) changes in the user's branches

- **bbl** `model/src`: changed: forward_step.F, packages_init_variables.F, the_main_loop.F
- **bbl_c68g** `model/src`: changed: forward_step.F, packages_init_variables.F, the_main_loop.F
- **wad** `model/src`: changed: apply_forcing.F, calc_3d_diffusivity.F, calc_r_star.F, calc_viscosity.F, external_forcing_surf.F, forward_step.F, initialise_varia.F, load_fields_driver.F, packages_boot.F, packages_check.F, packages_init_fixed.F, packages_init_variables.F, packages_readparms.F, packages_write_pickup.F, update_r_star.F, update_surf_dr.F
- **wad** `model/inc`: changed: PARAMS.h
- **wadcheckin** `model/src`: changed: calc_r_star.F, external_forcing_surf.F, forward_step.F, initialise_varia.F, packages_boot.F, packages_check.F, packages_init_fixed.F, packages_init_variables.F, packages_readparms.F, packages_write_pickup.F, update_r_star.F, update_surf_dr.F
- **wadcheckin** `model/inc`: changed: PARAMS.h
