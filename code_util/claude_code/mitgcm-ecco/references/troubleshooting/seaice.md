# Troubleshooting: pkg/seaice and thsice
Distilled from mitgcm-support (2003-2026) and MITgcm GitHub issues/PRs, core-developer answers only; names checked against upstream master (Oct 2026) and the source index. "Martin" = Martin Losch, "J-M" = Jean-Michel Campin. Month-thread URLs are http://mailman.mitgcm.org/pipermail/mitgcm-support/YYYY-Month/thread.html.

## Blow-ups, NaNs, runaway ice thickness

### "STOP in CALC_R_STAR" / r* error / surface cell dries, with HEFF of 10-50 m in a few cells
- Cause: ice piles up in corners/embayments (dynamical, rarely thermodynamic growth), loads the surface cell (z*/useRealFreshWaterFlux) until it runs dry. Usual culprits: replacement pressure (`SEAICEpressReplFac=1` default) lets stagnant ice compress without resistance; sea-ice OBCS in the middle of an active pack; shelfice.
- Fix: `SEAICEpressReplFac = 0.` in data.seaice (first thing Martin suggests); stable advection (`SEAICEadvScheme=77`) rather than `SEAICEdiffKh*`; keep `SEAICEadvSnow=.TRUE.` and `SEAICEuseFlooding=.TRUE.` (snow pile-up blew up a run once); `SEAICE_drag_south` / larger `SEAICE_drag` to push ice out of corners; make the top cell thick (5 m) if running NLFS=4 in plain z. Last resort: cap HEFF (unphysical). If shelfice is on, try `SHELFICEboundaryLayer=.FALSE.` (Dimitris saw runaway ice with it TRUE).
- `SEAICE_CAP_ICELOAD` only caps the ice WEIGHT on Eta, not ice thickness. It will not change the monitored h_eff max.
- Era: 2008-2025, still current. Doc has a "known issues" section (PR #472).
- Src: mitgcm-support 2025-August 'Too large of Sea Ice'; 2020-May 'SEAICE, very thick ice'; 2019-February 'changes in seaice default behaviour?'; https://github.com/MITgcm/MITgcm/pull/472

### Blow-up (MOM_IMPLICIT_R / extreme potential temperature) only with EVP; LSR runs fine
- Cause: classic EVP is noisy/unstable; alpha/beta of 20 far too small; the original EVP is not what the manual recommends.
- Fix: use mEVP (`SEAICEuseEVPrev`/`SEAICEuseEVP` with `SEAICE_evpAlpha = SEAICE_evpBeta = 500`, `SEAICEnEVPstarSteps` ~200) or aEVP (tune `SEAICEaEVPcoeff`; `SEAICEaEVPcstar=4` default; see manual section "More stable variants of EVP"). Set `monitorFreq` to 1 to find where it explodes. Martin's recommendation remains the default Picard/LSR solver. Near ice-free conditions EVP can divide by tiny denomU/V: use the new `SEAICE_evpAreaReg` (>0, default -1) (PR #929, merged Nov 2025; the CPP flag name in the PR title, SEAICE_EVP_REGULARIZE_DENOMUV, is not in master).
- Era: 2007-2025. `lab_sea.hb87` now tests aEVP (issue #407 / PR #416). Rule of thumb for deltaTevp: Dimitris used 60 s for a 40-km Arctic; 5-10 s at 4-km.
- Src: mitgcm-support 2021-January 'Crashes with EVP seaice'; 2018-June 'EVP subcycling'; 2007-February 'EVP stability'; https://github.com/MITgcm/MITgcm/pull/929

### NaNs / blow-up at high resolution (O(1 km)) in summer melt with LSR; smaller deltaT helps only temporarily
- Cause: Picard solver under-converged (defaults `SEAICEnonLinIterMax=2` for LSR, `LSR_ERROR`); tile-edge noise in shear/divergence.
- Fix: `SEAICEnonLinIterMax = 10`, `LSR_ERROR = 1.e-5` (or 1e-6). The second is cheap and improves parallel behaviour. Do NOT switch to EVP for stability. Martin runs this everywhere.
- Era: 2014-2021. Upstream default `LSR_ERROR` is now 1e-5 (was 2e-4); `SEAICEnonLinIterMax` default is still 2 for LSR (issue #171 discussion: kept for run-time).
- Src: mitgcm-support 2021-April 'seaice, LSR_ERROR, EVP'; 2018-June 'EVP subcycling'; https://github.com/MITgcm/MITgcm/issues/171

### Stripes / "scars" in ice area, shear, divergence along MPI tile edges
- Cause: LSR is parallelised by restricted additive Schwarz on tiles; unconverged linear solve leaves tile-edge effects. Not a bug.
- Fix: `LSR_ERROR = 1.e-6` (and/or more `SEAICEnonLinIterMax`); or JFNK/Krylov solvers (expensive); test by `SEAICEuseDYNAMICS=.FALSE.` that dynamics is the cause. Tiles below ~30x30 also scale badly.
- Era: 2014-2019, same cause reported by several users.
- Src: mitgcm-support 2014-July 'seance leakage at the cpu domain boundaries'; 2015-November 'seaice anomalous advection in doubly periodic domain'; 2016-September 'sea ice diagnostics'

### NaN in first steps with SEAICEuseDYNAMICS=.TRUE. and no ice present (old code)
- Cause: LSR inverts a singular matrix when there is no ice (zero main diagonal); also min drag `DWATN` was hard-wired.
- Fix: update the code (fixed long ago). Since checkpoint67s the minimum drag is the run-time `SEAICEdWatMin` (default 0.25; can be 0 with `SEAICE_waterDrag=0`, but freedrift still needs non-zero drag, `SEAICE_CHECK` catches it). Making `SEAICE_LSRrelaxU/V` tiny "stabilises" only by not solving the system; default is under-relaxation <1.
- Era: checkpoint65u-ish and older (2016); fixed upstream.
- Src: mitgcm-support 2019-August 'Digest Vol 194, Issue 26'; 2019-November 'Digest Vol 197, Issue 15'; 2022-December 'mitgcm seaice'

### Sudden crash in seaice_growth ("Invalid operation PROG=seaice_growth") only with ALLOW_SITRACER on a vector machine
- Cause: aggressive optimisation on SX-ACE after PR #348 (not a model bug).
- Fix: add seaice_growth.F to NOOPTFILES in the optfile.
- Era: Oct 2020, checkpoint67s; machine retired. Useful as a hint when a harmless diagnostic change breaks one routine: lower optimisation for that file.
- Src: mitgcm-support 2020-October 'Bug in seaice code since checkpoint67s?'

### Ice thickness jumps 50 m in one step next to a freshwater (EmP) patch
- Cause: very large precip on cold surface with no snow handling and no downward LW (atmosphere at 0 K); snow cannot be advected/flooded.
- Fix: `SEAICEadvSnow=.TRUE.`, `SEAICEuseFlooding=.TRUE.`; provide `lwdown`; check precip magnitude (0.08 Sv over 25000 km2 is 3e-6 m/s).
- Era: 2016.
- Src: mitgcm-support 2016-March 'SEAICE pkg, unstable HEFF'

### Super-cooled water from ice-shelf cavities makes 40 m of ice / "extreme potential temperature" at cavity mouth
- Cause: all frazil heat is converted to ice in one step (`SEAICE_frazilFrac=1` default).
- Fix: `SEAICE_frazilFrac = SEAICE_deltaTtherm/(3 days)` (or set `SEAICE_gamma_t_frz`); `SEAICE_mcPheeTaper`; cap HEFF as last resort. Checks: `SEAICE_frazilFrac` and `SEAICE_mcPheePiston` must lie in range (see below).
- Era: 2016-2019.
- Src: mitgcm-support 2016-November 'on the control of sea ice growth'; 2019-February 'Increasing sea-ice/atmosphere drag'

### "SEAICE_CHECK: SEAICE_mcPheePiston is out of bounds ... must lie within 0. and drF(1)/SEAICE_deltaTtherm"
- Cause: piston velocity (m/s) cannot exceed surface-layer thickness per thermodynamic step; hit when Nr/drF(1) or deltaT changed.
- Fix: choose `SEAICE_mcPheePiston` < `drF(1)/SEAICE_deltaTtherm`, or leave unset (derived from `SEAICE_availHeatFrac` if given). Message still present in seaice_check.F.
- Era: 2020 (1D_ocean_ice_column with 60 levels).
- Src: mitgcm-support 2020-October 'Questions about the seaice pkg (SEAICE_mcPheePiston)'; 2012-August 'Turbulent Ice-Ocean heat flux'

### Ocean+ice model blows up at high wind (40 m/s) at the seafloor next to the ice edge
- Cause: plain ocean CFL/time-step limit, not sea ice.
- Fix: halve deltaT. Chris Hill/Martin saw identical symptoms in the Southern Ocean.
- Era: 2007 (checkpoint57y_post). Same advice for "model blows after minutes" in 2007-July (use deltaT 1800 first).
- Src: mitgcm-support 2007-May 'problem in pkg/seaice'

### Water below the freezing point at depth in coarse ice-ocean runs (-2.4 C)
- Cause: either advection undershoot with centred 2nd-order (J-M, Jeff Scott) or unbalanced surface heat loss feeding convection (Martin). `useOldFreezing` resets T to -1.9 everywhere and adds heat.
- Fix: use a monotone advection scheme (33, 77, 7; 80); check net surface heat flux balance of the coupled system; do not use `useOldFreezing`.
- Era: 2008-2013.
- Src: mitgcm-support 2008-March 'Water at the bottom of NA below freezing point'; 2013-April 'ice package in coupled atmosphere ocean tutorial (cpl_aim+ocn)'

## Dynamics, solver and defaults

### Ice does not move at all / starts late at very high resolution (dx of metres to 100 m), LSR, free drift
- Cause: poor convergence of RAS-parallelised LSR; ice touching a tile edge is "held". Also happens for an isolated floe/melange with no initial ice motion.
- Fix: in SEAICE_OPTIONS.h `#define SEAICE_ALLOW_FREEDRIFT`; in data.seaice `LSR_mixIniGuess = 2` (or 4); or `SEAICEuseFREEDRIFT=.TRUE.` (the CPP flag only compiles it). Needs wind/water forcing; no forcing means no motion. NB: PR text calls it `LSR_minIniGuess` in places; the namelist name is `LSR_mixIniGuess`.
- Era: 2020-2021, issue #327 (open).
- Src: https://github.com/MITgcm/MITgcm/issues/327 ; mitgcm-support 2020-January 'Ice not Moving when Resolution is Very High'; 2021-August 'advection of ice-melange'

### Floating-point exception (div by zero) in seaice_lsr.F with SEAICEscaleSurfStress=T and LSR_mixIniGuess=1
- Cause: new defaults (PR #116) scale surface stress by AREA; with no ice, the free-drift guess divides by zero.
- Fix: upstream handled it (areaW/S regularisation, issue #148 closed); on old code set `SEAICEscaleSurfStress=.FALSE.` or `LSR_mixIniGuess=0`.
- Era: Aug 2018, checkpoint67-era; fixed upstream.
- Src: https://github.com/MITgcm/MITgcm/issues/148

### Results changed after updating: ice thicker / different drag / snow mass
- Cause: defaults changed in PR #116 (checkpoint 67 era, issue #107): `SEAICEscaleSurfStress` F->T, `SEAICEaddSnowMass` F->T, `SEAICE_useMultDimSnow` F->T, `SEAICE_drag` halved to 0.001, `SEAICE_waterDrag` 5.5e-3 (now scaled by rhoConstFresh), `SEAICEetaZmethod=3`. Before 2017 stresses were not scaled by concentration at all.
- Fix: set the old values explicitly in data.seaice to compare. Verification experiments `global_ocean.cs32x15.icedyn` override them to keep old results.
- Era: 2017-2019; still the reference list of "why did my ice change".
- Src: https://github.com/MITgcm/MITgcm/issues/107 ; mitgcm-support 2019-February 'changes in seaice default behaviour?'; 2017-August 'sea ice stresses multiplied by concentration or not?'

### pkg/thsice + SEAICE: snow mass missing in ice load (snowHeight not mapped to HSNOW)
- Cause: thsice `snowHeight` was never copied to seaice `HSNOW` although `SEAICEaddSnowMass=.TRUE.` default uses it in seaice_dynsolver.
- Fix: upstream PR #811 (Feb 2024) maps it; to reproduce old (wrong) results set `SEAICEaddSnowMass=.FALSE.`.
- Era: fixed upstream in PR #811 (2024).
- Src: https://github.com/MITgcm/MITgcm/pull/811

### Metric terms in sea-ice stress divergence changed results near the pole / on curvilinear grids
- Cause: PR #976 (Feb 2026) added missing tensor metric terms (k1*sigma12 - k2*sigma22 etc.). Param evolved: `SEAICEuseMetricTerms` (T) is replaced by `SEAICEselectMetricTerms` (default 2; 0 = off; forced 0 on Cartesian grids).
- Fix: to recover old results `SEAICEselectMetricTerms = 0` (or `SEAICEuseMetricTerms=.FALSE.`). A draft name `SEAICEuseExtraMetricTerms` was used in discussion but is NOT in master.
- Era: fixed upstream 2026; users on older checkpoints lack the terms.
- Src: https://github.com/MITgcm/MITgcm/pull/976

### Compile error "Expecting END DO statement" in seaice_solve4temp.F
- Cause: mismatched DO/IF nesting in one block (gfortran).
- Fix: PR #953 (Dec 2025) re-nests the loops. Locally swap the ENDIF/END DO order.
- Era: fixed upstream PR #953.
- Src: https://github.com/MITgcm/MITgcm/pull/953

### Compile error in seaice_evp.F (multiple definition of nEVPstepMax) when ALLOW_AUTODIFF_TAMC undefined; STOP with SEAICE_deltaTdyn > SEAICE_deltaTtherm + autodiff
- Cause: TAF-specific ifdefs wrapped non-TAF code; local TAUX/Y undefined when dynamics step is skipped.
- Fix: PR #731 (May 2023) switches to ALLOW_AUTODIFF where appropriate. Using `SEAICE_deltaTdyn > SEAICE_deltaTtherm` with pkg/autodiff is disallowed.
- Era: fixed upstream May 2023 (issues #730, #731).
- Src: https://github.com/MITgcm/MITgcm/issues/730

### "SEAICE_CHECK"/readparms error: deltaTdyn vs deltaTtherm consistency
- Cause: need `SEAICE_deltaTtherm = dTtracerLev(1)` and `SEAICE_deltaTdyn` an integer multiple >= deltaTtherm. Early versions (2014) had `.LT.` where `.GT.` was meant; today the check is `deltaTdyn .LT. deltaTtherm` or non-integer ratio.
- Fix: leave `SEAICE_deltaTdyn` at the default or choose an integer multiple of deltaTtherm. With a 1-day tracer step setting deltaTdyn=1200 makes ice dynamics run 72x too slowly (acceleration term small, so often tolerable, but conceptually wrong).
- Era: 2014; check still in seaice_readparms.F.
- Src: mitgcm-support 2014-February 'CS64 grid'

### Ice seems to ignore SEAICE_waterDrag / drag has no effect in idealised runs
- Cause: with `SEAICEuseDYNAMICS=.FALSE.` the drag coefficients are never computed (drag = 0 regardless of Cd); also min drag `SEAICEdWatMin`.
- Fix: keep dynamics on for drag; set `SEAICEdWatMin = 0` on checkpoint >=67s.
- Era: 2020 (Kalyan Shrestha thread).
- Src: mitgcm-support 2020-May 'Seaice water drag'; 2022-December 'mitgcm seaice'

### Ice velocity locked to a constant in a doubly periodic / uniform-forcing domain
- Cause: no divergence in forcing, AREA->1 gives huge ice strength; ice moves as a slab. Not a bug.
- Fix: add divergent wind, an obstacle (see offline_exf_seaice), or reduce `SEAICE_strength` / `SEAICE_waterDrag` to test.
- Era: 2016.
- Src: mitgcm-support 2016-February 'funky ice dynamics in doubly periodic domain'

### Sea ice on cubed-sphere corners / Cartesian grid / periodic channel blows up
- Cause: metric terms were missing (until PR #976) and corner tiles never tested; llc grids rotate corners onto land. Older code recomputed grid itself and did not run on Cartesian grids.
- Fix: use llc, or `SEAICE_strength=0` + EVP free drift as hack; test with `SEAICEuseDYNAMICS=.FALSE.`; use a current version.
- Era: 2004-2008, partly superseded.
- Src: mitgcm-support 2008-July 'ice layer modeling -> cube corners?'; 2004-November 'seaice'

## Package interaction (OBCS, shelfice, thsice, exf)

### Ice piles up or forms spurious velocity at open boundaries (SEAICE + OBCS); how to relax ice at boundaries
- Cause: OBCS+seaice is "far from complete": boundary values only for uice/vice/heff/hsnow/area, not ice strength `press`; Stevens BCs ignore ice.
- Fix: keep ice away from boundaries; `OBCS_SEAICE_AVOID_CONVERGENCE`, `OBCS_SEAICE_COMPUTE_UVICE` (Neumann du/dn=0, only implemented for OBCS_UVICE_OLD, port to seaice_apply_uvice.F), `ALLOW_OBCS_SEAICE_SPONGE` (verification `seaice_obcs.sponge`); `SEAICEpressReplFac=0`; RBCS has no seaice hook (extending is possible, restoring velocity is pointless). Dimitris' practical suggestion: a relaxation zone around the edge.
- Era: 2006-2022, still true.
- Src: mitgcm-support 2022-November 'RBCS for sea ice boundary conditions'; 2020-May 'SEAICE, very thick ice'; 2013-June 'OBCS for sea ice!'; 2024-November 'Restoring sea ice using RBCS and SEAICE'

### Periodic seaice but OBCS for the ocean
- Cause: obcs_*ice* routines overwrite the periodic ocean exchange.
- Fix: make obcs_adjust_uvice.F, obcs_apply_seaice.F, obcs_apply_uvice.F and obcs_seaice_sponge.F return immediately (EXIT / CPP flag); model is periodic by default.
- Era: 2018.
- Src: mitgcm-support 2018-February 'OBCS for ocean, but periodic for sea ice?'

### Thin ice-shelf (hFacC(k=1)>0) with seaice or thsice fails / shortwave warms ocean under shelf
- Cause: seaice masks assumed ice shelf fills the whole surface cell.
- Fix: PR #187 (thsice, shortwave) and #188 (seaice) prevent ice formation where there is shelf ice. Upstream since Jan 2019.
- Era: fixed upstream; old code needs these PRs.
- Src: https://github.com/MITgcm/MITgcm/pull/188 ; https://github.com/MITgcm/MITgcm/pull/187

### Large Eta / sea-ice loss next to a shelfice cavity (EVP/LSR regional, SHELFICE + SEAICE)
- Cause: inconsistent `SHELFICEloadAnomalyFile` vs T/S in cavity; Eta compensates, expels sea-ice outside; Eta is non-zero under shelves by design.
- Fix: study verification/isomip; recompute pressure load with the same T/S as initial conditions; see PR #251.
- Era: 2020.
- Src: mitgcm-support 2020-February 'Strange output fields/discontinuities'

### seaice + thsice together: initial ice thickness set in data.seaice has no effect
- Cause: thsice thermodynamic state overwrites seaice's.
- Fix: `thSIceThick_InitFile` in data.ice. Combination = seaice dynamics + thsice thermodynamics, just enable both in data.pkg (verification global_ocean.cs32x15/input.icedyn). Also `thsice` regularisations: PR #633 fixed `ustar` minimum (5e-3) and removed qicen/hlyr regularisation; `THSICE_REGULARIZE_CALC_THICKN` replaces `THICKNESS_AUTODIFF_REGULARIZE`.
- Era: 2011-2022; fixed upstream PR #633.
- Src: mitgcm-support 2020-December 'Questions about the thsice pkg'; 2011-February 'SEAICE and THSICE'; https://github.com/MITgcm/MITgcm/pull/633

### pkg/seaice cannot be used in coupled atmosphere (cpl_aim+ocn); thsice instead
- Cause: seaice requires exf and an ocean-based forcing; in cpl_aim+ocn sea ice is entirely done by the atmosphere component.
- Fix: use thsice (optionally seaice dynamics); see verification_other/offline_cheapaml `input.dyn` for cheapaml+seaice+thsice. Over ice, exf fluxes are modified, so use EXF_OPTIONS.h option (3): `ALLOW_ATM_TEMP`, `ALLOW_DOWNWARD_RADIATION`, `ALLOW_BULKFORMULAE`; with thsice, `EXF_READ_EVAP` stops in thsice_get_exf.F.
- Era: 2013-2018 (still true).
- Src: mitgcm-support 2015-April 'coupling issues'; 2017-October 'Prescribed atmospheric temperatures and freshwater flux with exf/seaice/thsice packages'

### NLFS=4 (nonlinear free surface) with plain z-coordinates + seaice: allowed? (surface cell can dry)
- Cause: top-cell thickness vs ice thickness (ice loads the surface).
- Fix: allowed; use top cell >= 5 m or z*; `sIceLoadFac = 0.` (in `data` PARM01, new) restores "levitating" ice that does not load the ocean.
- Era: Sept 2026.
- Src: mitgcm-support 2026-September 'NLFS=4 with z coordinates.. and sea ice'

### Imposing observed sea-ice (concentration from file) instead of computing it
- Cause: no direct option.
- Fix: `#define EXF_SEAICE_FRACTION` (EXF_OPTIONS.h) reads `exf_iceFraction` (areamask* params; relaxes AREA via d_AREAbyRLX / d_HEFFbyRLX), but is little used and you must adapt how it enters. Older idea: read AREA via exf and overwrite in seaice_model before seaice_growth. The old `SEAICE_ALLOW_AREA_RELAXATION` flag no longer exists.
- Era: 2012-2022.
- Src: mitgcm-support 2022-February 'imposing sea ice'; 2016-August 'forcing the model with offline sea-ice'

### Ice drift under a stationary ocean / stand-alone sea ice
- Cause: seaice cannot run alone.
- Fix: verification/offline_exf_seaice (1-layer, no dynamics, constant ocean).
- Era: 2015-2016.
- Src: mitgcm-support 2016-June 'Stand-alone sea ice model'; 2015-July 'oceanic forcing for sea ice package'

## Initialisation, restarts, parameter handling

### Initial ice area from HeffFile ignored / restart ignores SEAICE_initialHEFF
- Cause: `HeffFile` is thickness, not concentration; AREA is set to 1 wherever HEFF>0 unless AreaFile given. On restart with pickup the pickup_seaice values are used (nIter0 selects them automatically).
- Fix: provide both `HeffFile` and `AreaFile`; for restarts nothing else to do. HEFF = mean thickness x AREA (volume = HEFF*RAC), actual thickness = HEFF/AREA.
- Era: 2006-2016.
- Src: mitgcm-support 2006-September 'how to initialize seaice package'; 2009-June '(no subject)'; 2016-August 'sea ice volume conservation'

### Typo in data.seaice silently ignores everything after it
- Cause: old IOSTAT test `errIO .LT. 0`.
- Fix: update; test is `.NE. 0` in modern code. Check STDOUT parameter echo.
- Era: 2005, fixed.
- Src: mitgcm-support 2005-July 'parameter file read'

### Old SEAICE_OPTIONS.h / darwin code fails to compile after update ("COMMON block data object must not be an automatic", SEAICE.h)
- Cause: header split: SEAICE.h, SEAICE_SIZE.h, SEAICE_PARAMS.h.
- Fix: start from current pkg/seaice/SEAICE_OPTIONS.h and re-apply local edits; include SEAICE_SIZE.h wherever SEAICE.h is.
- Era: 2015, checkpoint65 era.
- Src: mitgcm-support 2015-July 'old & new sea ice model'

### Interpolating pickup_seaice to a coarser grid
- Fix: area-average fields that already carry AREA (HEFF, HSNOW, HSALT); average TICES weighted by AREA*RAC; area-average UICE/VICE with RAW/RAS weights.
- Era: 2011.
- Src: mitgcm-support 2011-August 'interpolate sea ice fields'

## Diagnostics and budgets

### Which diagnostic is stress on ocean under ice; heat/FW budgets that do not close
- Cause: several overlapping diagnostics.
- Fix: `oceTAUX/oceTAUY` = `SIfu/SIfv` = stress felt by the ocean (cell-averaged, includes ice damping); `EXFtaux/y` = stress before ice (on tracer points with bulk formulae). Atmosphere-ice stress: `SItaux/SItauy`. `oceQnet` = -`SIqnet` (positive down vs up); ice export: `ADVxHEFF/ADVyHEFF/ADVxSNOW/ADVySNOW` (volume per area, times rhoIce/rhoSnow). Thermodynamic growth = `SIdHbOCN+SIdHbATC+SIdHbATO+SIdHbFLO`. Total evaporation: `1000*EXFevap + SIfwSubl` (1000 = rhoConstFresh); atm FW into system: `-(SIatmFW - rhoFresh*(EXFpreci+EXFroff))`. Stress-related diagnostics (`SIsig1/2`, `SIpress`, `SIzeta`...) are filled at the end of seaice_dynsolver (fixed in 67s: 1st step was zero before). Internal stress divergence is not carried between steps: add diagnostics at the end of seaice_dynsolver.F.
- Era: 2013-2024; EXFevap contains only ocean evaporation since PR #203 (June 2019).
- Src: mitgcm-support 2016-September 'funky ice dynamics...' (oceTAUX list); 2024-August 'Calculating sea ice export'; 2019-July 'Evaporation in sea ice regions'; 2024-October 'Diagnosing internal ice stress with the seaice pkg'; 2024-August 'Diagnostic for Heat flux from ocean to sea ice'

### SEAICE_USE_GROWTH_ADX: heat/mass budget not closed, +25 C spike when snow floods
- Cause: seaice_growth_adx.F included flooding conversion (snow->ice, no energy) in the heat that goes to the ocean and sublimation in diagnostics SIatmFW.
- Fix: PR #703 (Feb 2023) and #721 (Apr 2023): flooding moved after QNET/EmPmR; budget closure for linear free surface; diagnostic `SIaaflux` replaced by `SIacflux`.
- Era: fixed upstream 2023; ECCO v4r5/r6-type configs on older code have it.
- Src: https://github.com/MITgcm/MITgcm/pull/721 ; https://github.com/MITgcm/MITgcm/pull/703

### Salt/FW from sea ice: SEAICE_saltFrac and seaice_salt0 do nothing
- Cause: `SEAICE_saltFrac` only matters with `SEAICE_VARIABLE_SALINITY` defined (HSALT advected, released on melt); `SEAICE_salt0` only when it is undefined. Sea-ice tracers (`ALLOW_SITRACER`, verification lab_sea input.salt_plume) are the more consistent route.
- Era: 2022.
- Src: mitgcm-support 2022-August 'seaice_saltFrac and SEAICE_VARIABLE_SALINITY'; 2017-February 'Ptracers in sea ice'

### Track or zero the freshwater from ice melt
- Fix: only place sea ice modifies FW flux is seaice_growth.F (`tmpscal1*convertHI2PRECIP*rhoConstFresh`, ~line 2357 in 2020-22 code); zero or save it there (sign <0 = melt) and pass to ptracers_forcing_surf.F.
- Era: 2020-2022.
- Src: mitgcm-support 2022-April 'Track freshwater flow from sea ice melting'; 2020-October 'Calculation of melting ice water in seaice pkg'
