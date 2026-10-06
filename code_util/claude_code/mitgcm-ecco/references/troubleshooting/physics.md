# Troubleshooting: physics (free surface, EOS/TEOS, shelfice, RBCS, OBCS/tides, drag/viscosity/Coriolis, advection, cg2d/cg3d)
Distilled from answered mitgcm-support threads (2003-2026) and MITgcm GitHub issues/PRs. Names were checked against the current-source index and `origin/master` (HEAD ae4c03af7, 2026-10-03); "Era" says when a bug was fixed upstream or a name changed. mitgcm-support month URLs are `http://mailman.mitgcm.org/pipermail/mitgcm-support/<YYYY-Month>/thread.html`.

## Solver convergence (cg2d / cg3d)

### cg2dTargetResWunit / cg3dTargetResWunit gives a tolerance that is too loose (off by sqrt(Nxy))
- Cause: the W-unit to cg2dTolerance conversion used globalArea^2 instead of rAc*globalArea. cg3d had a similar factor of sqrt(Nxy/Nz). Also, with implicDiv2DFlow<1 the W-unit tolerance should include that factor (not done).
- Fix: on master (PR #959; new GRID.h globals rAc_3dMean, n2dWetPts, n3dWetPts) the scaling is corrected. Old runs that used cg2dTargetResWunit must retune it (~sqrt(Nxy) smaller) to reproduce results. Or avoid W-units: set `cg2dTargetResidual` (normalised RHS; cg2dNormaliseRHS is simply cg2dTargetResWunit<=0). In adjoint runs prefer cg2dTargetResidual: the adjoint CG acts on fields that don't fit W scaling (JM-C/Ou Wang).
- Era: seen until 2025-12; fixed in PR #959 (2026-01). Checkpoints before that (incl. ECCO checkpoint68g/68o era) still have the wrong factor.
- Src: https://github.com/MITgcm/MITgcm/issues/956 ; https://github.com/MITgcm/MITgcm/pull/959 ; https://github.com/MITgcm/MITgcm/issues/860 (beta factor discussion)

### cg2d converges "too fast"/imprecisely after shifting initial SSH (cg2dNormaliseRHS=T)
- Cause: RHS normalisation = 1/max|rhs|, which is not tied to the initial residual. A large uniform SSH offset (e.g. pSurfInitFile + 100 m) inflates rhs, so the relative target is met in fewer iterations with a larger absolute error.
- Fix: use absolute criterion cg2dTargetResWunit (cg2dNormaliseRHS=F) when initial SSH has large offsets, or remove the mean offset. cg3dNormaliseRHS=F option was not coded when reported.
- Era: reported 2024-08, open.
- Src: https://github.com/MITgcm/MITgcm/issues/856

### cg2d: NaN / "STOP ... MOM_IMPLICIT_R (or GAD_IMPLICIT_R): error when solving 3-Diag problem" / "SOLUTION IS HEADING OUT OF BOUNDS ... exceeds allowed range (monSolutionMaxRange=1.E+03)" / "IEEE_INVALID_FLAG"
- Cause: model blew up (almost never an MPI/solver problem). Usual suspects, in thread order: CFL (advcfl_* > ~0.5, esp. vertical with thin cells, tides adding currents), initial shock, hFacMin/hFacMinDr silly (hFacMin must be <=1; ~0.1-0.3), bad bathymetry sign/units, input file endianness/precision (initial values ~1e32), a NaN in an exf forcing file (backtrace in exf_bulkformulae), unbalanced OBCS, explicit viscosity/diffusivity above stability (deltaT < dz^2/(8 viscAr); viscAhGridMax<=1).
- Fix: monitorFreq=deltaT (or small) and grep `advcfl_` and `dynstat_*` in STDOUT; reduce deltaT until initial adjustment passes, then raise; dump first steps; build with `-devel`/`-ieee` and read the backtrace; strip packages and re-add one by one; check STDERR (the "EEDIE: Only 0 threads completed" tail is just the stop being improper). monSolutionMaxRange (pkg/monitor, default 1e4 on master; digest shows 1e3 in 2023) is the range check that stops the run.
- Era: all years (2004-2025).
- Src: mitgcm-support 2025-April 'the floating-point exception issue' ; 2023-October 'Problem with tmin, tmax and monSolutionMaxRange' ; 2020-September 'STOP MOM_IMPLICIT_R: error when solving 3-Diag problem.' ; 2014-July 'Model based on exp4 experiment - stability problems' ; 2004-August 'solution is heading out of bounds' (all .../YYYY-Month/thread.html)

### cg2d hits cg2dMaxIters and does not converge (flat-bottom/idealised box, tides, rigid lid) then blows up
- Cause: poor conditioning of the 2-D problem for simple domains; also implicit free surface is far cheaper than rigid lid (rigid lid took >1500 iterations; implicitFreeSurface=T typically <~100 after the first steps for 1e-7).
- Fix: `cg2dUseMinResSol=1` in PARM02 (keeps the minimum-residual iterate instead of the last one; fixed a rigid-lid NH blow-up); tighten cg2dTargetResidual; try exactConserv=T if it is a convergence issue (it can change results when solver is poorly converged); prefer implicitFreeSurface=T.
- Era: 2007, 2011, 2013-2014 (parameter still exists).
- Src: mitgcm-support 2013-December 'Linear internal waves evolution under a free surface and a rigid lid.' ; 2007-April 'Help with CG2D settings' ; 2011-April 'Spurious internal buoyancy minima in long runs'

### cg2d is the scaling bottleneck (GLOBAL_SUM_TILE_RL dominates, SOLVE_FOR_PRESSURE fraction grows with nProcs)
- Cause: cg2d does 2 exchanges + 2-3 global sums per iteration; with GLOBAL_SUM_ORDER_TILES (bit-reproducible sum) the sum is costly. Odd/prime per-tile sizes (e.g. 139x283) also scale badly; halos dominate when tiles are small.
- Fix: (a) `#undef GLOBAL_SUM_ORDER_TILES` in eesupp/inc/CPP_EEOPTIONS.h (result then depends on tile/proc count); (b) single-reduction CG: `ALLOW_SRCG` (now #define by default in CPP_OPTIONS.h) + `useSRCGSolver=.TRUE.` in PARM03 (user saw ~15% saving); (c) relax cg2dTargetResidual/cg2dMaxIters if physics tolerates; (d) choose factorizable domain (e.g. 558x3400 as 18x17 or 18x20 tiles); (e) GLOBAL_SUM_SEND_RECV was *worse* on one OpenMPI cluster; (f) check cg2d_iters first.
- Era: 2007-2025.
- Src: https://github.com/MITgcm/MITgcm/issues/530 ; mitgcm-support 2025-February 'Reducing Runtime for High Resolution Model' ; 2009-July 'inefficient pressure solver'

### Non-hydrostatic run 10x slower than hydrostatic; cg3d dominates (SOLVE_FOR_PRESSURE ~ CG3D)
- Cause: cg3dMaxIters/cg3dTargetResidual=1e-13 with large Nx*Ny*Nr; NH is only worthwhile when dx ~ dz (<~100 m), not at 2-deg or eddy-resolving scales.
- Fix: nonHydrostatic=.FALSE. unless dx~dz; otherwise cap `cg3dMaxIters` (e.g. 40-100; keep > a couple of Nr) and loosen cg3dTargetResidual (1e-13 -> ~1e-8..1e-10). Don't ask the 3-D solver to converge on noisy velocity fields. cg3dMaxIters is ignored when nonHydrostatic=F.
- Era: 2013-2017.
- Src: mitgcm-support 2016-July 'Hi' thread (cg3d time) / 2013-September 'Runtime for non-hydrostatic model' ; 2017-August 'The criteria of cg3dMaxIters'

### Adjoint: tape/recompute problems and wrong cg2d adjoint with nonlinFreeSurf>2 or _RS=real*4
- Cause: (a) pre-PR #392 flow directives for self-adjoint cg2d were wrong for nonlinFreeSurf>2 (extra storage + new hand-written cg2d adjoint added); (b) CG2D_I_R_AD common block typed _RL instead of _RS -> floating overflow in update_cg2d_ad.f with cg2dFullAdjoint=T and _RS=real*4 (#526/#527); (c) pkg/rbcs has no level-N checkpointing -> endless recomputation; (d) viscFacInAd had no effect with vectorInvariantMomentum until #384 (needs `AUTODIFF_ALLOW_VISCFACADJ`, now defined in AUTODIFF_OPTIONS.h, and viscFacInFw).
- Fix: use a checkpoint after those PRs (2020-11 / 2021-08); for ECCO-style rStar adjoints make sure code_ad has current checkpoint*_directives.h; set cg2dTargetResidual rather than W-units (see above).
- Era: 2006-2021; fixed in PR #392, #527, #384.
- Src: https://github.com/MITgcm/MITgcm/pull/392 ; https://github.com/MITgcm/MITgcm/issues/526 ; https://github.com/MITgcm/MITgcm/pull/384 ; mitgcm-support 2006-July 'Hi' (rbcs adjoint)

## Free surface, volume and freshwater conservation

### "CONFIG_CHECK: nonlinFreeSurf not yet implemented in nonHydrostatic code" (and "nonHydrostatic NOT SAFE with non fully implicit barotropic solver")
- Cause: nonlinFreeSurf/rStar and the 3-D solver (nonHydrostatic or implicitIntGravWave) are mutually exclusive; the second message is the check on implicitNHPress*implicSurfPress*implicDiv2DFlow != 1 (Crank-Nicolson settings as in verification/internal_wave).
- Fix: either nonlinFreeSurf=0 (identical to #undef NONLIN_FRSURF) with exactConserv=T, or hydrostatic. For NH set implicSurfPress=implicDiv2DFlow=1 (comment out the CN lines). NH tide/internal-wave work: linear free surface + `selectNHfreeSurf` (JM-C, experimental), tune hFacMin/surface-layer thickness so the surface layer cannot drain.
- Era: 2005-2026, still a STOP in config_check.F.
- Src: mitgcm-support 2006-December 'nonlinear free surface' ; 2006-September 'OBC problem: spurious boundary jets with C-D coupling' thread (config_check) ; 2009-July '#undef nonlinear free surf and nonhydrostatic'

### useRealFreshWaterFlux with implicDiv2DFlow < 1 (or selectAddFluid=1) gives wrong volume
- Cause: cg2d_b misses the (1-implicDiv2DFlow)*EmP^(n-1) term unless exactConserv=T; addMass was added to cg2d_b/cg3d_b without the implicDiv2DFlow factor (double counting with exactConserv=T).
- Fix: use exactConserv=T whenever implicDiv2DFlow<1 (config_check now STOPs otherwise, "implicDiv2DFlow < 1 requires exactConserv= T"). Fix for addMass in PR #869.
- Era: 2024-2025.
- Src: https://github.com/MITgcm/MITgcm/issues/860

### addMassFile / ALLOW_ADDFLUID shows no effect, or surface mass source "does nothing"
- Cause: ALLOW_ADDFLUID + selectAddFluid=1 is for sources inside the water column. For surface-only mass input it is the wrong tool. Momentum is not changed (source assumed to carry local velocity); temp_addMass/salt_addMass default UNSET -> neutral (local T,S).
- Fix: for surface input use EmPmR/EXF precip with `useRealFreshWaterFlux=.TRUE.`, `rigidLid=.FALSE.`; for add-mass set selectAddFluid=1, addMassFile (units kg/s per cell), and temp_addMass/salt_addMass if the source has a fixed T/S (energy reference: fresh water at freezing point). Note the implicDiv2DFlow issue above.
- Era: 2010-2025.
- Src: mitgcm-support 2025-April 'About "addmass" in MITgcm' ; 2010-August 'Hi' thread (mass flux/runoff) ; 2011-August 'ALLOW_ADDFLUID option and U, V momentum equation'

### Sea level / mean salinity drifts in a regional or lake domain (SSH decreasing, dynstat_eta_mean drifting)
- Cause: net E-P-R not zero, unbalanced OBCS volume transport, or SSS restoring adding volume with useRealFreshWaterFlux=T.
- Fix: look at dynstat_eta_mean and dynstat_salt_min in STDOUT; set `balanceEmPmR=.TRUE.` (or tune a multiplicative factor on precip offline); `useOBCSbalance=.TRUE.` with OBCS_balanceFacE/W/N/S (nonlinear FS: balancing is applied before tides are added); balance inflow/outflow offline for NH runs; for very small domains make the boundary far away/sponged. Restoring salinity alone cannot fix an EmPmR imbalance.
- Era: 2013-2025.
- Src: mitgcm-support 2013-March 'Underestimated ETA values for Global Ocean!' ; 2016-April (Bay of Bengal salinity loss, 'Hi' thread) ; 2018-August 'useOBCSbalance not working with RBCS for salinity' ; 2014-November 'Running off of the seawater'

### Very negative SSH (-20 m...), "WARNING: hFac < hFacInf", "STOP in CALC_R_STAR : too SMALL rStarFac", "ABNORMAL END: S/R CALC_R_STAR"
- Cause: nonlinear free surface + surface layer too thin vs |eta| (hFacInf<0.1 does not rescue it), or sea ice piling up in inlets and loading the surface (useRealFreshWaterFlux=T: ice is "submerged"), or large EmPmR.
- Fix: `select_rStar=2, nonlinFreeSurf=4, hFacInf=0.2, hFacSup=2.0` (only nonlinFreeSurf=4 is the complete scheme); surface layer >> expected |eta| (wave height < ~20% of drF(1)); drF(k+1)/drF(k) < 1.4; for sea ice use `SEAICE_no_slip=.FALSE.` and/or `#define SEAICE_CAP_ICELOAD`; smooth/fill the coast where ice accumulates; balanceEmPmR. Setting `linFSConserveTr=T` helps when diagnosing tracer non-conservation in the surface cell (linear FS + realFreshWater). If you only want a linear FS: drop select_rStar/nonlinFreeSurf/hFacInf/hFacSup.
- Era: 2006-2022.
- Src: mitgcm-support 2020-May 'SSH problems with nonlinFreeSurf and Star' ; 2014-July 'Hi' (SSH < -20 m thread) ; 2012-September 'Non-linear free-surface and vertical resolution' ; 2017-March 'cheapAML+global_ocean.90x40x15 - salty ocean'

### Surface-cell tracer not conserved / "layers" surface term wrong (linear FS + useRealFreshWaterFlux)
- Cause: without linFSConserveTr the surface cell loses/gains tracer when eta changes (w(k=1)!=0); pkg/layers called layers_wsurf_tr regardless of LinFSConserveTr. linFSConserveTr affects only T,S (not ptracers).
- Fix: `linFSConserveTr=.TRUE.` for linear FS + real FW flux; pkg/layers bug fixed in PR #988 (2026).
- Era: 2010; PR #988 merged 2026-05.
- Src: mitgcm-support 2010-June 'Mass balance... again...' ; https://github.com/MITgcm/MITgcm/pull/988

### Implicit vertical advection with limiter schemes crashes in CALC_R_STAR (e.g. VertAdvScheme=33 + ImplVertAdv=T on llc90)
- Cause: all schemes that are nonlinear in the tracer (every limiter scheme) behave badly with the current implicit vertical advection; undocumented. Also tempImplVertAdv + 33 was linked to blow-ups in a rigid-lid test.
- Fix: use tempVertAdvScheme=3 (or linear) with tempImplVertAdv=T, or explicit vertical advection with limiter schemes. (Warning/doc requested.)
- Era: 2013, 2020 (open issue #331).
- Src: https://github.com/MITgcm/MITgcm/issues/331 ; mitgcm-support 2013-December 'Linear internal waves...'

### Eta below column depth / dry-looking cells at ebb tide
- Cause: no wetting/drying in the standard model; with linear FS a wet column stays wet with all Nr levels; eta can drop below the top level or below the depth.
- Fix: none in core (use nonlinFreeSurf=4/rStar for large amplitude; for NH only linear FS is allowed). Wetting-drying is a package job.
- Era: 2018.
- Src: mitgcm-support 2018-January 'Simple problems about the output of salinity and eta'

### Tides come out damped / weak surface gravity waves with implicitFreeSurface
- Cause: backward-in-time (fully implicit) barotropic solver damps external waves, tide phase/amplitude wrong, particularly at long deltaT.
- Fix: Crank-Nicolson barotropic: `implicSurfPress=0.5, implicDiv2DFlow=0.5` (as in verification/internal_wave; 0.6 is the dissipative compromise) with exactConserv=T; not allowed with NH unless both=1 (see above). Otherwise reduce deltaT, prefer rStar if min depth >~20 m.
- Era: 2007-2021.
- Src: mitgcm-support 2007-April 'Hi' (damped 2 Hz waves) ; 2017-June 'Temperature Advection Scheme : Internal Gravity Wave Simulation' ; 2009-June 'Surface gravity waves?'

## Equation of state, pressure, TEOS-10

### Which pressure goes into the EOS (JMD95Z vs MDJWF vs TEOS10): selectP_inEOS_Zc
- Cause: JMD95Z/UNESCO use p=-rhoConst*g*z; JMD95P/MDJWF use hydrostatic pressure (lagged one step). Boussinesq energy conservation is lost if full pressure enters EOS, so JMD95Z is the safe default for incompressible Boussinesq (Paola Cessi).
- Fix: set `selectP_inEOS_Zc` in PARM01: 0 = -g*rhoConst*z, 1 = pRef integral, 2 = hydrostatic dynamical pressure, 3 = full (Hyd+NH, requires nonHydrostatic). Defaults follow eosType (MDJWF -> 2, JMD95Z -> 0); value printed in STDOUT. eosType='MDJWF' + selectP_inEOS_Zc=0 works. Switching a pickup from JMD95Z to MDJWF/JMD95P needs phiHyd added to the pickup (and fldList in .meta). N2 offline: use diagnostic DRHODR (N2=-g/rho0*DRHODR), not density at depth. Don't use UNESCO (takes in-situ T; model carries theta).
- Era: 2005-2024.
- Src: mitgcm-support 2024-March 'Nonlinear equations of state' ; 2022-September 'Diagnostic calculation of buoyancy frequency with nonlinear equation of state' ; 2013-January 'restart from pickup when using MDJWF eos'

### eosType='TEOS10': FIND_ALPHA/FIND_BETA derivative errors; surface heat flux bias from theta=CT
- Cause: THETA/SALT are interpreted as conservative T/absolute S, so bulk formulae saw CT instead of in-situ SST (PR #812). Separately, find_alpha.F had wrong terms in the KPP-relevant expansion: line `sa*2.*teos(16)` should be `sa*teos(16)`, `teos(42)` should be `sa*teos(42)`, FIND_BETA `p*teos(46)` -> `p*ct*teos(46)`; GSW constants lacked `_d 0` (precision).
- Fix: on master PR #812 (new surface-level temperature for exf/bulk formula; sstExtrapol stays in pkg/exf) and PR #1027 (derivative + `_d 0` fixes, 2026-08). Only global_ocean.cs32x15.in_p tests TEOS10 (not KPP). Pre-fix checkpoints: do not use TEOS10 with KPP or with bulk formulae. TEOS-10 still treats salt as absolute salinity, no preformed-salinity switch.
- Era: 2015-2019 "outdated/incomplete" (issue #115); fixed 2024-06 (#812) and 2026-08 (#1018/#1027).
- Src: https://github.com/MITgcm/MITgcm/pull/812 ; https://github.com/MITgcm/MITgcm/issues/1018 ; https://github.com/MITgcm/MITgcm/pull/1027 ; mitgcm-support 2019-July 'Hi' (TEOS-10 status) ; https://github.com/MITgcm/MITgcm/issues/115

### Atmospheric loading (pLoad / apressure) enters the EOS as full pressure
- Cause: with selectP_inEOS_Zc>=2 the EOS needs sea pressure (absolute - 1 atm); pLoad held full atmospheric pressure.
- Fix: PR #333 (2020-04) adds `surf_pRef` (default 1 atm=101325 Pa); exf subtracts it when filling pLoad; apressure stays full for air density. Old checkpoints: define pLoad as anomaly yourself.
- Era: reported 2019-01; fixed 2020-04.
- Src: https://github.com/MITgcm/MITgcm/issues/199

### NaN/blow-up with extreme T/S initial conditions (hypersaline lake, paleo 130 psu, huge T/S range)
- Cause: JMD95/MDJWF/TEOS10 are fits for ocean ranges (errors ~1% at S=70); sharp S gradients give huge density gradients.
- Fix: smooth initial fields, short initial deltaT, try eosType='LINEAR' first; for salinities well outside the range implement an EOS (model/src/find_rho.F, find_alpha.F, find_beta.F; add parameter in PARAMS.h + set_defaults.F + ini_parms.F namelist + config_summary.F). Fill land points with finite values (not NaN) in input files.
- Era: 2015-2023.
- Src: mitgcm-support 2023-February (T,S initial NaN thread) ; 2015-December 'Equation of state for extreme paleo scenarios' ; 2021-December 'How to add a new equation of state in mitgcm'

### Linear EOS: RHOAnoma zero / weird, "dRho", tRef semantics
- Cause: with eosType='LINEAR' density anomaly = rhoNil*(-tAlpha*(T-tRef)+sBeta*(S-sRef)); tRef(k)/sRef(k) also initialise T,S. If initial T == tRef the anomaly is 0. rhoNil (EOS) and rhoConst (momentum scaling) are different. Constant tAlpha mismatch: do not use in-situ densities to build tRef (N2 too large, blow-up).
- Fix: initialise T via hydrogThetaFile and keep tRef constant if you want a gradient to show; with sBeta=0 and saltStepping=.FALSE. use T as the buoyancy variable; check rhoNil/rhoConst in STDOUT; INCLUDE_PHIHYD_CALCULATION_CODE must stay #defined (CPP_OPTIONS.h) or no hydrostatic pressure is computed and stratification has no dynamical effect.
- Era: 2004-2013.
- Src: mitgcm-support 2008-February 'Stratification simulations' ; 2014-January 'MITgcm-support Digest, Vol 127, Issue 7' ; 2005-June (rhoNil vs rhoConst thread) ; 2004-June 'Hi' (tRef from rho)

## Ice shelves (pkg/shelfice, frazil, GM/Redi, mixing)

### shelfice_thermodynamics.F NaN: freshWaterFlux = shiTransCoeffS*(1 - sLoc/saltFreeze) with saltFreeze = 0
- Cause: sLoc (or the quadratic's saltFreeze) = 0 under cavities (initial melt undershoot) -> 0/0. Also a minor bug for SHELFICEadvDiffHeatFlux=T in the 3-eq solve.
- Fix: PR #968 (2026-01) diagnoses freshWaterFlux from the heat balance by default (as pkg/icefront and steep_icecavity); old salt-balance form via `#define SHI_SALTBAL_FWFLX` in SHELFICE_OPTIONS.h (now division-by-zero guarded). Locally: `IF (saltFreeze .EQ. 0. _d 0)` branch. If salinity goes to 0 in a cavity something else (overshoot) is wrong.
- Era: llc4320 2026-01; fixed upstream in PR #968.
- Src: https://github.com/MITgcm/MITgcm/pull/968 ; mitgcm-support 2026-January 'shelfice_thermodynamics.F crash'

### Supercooled water (theta ~ -100 degC) in cavity cells, sea ice growing next to cavity
- Cause: SHELFICEboundaryLayer=.TRUE. (with frazil on) let water supercool in isolated cells; with SHELFICEuseGammaFrict=T u* floors at 1e-3 m/s (hard-wired) so exchange never fully vanishes. pkg/frazil only shuffles heat and does not fix it.
- Fix: `SHELFICEboundaryLayer=.FALSE.` and instead `pCellMix_select=22` (enhanced mixing at bottom and under ice; `20` = under ice only; 2 = bottom only) with pCellMix_viscAr/diffKr; drop pkg/frazil; fill 1-cell holes in ice-shelf topography; fix bathymetry/draft errors if pCellMix_select=22 crashes. Dimitris asked for a warning on the BL option.
- Era: 2025-08 (llc1080/llc4320 cavities); no upstream change.
- Src: mitgcm-support 2025-August 'Supercooled waters in ice shelf cavities'

### Ice-shelf masking wrong where draft is thinner than top cell / Ro_surf ~ 1e-13 (kTopC = kSurfC in open ocean); remeshing cannot merge thin ice
- Cause: kTopC derived from Ro_surf after hFacMin/hFacMinDr lopping; Ro_surf can be O(1e-13) negative; shelfice_thermodynamics tests pLoc (R_shelfIce) not hFacC; remeshing of very thin ice.
- Fix: ensure draft >= ~0.5*drF(1) (set hFacMinDr >= drF(1) so thin ice is removed consistently), fill/clip drafts offline; use `SHI_update_kTopC` (PR #465/#379 fix, kTopC/shelficeMass stored for TAF) when remeshing; checkpoint after PR #312 for the obcs_check mismatch ("S/R OBCS_CHECK: Inside Mask and OB locations disagree", hFacMin=0.1 at ice/bedrock junction; workaround hFacMin=0.11-0.2).
- Era: 2018-2023; #99, #324 (PR #312), #379 fixed; #798 (Ro_surf round-off) opened 2023-11.
- Src: https://github.com/MITgcm/MITgcm/issues/99 ; https://github.com/MITgcm/MITgcm/issues/324 ; https://github.com/MITgcm/MITgcm/issues/379 ; https://github.com/MITgcm/MITgcm/issues/798

### pkg/gmredi under ice shelves loses/gains tracer; 'linear' GM taper extrapolates through ice; KPP/GGL90 under shelves
- Cause: gmredi_calc_psi_b masked psiX/psiY wrongly at dry-top points; gmredi_x/ytransport applied Kuz/Kvz wrongly when GM and Redi kappa differ. GGL90 had no shelfice logic; KPP always uses surface flux at k=1.
- Fix: PR #593 (2022-02): corrected masking (skew-flux form with same kappa for GM and Redi is unaffected). Use GM_taper_scheme other than 'linear' under cavities ('gkw91'/'dm95'); GM advective form can blow up under ice. GGL90 fixed for shelfice in PR #597 (2022); pkg/kpp still unsupported under ice shelves (issue #588 open).
- Era: 2016-2023.
- Src: https://github.com/MITgcm/MITgcm/issues/591 ; https://github.com/MITgcm/MITgcm/issues/588 ; mitgcm-support 2016-March 'Hi' (shelfice and GMREDI)

### Turn off / prescribe ice-shelf thermodynamics; remove conduction; prescribe melt rate
- Cause: no built-in prescribed-melt option.
- Fix: `SHELFICEheatTransCoeff=0., SHELFICEsaltTransCoeff=0.` (with SHELFICEuseGammaFrict=T set the gamma coefficients to 0) turns off heat/FW fluxes (cavity -> sill/obstacle); `SHELFICEconserve=.FALSE.` removes the FW-associated heat flux; `SHELFICEkappa=0.` (default SHELFICEadvDiffHeatFlux=F) removes conductive heat flux; for prescribed melt use `useISOMIPTD=.TRUE.` and set shelfIceFreshWaterFlux in the ISOMIP block of shelfice_thermodynamics.F. Holland-Jenkins 99 needs `#define SHI_ALLOW_GAMMAFRICT` + `SHELFICEuseGammaFrict=.TRUE.` (etaStar hard-set to 1).
- Era: 2014-2021.
- Src: mitgcm-support 2017-August 'on the shutdown of the thermodynamic effect of an ice shelf' ; 2021-November 'ISOMIP melt rate settings' ; 2019-May 'Remove Conductive Heat Flux in Shelfice package' ; 2014-June 'using H&J '99 param for shelfice'

### Ice-shelf cavity blows up at start / pressure at ice base (SHELFICEloadAnomalyFile, rhoConst/rhoNil)
- Cause: initial pressure load inconsistent with ambient density, large initial adjustment; pLoc for melt used R_shelfIce (inaccurate).
- Fix: compute load anomaly as in verification/isomip/input/gendata.m, or set rhoConst=rhoNil to the cavity-water density and let pkg/shelfice initialise load; short deltaT during spin-up; PR #456 (2021-04) uses shelficeMass*gravity for pLoc (changes all shelfice results). Dig out the seabed under floating ice at the grounding zone so every floating column has ocean beneath (Holland); keep minimum water-column thickness.
- Era: 2014-2021.
- Src: mitgcm-support 2014-October 'Calculation of the density in the SHELFICEloadAnomalyFile in data.shelfice' ; 2020-November 'Treatment of Antarctic grounding zone regions' ; https://github.com/MITgcm/MITgcm/pull/456

### SHI_withBL_realFWFlux / SHELFICEboundaryLayer with autodiff -> STOP
- Cause: a STOP (added 2015) when AD compiled with real-FW flux + BL; no technical reason found, adjoint ran fine after removing it (timothyas).
- Fix: remove the STOP in shelfice (needs PR); SHELFICEMassStepping step_icemass block should not be dropped under ALLOW_AUTODIFF. Remeshing itself is not adjointed.
- Era: 2020-2021.
- Src: https://github.com/MITgcm/MITgcm/issues/349

### pkg/frazil: not coupled to shelfice; freezing point differs from shelfice/seaice
- Cause: frazil hard-codes its own freezing-point coefficients; it only adds a heat sink + surface source (hack), no accretion at the shelf.
- Fix: do not expect accretion; prefer removing frazil (see supercooling entry); hard-code consistent coefficients (issue #131 open).
- Era: 2018-2025.
- Src: https://github.com/MITgcm/MITgcm/issues/131 ; mitgcm-support 2023-December 'Question about FRAZIL package' ; 2022-June 'Question on SHELFICE: ice accretion under ice shelf'

### Ice-shelf diagnostics confusion (SHIForcS vs SHIfwFlx, tendencies not where flux applied)
- Cause: SHIForcS/SHIfwFlx/SHIhtFlx give the forcing, but with SHELFICEboundaryLayer it is spread over cells kTopC and kTopC+1 (unequal denominators); MNC useMissingValue masked shelfice diags wrongly (2018).
- Fix: use SHIForcS (salt, g/m2/s) / SHIfwFlx (kg/m2/s) for budgets; integrate over ice area for melt rate; for where it lands save gT_Forc/gS_Forc (3-D, since 2014).
- Era: 2009-2018.
- Src: mitgcm-support 2017-August 'on the salinity flux at the ice shelf-ocean interface' ; 2010-May 'shelficeForcingS addition to boundary layer' ; 2018-May 'Problem with MNC masking of ice shelf diagnostics'

## RBCS (restoring) checklist

### RBCS appears to do nothing / "pkg/rbcs was not compiled (ALLOW_RBCS undef)" / namelist error (exit 19)
- Cause: usual list: package not in packages.conf (or not recompiled), useRBCS missing in data.pkg, useRBCtemp/useRBCsalt/useRBCuVel not set, mask not 3-D or all zero, relax file not 3-D, private external_forcing.F in code/ overriding the package, #define DISABLE_RBCS_MOM (RBCS_OPTIONS.h, default #undef) for velocities, number of masks > `maskLEN` in RBCS_SIZE.h (must change RBCS_SIZE.h, not old RBCS.h), retired rbcsIniter, or tauRelax too long.
- Fix: edit packages.conf then `make CLEAN`/genmake2 again; print-debug in rbcs_add_tendency.F; use diagnostics Um_Ext/Vm_Ext (RBCS part, k>1) and gT_Forc/gS_Forc (3-D); debugLevel>=3 writes the masks; the loops 0..sNx+1 in rbcs_add_tendency are intentionally not the full halo.
- Era: 2010-2025; rbcsIniter retired ~2011; RBCS.h split into RBCS_SIZE/PARAMS/FIELDS.
- Src: mitgcm-support 2015-August 'Hi' (ALLOW_RBCS) ; 2014-June 'Namelist Error with RBCS and multiple tracers' ; 2010-July 'Rotating tank in Cartesian' ; 2025-June 'Question about RBCS gridding in rbcs_add_tendency.F' ; 2014-January 'RBCS relax u velocities problem'

### RBCS with tauRelax < deltaT: noise or no restoring
- Cause: relaxation is explicit: deltaT/tauRelax<2 for stability (<1 with AB2), <1 (<1/2 with AB-2) for non-oscillatory; tauRelax=deltaT replaces the field completely.
- Fix: tauRelax >= deltaT (>= 2 deltaT with AB2); implicit relaxation not available. Applies to tauThetaClimRelax/climsstTauRelax too.
- Era: 2014-2018.
- Src: mitgcm-support 2018-March 'RBCS relaxation time scale limits' ; 2014-February 'Adding constant submerged and surface ice.'

### RBCS time axis: forcing offset, single-time files, cycles, error "startTime before rbcsForcingOffset+0.5*rbcsForcingPeriod"
- Cause: default rbcsForcingOffset=0 assumes records are averages labelled at end of the period; snapshot files need an offset of half the period; cyclic data in one file needs 3-D files; linear interpolation between records applies as for other forcing.
- Fix: for snapshots every P seconds set `rbcsForcingOffset = P/2` (e.g. 900 for 1800 s); `rbcsSingleTimeFiles=.TRUE.` uses file suffix = iteration number of records (independent of deltaT); `debugLevel>=debLevB` prints RBCS_FIELDS_LOAD record numbers; to stop forcing mid-run, restart with pickupStrictlyMatch=.FALSE. and no RBCS file. rbcsVanishingFac/rbcsVanishingTime ramps forcing down.
- Era: 2011-2022.
- Src: mitgcm-support 2011-June 'Hi' (RBCS offset) ; 2021-June 'RBCS forcing period spec' ; 2022-February 'RBCS scheme query' ; 2018-June (RBCS stop forcing thread)

### RBCS velocities blow up / asymmetry when fields are on the wrong grid
- Cause: relaxU/V files and masks must be on U/V (C-grid) points, T/S on centres; criteria on theta evaluated at T points break symmetry for U/V.
- Fix: put each field on its own grid; in a conditional use the average of the two adjacent theta; use separate U/V masks (not the T mask); relax baroclinic internal tides in an RBCS sponge rather than OBCS when lateral structure matters.
- Era: 2014-2023.
- Src: mitgcm-support 2017-July 'RBCS, OBCS and IC grid' ; 2023-October 'Asymmetric surface temperature under symmetric forcing when rbcs package is used' ; 2014-January 'RBCS relax u velocities problem'

### Model much slower with RBCS (3-D relax fields)
- Cause: reading 3-D forcing every record (I/O), not compute.
- Fix: useSingleCpuIO=.TRUE., readBinaryPrec=32 inputs, check TIMER sections (LOAD_FIELDS_DRIVER / RBCS_FIELDS_LOAD).
- Era: 2018.
- Src: mitgcm-support 2018-January 'Hi' (RBCS slowdown)

## OBCS, tides and sponges

### SSH does not follow the parent model / drifts after enabling tides (useOBCSbalance, useOBCStides)
- Cause: `useOBCSbalance=.TRUE.` removes all net transport except prescribed tides, killing low-frequency SSH from the parent; without it unbalanced inflow makes eta drift; there is no SSH nudging (OBCS eta is not used); do not double count tides (useOBCStides plus tides inside useOBCSprescribe files); balance is applied *before* tides are added.
- Fix: balance the boundary velocities offline over a few days (keep tides separate), or use OBCS_balanceFac* selectively (only one boundary adjusted); `useOBCSbalance=F` for low-frequency parent SSH; sea-level nudging would need RBCS/OBCS code. CPP flag is ALLOW_OBCS_BALANCE.
- Era: 2018-2025.
- Src: mitgcm-support 2025-July 'Question About Sea Surface Height (Eta) Behavior in MITgcm' ; 2019-December 'Volume conservation with open boundaries and tides' ; 2022-August 'model unstability caused by increasing sea surface elevation(eta)' ; 2021-November '[EXTERNAL] Warm start and cold start yielding same results'

### Ways to force tides; units; division by zero in obcs_add_tides
- Cause/Fix: regional: `useOBCStides=.TRUE.` + tidalPeriod(s), OB*amFile/OB*phFile (amplitude m/s normal velocity, phase in seconds, not hours), done in pkg/obcs/obcs_add_tides.F. Large domains also need tidal potential: exf `tidePot` (`#define EXF_ALLOW_TIDES`, 2-D geopotential files, tidePotFile/tidePotStartdate*) or Oliver Jahn's SPICE-based potential; LLC4320 tides were prescribed through apressure files. Tides in a periodic channel need body force or tidePot. PR #602: tidalPeriod=0 for unused components divided by zero (fixed 2022-02); "Non-existing record number" = too-short OB file or wrong precision.
- Era: 2014-2025.
- Src: mitgcm-support 2024-December 'How to include tides by the atmospheric pressure like the tidal forcing in LLC4320' ; 2021-December 'Tides in channel model' ; 2017-June 'Problem with tides in the OBCS' ; 2021-March 'A problem about tides data' ; https://github.com/MITgcm/MITgcm/pull/602

### Spurious boundary jets / noise at open boundaries (CD scheme, Orlanski, prescribed-only, sponge)
- Cause: OBCS is not implemented with the CD scheme (obcs_check STOP "OBCS not yet implemented in CD-Scheme"); prescribed OBCS reflect everything not exactly matching; Orlanski is good for a single mode; linear FS leaves eta=0 at boundary; U-only sponge toward barotropic flow gives artefacts; Orlanski + nonlinear free surface was not implemented in 2009-2012 threads (re-check obcs_check before relying on it).
- Fix: `useCDscheme=.FALSE.` (or tauCD=deltaTmom is a no-op); use Stevens BCs (useStevensNorth etc.) and/or a sponge: `useOBCSsponge=.TRUE.`, `spongeThickness` (grid cells), `Urelaxobcsinner`/`Urelaxobcsbound` (and V/T/S) -- sponge T,S toward a profile (not constant tRef) and confirm the file path; large telescoped reservoir beyond the control region; for internal tides use an RBCS sponge with the nudged wave. Make boundary far from interest; balance transport (see above).
- Era: 2006-2026 (issue #1031 open 2026-08).
- Src: mitgcm-support 2006-September 'OBC problem: spurious boundary jets with C-D coupling' ; 2022-March 'Issue with surface ocean currents' ; 2022-September 'Can not run a regional model for long time' ; https://github.com/MITgcm/MITgcm/issues/1031

### Model blows up / strange fresh river inflow through OBCS (negative salinity, NaN in cg2d)
- Cause: advection undershoot at the OB when OB salinity is lower than interior, with big vertical CFL near surface (2.5 m cells).
- Fix: shrink deltaT (3 s fixed a fjord case), check advcfl_*, use flux-limited scheme, ramp river up, consider RBCS/sponge, useRealFreshWaterFlux.
- Era: 2012.
- Src: mitgcm-support 2012-December 'Freshwater OBCS problem'

### Unexpected periodicity: same tide at E/N and W/S, wrap-around (MITgcm is periodic by default)
- Cause: tile edges wrap unless land blocks them; "same tide at E/N and W/S" means periodic.
- Fix: put a closed (depth 0) row/column at the edges, or `notUsingXPeriodicity/notUsingYPeriodicity` in eedata, or ALWAYS_PREVENT_X/Y_PERIODICITY in CPP_EEOPTIONS.h.
- Era: 2006-2023.
- Src: mitgcm-support 2023-May 'Problem with Obcs when simulating tide' ; 2011-March 'How to simulate the geostrophic current?' ; 2006-December 'nonlinear free surface'

## Drag, viscosity, Coriolis

### Bottom drag confusion: no_slip_bottom vs bottomDragLinear/Quadratic vs sideDragFactor; too-strong tidal channel currents
- Cause: no_slip_bottom sets the BC of the vertical viscosity operator (friction ~ viscAr / dz); bottomDragLinear/Quadratic add a separate stress (independent of viscAr); all three add together. Quadratic drag is a stress divided by bottom cell thickness (JM-C). sideDragFactor default 2 = no slip, 1 = half slip. no_slip_bottom=T with no drag still gives linear-like drag.
- Fix: for narrow/shallow channels combine no_slip_sides=T, bottomDragQuadratic ~1e-3-2.5e-3 (or `zRoughBot=0.01` log-law drag, PR #574), viscAhGridMax=0.5 (max 1), viscAhGrid/viscA4Grid ~0.01 instead of fixed viscAh, deepen/widen unresolved channels. Spatially varying drag: pkg/ctrl 2-D map, PR #574 (2-D quadratic coefficient), or Klymak's var_bot_drag. UBotDrag is inside Um_Diss only with selectImplicitDrag=0.
- Era: 2008-2024.
- Src: mitgcm-support 2019-June 'Bottom drag setup question' ; 2024-June 'bottomDragQuadratic magnitude' ; 2024-May 'Problem: Extreme Tidal Currents in the Simulation at Narrow Channels.' ; 2016-June 'Side Drag in a narrow canal' ; https://github.com/MITgcm/MITgcm/issues/423

### Grid-scale noise / waves near boundaries on a lat-lon beta-plane (tutorial_baroclinic_gyre)
- Cause: dx=cos(lat)*delta (spacing in degrees) but dy fixed -> non-square cells; unresolved gravity waves from boundaries; CD scheme masks problems.
- Fix: `selectCoriScheme=3` (Jamart & Ozer; noise confined to western boundary) with viscAh ~20e3 / viscA4Grid~0.01; turn off CD scheme; make cells near-square (delY from cos(lat) recurrence); use useSingleCpuIO (not globalFiles); synchronous time stepping did not help.
- Era: 2025-04.
- Src: mitgcm-support 2025-April 'Irregular Velocity Fields Near Basin Boundaries in MITgcm'

### Viscosity/biharmonic limits and cosPower on isotropic grids
- Cause: stability limits for viscA4 differ by cosPower/area-scaling options (older data point for 1/4 deg ECCO2 grid).
- Fix: with `#define COSINEMETH_III` and `ISOTROPIC_COS_SCALING` (CPP_OPTIONS.h or MOM_FLUXFORM_OPTIONS.h) and cosPower=2 viscA4 can go ~30x higher (cosPower 3-4 for larger). Use dimensionless viscAhGrid/viscA4Grid + viscAhGridMax/viscA4GridMax instead.
- Era: 2003-2019.
- Src: mitgcm-support 2003-December '[Fwd: RE: grid-scale noise in 1/4-deg model]' ; 2019-January 'non-hydrostatic pressure and KPP' (viscAhGrid)

### vectorInvariantMomentum + highOrderVorticity reverses flow over topography (ACC channel)
- Cause: 4th-order vorticity (highOrderVorticity) in mom_vecinv mishandles masks/partial cells at lateral boundaries near bumps; positive feedback creates large meridional currents (JM-C could reproduce). Flux-form has no analogue.
- Fix: do not set highOrderVorticity with topography (or use flux form); slightly more energy at high wavenumbers anyway. useAreaViscLength, KPP and viscA4GridMax masked the blow-up.
- Era: 2017-02 (checkpoint ~65x).
- Src: mitgcm-support 2017-February "Effect of 'vectorInvariantMomentum' flag on first order circulation" ; 2006-April 'Higher order advection scheme for momentum'

### Momentum budget will not close (TOTUTEND vs terms)
- Cause: missing terms or units. TOTUTEND is per day (divide by 86400). Close it as TOTUTEND/86400 = Um_Advec (includes Coriolis if no CD scheme) + Um_dPhiX (older name Um_dPHdx) + Um_Diss (explicit, incl. bottom drag if selectImplicitDrag=0) + Um_ImplD (implicit viscosity; older: VISrI_Um differences) + Um_Ext + AB_gU, with the surface pressure gradient already inside Um_dPhiX (older code needed -g*d(ETAN)/dx). With atmospheric/sea-ice loading add -grad(pLoad) with rhoConst (not rhoConstFresh). In NH, PHI_SURF + PHI_NH, not g*ETAN.
- Fix: 64-bit output; AB_gU/AB_gV diagnostics had a typo before timestep.F rev 1.55 (Dec 2013); VISrE/I are fluxes (multiply by area to compare; tendency = difference/volume); vorticity budget: take curl of closed momentum budget (stencil that makes pxy-pyx=0). Extra diagnostics `ALLOW_MOM_TEND_EXTRA_DIAGS` (PR #817) only if compiled.
- Era: 2010-2024.
- Src: mitgcm-support 2010-December 'Hi' (Shroyer) ; 2013-December 'MITgcm-support Digest, Vol 126, Issue 13' ; 2014-February 'Hi' (Piecuch) ; 2017-August 'non-hydrostatic momentum budget' ; https://github.com/MITgcm/MITgcm/issues/820

### Pressure diagnostics: PHIHYD, PHL, PH, PHRefC, z* differences
- Cause: totPhiHyd (PHIHYD) = phiHyd anomaly + Bo_surf*etaN + phi0surf (pressure/rhoConst, no reference part); with z* the grid-cell centre depth is not uniform so PHIHYD != z-run PHIHYD; PHL = bottom hydrostatic anomaly; PHRefC = g*z-type reference, not full.
- Fix: total p = rhoConst*(PHIHYD) + g*rhoConst*|z| (use rC; for z* compare against diagnostic `PHIHYDcR` at fixed depth rC); PNH separate and never part of PHIHYD; conversion uses rhoConst (rhoNil is the EOS density).
- Era: 2007-2022.
- Src: mitgcm-support 2022-October 'the implicit pressure components of PHIHYD?' ; 2017-May (PH/PHL/PHRefC thread, Klymak gist) ; 2018-November 'Pressure field in MITgcm experiment'

## Advection schemes and noise

### Instability/very different circulation with DST schemes (30, 33, 77) unless staggerTimeStep
- Cause: these schemes use forward time stepping (no Adams-Bashforth) so internal waves are unstable unless T,S are staggered in time; with synchronous steps stratification erodes via spurious w.
- Fix: `staggerTimeStep=.TRUE.` whenever tempAdvScheme/saltAdvScheme are 30/33/77 (and other schemes without AB-2). tempAdvScheme=7 (OS7MP) is the ECCO/Menemenlis default, far less diffusive than 33; scheme 33 can be too diffusive (erodes thermocline/halocline) at 18 km. Scheme 2 is dispersive; scheme 3/33 diffusive. staggerTimeStep does not change ptracer numerics directly.
- Era: 2003-2021.
- Src: mitgcm-support 2003-September 'problem higher order advections schemes' ; 2005-April 'staggerTimeStep and advection scheme' ; 2010-February (negative salinity/advection thread) ; 2004-March 'Temperature Advection schemes'

### Negative salinity / overshoots (T=49 degC, S=57 at bottom) with centered advection
- Cause: scheme 2 and 4 generate spurious extrema.
- Fix: use flux-limited (33, 77, 7) with staggerTimeStep=T; reduce explicit horizontal diffusion; check forcing (runoff) too. Scheme 81 (Prather SOM with limiter) preserves positivity only (negative T breaks it; Prather 1986); scheme 80/81 needs `GAD_ALLOW_TS_SOM_ADV`; restart from non-SOM pickup needs pickupStrictlyMatch=.FALSE..
- Era: 2010-2017.
- Src: mitgcm-support 2010-July 'Abnormal high and low temperature in model' ; 2017-November 'unexpected behavior of advection scheme 81' ; 2008-September 'Prather advection scheme' ; 2021-February 'Spurious mixing with internal tides'

### Grid-scale noise in non-hydrostatic wake/internal wave runs
- Cause: numerical dispersion of the central/DST scheme (33 in the 2026 case).
- Fix: switch to a more diffusive scheme (77) or 7, accept dissipation; add viscosity; reduce hFacMin to 0.1-0.01 removed slope noise in one 2-D tide run; Smagorinsky/Leith only act on horizontal viscosity.
- Era: 2019-2026.
- Src: mitgcm-support 2026-June 'numerical instability in a non-hydrostatic simulation' ; 2019-December 'help-noise with the internal wave simulation' ; 2019-July (Yangxin He slope noise thread in 2019-December)

### Non-hydrostatic + KPP warnings; 3-D Smagorinsky coefficients
- Cause: KPP is a mixed-layer parameterisation for km-scale grids; config_check warns "Implicit viscosity applies to provisional u,vVel" in NH; useSmag3D exists but with smag3D_coeff default 1e-2, smag3D_diffCoeff=0.
- Fix: for NH LES-like runs use `useSmag3D=.TRUE.` (compile flag ALLOW_SMAG_3D), no KPP; ivdc_kappa only if convection unresolved (ivdc 1000 too large); with GGL90 keep convective adjustment (cAdjFreq/ivdc_kappa).
- Era: 2016-2023.
- Src: mitgcm-support 2019-January 'non-hydrostatic pressure and KPP' ; 2016-August 'nonhydrostatic and kpp' ; 2023-December 'Vertical diffusion for non-hydrostatic ocean model with a stretched grid' ; 2013-November 'question about convection in nonhydrostatic and hydrostatic'

## Sea ice dynamics that interact with ocean physics (brief)

### Stripes/melt at tile (CPU) boundaries and immobile ice at 20 m grids with LSR
- Cause: LSR parallel solver accuracy (LSR_ERROR=2e-4 default is too loose), and local residuals at very high resolution.
- Fix: `LSR_ERROR=1e-5..1e-6`, `SEAICEnonLinIterMax=10` (default 2 barely converges), or use JFNK/EVP; `LSR_mixIniGuess>=2` fixed immobile floes in free drift; switch off SEAICE dynamics (SEAICEuseDYNAMICS=F) to test.
- Era: 2014-2021.
- Src: https://github.com/MITgcm/MITgcm/issues/171 ; https://github.com/MITgcm/MITgcm/issues/327 ; mitgcm-support 2014-July 'seance leakage at the cpu domain boundaries'

### EVP solver unstable in near ice-free cells (denomU/V ~ 0)
- Cause: division by near-zero denomU/V.
- Fix: PR #929: `SEAICE_evpAreaReg` (RL, run-time; regularises ice mass and area in denomU/V). Final name differs from proposals (SEAICEuseEVPreg / SEAICEevpRegDenomUV never merged). Also SEAICE_waterDrag silently changed units (kg/m3 factor 1000) around 2018 -> old data.seaice blows up (issue #107).
- Era: 2018, 2025-10.
- Src: https://github.com/MITgcm/MITgcm/pull/929 ; mitgcm-support 2018-November 'MOM_IMPLICIT_R: error when solving 3-Diag problem.'

## Build / run gotchas touching physics

### Segfault in adjoint or big runs: "Caught signal 11 (Segmentation fault: address not mapped...)"
- Cause: stack limit.
- Fix: `ulimit -s unlimited` (bash) / `limit stacksize unlimited` (csh) in the job script.
- Era: 2024-12 (Discover).
- Src: mitgcm-support 2024-December 'Running on Discover'

### Output "zero every 4th level" / reading bathy gives zeros / MON stats all zero except time and vorticity
- Cause: ifort missing `-assume byterecl` with -DWORDLENGTH=4 (zero planes); bathymetry must be negative (or zero); hFacMin must be <=1 (hFacMin=2.5 => empty domain, all-zero hFacC); OLx,OLy=1 too small for the scheme (NaNs, use 3 for 33/77).
- Fix: use standard optfile (linux_amd64_ifort+impi etc.), `genmake2 -mpi` (no manual -DALLOW_USE_MPI), check hFacC.data, bathy sign, SIZE.h overlaps.
- Era: 2007-2019.
- Src: mitgcm-support 2018-November 'No output every 4th depth' ; 2014-June 'Monitor statistics are all zeros except time and vorticity' ; 2007-February 'cg2d_init_res=NaN, cg2d_res=NaN'
