# Troubleshooting: instability, blow-ups, NaNs, CFL, time step
Distilled from answered mitgcm-support threads (2003-2026) and MITgcm GitHub issues/PRs. Entries merged by symptom; names checked against current origin/master (Oct 2026) unless Era says otherwise. "Community" = answered by non-core users (Klymak, Mazloff, Dustin, etc.), still consistent with core advice.
Thread URLs are month indexes: http://mailman.mitgcm.org/pipermail/mitgcm-support/YYYY-Month/thread.html

## How to diagnose a blow-up
Checklist distilled from repeated advice (Losch, Menemenlis, Campin, Doddridge, Klymak):
1. Rerun with `monitorFreq` = `deltaT` (or smaller; data, PARM03) so %MON stats are written every step; a larger value hides the step where it died. Production: 20-50*deltaT. Also `debugLevel` default or >=1 so `cg2d: Sum(rhs),rhsMax` is printed every step (`grep -c 'cg2d:' STDOUT.0000` = steps done). `debugLevel=-1` suppresses warnings (e.g. GAD_CHECK).
2. `grep advcfl STDOUT.0000`: `advcfl_uvel_max`, `advcfl_vvel_max`, `advcfl_wvel_max`, `advcfl_W_hf_max` (vertical CFL including hFac; the one that goes first with thin partial cells). Keep < 0.5 (many users aim < 0.1-0.2 for W). Find which one rises exponentially a few steps before the first NaN.
3. Read the END of both STDOUT.* and STDERR.*: the message is often split between them; STDOUT can be buffered, so the last line is misleading. Some ranks flush more than others; check all STDOUT.xxxx tails. `debugMode=.TRUE.` (eedata) flushes but is very slow (use `debugLevel=2` with it).
4. Locate where/when it starts: write snapshots every step (`dumpFreq` < `deltaT`, or pkg/diagnostics with `frequency=-deltaT`) and look at the first cell with extreme value (cube corners, narrow bays/inlets, steep topography, sponge edges, open boundaries, ice edge, thin hFac cells). Failing cell indices are printed by CALC_R_STAR and by exf range checks.
5. Check dynstat_eta_min/max, theta/salt min/max in %MON for the first drift (eta_mean drifting = unbalanced OBCS or unbalanced fresh water; theta min < freezing/-2 = flux limiter masking a problem).
6. Simplify: switch off packages in data.pkg (OBCS, exf, seaice, shelfice, KPP/GM) and re-add one at a time; run with initial T/S only and no forcing; run without flux limiters (e.g. tempAdvScheme=2/4) to see the underlying problem; try smaller deltaT to separate CFL from setup error.
7. Sanity check inputs before digging into code: NaN in T/S/forcing/bathymetry files (land filled with NaN or -9999), bathymetry sign (negative depth), correct record counts/precision (`readBinaryPrec`), forcing units/sign, wind/heat flux range (exf `useExfCheckRange`).
8. If it blows up only after long runs, restarting from the last pickup often does not reproduce it (round-off chaos); that does not mean it is fine (see restart.md).
9. Trap NaNs early: build with `genmake2 -ieee -devel` + FP-trap flags (`-ffpe-trap=invalid,zero,overflow -fbacktrace` for gfortran) and keep the .f files to read the traceback line (see "NaNs continue" below).
10. Compare with the nearest verification experiment (testreport) built with same optfile: if that blows up on your machine it is a compiler/MPI/optimization problem, not the setup.

## CFL, time step, viscosity/diffusion limits
### SOLUTION IS HEADING OUT OF BOUNDS: tMin,tMax= ...  (ABNORMAL END: S/R MON_SOLUTION, stops due to EXTREME Pot.Temp)
- Cause: pkg/monitor mon_solution.F stops the run when theta range exceeds `monSolutionMaxRange` (1.e3): symptom of a model that already exploded (CFL violation, strong initial adjustment, bad input), not a parallel-environment error. The "S/R EEDIE: Only 0 threads have completed" line printed next to it is a harmless artifact of the abrupt stop.
- Fix: monitorFreq=deltaT, inspect advcfl_*; reduce deltaT; check viscosity/diffusion CFL; check where it blows up (see checklist). Reduce deltaT first for a few days after cold start, then raise.
- Era: 2004-2021 (Martin Losch repeatedly; same text in checkpoint6x-69). Message and parameter still exist.
- Src: mitgcm-support 2006-August 'SOLUTION IS HEADING OUT OF BOUNDS: tMin,tMax= ?' http://mailman.mitgcm.org/pipermail/mitgcm-support/2006-August/thread.html ; 2016-December 'Unstability in a Model (SOLUTION IS HEADING OUT OF BOUNDS)' http://mailman.mitgcm.org/pipermail/mitgcm-support/2016-December/thread.html ; 2015-March 'how many timesteps?' (Doddridge: message + "STOPPING CALCULATION at Iter=" in STDERR) http://mailman.mitgcm.org/pipermail/mitgcm-support/2015-March/thread.html

### Model fine for N steps then NaN / "tMin,tMax" large, after changing resolution or grid
- Cause: parameters (deltaT, viscAh, diffKh) of a coarse tutorial reused at higher resolution; advective CFL u*dt/dx < ~0.5 (use the smallest cell: near poles dx = R*dphi*cos(lat)), plus viscous/diffusive limits S_l = 4*Ah*dt/dx^2 < ~0.3 and 4*Kh*dt/dx^2 < 1 (doc: tutorial barotropic/global_ocean stability section; Mazloff: guidelines, tuned by trial). Diffusion example: Kh=4e4 at 1/20 deg dt=360 s gave S=1.86 -> NaN; Kh=1e4 marginal.
- Fix: reduce deltaT (by ~the refinement factor), or use scaled viscosity `viscAhGrid` (< 1; start ~0.1), `viscA4Grid` (~0.1), or Leith/Smagorinsky; set viscAh=0 when using Grid-scaled options. Vertical: dt < dz^2/(8*viscAr) if viscosity is explicit (implicitViscosity=.TRUE. removes the limit).
- Era: all eras. 
- Src: mitgcm-support 2012-June 'Problem with Global Ocean 1deg run' (parameters of coarse model) http://mailman.mitgcm.org/pipermail/mitgcm-support/2012-June/thread.html ; 2013-June 'High diffusion offline runs blow up' http://mailman.mitgcm.org/pipermail/mitgcm-support/2013-June/thread.html ; 2013-July 'model crashing & time step issues' http://mailman.mitgcm.org/pipermail/mitgcm-support/2013-July/thread.html ; 2014-July 'Model based on exp4 experiment - stability problems' (Campin: dt < dz^2/(8 viscAr)) http://mailman.mitgcm.org/pipermail/mitgcm-support/2014-July/thread.html ; 2016-March 'CFL condition for a Spherical Case' http://mailman.mitgcm.org/pipermail/mitgcm-support/2016-March/thread.html

### viscAhGrid / viscAh details; numerical noise although advective CFL << 1 (e.g. 0.02)
- Cause: advective CFL is only one criterion. mom_calc_visc.F: harmonic visc = viscAh + 0.25*L^2*viscAhGrid/deltaT, so viscAh and viscAhGrid ADD (a 2 deg viscAh=1e4 plus viscAhGrid at 1/6 deg overshoots). viscAhGrid effective viscosity grows as deltaT shrinks, so changing deltaT is not a clean stability test. Too small viscAhGrid gives grid-scale noise or blow-up (Menemenlis).
- Fix: viscAhGrid ~0.1 (<1 guarantees the viscous CFL), optional viscA4Grid=0.1 to remove C-grid noise; smaller `ivdc_kappa` (1-10) if noise is vertical; Leith (`useFullLeith=.TRUE.`, viscC4Leith/viscC4Leithd ~1-2, `viscAhGridMax`); one user (Orlanski-open domain, 2019) switched viscAhGrid -> `viscC2smag` successfully.
- Era: 2005-2019. Parameters exist today.
- Src: mitgcm-support 2018-December 'numerical noise with CFL<0.02' http://mailman.mitgcm.org/pipermail/mitgcm-support/2018-December/thread.html ; 2007-May 'ploblem with a variable-grid global ocean circulation model' (viscAh + viscAhGrid sum, formula quoted) http://mailman.mitgcm.org/pipermail/mitgcm-support/2007-May/thread.html ; 2014-March 'MITgcm-support Digest, Vol 128, Issue 38' (grid noise, blow-up) http://mailman.mitgcm.org/pipermail/mitgcm-support/2014-March/thread.html ; 2019-July 'SSH blow-up at Orlanski BC with small deltaT' http://mailman.mitgcm.org/pipermail/mitgcm-support/2019-July/thread.html

### Blow-up with thin partial cells (hFacMin small, hFacMinDr small): vertical CFL
- Cause: thinnest cell = MAX(hFacMin, MIN(hFacMinDr*recip_drF(k), 1))*drF; w < dz_min/dt, so a tiny `hFacMin` (e.g. 0.001) or `hFacMinDr` makes the vertical CFL impossible; noise often starts in bottom partial cells near steep topography. `advcfl_W_hf_max` includes hFac, `advcfl_wvel_max` does not (identical values mean no thin cells active).
- Fix: raise hFacMin (>=0.1-0.2) / hFacMinDr, reduce deltaT, add no-slip or quadratic drag, check bathymetry deeper than sum(delR) (reported blow-up with Prather scheme); test linear free surface when debugging.
- Era: 2011-2018 (Klymak, Losch). No upstream change.
- Src: mitgcm-support 2018-December 'Tracer instabilities when I reduce hFacMinDr' (Losch formula, unresolved by user) http://mailman.mitgcm.org/pipermail/mitgcm-support/2018-December/thread.html ; 2017-June 'nTimeStep' (Klymak, hFacMin arithmetic) http://mailman.mitgcm.org/pipermail/mitgcm-support/2017-June/thread.html ; 2012-May 'Strips on the surface' (hFacMin=0.001, Losch: use >=0.1) http://mailman.mitgcm.org/pipermail/mitgcm-support/2012-May/thread.html

### Blows up in first hours/days after starting from climatology / interpolated / perturbed initial conditions
- Cause: initial shock: dense-over-light water (vertical instabilities, not horizontal smoothness), strong horizontal density gradient with no balanced flow (esp. next to open boundaries), perturbed forcing, or velocity fields interpolated inconsistently from another grid (advcfl_wvel already ~0.7 at iteration 0).
- Fix: short deltaT for the first day(s) with monitorFreq=deltaT until advcfl_w falls below 0.5, then increase gradually (Menemenlis llc1080 spin-up: 30 s x1 d, 90 s x9 d, 120 s x10 d, 180 s x25 d, then 240 s); keep adjusting deltaT in restarts (set via nIter0/pickup). Remove static instability in the IC (convective-adjust IC offline); `ivdc_kappa`=1-10 helps first steps. When regridding a state (1/12 -> 1/6 deg) bin-average rather than fill land points with defaults; never fill land/dry cells with a constant (Mazloff, Shevchenko 2020 fix); use nearest-neighbour or density-aware interpolation to avoid unstable stratification.
- Era: 2013-2024.
- Src: mitgcm-support 2014-February 'is it feasible to use tRef, sRef as initial conditions' (Menemenlis, Spall) http://mailman.mitgcm.org/pipermail/mitgcm-support/2014-February/thread.html ; 2014-May 'run speed unaffected by timestep' (llc1080 schedule) http://mailman.mitgcm.org/pipermail/mitgcm-support/2014-May/thread.html ; 2014-July 'model crashing, cfl parameter' http://mailman.mitgcm.org/pipermail/mitgcm-support/2014-July/thread.html ; 2020-March 'Solution blows up' http://mailman.mitgcm.org/pipermail/mitgcm-support/2020-March/thread.html ; 2025-April 'the floating-point exception issue' (Menemenlis: CFL < 0.5, shorten deltat after grid change) http://mailman.mitgcm.org/pipermail/mitgcm-support/2025-April/thread.html

### NaNs/blow-up caused by tracer reference profile ordering or EOS parameters (tRef reversed, tAlpha large, extreme salinity)
- Cause: `tRef`/`sRef` are listed from the surface downward; a warm-at-depth profile is statically unstable and blows up in a few steps. Larger `tAlpha` makes N^2 larger and can reduce stable dt. Very high salinity (hypersaline lakes, 30-80+) is outside the EOS fit range: large density gradients -> NaN at step 2.
- Fix: check tRef ordering; shorter deltaT for larger tAlpha; for exotic salinity try eosType='LINEAR' first, smooth the initial fields; Dead-Sea type problems need own EOS (find_rho.F). Energy-consistent Boussinesq pressure choice: `selectP_inEOS_Zc` 0..3 (defaults: JMD95Z->0, MDJWF->2; reported in config summary).
- Era: 2021-2024.
- Src: mitgcm-support 2021-April 'internal wave tutorial - tRef' http://mailman.mitgcm.org/pipermail/mitgcm-support/2021-April/thread.html ; 2021-December 'How to select the thermal expansion coefficient' http://mailman.mitgcm.org/pipermail/mitgcm-support/2021-December/thread.html ; 2023-February 'Initial conditions ... higher initial values for temperature and salinity, model yields nan' http://mailman.mitgcm.org/pipermail/mitgcm-support/2023-February/thread.html ; 2024-March 'Nonlinear equations of state' http://mailman.mitgcm.org/pipermail/mitgcm-support/2024-March/thread.html

### STOP ABNORMAL END: S/R INI_THETA / INI_SALT "found N wet grid-pts with theta=0 identically"
- Cause: `hydrogThetaFile` (PARM05) OVERWRITES tRef; a missed point in your interpolated field is exactly 0 (or NaN). Same check exists for salt (ini_salt.F).
- Fix: fill all wet points; for genuine zeros set `checkIniTemp=.FALSE.` / `checkIniSalt=.FALSE.` (PARM05). To use tRef/sRef alone, comment out hydrogThetaFile/hydrogSaltFile. Fresh (salt=0) ocean points in a salty IC also seed blow-ups (INI_SALT warning).
- Era: 2007-2017; checks and flags still present (ini_theta.F, ini_salt.F).
- Src: mitgcm-support 2017-January 'ABNORMAL END: S/R INI_THETA (when changing geometry)' http://mailman.mitgcm.org/pipermail/mitgcm-support/2017-January/thread.html ; 2008-June 'abnormal end' http://mailman.mitgcm.org/pipermail/mitgcm-support/2008-June/thread.html ; 2011-September 'Abnormal end' (salt=0 points) http://mailman.mitgcm.org/pipermail/mitgcm-support/2011-September/thread.html

### Immediate NaN / "Initial hydrostatic pressure did not converge after 15 steps" / "Iteration 1, RMS-difference = NaN" / seaice-dynamics NaN after one step
- Cause: NaN (or -999/-9999 fill) in input T/S, wind or bathymetry files (including land-masked atmospheric forcing: exf averages winds to u/v points so masking the forcing leads to NaN near coasts). With seaice, NaNs enter at the start of the step via the ice solver and are blamed on LSR. Another cause: OLx/OLy too small.
- Fix: replace NaN by 0 (or valid values) in every binary input; define atmospheric forcing over the entire domain; test `debugLevel=5`, `dumpFreq=1`; check OLx,OLy >= 2-3 in SIZE.h. Increasing the loop count in ini_pressure.F does not help.
- Era: 2008-2020.
- Src: mitgcm-support 2020-February 'STOP ABNORMAL END: S/R INI_PRESSURE' http://mailman.mitgcm.org/pipermail/mitgcm-support/2020-February/thread.html ; 2019-August 'Sea ice dynamics causing model to output NaN's after changing grid configuration' (Losch: NaN in wind forcing; debugLevel=5) http://mailman.mitgcm.org/pipermail/mitgcm-support/2019-August/thread.html ; 2018-August '(no subject)' (Klymak: fill file with valid values) http://mailman.mitgcm.org/pipermail/mitgcm-support/2018-August/thread.html ; 2007-February 'cg2d_init_res=NaN, cg2d_res=NaN' (cause OLx=OLy=1; fixed with 3) http://mailman.mitgcm.org/pipermail/mitgcm-support/2007-February/thread.html

### Uniform/zero-forcing run still unstable: central advection without diffusion (default tempAdvScheme=2)
- Cause: default centered 2nd-order tracer advection needs explicit diffusion; with diffKh=0 it is unstable (lakes, small domains, 1-m vertical grid). Vertical CFL then limits dt (dt ~ 1 s at 1-m dz).
- Fix: use a scheme with numerical diffusion (`tempAdvScheme` 3, 30, 33, 7; see GAD.h), or set `diffKhT`>0 and nonzero `diffKrT`; `saltStepping=.FALSE.` for freshwater lakes; explicit vertical diffusion + vertical walls reported unstable (removing diffKzT cured a 2-D case).
- Era: 2011-2024.
- Src: mitgcm-support 2024-August 'Seeking Assistance with Numerical Instability in MITgcm Simulation of a Plateau Lake' http://mailman.mitgcm.org/pipermail/mitgcm-support/2024-August/thread.html ; 2011-August 'Model crashing due to vertical velocities' http://mailman.mitgcm.org/pipermail/mitgcm-support/2011-August/thread.html

## Advection scheme / time-stepping scheme
### ** WARNING ** GAD_CHECK: potentially unstable time-stepping (Internal Wave) / need "staggerTimeStep=.TRUE." in "data", nml PARM01
- Cause: DST/flux-limited multi-dimensional schemes (30, 33, 77, 7, ...) do not use Adams-Bashforth for theta/salt, which leaves the internal-wave mode explicit; unless staggerTimeStep is on, the stratification can overturn during initial adjustment and blow up (Campin/Adcroft since 2003).
- Fix: `staggerTimeStep=.TRUE.` (PARM01) with any of these schemes. Warning is written to STDERR; `debugLevel=-1` can hide similar messages. Still checked in pkg/generic_advdiff/gad_check.F.
- Era: 2003-2018 (still current).
- Src: mitgcm-support 2015-March 'stagger time step' http://mailman.mitgcm.org/pipermail/mitgcm-support/2015-March/thread.html ; 2015-October 'Potentially unstable time-stepping' http://mailman.mitgcm.org/pipermail/mitgcm-support/2015-October/thread.html ; 2005-April 'staggerTimeStep and advection scheme' (explicit rule) http://mailman.mitgcm.org/pipermail/mitgcm-support/2005-April/thread.html ; 2003-September 'problem higher order advections schemes' http://mailman.mitgcm.org/pipermail/mitgcm-support/2003-September/thread.html

### Instability/odd 2-D or poles runs fixed by Adams-Bashforth parameter
- Cause: abEps = 0.01 too weak to damp computational mode (pole rebound in 2-D lat-depth, internal waves).
- Fix: `abEps=0.1` (PARM03); typical production value in ECCO-type configs.
- Era: 2023 (Flynn Ames, Enceladus 2-D).
- Src: mitgcm-support 2023-March 'Anomalous rebound circulations at poles in 2D setup leading to numerical instability' http://mailman.mitgcm.org/pipermail/mitgcm-support/2023-March/thread.html

### Blow-up / NaN with higher-order scheme and tiny overlap (OLx, OLy too small; sNy=1 2-D)
- Cause: advection schemes need overlap width: 33/30/77 need OLx,OLy >= 2-3 (cube grids 2x), scheme 7 needs 4, PTRACERS_advScheme=33 needs >=2 (+1 with GMadvectiveForm+Visbeck); OLx=OLy=0/1 blows up immediately. In 2-D (Ny=1) slabs keep OLy>=2 (>=4 for scheme 7); channel with Ny=1 needs a wall in x or it is doubly periodic.
- Fix: set OLx, OLy in SIZE.h accordingly and rebuild (`make CLEAN` after header edit). Not always trapped at run time.
- Era: 2005-2010.
- Src: mitgcm-support 2008-September 'problem running tutorial_global_oce_latlon' (OLx 2->4 fixed NaN) http://mailman.mitgcm.org/pipermail/mitgcm-support/2008-September/thread.html ; 2010-July '2D setup' (Campin) http://mailman.mitgcm.org/pipermail/mitgcm-support/2010-July/thread.html ; 2008-May 'Re: MITgcm files' (overlaps 0 blow up) http://mailman.mitgcm.org/pipermail/mitgcm-support/2008-May/thread.html

### Quasi-hydrostatic + staggerTimeStep fatal instability (large QH terms, weak stratification)
- Cause: QH contributions (Coriolis-prime, NH metric terms) were advanced with unstable forward Euler under staggerTimeStep.
- Fix: fixed upstream in PR #433 (2021-02), new CPP `ALLOW_QHYD_STAGGER_TS` in CPP_OPTIONS.h; adds a field to pickup (restart from old pickup possible with pickupStrictlyMatch=.FALSE., not perfect).
- Era: before checkpoint68-ish (Feb 2021).
- Src: https://github.com/MITgcm/MITgcm/pull/433

### Instability after updating code past checkpoint64e with SOM advection (schemes 80/81)
- Cause: in gad_som_advect.F `noFlowAcrossSurf = rigidLid .OR. nonlinFreeSurf.GE.1 .OR. select_rStar.NE.0` evaluates TRUE even with linear free surface (user's finding).
- Fix: user workaround: set noFlowAcrossSurf=.FALSE. in a local copy; expression is still present in current gad_som_advect.F / gad_som_adv_r.F.
- Era: 2018 (checkpoint64e+).
- Src: mitgcm-support 2018-August 'SOM instability after checkpoint 64e' http://mailman.mitgcm.org/pipermail/mitgcm-support/2018-August/thread.html

### Hydrostatic grid-scale W noise ("columns of alternating w") over slopes/mixed layer
- Cause: hydrostatic convection unresolved; GGL90 or no convection scheme with `ivdc_kappa` and `cAdjFreq` default 0.
- Fix: with GGL90 still enable convective adjustment `ivdc_kappa` (>0; Heimbach), and optionally `#define ALLOW_GGL90_SMOOTH` (GGL90_OPTIONS.h); KPP does not need it. Accept some noise at coarse resolution.
- Era: 2014 (names still exist).
- Src: mitgcm-support 2014-May 'Help me get rid of my noise!' http://mailman.mitgcm.org/pipermail/mitgcm-support/2014-May/thread.html

### Surface-wave / implicit free surface damping or rigid lid blow-ups
- Cause: implicit free surface damps external waves; rigidLid=.TRUE. case blew up as dx grew (eta_mean drifted before CFL rose); `implicitFreeSurface=.TRUE.` was stable. NH + non-fully-implicit barotropic solver is flagged "NOT SAFE" by config_check.
- Fix: for accurate surface waves reduce deltaT or Crank-Nicolson `implicSurfPress=0.5`, `implicDiv2DFlow=0.5` (with hydrostatic; NH config_check message exists); otherwise prefer implicit free surface to rigid lid.
- Era: 2007-2017 (config_check NOT SAFE message still present).
- Src: mitgcm-support 2017-April 'surface treatment and resolution dependent instability?' http://mailman.mitgcm.org/pipermail/mitgcm-support/2017-April/thread.html ; 2009-June 'Surface gravity waves?' http://mailman.mitgcm.org/pipermail/mitgcm-support/2009-June/thread.html ; 2007-April 'wave damping' http://mailman.mitgcm.org/pipermail/mitgcm-support/2007-April/thread.html

### NH: NaN after minutes with sharp density step / cg3d not converging
- Cause: cg3dMaxIters too small (20 with Nr 40+), `cg2dUseMinResSol` default 0; initial profile with delta-function pycnocline and extreme vertical resolution; hydrostatic OBCS inflow unbalanced.
- Fix: `cg3dMaxIters` >= few*Nr (Campin), try `cg2dUseMinResSol=1`; smooth initial density; re-evaluate setup (Campin 2014 thread).
- Era: 2010-2014.
- Src: mitgcm-support 2010-June 'A very strange problem' http://mailman.mitgcm.org/pipermail/mitgcm-support/2010-June/thread.html ; 2014-April 'Weird numerical noise leads to the model blowing up' http://mailman.mitgcm.org/pipermail/mitgcm-support/2014-April/thread.html

## Time-step parameters (deltaT, deltaTmom, deltaTtracer, deltaTClock)
### Setting deltaTmom/deltaTtracer/deltaTfreesurf/deltaTClock separately: no speed-up, "wrong" results, ice not forming, Kelvin waves too slow
- Cause: they do NOT mean sub-stepping. Number of steps = (endTime-startTime)/deltaTClock (deltaTClock defaults to deltaTtracer if only deltaTmom/Tracer set; deltaTfreesurf defaults to deltaTmom unless set), and each step advances momentum by deltaTmom only. This is Bryan (1984) tracer acceleration: distorts physics (slowed waves, damped seasonal cycle); for spin-up only. sea-ice uses dTtracer as default step; very small deltaTtracer -> nothing develops.
- Fix: for normal runs set only `deltaT` (copied to all). If accelerating: deltaTmom=1200, deltaTtracer=deltaTfreesurf=deltaTClock=larger (global_ocean.cs32x15 example), check STDOUT config summary for actual values; finish with synchronous steps (Losch: asynchronous until drift negligible, then ~30 y synchronous).
- Era: 2003-2024.
- Src: mitgcm-support 2024-September 'Inquiry about deltaTtracer, deltaTmom, and deltaTClock in MITgcm' http://mailman.mitgcm.org/pipermail/mitgcm-support/2024-September/thread.html ; 2016-October 'Varying deltaTtracer in a barotropic model' http://mailman.mitgcm.org/pipermail/mitgcm-support/2016-October/thread.html ; 2005-November 'tracer's time step' http://mailman.mitgcm.org/pipermail/mitgcm-support/2005-November/thread.html ; 2018-December 'Time stepping for momentum and tracers' http://mailman.mitgcm.org/pipermail/mitgcm-support/2018-December/thread.html ; 2018-June 'Smaller deltaTmom speeds up the model?!' (cg2dTargetResidual 1e-13 only for adjoint; 1e-7/-8 fine) http://mailman.mitgcm.org/pipermail/mitgcm-support/2018-June/thread.html ; 2008-November 'The global flow' http://mailman.mitgcm.org/pipermail/mitgcm-support/2008-November/thread.html

### Offline tracer run unstable: deltaTtracer=86400 with daily flow
- Cause: offline `deltaToffline` only selects velocity file names/times; `deltaTtracer` is the advection step and must satisfy CFL.
- Fix: deltaTtracer = small enough (CFL << 1), keep offlineForcingPeriod/Cycle consistent; check files by iteration numbers.
- Era: 2019.
- Src: mitgcm-support 2019-June 'nondivergent flow advecting a tracer' http://mailman.mitgcm.org/pipermail/mitgcm-support/2019-June/thread.html

### cal_Set: The time step is less than a second. / contains fractions of a second.
- Cause: pkg/cal (hence exf) requires integer-second timestep.
- Fix: use integer seconds, or use pkg/exf without pkg/cal (verification offline_exf_seaice) which allows deltaT < 1 s; or pkg/bulk_force. Message still in cal_printerror.F.
- Era: 2009-2019.
- Src: mitgcm-support 2019-October 'Sub-1s Timestep with EXF' http://mailman.mitgcm.org/pipermail/mitgcm-support/2019-October/thread.html ; 2012-July 'using time step less than a second ??' http://mailman.mitgcm.org/pipermail/mitgcm-support/2012-July/thread.html

### Diagnostics/dump output never written at the requested time
- Cause: dumpFreq / diagnostics frequency not a multiple of deltaT (e.g. 400 s with dt=60); last-iteration snapshot diagnostics are not written (dumpAtLast only for averages).
- Fix: choose multiples; run one extra iteration to get the final snapshot.
- Era: 2007-2009.
- Src: mitgcm-support 2009-February 'mnc output problem-multiple time step' http://mailman.mitgcm.org/pipermail/mitgcm-support/2009-February/thread.html ; 2008-August 'No diagnostics output on last timestep' http://mailman.mitgcm.org/pipermail/mitgcm-support/2008-August/thread.html

## Surface forcing (EXF), heat/salt/fresh-water
### ABNORMAL END: S/R EXF_CHECK_RANGE
- Cause: an exf input field outside plausible range (e.g. wind > 100 m/s). Its message goes to STDOUT (not STDERR).
- Fix: inspect the field (units/sign/scale via exf_inscal_*, missing-value fill); only as last resort `useExfCheckRange=.FALSE.` (data.exf) which also removes protection against absurd heat fluxes. Typical sane ranges: hflux -250..600 W/m2, swflux -350..0 (exf_fields.h comments).
- Era: 2010-2026.
- Src: mitgcm-support 2022-March 'ABNORMAL END: S/R EXF_CHECK_RANGE' http://mailman.mitgcm.org/pipermail/mitgcm-support/2022-March/thread.html ; 2026-March 'LabSea model collapse' http://mailman.mitgcm.org/pipermail/mitgcm-support/2026-March/thread.html ; 2019-August 'Sea ice dynamics causing model to output NaN's ...' (useExfCheckRange=.FALSE. masks ridiculous fluxes)

### SST 30C+ ("heat trapping") then SURF_ADJUSTMENT / cg2d divergence in summer
- Cause: unrealistic Qnet/Qsw (sign error, units, shortwave in wrong channel, missing downward radiation or humidity files required by EXF_OPTIONS.h, global-mean flux -4800 W/m2, global qnet min/max -18/+32 kW/m2).
- Fix: check hflux/swflux sign convention (positive = cooling) and magnitude; keep useExfCheckRange=.TRUE.; fix forcing; check forcing_qnet_mean/min/max in %MON.
- Era: 2012-2026.
- Src: mitgcm-support 2026-March 'LabSea model collapse' (Menemenlis) http://mailman.mitgcm.org/pipermail/mitgcm-support/2026-March/thread.html ; 2012-February 'NaN output' (Campin: qnet -4800 W/m2) http://mailman.mitgcm.org/pipermail/mitgcm-support/2012-February/thread.html

### NaN with cheapaml / negative specific humidity / negative relative humidity
- Cause: slightly negative humidity (reanalysis overshoot, CORE-corrected RH) under SQRT in the long-wave bulk code; cheapAML also needs cheapamlYperiodic=.FALSE. on non-periodic lat-lon grids.
- Fix: clip humidity >= 0 in preprocessing; in current cheapaml.F the term uses SQRT(ABS(qair)) (fixed upstream after 2016); set cheapaml_ntim ~ O(10-50), use `FluxFormula='COARE3'` or provide downward long-wave.
- Era: 2010, 2016 (cheapaml fix upstream).
- Src: mitgcm-support 2016-November "cheapAML+global_ocean.90x40x15: it's all NaNs" http://mailman.mitgcm.org/pipermail/mitgcm-support/2016-November/thread.html ; 2010-March 'Relative Humidity' http://mailman.mitgcm.org/pipermail/mitgcm-support/2010-March/thread.html

### exf_interp_uv produces |vector| > 1e6 near the pole
- Cause: bicubic vector interpolation (interpMethod 12/22) with last data latitude extremely close to 90N (e.g. 89.999): denominator in LAGRAN -> 0; scalars use bilinear (1) so are fine.
- Fix: keep last input latitude well away from 90N; or use bilinear for vector fields; (suggestion to add a warning was not confirmed implemented).
- Era: 2017.
- Src: mitgcm-support 2017-May 'numerical issues with exf_interp_uv?' http://mailman.mitgcm.org/pipermail/mitgcm-support/2017-May/thread.html

### Crash with forcing/BC file read: "Non-existing record number" / read past end (looks like a blow-up)
- Cause: forcing file too short: model linearly interpolates between 2 consecutive records, so N days of daily forcing starting at 00:00 need N+1 records; same for OBCS files (MDS_READ_SEC_YZ ABNORMAL END). Changing externForcingPeriod without changing Cycle (or files) also triggers it.
- Fix: add records/shift startdate; use exf with per-field period; do not set periodicExternalForcing/externForcing* when using pkg/exf (they apply to all non-exf fields incl. OBCS).
- Era: 2012-2020.
- Src: mitgcm-support 2020-December 'Getting Backtrace error' http://mailman.mitgcm.org/pipermail/mitgcm-support/2020-December/thread.html ; 2020-March 'jobs died suddenly' http://mailman.mitgcm.org/pipermail/mitgcm-support/2020-March/thread.html ; 2016-September 'About pickup' (Voelker: Period/Cycle) http://mailman.mitgcm.org/pipermail/mitgcm-support/2016-September/thread.html ; 2012-November 'add the time dependent wind by zonalWindFile' http://mailman.mitgcm.org/pipermail/mitgcm-support/2012-November/thread.html

## Open boundaries (OBCS)
### Eta/sea level drifts then NaN: unbalanced OBC inflow/outflow (SSH growth, "unrealistic sea surface height growth")
- Cause: net volume through open boundaries not zero (Boussinesq: instantaneous action at a distance, barotropic waves); non-hydrostatic runs cannot use nonlinear free surface to absorb it.
- Fix: balance OBC transports offline (account for hFac geometry, store as 64-bit: `exf_iprec_obcs=64` or readBinaryPrec=64); or `useOBCSbalance=.TRUE.` (OBCS_balanceFacE etc., crude hack that removes tide signal). Check eta_mean in %MON.
- Era: 2016-2019.
- Src: mitgcm-support 2016-March 'instability at Eta' http://mailman.mitgcm.org/pipermail/mitgcm-support/2016-March/thread.html ; 2017-December 'Water balance problem' http://mailman.mitgcm.org/pipermail/mitgcm-support/2017-December/thread.html ; 2019-October 'Unrealistic sea surface height growth' http://mailman.mitgcm.org/pipermail/mitgcm-support/2019-October/thread.html ; 2018-June 'Eta calculation on different grid sizes' http://mailman.mitgcm.org/pipermail/mitgcm-support/2018-June/thread.html

### Spurious boundary jets / blow-up with OBCS and useCDscheme=.TRUE.
- Cause: CD-scheme + OBCS not implemented (D-grid velocities along boundaries).
- Fix: disable CD scheme; since Sept 2006 obcs_check.F stops with "OBCS not yet implemented in CD-Scheme (useCDscheme=T)". Use viscosity, `viscA4`, staggerTimeStep instead.
- Era: 2006 (check still in obcs_check.F).
- Src: mitgcm-support 2006-September 'OBC problem: spurious boundary jets with C-D coupling' http://mailman.mitgcm.org/pipermail/mitgcm-support/2006-September/thread.html

### Blow-up at sponge layer / Orlanski boundary / inconsistent lateral boundaries, "horizontal strips in currents"
- Cause: sponge restores toward OBNt/OBNs/OBNu/OBNv values that differ from interior (e.g. no-forcing case still moves); steep resolution telescoping; inconsistent OB files (e.g. wrong boundary precision, transposed z-by-x arrays); Orlanski net inflow drifts sea level by hundreds of m before crash; insufficient viscosity.
- Fix: make boundary data match the interior state at t=0; compare against the parent solution animations; smooth telescoping; turn off Orlanski as a test, use sponge (obcs or pkg/rbcs); check OBC array shapes (Nz x Nx vs Nx x Nz) and OB_Ieast/OB_Jnorth indices (-1 = last point).
- Era: 2004-2023.
- Src: mitgcm-support 2012-February 'instability issue with sponge' http://mailman.mitgcm.org/pipermail/mitgcm-support/2012-February/thread.html ; 2023-June 'Query regarding horizontal strips in currents' http://mailman.mitgcm.org/pipermail/mitgcm-support/2023-June/thread.html ; 2004-August 'Re: Open BCs' http://mailman.mitgcm.org/pipermail/mitgcm-support/2004-August/thread.html ; 2014-July 'Model based on exp4 experiment' (Dustin: transposed OBCS arrays) http://mailman.mitgcm.org/pipermail/mitgcm-support/2014-July/thread.html ; 2020-March 'suddenly producing NaNs' http://mailman.mitgcm.org/pipermail/mitgcm-support/2020-March/thread.html

## Sea ice
### STOP in CALC_R_STAR : too SMALL rStarFac[C,W,S] ! / ABNORMAL END: S/R CALC_R_STAR  (very thick ice, big eta excursions)
- Cause: with nonlinear free surface (select_rStar=2, nonlinFreeSurf=4) the surface cell thickness ratio rStarFac fell below hFacInf: model exploded, OR ice (+snow) load so heavy (tens-hundreds of m, usually in embayments/inlets where ice cannot be advected away, LGM runs, shelfice) that eta drops by hundreds of m, OR sea level drained through open boundaries. Reducing deltaT only helps in the first case. Failing cell indices are printed ("fail at i,j=...").
- Fix: read STDERR too; check ice thickness/eta; `SEAICE_no_slip=.FALSE.` (PARM, default FALSE), `#define SEAICE_CAP_ICELOAD` (only limits the load on eta, NOT ice thickness), smooth narrow inlets; if shelfice+seaice try SHELFICEboundaryLayer=.FALSE. (thick ice outside cavities); if r* not needed `select_rStar=0`; lower `hFacInf`; reduce deltaT if CFL. In AD runs it may come from tape recomputation (see Adjoint).
- Era: 2011-2025.
- Src: mitgcm-support 2015-June 'R_STAR issue' (Losch) http://mailman.mitgcm.org/pipermail/mitgcm-support/2015-June/thread.html ; 2014-July 'low eta_min' http://mailman.mitgcm.org/pipermail/mitgcm-support/2014-July/thread.html ; 2014-January 'Warning and time stepping issue' http://mailman.mitgcm.org/pipermail/mitgcm-support/2014-January/thread.html ; 2025-August 'Too large of Sea Ice' http://mailman.mitgcm.org/pipermail/mitgcm-support/2025-August/thread.html ; 2019-September 'MITgcm-support Digest, Vol 195, Issue 14' (ice velocities hundreds m/s) http://mailman.mitgcm.org/pipermail/mitgcm-support/2019-September/thread.html

### MOM_IMPLICIT_R: error when solving 3-Diag problem. (also after sea-ice crash)
- Cause: just the first routine to notice an exploded state (CFL/viscosity), here with sea-ice EVP velocities exploding; it STOPs only the rank that notices, which can hang MPI jobs (PR/issue #439).
- Fix: monitorFreq=1, find first blowing field; see EVP entry below. Hang: grep `ABNORMAL END`/`Execution ended Normally` in STDOUT.0000 to detect failure in job scripts.
- Era: 2019-2021 (#439 open as of 2026).
- Src: mitgcm-support 2021-January 'Crashes with EVP seaice' http://mailman.mitgcm.org/pipermail/mitgcm-support/2021-January/thread.html ; https://github.com/MITgcm/MITgcm/issues/439 ; https://github.com/MITgcm/MITgcm/issues/197

### EVP sea-ice solver explodes (ice velocities blow up at a point, no LSR problem)
- Cause: standard EVP with small alpha/beta (<=20) not converged; very smooth ice fields at 5 km are another sign.
- Fix: use mEVP with `SEAICE_evpAlpha=SEAICE_evpBeta` ~ 500 (>= 100) or aEVP: `SEAICEuseEVPstar=.TRUE.`, `SEAICEuseEVPrev=.TRUE.`, `SEAICEaEVPcoeff=0.5` (tuning), `SEAICEnEVPstarSteps=200-500` (cost +15-50%); or default LSR/Picard solver; for VP `SEAICEnonLinIterMax=10`. Docs section "more stable variants of EVP".
- Era: 2016-2021.
- Src: mitgcm-support 2021-January 'Crashes with EVP seaice' http://mailman.mitgcm.org/pipermail/mitgcm-support/2021-January/thread.html ; 2016-March 'SEAICE pkg, unstable HEFF' http://mailman.mitgcm.org/pipermail/mitgcm-support/2016-March/thread.html

### Runaway ice thickness (HEFF +50 m in one step), large precipitation / snow
- Cause: too much freshwater/snow on ice with snow advection and flooding turned off (`SEAICEadvSnow=.FALSE.`, `SEAICEuseFlooding=.FALSE.`) and no downward long-wave (surface loses heat to a 0 K atmosphere); implausible precipitation rate.
- Fix: enable (default TRUE) SEAICEadvSnow/SEAICEuseFlooding, supply downward LW, check precip (1.5e-6 m/s over 5e4 km2 unstable vs 1e-6 stable); SEAICE_deltaTdyn must be consistent with deltaT (the 2014 consistency check bug IF(SEAICE_deltaTdyn .LT. SEAICE_deltaTtherm) was corrected by Menemenlis).
- Era: 2014-2016.
- Src: mitgcm-support 2016-March 'SEAICE pkg, unstable HEFF' http://mailman.mitgcm.org/pipermail/mitgcm-support/2016-March/thread.html ; 2014-February 'CS64 grid' http://mailman.mitgcm.org/pipermail/mitgcm-support/2014-February/thread.html

### Huge Qnet (1e4 W/m2) where supercooled ice-shelf water surfaces
- Cause: all negative heat content below freezing point turned into ice in one step (frazil).
- Fix: `SEAICE_frazilFrac` = SEAICE_deltaTtherm/(3 days) (default 1); related SEAICE_gamma_t_frz, SEAICE_mcPheeTaper (SEAICE_PARAMS.h).
- Era: 2016.
- Src: mitgcm-support 2016-November 'on the control of sea ice growth' http://mailman.mitgcm.org/pipermail/mitgcm-support/2016-November/thread.html

### Sea ice dynamics gives NaN with no ice present (checkpoint65u LSR)
- Cause: older LSR code inverted a singular matrix when no ice anywhere; fixed in later code.
- Fix: update MITgcm (Losch declined to debug old code).
- Era: 2019 report using checkpoint65u (2016).
- Src: mitgcm-support 2019-August 'MITgcm-support Digest, Vol 194, Issue 26' http://mailman.mitgcm.org/pipermail/mitgcm-support/2019-August/thread.html

### Huge heat release to ocean (+25 C) when snow floods, SEAICE_USE_GROWTH_ADX
- Cause: seaice_growth_adx.F counted the latent heat of snow->ice flooding in QNET (flooding must be energy neutral); budget not closed.
- Fix: fixed upstream in PR #721 (2023-04) (flooding term removed from QNET; budget closure for linear free surface); results change for anyone who used flooding with ADX. Related diagnostic fix PR #703 (SIatmFW sublimation, new SIeprflx).
- Era: bug present before 2019; fixed April 2023 (checkpoint68s+ era).
- Src: https://github.com/MITgcm/MITgcm/pull/721 ; https://github.com/MITgcm/MITgcm/pull/703

### Sea-ice stress divergence missing metric terms (curvilinear grids)
- Cause/Fix: PR #976 (2026) adds the missing extra metric terms; controlled by `SEAICEselectMetricTerms` (SEAICEuseMetricTerms is partially retired); results change slightly on llc/cs grids (not a stability bug, but changes ice solution).
- Era: 2026 (PR closed Apr 2026; check tag-index).
- Src: https://github.com/MITgcm/MITgcm/pull/976

## Ice shelf (pkg/shelfice, icefront)
### shelfice_thermodynamics NaN: freshWaterFlux = NaN when saltFreeze = 0 (sLoc = 0), e.g. llc4320 cavity init
- Cause: 1 - sLoc/saltFreeze with sLoc=0 and saltFreeze=0 (salinity went to zero, e.g. undershoot at excessive initial melting).
- Fix: fixed upstream (PR #968, Jan 2026): code now skips the division when saltFreeze==0 (SHI_SALTBAL_FWFLX branch) and the default heat-balance formulation handles identical salinities; sLoc is already MAX(salt,0). Ian Fenty: zero ocean salinity signals an upstream problem (overshoot).
- Era: reported Jan 2026 (checkpoint69x); fixed in master.
- Src: mitgcm-support 2026-January 'shelfice_thermodynamics.F crash' http://mailman.mitgcm.org/pipermail/mitgcm-support/2026-January/thread.html ; https://github.com/MITgcm/MITgcm/pull/968

### Blow-up within minutes after enabling useSHELFICE (regional, OBCS+seaice)
- Cause: inconsistent `SHELFICEloadAnomalyFile` (pressure load anomaly at ice base) -> large initial free-surface adjustment -> CFL; or setting SHELFICElatentHeat/HeatCapacity=9.e99 (icefront too) to "freeze" the ice.
- Fix: recompute pload from the initial T/S (see verification/isomip/input/gendata.m or icefront 2D_example m-files); trick: set `rhoConst` ~ mean cavity density so shelfice_init_fixed.F computes a reasonable load; test with transfer coefficients (gamma_T/S in data.shelfice) = 0; shelfice/icefront geometry is fixed in time; try cAdjFreq=-1 instead of KPP while debugging (KPP works under ice shelves but only the interior-mixing part).
- Era: 2014 (still valid; names exist).
- Src: mitgcm-support 2014-September 'instability with shelf ice package' http://mailman.mitgcm.org/pipermail/mitgcm-support/2014-September/thread.html ; 2014-February 'Icefront returnin NaN based on topology.' http://mailman.mitgcm.org/pipermail/mitgcm-support/2014-February/thread.html

## Adjoint / TAF
### Adjoint run dies with CALC_R_STAR "too SMALL rStarFac" in tape (forward) computations although forward run is fine
- Cause: (a) tape levels stored in single precision: `doSinglePrecTapelev` (data.ctrl) default in ECCO; (b) with nonlinFreeSurf>2 (r*) TAF did not recompute update_cg2d in inner tapes (cg2d.flow says self-adjoint), so the inner tape forward steps drift (issue #391).
- Fix: (a) `doSinglePrecTapelev=.FALSE.` in data.ctrl to test; (b) fixed upstream in PR #392 (Dec 2020, fake adjoint of cg2d) - update code. Also use `#define ALLOW_AUTODIFF_WHTAPEIO`.
- Era: 2020-2021 (fixed after checkpoint67r-ish).
- Src: mitgcm-support 2020-August 'mitgcmuv_ad explodes during tape computations' http://mailman.mitgcm.org/pipermail/mitgcm-support/2020-August/thread.html ; https://github.com/MITgcm/MITgcm/issues/391

## Input/setup errors that look like instabilities
### Crash "forrtl: error (65): floating invalid" right after reading an edited pickup / NaN at step 1-2 with new input
- Cause: invalid values written into your hand-edited pickup/IC; or explode quickly; see restart.md for pickup editing.
- Fix: run with traceback + debugLevel=4; verify the file with rdmds (precision real*8, ieee-be).
- Era: 2016.
- Src: mitgcm-support 2016-August 'floating invalid with new pickup' http://mailman.mitgcm.org/pipermail/mitgcm-support/2016-August/thread.html

### Segmentation fault / NaN: SIZE.h inconsistent with data
- Cause: Nx = sNx*nSx*nPx must equal the number of delX values; Ny likewise (a SIZE.h with Ny=40 for an 80-row domain; cube-sphere tile size must divide the face size: ABNORMAL END: S/R W2_SET_MAP_TILES; "S/R EEBOOT_MINIMAL: No. of processes not equal to nPx*nPy"); tile sizes with large prime factors give poor decomposition.
- Fix: edit SIZE.h consistently, rebuild; for cs510 use tile sizes dividing 510 (nPy need not be 1, Campin); use dxSpacing/dySpacing for constant spacing.
- Era: 2006-2025.
- Src: mitgcm-support 2015-March 'crash with a new processor / grid size setup' http://mailman.mitgcm.org/pipermail/mitgcm-support/2015-March/thread.html ; 2017-January 'ABNORMAL END: S/R LOAD_GRID_SPACING' http://mailman.mitgcm.org/pipermail/mitgcm-support/2017-January/thread.html ; 2006-November 'What is the demand for bathyFile ?' http://mailman.mitgcm.org/pipermail/mitgcm-support/2006-November/thread.html ; 2025-February 'Reducing Runtime for High Resolution Model' http://mailman.mitgcm.org/pipermail/mitgcm-support/2025-February/thread.html

### Zero ocean velocity/odd values: bathymetry sign, grid latitude
- Cause: `bathyFile` must contain negative depths (zeros -> land); `phiMin`/ygOrigin is the SOUTHERN edge (a grid starting at -50 but meant -80 shifts Coriolis, can blow up at poles; ygOrigin near 90N tiny dx).
- Fix: check hFacC in output; check YC/fCori in grid files; polar caps need a fine dt.
- Era: 2006-2018.
- Src: mitgcm-support 2018-November 'No output every 4th depth' http://mailman.mitgcm.org/pipermail/mitgcm-support/2018-November/thread.html ; 2006-September 'Re: ...' (Losch phiMin) http://mailman.mitgcm.org/pipermail/mitgcm-support/2006-September/thread.html ; 2016-April 'Troubleshooting Gridding for Near N-Pole Simulation' http://mailman.mitgcm.org/pipermail/mitgcm-support/2016-April/thread.html

## Packages
### aim_v23 / SPEEDY: "AIM_DYN2AIM: Temp out of range 100 400" / "SUFLUX_POST: TS out of range 100 400" after centuries
- Cause: rare storm at cube corner makes the atmosphere locally unstable.
- Fix: restart from last pickup after slightly changing a filter parameter, e.g. `Shap_Trtau` 5400 -> 5300 (data.shap) (Ferreira/Rose), or smaller dt; may recur.
- Era: 2015-2017.
- Src: mitgcm-support 2015-May 'SPEEDY out of range surface air temp' http://mailman.mitgcm.org/pipermail/mitgcm-support/2015-May/thread.html ; 2017-May 'very warm climates and aim_v23' http://mailman.mitgcm.org/pipermail/mitgcm-support/2017-May/thread.html

### NaN in DIC/GCHEM tracers
- Cause: pH solver fed with wild T or S, or huge initial atmosphere-ocean pCO2 disequilibrium with interactive atmospheric CO2 (dic_int1=3) at start-up.
- Fix: check T/S first; `#define DIC_BOUNDS` (DIC_OPTIONS.h; a hack to clip inputs); increase initial atm pCO2 box.
- Era: 2011-2015.
- Src: mitgcm-support 2015-August 'NaNs in the DIC package' http://mailman.mitgcm.org/pipermail/mitgcm-support/2015-August/thread.html ; 2011-June 'NaNs in dic' http://mailman.mitgcm.org/pipermail/mitgcm-support/2011-June/thread.html

### Passive homogeneous salinity blows up after very long run
- Cause: roundoff drift in gS of a salinity that is passive/constant.
- Fix: `saltStepping=.FALSE.`.
- Era: 2016.
- Src: mitgcm-support 2016-June 'passive constant homogeneous salinity blows up' http://mailman.mitgcm.org/pipermail/mitgcm-support/2016-June/thread.html

### Floats (pkg/flt) segfault when model near blow-up
- Cause: huge velocities give huge float indices in FLT_TRILINEAR.
- Fix: fix the model instability first; (no bounds check added upstream as far as the thread says).
- Era: 2012.
- Src: mitgcm-support 2012-June 'FLT package' http://mailman.mitgcm.org/pipermail/mitgcm-support/2012-June/thread.html

## Platform, compiler, NaNs continue
### Model keeps running writing NaN after blow-up (burning CPU)
- Cause: compiler/platform does not trap FP exceptions; monitor check only fires if monitorFreq hits before NaN.
- Fix: build with `genmake2 -ieee -devel` and your optfile's FP-trap flags (gfortran `-ffpe-trap=invalid,zero,overflow -fbacktrace`; first NaN -> backtrace; keep .f files); `monitorFreq` small during tests; old hacks: STOP if rhsMax > ~1e8-1e10 after the cg2d print in cg2d.F (Mazloff, Brostrom), or dump fields in mon_solution.F on crash (Holland/Bruneau); a shell loop with `tail|grep NaN` also works.
- Era: 2006-2017.
- Src: mitgcm-support 2017-November 'NaN - stop running?' http://mailman.mitgcm.org/pipermail/mitgcm-support/2017-November/thread.html ; 2011-July 'model blew up but not terminated' http://mailman.mitgcm.org/pipermail/mitgcm-support/2011-July/thread.html ; 2015-March 'MITgcm-support Digest, Vol 141, Issue 6' http://mailman.mitgcm.org/pipermail/mitgcm-support/2015-March/thread.html ; 2007-March 'NaNQ' http://mailman.mitgcm.org/pipermail/mitgcm-support/2007-March/thread.html

### Note: IEEE_INVALID_FLAG / IEEE_OVERFLOW_FLAG printed at "STOP NORMAL END"
- Cause: floating-point exceptions were signalled somewhere (NaN/overflow) but the run did not stop; with CFL fine in the monitor lines it can still be vertical/other.
- Fix: find where (FP trap flags + debugger, see above) before changing dt; Menemenlis: look at the CFL numbers first.
- Era: 2017, 2025.
- Src: mitgcm-support 2017-November 'IEEE_INVALID_FLAG and IEEE_OVERFLOW_FLAG' http://mailman.mitgcm.org/pipermail/mitgcm-support/2017-November/thread.html ; 2025-April 'the floating-point exception issue' http://mailman.mitgcm.org/pipermail/mitgcm-support/2025-April/thread.html

### Blow-up depends on compiler flags / optimization / tile count (non-reproducible)
- Cause: marginally stable setups are sensitive to round-off; compiler optimization with fused ops or higher internal precision; odd per-tile sizes reported to trigger optimizer noise (Klymak 2025, unconfirmed). Symmetric flat-bottom setups can also appear "stable" at -O0 only because no instability seeds (add small noise).
- Fix: test with `-ieee`/-O0; run verification via `testreport` with the same optfile; increase viscosity slightly or reduce dt; MPI vs non-MPI comparison of a simple experiment.
- Era: 2008-2025.
- Src: mitgcm-support 2018-June 'Reproducibility of blowup' http://mailman.mitgcm.org/pipermail/mitgcm-support/2018-June/thread.html ; 2015-January 'Baroclinic instability with MPI run' http://mailman.mitgcm.org/pipermail/mitgcm-support/2015-January/thread.html ; 2025-February 'Reducing Runtime for High Resolution Model' http://mailman.mitgcm.org/pipermail/mitgcm-support/2025-February/thread.html ; 2008-September 'problem running tutorial_global_oce_latlon' http://mailman.mitgcm.org/pipermail/mitgcm-support/2008-September/thread.html

### Build/run: "relocation truncated to fit: R_X86_64_PC32 ... COMMON" / gfortran "-mcmodel=medium" unrecognized (arm64)
- Cause: static arrays > 2 GB per process (large tile): need medium/large code model; on aarch64 gcc does not support -mcmodel=medium.
- Fix: x86: add `-mcmodel=medium` (and `-fPIC`, NetCDF built compatibly) or use more MPI processes (smaller tiles); aarch64/arm64 (Mac M-series, Docker): genmake2 maps aarch64 to arm64 optfiles (commit 105ccd98, 2024-06, issue #847) - use `linux_arm64_*`/darwin arm64 optfile, not linux_amd64_gfortran ('-mcmodel=small' as a stop-gap).
- Era: 2006-2024 (#847 checkpoint68y).
- Src: https://github.com/MITgcm/MITgcm/issues/847 ; mitgcm-support 2006-September 'Bug R_X86_64_PC32' http://mailman.mitgcm.org/pipermail/mitgcm-support/2006-September/thread.html ; 2010-April "problem of 'relocation truncated to fit'" http://mailman.mitgcm.org/pipermail/mitgcm-support/2010-April/thread.html

### Job hangs after a rank STOPs (MPI deadlock) / run "dies" but SLURM keeps burning hours
- Cause: individual STOP from one rank (MOM_IMPLICIT_R, MDS_READ_*, CALC_R_STAR) without MPI_ABORT; only ALL_PROC_DIE is MPI-safe; also disk/network glitches on one node.
- Fix: scripts: detect `ABNORMAL END` in STDOUT/STDERR or missing "Execution ended Normally"; `writePickupAtEnd` or pChkptFreq so a pickup exists at the end; add watchdog/timeouts; issue #439 open (idea: EMERGENCY_STOP calling MPI_ABORT).
- Era: 2019-2021 (#197 closed 2023, #439 open).
- Src: https://github.com/MITgcm/MITgcm/issues/197 ; https://github.com/MITgcm/MITgcm/issues/439 ; mitgcm-support 2020-March 'jobs died suddenly' http://mailman.mitgcm.org/pipermail/mitgcm-support/2020-March/thread.html

### Stderr/stdout files: STOP/ABNORMAL END but STDERR empty / "Cannot open file '#dev/null'" with SINGLE_DISK_IO
- Cause: `#define SINGLE_DISK_IO` (CPP_EEOPTIONS.h) makes ranks > 0 open /dev/null (problem on ARCHER2 with gfortran) and drops all messages from ranks != 0; (warning in the file).
- Fix: keep SINGLE_DISK_IO undefined unless the set-up is proven; debugMode=.FALSE. in eedata also reduces STDOUT volume.
- Era: 2014, 2021.
- Src: mitgcm-support 2021-March 'Crashing at run-time because /dev/null is not available' http://mailman.mitgcm.org/pipermail/mitgcm-support/2021-March/thread.html ; 2014-May 'reducing the size of STDOUT files etc' http://mailman.mitgcm.org/pipermail/mitgcm-support/2014-May/thread.html
