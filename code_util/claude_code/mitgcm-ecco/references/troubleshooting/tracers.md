# Troubleshooting: tracers, BGC, mixing and viscosity
Scope: pkg/ptracers, gchem, dic, bling, darwin (as seen on the list), RBCS/OBCS-for-tracers, offline; GMRedi, KPP, GGL90, other vertical mixing, horizontal viscosity (Leith/Smag/viscAhGrid), convection. Distilled from mitgcm-support 2003-2026 + MITgcm GitHub issues/PRs; names checked against upstream origin/master (checkpoint ~2026-10) and darwin3. "Src" thread URLs are the month index at mailman.mitgcm.org/pipermail/mitgcm-support/<YYYY-Month>/thread.html.

## ptracers set-up, restart, namelists

### "invalid reference to variable in NAMELIST input" / "namelist not terminated with / or &end" / "End of file" in ptracers_readparms or any data.* read
- Cause: (a) a stray "&" mid-namelist terminates the read (people put "&" between tracers); (b) PTRACERS_num in PTRACERS_SIZE.h smaller than the indices used in data.ptracers (default 1); (c) compiler wants "/" not "&" as terminator; (d) float typos like `1e-2.,`.
- Fix: one `&PTRACERS_PARM01 ... &` block only; set PTRACERS_num (PTRACERS_SIZE.h, via -mods code dir) >= PTRACERS_numInUse and rebuild; for "/" compilers add `-DNML_TERMINATOR` to DEFINES in the optfile (exists in tools/build_options/darwin_amd64_gfortran; handled in eesupp/src/nml_set_terminator.F). Minimal debug: strip `#` lines from data.* to mimic the temp file the model reads.
- Era: all years (2003-2026). data.kpp used to swallow bad values silently (IOSTAT read); in current kpp_readparms.F it is a plain READ again, so it now stops.
- Src: mitgcm-support 2018-February 'ptracer error' http://mailman.mitgcm.org/pipermail/mitgcm-support/2018-February/thread.html ; 2017-October 'Digest Vol 172 Issue 9' http://mailman.mitgcm.org/pipermail/mitgcm-support/2017-October/thread.html ; 2014-February 'Fortran problem with MITgcm' (NML_TERMINATOR, JMC) http://mailman.mitgcm.org/pipermail/mitgcm-support/2014-February/thread.html ; 2017-October 'KPP parameters not read?'

### PTRACERS_Iter0 > nIter0 stops run ("wrong setting of PTRACERS_Iter0 ... nIter0 < PTRACERS_Iter0 < nEndIter not supported"), or tracers are all zero
- Cause: PTRACERS_Iter0 is the iteration at which tracers are *initialised*; starting them later in a running integration was never supported (ptracers_check.F now stops).
- Fix: (1) per-tracer `PTRACERS_startStepFwd(n)` = time in s at which tracer n starts being stepped (saves CPU, tracer stays at initial value); (2) split run in two and set PTRACERS_Iter0 = nIter0 for segment 2 (reads PTRACERS_initialFile, no pickup_ptracers needed); (3) PTRACERS_resetFreq/PTRACERS_resetPhase re-initialise tracers periodically (does NOT skip stepping). PTRACERS_Iter0 < nIter0 requires pickup_ptracers.* (cfc_example/input shows physics-from-pickup, tracers-from-file).
- Era: 2008-2025; check stop present in current ptracers_check.F.
- Src: https://github.com/MITgcm/MITgcm/issues/911 ; mitgcm-support 2008-November 'PTRACERS Initialization' ; 2020-November 'Restart simulation with biogeochemical packages' ; 2010-March 'ptracer start at advanced time'

### Restart with pickupSuff='ckptA' or after enabling bgc: DIC/ptracers NaN or zeros at domain edges
- Cause: nIter0 (or startTime) not set when using pickupSuff; tracers restarted without matching pickup_ptracers/pickup_dic; older versions: missing pickup_dic/cpl fields.
- Fix: always set nIter0 explicitly (pickupSuff only picks the file name); to restart physics only and start bgc fresh use PTRACERS_Iter0 = nIter0 plus PTRACERS_initialFile (or PTRACERS_ref profile; all-zero if neither); pkg/dic pickup now read with meta check, `pickupStrictlyMatch=.FALSE.` to tolerate missing fields.
- Era: 2010-2023 (PR #757, 2023: dic_read_pickup.F meta handling).
- Src: mitgcm-support 2014-August 'Forcing files in data.dic (biogeochemistry)' http://mailman.mitgcm.org/pipermail/mitgcm-support/2014-August/thread.html ; https://github.com/MITgcm/MITgcm/pull/757

### All passive tracers behave like salinity / tracer concentrations track S
- Cause: ptracers defaults copy salt: PTRACERS_advScheme/diffKh/diffK4/diffKr default to saltAdvScheme/diffKhS/diffK4S/diffKrS (ptracers_readparms.F), vertical mixing (KPP/GM/BL79) uses the salt kappa, and old ptracers_forcing_surf.F example applied salt forcing to tracers.
- Fix: set PTRACERS_diffKh/diffK4/diffKr/advScheme explicitly in data.ptracers; zero or replace surfaceForcingPTr in ptracers_forcing_surf.F; PTRACERS_EvPrRn(n) = tracer concentration of E-P-R water: 0. lets rain dilute/evaporation concentrate (DIC, ALK), leave unset to keep conserved (PO4/DOP); vertical diffusivity of ptracers cannot differ from S under KPP except the background part.
- Era: 2003-2020 (behaviour unchanged).
- Src: mitgcm-support 2004-September 'New user with problems with ptracers' ; 2003-September 'problems with ptracer package' ; 2009-September 'Bryan-Lewis vertical diffusivity for passive tracers' ; 2020-November (Lauderdale on EvPrRn)

### How to add a source/sink, surface flux, runoff tracer, melt-water tracer
- Cause: pkg/ptracers only advects/mixes; any forcing is user code.
- Fix: surface fluxes -> `ptracers_forcing_surf.F` (surfaceForcingPTr, applied in ptracers_apply_forcing even under shelfice); 3-D sources -> `ptracers_apply_forcing.F` (called per k, add to gPtracer); with gchem -> `gchem_calc_tendency.F`; runoff tag: add runoff*rA*drF*hFac to surfaceForcingPTr at end of routine and PTRACERS_EvPrRn(n)=0.; local->global indices: ig = myXGlobalLo-1+(bi-1)*sNx+i. Age tracer example: verification/tutorial_global_oce_latlon/code (relax surface to 0, source 1 below; better use RBCS for the relaxation so tau is a parameter). Sea-ice uptake: `#define ALLOW_SITRACER` (SEAICE_OPTIONS.h), edit seaice_tracer_phys.F; when including SEAICE.h in ptracers code also include SEAICE_SIZE.h and SEAICE_PARAMS.h or you get "COMMON block data object must not be an automatic object".
- Era: 2009-2020, still current.
- Src: mitgcm-support 2020-July 'ptracers and shelfice package integration' ; 2009-September 'tracers' (runoff) ; 2017-December 'time continuous source for ptracers' ; 2020-April 'use seaice variables in ptracers' ; 2017-February 'Ptracers in sea ice' ; 2017-May 'ideal age tracer'

### Turn off advection (or change vertical scheme) for a single passive tracer
- Cause: no per-tracer vertical scheme; tempAdvection/saltAdvection switches don't exist for ptracers.
- Fix: `PTRACERS_advScheme(n)=0` disables advection (forcing + diffusion remain; undocumented, JMC). Vertical scheme = PTRACERS_advScheme (no PTRACERS_vertAdvScheme; add one in ptracers_integrate.F if needed).
- Era: 2017 (advScheme=0), 2020.
- Src: mitgcm-support 2017-November 'switch to turn off advection' ; 2020-October 'Vertical Advection Scheme for PTRACERS'

### RBCS relaxation of tracers: "Non-existing record number", namelist error, tracer never reaches target
- Cause: maskLEN in RBCS_SIZE.h (not RBCS.h, gone) too small: tracer iTr uses relaxMaskFile(2+iTr) (T=1, S=2); mask/target file records must match rbcsForcingPeriod/Cycle; with tauRelaxPTR > deltaT the tracer only approaches the target (T_new = T_r only if tau = deltaT).
- Fix: raise maskLEN >= 2+PTRACERS_numInUse and rebuild (or leave relaxMaskFile(3..) unset to share one mask); use rbcsForcingPeriod=0. for time-constant fields; tauRelaxPTR(n)=deltaT for a hard source (point source with tiny tau can go unstable: smooth the mask); rbcsIniter retired. Current rbcs_readparms.F prints "Increase maskLEN (in RBCS_SIZE.h) and recompile".
- Era: 2014-2023; clean STOP added Nov 2017.
- Src: mitgcm-support 2017-November 'Adding multiple tracers to MITgcm' ; 2014-June 'Namelist Error with RBCS and multiple tracers' ; 2023-December 'Trying to relax ptracer using RBCS' ; 2017-December 'time continuous source...'

### Open boundaries and ptracers: tracer explodes (1e38, negative), accumulates at boundary, or zeros leak in
- Cause: (a) no OB file for ptracers -> default is zero-gradient (Neumann), not "zero inflow"; (b) Orlanski + ptracers: obcs_check.F stops ("useOrlanski* OBC not yet implemented for pTracers"); commenting it out lets tracer pile up at the OB; (c) open boundary point next to land: old obcs_apply_ptracer put 0 in wet cells (-> NaN in carbon chemistry).
- Fix: prescribe `OBNptrFile(n)/OBSptrFile/OBEptrFile/OBWptrFile(n)` (verification/exp4, so_box_biogeo); `OBCS_u1_adv_Tr(n)=1` forces 1st-order upwind at outflow; or restore with RBCS; `OBCSfixTopo=.TRUE.` (data.obcs) removes topography gradients across OB; use useOrlanski* per boundary so the ptracer boundary is non-Orlanski. OBCS sponge does not act on ptracers: use pkg/rbcs.
- Era: 2004-2022; Orlanski+ptracer check still present in master.
- Src: mitgcm-support 2019-January 'OBCS Tracer Errors' ; 2019-April 'Behavior of Ptracers at Orlanski boundary' ; 2012-October 'Way around ptracers & Orlanski OBCS?' ; 2008-April 'obcs_apply_ptracer' ; 2018-July 'OBCS sponge for ptracers' ; 2017-April 'ptracer and default boundary conditions'

### Negative tracer concentrations after advection (undershoot)
- Cause: linear/centered schemes (default 2) and even limited schemes give tiny (1e-8 relative) undershoots; sharp initial patches worsen it.
- Fix: flux-limited schemes (33, 77) or 7 (needs OLx=OLy=4), smooth initial patch, small Kh; `PTRACERS_stayPositive(n)=.TRUE.` (Smolarkiewicz hack, in PTRACERS_PARAMS.h today); DIC_NO_NEG / BLING_NO_NEG reset negatives but break conservation. Scheme choice trades spurious diffusion vs extrema (7, 80/81 Prather as compromise).
- Era: 2014-2021.
- Src: mitgcm-support 2014-December 'Negative concentrations in pTracer' ; 2021-February 'Negative Passive Tracer Concentration' ; 2021-February '[EXTERNAL] Spurious mixing with internal tides' ; 2015-August 'pkg/dic options ?'

### SIGFPE / blow-up in gad_fluxlimit_* (Rj ~ 1e-307, "Rjp,Rj" printed)
- Cause: Cr = Rjm/Rj with tiny non-zero Rj (seen on sea ice at high res).
- Fix: fixed upstream (issue #459; current gad_fluxlimit_adv_x.F tests `ABS(Rj)*CrMax .LE. ABS(Cr)` rather than a fixed 1e-20 threshold).
- Era: 2021 (April), fixed.
- Src: https://github.com/MITgcm/MITgcm/issues/459

### Tracer content drifts / surface-cell budget does not close (linear free surface)
- Cause: with linFSConserveTr false and a linear free surface, w(k=1)!=0 changes surface-cell tracer; ptracers do not get the T/S correction unless asked; fresh-water dilution of ptracers handled by PTRACERS_EvPrRn/PTRACERS_ref.
- Fix: `PTRACERS_linFSConserve(n)=.TRUE.` (data.ptracers), `linFSConserveTr=.TRUE.` (data, T/S), `darwin_linFSConserve=.TRUE.` (data.darwin) - the three are separate switches, set per package in use; `exactConserv=.TRUE.`; best: r* (`#define NONLIN_FRSURF`, `nonlinFreeSurf=4, select_rStar=2`, hFacInf/hFacSup) gives near machine conservation for DIC (Munday). PTRACERS_addSrelax2EmP=.TRUE. adds S-restoring freshwater to tracer dilution (check global DIC/ALK/PO4 afterwards: Voelker+Munday saw drift, unresolved).
- Era: 2010-2024.
- Src: mitgcm-support 2024-March 'Closing DIC budget' http://mailman.mitgcm.org/pipermail/mitgcm-support/2024-March/thread.html ; 2010-June 'Mass balance... again...' ; 2010-August 'Salinity restoring and tracers' ; 2011-June 'dic pkg and PTRACERS_addSrelax2EmP'

### Closing a tracer budget from diagnostics (ADV*, DF*, KPPg_*, TOTTTEND) does not balance
- Cause: need all pieces: ADVx/y/r_*, DFxE/DFyE/DFrE/DFrI_*, KPPg_* (non-local KPP flux, diag name KPPg_SLT/KPPg_TH/KPPgTr##; KPPghatK is only the shape function), surface forcing; fluxes already include cell area (m^3/s) so divide by hFac*rA*drF; AB-based advection schemes (2,3,4) are not stored the same way as multi-dim schemes (30,33,77: straightforward); time-average vs snapshot timePhase (default snapshot phase = frequency/2: set timePhase=0 to align). NLFS: d(hS)/dt not available directly (issue #1010).
- Fix: use diagnostics TOTTTEND/TOTSTEND for dT/dt, add KPP non-local flux, snapshot with identical timePhase; for NLFS save snapshots of hFac*S and difference offline (Losch); read doc/Heat_Salt_Budget_MITgcm.pdf.
- Era: 2007-2026 (open design discussion in #1010).
- Src: mitgcm-support 2007-September 'balancing the heat equation' ; 2010-April 'KPP and Heat Budget' ; 2017-September 'On the exact calculation of ADVx/y/r_TH and DFrx/y/rE_TH' ; https://github.com/MITgcm/MITgcm/issues/1010

### Which flux diagnostics hold which tracer-mixing contribution (DFrE vs DFrI, GM, KPP)
- Cause: with useGMRedi the Redi/GM terms proportional to horizontal gradients (Kwx, Kwy) are always explicit -> DFrE_* even if implicitDiffusion=.TRUE.; the Kwz part is in DFrI_* (implicit) or DFrE_* (explicit). Skew-flux GM (GM_AdvForm=.FALSE., default) lands in DFxE/DFyE/DFrE; advective-form GM is inside ADVx/y/r. With KPP and no GMRedi DFrE_* is unfilled ("has not been filled ... write ZEROS").
- Fix: sum DFrE+DFrI (+KPPg_); GM_PsiX/Y, GM_Kuz etc. are zero in skew-flux mode. GM_VisbK is only the Visbeck K diagnostic.
- Era: 2007-2019, behaviour stable.
- Src: mitgcm-support 2007-March 'advective and diffusive diagnostics' (JMC 2017 reply) ; 2013-November 'diagnosing gmredi fluxes for tracers' ; 2019-February 'GM-Redi question'

### Slow runs with many ptracers (BLOCKING_EXCHANGES dominates), 51 tracers 13x slowdown
- Cause: PTRACERS_FIELDS_BLOCKING_EXCH exchanges each tracer separately (n MPI exchanges of 3-D fields).
- Fix: copy tracers into one 5-D array (k = nTr*Nr) and exchange once (Losch recipe, needed exch hack); LONGSTEP package (deltaT*6) for the tracer step; consider offline mode.
- Era: 2010-2014 reports; check current ptracers_fields_blocking_exch.F before hacking.
- Src: mitgcm-support 2010-March 'BLOCKING_EXCHANGES slowdown when using pkg/ptracers' ; 2014-September 'Model optimization when using many ptracers'

### Segmentation fault at start/pickup write/first step with more levels or Darwin
- Cause: usually stack/static memory: executable bigger than limits; occasionally mismatched subroutine args under CPP.
- Fix: `ulimit -s unlimited` in the job script; ifort/pgi `-mcmodel=medium` in FFLAGS (not CPP); check `size mitgcmuv`; rebuild with `genmake2 -devel` (bounds/fpe traps, `debugMode=.TRUE.` in eedata flushes STDOUT); recompile after any PTRACERS_SIZE.h/SIZE.h edit.
- Era: 2005-2021.
- Src: mitgcm-support 2018-May 'segmentation fault' ; 2017-December 'Failure to create pickup file for DIC' ; 2006-August 'relocation truncated to fit: R_386_32' ; 2013-August 'Mysterious initialisation problem'

### Offline tracer runs: set-up and gotchas
- Cause/Fix: `useOFFLINE=.TRUE.` in data.pkg (not just packages.conf) or fields are never read (zeros); fields are assumed centred in each forcing period, so files named 000..720 are interpolated Dec->Jan; shift with `offlineOffsetIter` (data.off); `deltaToffline` only names/selects files, advection uses deltaTtracer (check CFL, advcfl_* in STDOUT); offline cannot use GGL90 (no diffusivity input), KPP needs the KPP files; KPP_OUTPUT is now called for offline (kpp diagnostics KPPghatK used to be unfilled); online-vs-offline differences at 1e-8 come from time interpolation. Cost warning: offline cfc/bgc is often not faster (read cost).
- Era: 2011-2026; Menemenlis 2026: tutorial_cfc_offline.
- Src: mitgcm-support 2011-June 'Offline tracer advection' ; 2026-June 'Sequential Physics and Biogeochemistry Simulations' http://mailman.mitgcm.org/pipermail/mitgcm-support/2026-June/thread.html ; 2019-June 'nondivergent flow advecting a tracer' ; 2016-August 'KPP diffusion of ptracers' ; 2018-January 'offline KPP ptracers implicit diffusion'

### deltaTtracer / deltaTmom / deltaTClock: raising deltaTtracer alone gives no speed-up, ice does not form, tendencies tiny
- Cause: these do not mean "sub-steps": number of steps = (endTime-startTime)/deltaTClock; setting deltaTmom+deltaTtracer defaults deltaTClock to deltaTmom; asynchronous stepping is Bryan tracer acceleration (physics distortion, spin-up only; seaice defaults to deltaTtracer).
- Fix: use only deltaT for normal runs; for accelerated spin-up set deltaTmom, deltaTtracer=deltaTfreesurf=deltaTClock; finish with synchronous steps (30 y after tracer-accelerated spin-up).
- Era: 2004-2024.
- Src: mitgcm-support 2024-September 'Inquiry about deltaTtracer, deltaTmom, and deltaTClock' ; 2005-November "tracer's time step" ; 2018-December 'Time stepping for momentum and tracers' ; 2016-October 'Varying deltaTtracer in a barotropic model'

## gchem / DIC / BLING / Darwin

### GCHEM_CHECK: "Number of GCHEM tracers gchem_Tracer_num = N exceeds number of pTr: PTRACERS_numInUse"
- Cause: gchem_Tracer_num is set by each bgc package (gchem_tr_register.F); ptracers array must hold at least that many.
- Fix: PTRACERS_num (PTRACERS_SIZE.h) and PTRACERS_numInUse >= required; to drop a tracer change the package option first (e.g. `#undef ALLOW_FE`/`ALLOW_O2` in DIC_OPTIONS.h) and then reduce numInUse; more ptracers than gchem uses is fine.
- Era: 2021 (message present today).
- Src: mitgcm-support 2021-July 'A question about gchem pkg' http://mailman.mitgcm.org/pipermail/mitgcm-support/2021-July/thread.html

### "useDIC and useCFC cannot both be .TRUE." / running two bgc packages together
- Cause: shared tracer/surface-flux arrays in gchem.
- Fix: upstream gchem_check.F no longer lists DIC+CFC as exclusive (cfc has its own branch), but BLING+DIC, BLING+DARWIN, DARWIN+DIC still stop; tracer order = package registration order.
- Era: 2015 report; fixed upstream (verify with your checkpoint).
- Src: mitgcm-support 2015-June 'Running the DIC and CFC under the GCHEM' http://mailman.mitgcm.org/pipermail/mitgcm-support/2015-June/thread.html

### BLING with `#define USE_SIBLING`: SIGSEGV in bling_bio_nitrogen_ (G_SI)
- Cause: call to BLING_BIO_NITROGEN in bling_main.F was missing the G_SI argument; PTRACERS_num must be 9 (10 with ADVECT_PHYTO; biomass is a later tracer, order matters in data.ptracers).
- Fix: fixed in PR #533 (Sep 2021); else add `#ifdef USE_SIBLING O G_SI #endif` in the call and rebuild clean.
- Era: 2018 (67j) - 2021; fixed upstream.
- Src: https://github.com/MITgcm/MITgcm/pull/533 ; mitgcm-support 2021-September 'Error related to BLING' http://mailman.mitgcm.org/pipermail/mitgcm-support/2021-September/thread.html ; 2022-March 'Issue with nitrate and total phytoplankton biomass in BLING'

### BGC forcing files (wind, ice, atmospheric pCO2, silica, iron, PAR) not read / wrong units / "Cannot match namelist object name apco2startdate1" / "attempt to access non-existent record"
- Cause: forcing filenames live in the package namelist (data.dic: DIC_windFile, DIC_iceFile, DIC_silicaFile, DIC_parFile...; data.bling: BLING_windFile, BLING_ironFile...), not data.exf; windfile is wind *speed* (not u,v); with seaice/thsice the ice file is overwritten in dic_fields_update.F; EXF wspeed used only by BLING with USE_BLING_V1; READ_PAR cannot be combined with USE_QSW; exf-style BGC forcing via pkg/gchem is still an open PR (#876, issues #950/#138).
- Fix: put exf-style apco2 parameters in data.bling (apco2file etc.); DIC_forcingCycle/Period (or main externForcing*) control record timing; for wind from exf in DIC add `#include "EXF_FIELDS.h"` and overwrite `wind` (Losch). Darwin nutrient flux files are volume sources in mmol/m^3/s.
- Era: 2012-2025; USE_EXFCO2 no longer exists in BLING_OPTIONS.h.
- Src: mitgcm-support 2021-June 'Using wind from exf in DIC package' ; 2023-August 'windspeed calculation in dic package' ; 2012-June 'Biogeochemical package with #define READ_PAR' ; https://github.com/MITgcm/MITgcm/pull/876 ; https://github.com/MITgcm/MITgcm/issues/950 ; mitgcm-support 2014-September 'Darwin nutrient flux file units'

### DIC / BLING NaN in pH solver, DIC -> NaN after restart or at river runoff points
- Cause: pH/pCO2 solver fed out-of-range T/S/DIC/ALK (low-salinity runoff dilutes ALK but prescribed silica stays high -> carbonate alkalinity <0; huge initial outgassing with interactive atmosphere dic_int1=3; OBCS corner problems).
- Fix: check T,S first; `#define DIC_BOUNDS` (hack clamping solver inputs); `DIC_NO_NEG` (reset negatives); BLING has caps on silicate alk (<=20% of TA) and carbonate alk (>=10% TA) in bling_carbon_chem.F; give runoff high ALK via PTRACERS_ref or reduce silica; `CARBONCHEM_SOLVESAPHE` + `selectPHsolver` (BLING version had bugs, see next).
- Era: 2011-2023.
- Src: mitgcm-support 2015-August 'NaNs in the DIC package' ; 2023-September 'Query regarding BLING' http://mailman.mitgcm.org/pipermail/mitgcm-support/2023-September/thread.html ; 2011-June 'NaNs in dic'

### Carbonate system bugs seen in old checkouts (check your checkpoint)
- Cause/Fix: (a) Ksp calcite/aragonite used ln(T) instead of log10 (values ~1e224), omegaC mixed mol/m3 vs mol/kg (dic issue #513/#515 -> PR #514; BLING #521 -> PR #776, Oct 2023); (b) CARBONCHEM_TOTALPHSCALE filled ak2 from ak1, pH conversions in DIC_COEFFS_SURF/DEEP wrong (#269 -> PR #281, Nov 2019); (c) out-of-bounds hFacC(i,j,Nr+1) in car_flux_omega_top (PR #176, 2018); (d) same typo fixed in pkg/dic but not bling carbon_chem (2025, motivation for common gchem carbon-chem, issue #138); (e) gsm_s uninitialised (2009, fixed).
- Era: fixed upstream as listed; calcite-saturation option is now `DIC_CALCITE_SAT` (was CAR_DISS, expensive, evaluated infrequently), 3-D silica read supported since PR #620/#683.
- Src: https://github.com/MITgcm/MITgcm/issues/521 ; https://github.com/MITgcm/MITgcm/pull/281 ; https://github.com/MITgcm/MITgcm/pull/176 ; https://github.com/MITgcm/MITgcm/issues/138 ; mitgcm-support 2021-August 'solubility product of calcite and aragonite' ; 2021-December 'silica fields in DIC package'

### pkg/dic forcing/initial-field loading bugs
- Cause: `DIC_forcingCycle=0` skipped loading constant-in-time forcing files; initial pH used silica file with wrong timing (main model externForcingCycle instead of DIC_forcingCycle) -> wrong record or too few records; READ_PAR failed to compile (2012).
- Fix: fixed in PR #757 (Aug 2023, plus dic_ini_forcing.F restructure); update or patch dic_init_*.F.
- Era: before 2023-08.
- Src: https://github.com/MITgcm/MITgcm/pull/757 ; mitgcm-support 2012-June 'Using biogeochemical package with #define READ_PAR'

### Closing a DIC budget: DICBIOA/DICCARB don't match tendency
- Cause: DICBIOA is PO4 uptake (mol P); biological DIC term = Rcp*(-DICBIOA + DICPFLUX + DICRDOP) plus DICCARB; surface cell also needs air-sea flux DICTFLX (and virtual flux/dilution with linear free surface).
- Fix: multiply by Rcp (117), add PFLUX/RDOP (PFLUX has no diagnostic, may need code), use ForcTr01/gchem tendencies; see Lauderdale notebook; r* helps closure.
- Era: 2024.
- Src: mitgcm-support 2024-March 'Closing DIC budget' http://mailman.mitgcm.org/pipermail/mitgcm-support/2024-March/thread.html

### Fresh-water dilution of DIC/ALK/PO4 and PTRACERS_ref confusion
- Cause: with r*/nonlinear free surface dilution is automatic; with linear free surface PTRACERS_EvPrRn(n)=0. gives dilution (DIC, ALK), unset conserves (PO4/DOP); `ALLOW_OLD_VIRTUALFLUX` is the old virtual-flux alternative.
- Fix: choose one mechanism; unset PTRACERS_initialFile -> tracer = PTRACERS_ref profile.
- Era: 2015-2023.
- Src: mitgcm-support 2015-August 'pkg/dic options ?' ; 2020-November 'Restart simulation with biogeochemical packages'

### Darwin3 (ecco_darwin) compile errors "DARWIN_partScav / DARWIN_minFeLoss / ironSedFlux has no type"; radtrans
- Cause: ecco_darwin setup code (code_darwin) out of sync with darwin3 pkg version.
- Fix: update both repos together (ecco_darwin v05 and darwin3 master); radtrans/darwin3 to be added to the main MITgcm tree (issue #994, radtrans first).
- Era: 2022, 2026.
- Src: mitgcm-support 2022-March 'A question about darwin package' http://mailman.mitgcm.org/pipermail/mitgcm-support/2022-March/thread.html ; https://github.com/MITgcm/MITgcm/issues/994

### Delaying biogeochemistry until physics is spun up
- Cause: no run-time switch in gchem to skip bgc for first N steps.
- Fix: diagnostics start later with `timePhase(n)`; for bgc either restart as above (PTRACERS_Iter0=nIter0) or guard bling/darwin call in gchem_forcing_sep.F with a time test (gchem_int*/gchem_rl* placeholders); ptracers also still step unless PTRACERS_startStepFwd(n) is set.
- Era: 2019.
- Src: mitgcm-support 2019-May 'Delaying diagnostics, gchem'

## GMRedi

### Only Redi (no GM), or only GM; per-tracer GM
- Cause: skew-flux form (default) makes GM part of the Redi tensor; many diagnostics are zero.
- Fix: Redi only: `GM_background_K=0.`, `GM_isopycK=<K_redi>` (needs `#define GM_EXTRA_DIAGONAL`, default); if only one is set the other copies it; GM_isopycK != GM_background_K not possible with Visbeck variable K. diffKhT/S should normally be 0 with GMRedi (allowed, but adds horizontal diffusion; tutorial doc had 1e3 wrongly). GM for ptracers cannot be switched independently of T/S in 2004/2018 answers (PTRACERS_useGMRedi/PTRACERS_useKPP params exist today in PTRACERS_PARAMS.h; test before relying); hack on tracer identity in gmredi_calc_diff/x,ytransport. Eddy-resolving: turn GM off, Redi optional.
- Era: 2004-2019.
- Src: mitgcm-support 2007-October 'gmredi without gm' ; 2019-February 'GM-Redi question' ; 2011-April 'Temperature drift in tutorial (tutorial_global_oce_latlon)' ; 2018-November 'Different GMRedi coefficients for ptracers' ; 2004-July 'turning off gent-mcwilliams'

### 3-D (or 1-D x 2-D) GM/Redi diffusivity from files
- Cause: 3-D K was only available via pkg/ctrl (KapGM, KapRedi).
- Fix: CPP `GM_READ_K3D_REDI` / `GM_READ_K3D_GM` (GMREDI_OPTIONS.h; replaces old ALLOW_KAPGM_3DFILE/KAPREDI names) with `GM_isopycK3dFile`/`GM_background_K3dFile`; or without ctrl: K = GM_background_K*GM_bolFac1d(k)*GM_bolFac2d(i,j) via GM_bol1dFile/GM_bol2dFile/GM_iso1dFile/GM_iso2dFile. PR #561: crash if files set but useGMRedi=.FALSE. (fixed). Bates scheme renamed GM_BATES_K3D / GM_useBatesK3D; GEOMETRIC = GM_GEOM_VARIABLE_K.
- Era: 2018-2022.
- Src: mitgcm-support 2018-October 'GMRedi ptracers Diagnostics' ; https://github.com/MITgcm/MITgcm/issues/566 ; https://github.com/MITgcm/MITgcm/pull/561

### GMRedi under ice shelves / dry cells at top (shelfice): strong spurious mixing, tracer non-conservation, blow-up with advective GM
- Cause: pkg/gmredi assumed top at k=1; missing masks in gmredi_calc_psi_b.F, x/ytransport dTdz, wrong tapering ('linear' extrapolates psi to surface).
- Fix: PR #593 (Jan 2022) fixes masking; gmredi_check.F now stops for tapering schemes/Visbeck/mixed-layer schemes not yet fixed (use 'gkw91'/'ldd97' style options, not 'linear'/'fm07'); skew-flux with K_GM=K_Redi avoids the worst bug (a); advective GM under ice still unstable for some users; pkg/kpp has no ice-shelf support (issue #588), ggl90 fixed by PR #597.
- Era: 2016-2023.
- Src: https://github.com/MITgcm/MITgcm/pull/593 ; https://github.com/MITgcm/MITgcm/issues/591 ; https://github.com/MITgcm/MITgcm/issues/588 ; mitgcm-support 2016-March 'shelfice and GMREDI'

### GMREDI_WITH_STABLE_ADJOINT / GM_taper_scheme='stableGmAdjTap' (ECCO v4r5): wrong TAF keys in gmredi_slope_limit.F
- Cause: key computation ignored the k index and the 3 calls per level; tapes held slope of level 1 of last call.
- Fix: PR #686 (Jan 2023) uses local tapes; no change seen in verification or ice-shelf runs (gradients machine-precision identical). Related: `useGMRediInAdMode/useKPPinAdMode=.FALSE.` drop their contribution from gradients only (forward unaffected); GGL90 has stable AD.
- Era: 2022-2023.
- Src: https://github.com/MITgcm/MITgcm/issues/668 ; https://github.com/MITgcm/MITgcm/pull/686 ; mitgcm-support 2024-November 'Objective Function for Ptracers' (adjoint flags)

### Results changed in 2026: GMRedi Bates/GEOMETRIC/QGLeith fixes (PRs #1024, #1025 from issue #1012)
- Cause: signed cap on Bates wave speed in gmredi_calc_bates_k.F wrong; missing constant factor c_rosY in gmredi_calc_geom.F; wrong grid metric in gmredi_calc_qgleith.F and mom_calc_visc.F (matters on curvilinear/llc grids).
- Fix: update to a checkpoint after Aug 2026; verification output only changed for global_ocean.gm_k3d and gm_res; Leith QG: also config_check stops if viscC2LeithQG != 0 without `ALLOW_LEITH_QG` (MOM_COMMON_OPTIONS.h).
- Era: fixed Aug 2026.
- Src: https://github.com/MITgcm/MITgcm/pull/1024 ; https://github.com/MITgcm/MITgcm/pull/1025 ; https://github.com/MITgcm/MITgcm/issues/1012

## KPP

### "*** ERROR *** Some form of convection has been enabled" (KPP_CHECK)
- Cause: KPP contains its own interior convection; cAdjFreq or ivdc_kappa non-zero.
- Fix: leave cAdjFreq=0 and ivdc_kappa=0 (defaults) with pkg/kpp (not required for GGL90; MY82 same stop). KPP also needs implicitDiffusion/implicitViscosity=.TRUE. and uses diffKr*/viscAr as background.
- Era: 2005-2025 (still in kpp_check.F).
- Src: mitgcm-support 2025-June 'KPPhbl output' http://mailman.mitgcm.org/pipermail/mitgcm-support/2025-June/thread.html ; 2009-November 'KPP parameter' (Ferreira) ; 2005-January 'Convective Adjustment Clarification'

### KPP stops/needs SHORTWAVE_HEATING; "KPP cannot work in the dark"
- Cause: old kpp_check.F stop (2000) because swfrac calls were not inside #ifdef SHORTWAVE_HEATING; zero Qsw itself is fine (daily cycles).
- Fix: fixed upstream with PR #750 (selectPenetratingSW run-time option; kpp_check.F no longer has the stop); older code: `#define SHORTWAVE_HEATING` in CPP_OPTIONS.h and give Qsw separate from Qnet (Qnet net upward incl. SW; use surfQswFile with Qres or exf swdown).
- Era: 2003-2025; fixed 2025.
- Src: https://github.com/MITgcm/MITgcm/issues/913 ; mitgcm-support 2007-March 'KPP & OB file' ; 2003-October 'more KPP blues, maybe related?'

### KPP checkerboard / grid-scale noise in KPPviscAz/KPPdiffKzT, vertical banding at fine dz
- Cause: positive feedback (stratification <-> mixing) amplifying 2dx noise; very small dz (<~1 m) and shear/dbloc thresholds; unmatched boundary-layer vs interior diffusivities.
- Fix: keep `KPP_SMOOTH_SHSQ` + `KPP_SMOOTH_DBLOC` (defaults); try `KPP_SMOOTH_DVSQ` (high lat), `KPP_SMOOTH_DENS`, `KPP_SMOOTH_VISC/DIFF`, `ALLOW_KPP_VERTICALLY_SMOOTH` with large num_v_smooth_Ri/Rm; PR #346 (2020) adds `KPP_DO_NOT_MATCH_DIFFUSIVITIES`, `KPP_DO_NOT_MATCH_DERIVATIVES`; reduce difm0/difs0/dift0 by 10x (shear instability noise); consider GGL90 (no noise in Walberg test case). Sub-1 m dz / non-hydrostatic: KPP not meant for it.
- Era: 2003-2020.
- Src: mitgcm-support 2020-April '[EXTERNAL] Banding, Checkerboarding in KPP Viscosities' http://mailman.mitgcm.org/pipermail/mitgcm-support/2020-April/thread.html ; 2003-October 'more KPP blues, maybe related?' ; 2006-March 'noise in high resolution run' ; https://github.com/MITgcm/MITgcm/pull/346 (cited in thread)

### KPP_ESTIMATE_UREF does not compile; surface-velocity reference
- Cause: dbloc/work1 not passed to kpp_forcing_surf.F (unused since 2007); zref logs with rF<0.
- Fix: PR #340 (Apr 2020) makes it compile (zF=abs(rF)) but "does not mean it does what it is supposed to"; default uses top-cell velocity (resolution dependent, noisy at high res). Option kept for WW3 coupling interest.
- Era: 2007-2020.
- Src: https://github.com/MITgcm/MITgcm/issues/336 ; https://github.com/MITgcm/MITgcm/pull/340 ; mitgcm-support 2020-April 'Banding, Checkerboarding...'

### KPP diffusivity maxima at the bottom (shear mixing), KPPhbl stuck at minKPPhbl or 5 m, KPP in unforced runs
- Cause: neutral stratification + any shear gives large Ri-based mixing (Riinfty=0.7); minKPPhbl = 0.5*drF(1); weakly stratified top layers legitimately give deep hbl; offline runs have no forcing/shear so hbl = minimum.
- Fix: test with `#define EXCLUDE_KPP_SHEAR_MIX` or difm0=difs0=dift0=0.; use diagnostic MXLDEPTH (calc_oce_mxlayer.F, hMixCriteria) as independent MLD; KPPhbl in diagnostics (see vermix/input/data.diagnostics); no bottom boundary layer in MITgcm KPP; KPPdiffKz* diagnostics include background (kappa in KPP convection >1 m2/s is normal); KPP_ghatUseTotalDiffus adds GM-Redi kappa to the non-local term.
- Era: 2003-2025.
- Src: mitgcm-support 2005-November 'KPP diffusivities' ; 2004-January 'I love KPP' ; 2025-June 'KPPhbl output' ; 2018-June 'Student with a question/issue with the KPP package' ; 2019-September 'Diagnosed KPP diffusivities' ; 2015-July 'the calculation of the MXLDEPTH'

### KPP and ice shelves / non-hydrostatic / inverse stratification
- Cause: KPP surface fluxes always applied at k=1 (no forcing under ice shelf, issue #588); config_check warns "Implicit viscosity applies to provisional u,vVel" for nonhydrostatic+KPP; KPP is a mixed-layer scheme (resolving convection explicitly makes it redundant); fresh-water (density maximum) cases are handled via the EOS alpha, not a hard-coded linear EOS.
- Fix: use GGL90 (shelfice OK since PR #597) or none; turn KPP off in non-hydrostatic LES; hardcoded-linear-EOS worry was a comment only.
- Era: 2016-2023 (ice shelf: still open for KPP).
- Src: https://github.com/MITgcm/MITgcm/issues/588 ; mitgcm-support 2019-January 'non-hydrostatic pressure and KPP' ; 2020-January 'KPP with inverse stratification' ; 2016-August 'nonhydrostatic and kpp'

## GGL90, other vertical mixing, convection

### GGL90: results changed / bugs fixed (check your checkpoint)
- Cause/Fix: PR #714 (Jun 2023) matrix coefficients now always scaled by recip_hFacI (CPP `GGL90_MISSING_HFAC_BUG` restores old, undef default) + IDEMIX Ri/Prandtl bug + IDEMIX cvmix version `GGL90_IDEMIX_CVMIX_VERSION`; PR #597 ice-shelf masking; #511 (2021) IDEMIX restart needs EXCH of GGL90TKE/IDEMIX_E; #755 Langmuir (`ALLOW_GGL90_LANGMUIR`) p-coordinate loops wrong; 2026 PR #1015: loop limits in ggl90_calc.F (i loop used jMin/jMax; only matters for non-square tiles under ice shelf) + obcs_apply_r_star.F OBN/OBSeta indices; `ALLOW_GGL90_HORIZDIFF` path had a missing line (PR #690).
- Era: 2021-2026, fixed upstream.
- Src: https://github.com/MITgcm/MITgcm/pull/714 ; https://github.com/MITgcm/MITgcm/issues/755 ; https://github.com/MITgcm/MITgcm/pull/1015 ; https://github.com/MITgcm/MITgcm/issues/511

### GGL90 too much deep mixing (Antarctic) / grid-scale velocity patterns on llc90; mixing-efficiency default
- Cause: mxlMaxFlag choice and settings copied from verification (idemix) experiments; default GGL90ck=0.1 corresponds to gamma=0.286 not 0.2 (c_k = 0.5*gamma*c_eps*Prt, c_eps=0.7).
- Fix: use ECCO v4r4/r5 data.ggl90 as reference (mxlMaxFlag=2); for gamma=0.2 set `GGL90ck=0.07` (issue #169 still open, default unchanged); ggl90 has no MLD diagnostic like KPPhbl (Blanke-Delecluse local mixing length, no l_u/l_d).
- Era: 2018-2021.
- Src: https://github.com/MITgcm/MITgcm/issues/169 ; mitgcm-support 2021-December 'best practice for ggl90' http://mailman.mitgcm.org/pipermail/mitgcm-support/2021-December/thread.html

### Weak/no effect of RiLimit (pp81) / RiMax (my82); MY82 stops; MY82 units
- Cause: those constants are hard-coded, not in namelists; MY82 is unmaintained 2nd-order (not 2.5) scheme; MY82 stops if cAdjFreq/ivdc_kappa set.
- Fix: add RiLimit to pp81_readparms.F / RiMax to my82_readparms.F; pp81 strength via PPnu0; MY82 MYviscAr/MYdiffKr are at W-points (diagnostic code 'L'); possible units bug in my82_calc.F (tke vs tkel) unresolved.
- Era: 2005-2024.
- Src: mitgcm-support 2023-December 'MITgcm: Mixing in Vermix Experiment' ; 2024-February 'Mellor-Yamada in MITgcm' ; 2020-April 'Regional high-res model configuration'

### Convective adjustment: cAdjFreq vs ivdc_kappa; AD bug in convective_adjustment.F
- Cause: three schemes: cAdjFreq (>0 every n s, <0 every step) swaps densities; ivdc_kappa (with cAdjFreq=0, needs implicitDiffusion/implicitViscosity) uses high diffusivity; KPP/GGL90 do it internally. Old: store directives for theta/salt at k-1 wrong in AD (PR #457, 2021).
- Fix: ivdc_kappa: Losch thought 1000 too large (value is a tuning choice); unforced test: salinity restoring with tauSaltClimRelax but no SSS file restores salt to 0 and stops convection; CONVADJ diagnostic shows counter.
- Era: 2005-2021.
- Src: mitgcm-support 2005-January 'Convective Adjustment Clarification' ; 2003-October 'convection problem' ; https://github.com/MITgcm/MITgcm/pull/457

### KL10 spurious huge viscosity with non-linear EOS
- Cause: FIND_RHO_SCALAR called with totPhiHyd not pressure (kl10_calc.F ~l.93), in-situ density overturn.
- Fix: issue #647 (2022): proposed use SigmaR (iso-neutral gradient) to build density profile; check status before use (linear EOS ok).
- Era: 2022-2023.
- Src: https://github.com/MITgcm/MITgcm/issues/647 ; mitgcm-support 2022-August 'Help with internal wave breaking and shoaling parameterization package'

### "CONFIG_CHECK: diffKrFile is set but never used" / NaN with a 3-D diffKr file / changing diffKrT does nothing
- Cause: `ALLOW_3D_DIFFKR` not defined in the *build-directory* CPP_OPTIONS.h (stale build); large local K violates dt < dz^2/K (esp. partial cells); with diffKrFile present total K comes from file (+KPP/GM) so diffKrT/S edits are masked.
- Fix: `#define ALLOW_3D_DIFFKR`, clean build dir, genmake2/make depend/make; check CFL; file applies to T, S and (by default) ptracers; `diffKrFile` is read in ini_mixing.F; no 3-D viscAr (only viscArNr profile; modify calc_viscosity.F using diffKr); no 3-D horizontal diffusivity (modify gad_calc_rhs.F/gad_diff_x,y.F); ALLOW_3D_VISCAH/VISCA4 give 3-D horizontal viscosity.
- Era: 2007-2025.
- Src: mitgcm-support 2025-February 'Specifying a space dependent vertical diffusivity coefficient for temperature' http://mailman.mitgcm.org/pipermail/mitgcm-support/2025-February/thread.html ; 2017-October 'Adding 3D diffusivity field to model input' ; 2019-December 'how to change the vertical mixing strength' ; 2018-March 'Prescribe profile of vertical eddy viscosity' ; 2010-January 'how to specify the vertical and horizontal diffusivity'

## Horizontal viscosity, noise, stability

### Choosing viscAh / viscA4 and viscAhGrid; "SOLUTION IS HEADING OUT OF BOUNDS" with viscosity CFL
- Cause: stability limited by grid Re/Munk layer and viscosity CFL; viscAh and viscAhGrid are additive.
- Fix: use non-dimensional `viscAhGrid` (<1; viscAh = 0.25*L^2*viscAhGrid/deltaT, so viscAhGrid=4*deltaT*viscAh/L^2) and/or `viscA4Grid` (viscA4 ~ 0.25*0.125*L^4*viscA4Grid/deltaT; ~0.01), cap with viscAhGridMax<=1 (Leith/Smag too); Munk layer (Ah/beta)^(1/3) >= 2-3 grid cells; too-small viscAhGrid -> grid noise, smaller allowable dt, blow-up; dz(k+1)/dz(k) <= 1.4; no_slip_bottom plus bottomDrag* double-counts drag; header of pkg/mom_common/mom_calc_visc.F lists recommended values (viscC2Leith/D 1-3, viscC4Leith 1-3, viscC4LeithD 1.5-3, viscC2smag 2.2-4 or 0.2-0.9, ...).
- Era: 2004-2025.
- Src: mitgcm-support 2020-October 'Setting up viscosities for the model' ; 2023-June 'Query regarding horizontal strips in currents' ; 2014-February 'Viscosity parameters' ; 2007-May 'ploblem with a variable-grid global ocean circulation model' ; 2009-July 'Minimum Viscosities'

### Striping/noise along basin boundaries on spherical grids (tutorial_baroclinic_gyre noisy)
- Cause: dx shrinks as cos(lat) while dy fixed -> non-square cells; unresolved gravity waves from boundaries; Coriolis scheme.
- Fix: scale delY by cos(lat) (isotropic cells) + large viscA4Grid; `selectCoriScheme=3` (Jamart-Ozer, JMC) reduces noise to western boundary; drop CD scheme (useCDscheme false) and use small biharmonic instead; strip forcing/IC to bisect; check forcing array orientation (transposed wind field).
- Era: 2025.
- Src: mitgcm-support 2025-April 'Irregular Velocity Fields Near Basin Boundaries in MITgcm' http://mailman.mitgcm.org/pipermail/mitgcm-support/2025-April/thread.html

### Biharmonic viscosity: ISOTROPIC_COS_SCALING/cosPower, L4rdt change, noise at 1/4-1/6 deg
- Cause: ISOTROPIC_COS_SCALING (retired; viscAhGrid/viscA4Grid supersede) did not apply to momentum if set in GAD_OPTIONS.h only, and COSINEMETH_III conflicts with it; 2009 L4rdt fix lowered effective A4 by ~30-37% so old viscA4Grid/viscA4GridMax values needed raising.
- Fix: use viscA4Grid/viscA4GridMax (CS510 values: viscC4Leith=1.5, viscC4Leithd=1.5, viscA4GridMax=0.5; useFullLeith etc.); at <=1 km hydrostatic gets grid-cell storms -> nonhydrostatic.
- Era: 2003-2009 (cos scaling), still relevant for comparing old runs.
- Src: mitgcm-support 2006-April 'noise in high resolution run' ; 2003-December '[Fwd: RE: grid-scale noise in 1/4-deg model]' ; 2009-August 'L4rdt' ; 2018-April 'cosPower with GMRedi and Leith viscosity'

### Leith / Smagorinsky usage pitfalls
- Cause: Leith/Smag only act on horizontal momentum viscosity (not tracers, not vertical viscosity); Leith for 2-D/QG turbulence, Smag for 3-D (<= deformation radius scales); Leith+Smag together over-dissipate (Losch advice); Leith originally needed vectorInvariantMomentum (old behaviour).
- Fix: Leith: viscC2Leith (+viscC2LeithD to catch divergence) or biharmonic viscC4Leith/D, with viscAhGridMax/viscA4GridMax; QG Leith: `ALLOW_LEITH_QG` + `viscC2LeithQG` (config stop otherwise, PR #320); `useSmag3D` + `smag3D_coeff` (tutorial_deep_convection/input.smag3d; today also `smag3D_diffCoeff` for tracer diffusivity); known issues: 3-D Smag with no-slip BCs, Um_Diss/Vm_Diss may not include Smag3D (2014), no Germano dynamic coefficient; use Leith D term at <deformation radius; offline Kh must be uniform.
- Era: 2003-2020.
- Src: mitgcm-support 2005-July 'viscAh & viscAz!' ; 2018-April 'Choosing Leith biharmonic co-efficient?' ; 2019-May 'vertical smagorinsky viscosity' ; 2014-December 'Smagorinsky 3D tendency diagnostics' ; 2015-March 'Smagorinsky viscosity' ; https://github.com/MITgcm/MITgcm/pull/320

### Dissipation diagnostics: Um_Diss etc. don't balance, units m^4/s^2
- Cause: Um_Diss/Vm_Diss (guDissip) hold only explicit terms; VISrE_Um/VISrI_Um and DFrE_*/DFrI_* are fluxes including cell area.
- Fix: add implicit part from VISrI_*; tendency = (flux[k]-flux[k+1])/(rA*drF*hFac); viscAh/viscAr are kinematic (m^2/s).
- Era: 2019.
- Src: mitgcm-support 2019-October 'dissipation and diffusion rates' ; 2019-September 'MITgcm viscocity'

### Time-dependent or state-dependent viscosity/diffusivity (custom code)
- Cause: mom_calc_visc.F lacks myTime/myIter in its argument list; 3-D vertical coefficients are computed per step in calc_viscosity.F / calc_3d_diffusivity.F.
- Fix: extend the formal parameter list (mom_fluxform.F, mom_vecinv.F, mom_calc_visc.F); box tests with XC/YC only inside 1<=i<=sNx,1<=j<=sNy; start vertical schemes from pp81.
- Era: 2022-2025.
- Src: mitgcm-support 2025-June 'use the variable myIter' ; 2022-March 'How to make the viscosity and diffusivity change with temperature and pressure' ; 2021-December 'How to add a new viscosity scheme in MITgcm?'

### Spurious mixing / variance destruction with flux-limited schemes; choosing advection for tracers
- Cause: any scheme either conserves variance (centered, noisy) or extrema (limited, diffusive); 33/77 give "spurious" variance sink even for non-divergent flow.
- Fix: 7 (OS7MP) or 80/81 Prather (extra pickup_somTRAC*, Prather moments, ignore when editing pickups); do external variance/energy budgets and treat residual as numerical; synchronous staggerTimeStep with 33.
- Era: 2021.
- Src: mitgcm-support 2021-January 'Spurious mixing with internal tides' ; 2021-February '[EXTERNAL] Spurious mixing with internal tides' ; 2017-May 'pickup file for ptracer'

### Mixed layer / vertical-mixing background choices (what to set under KPP)
- Cause: background diffKr/viscAr still matter under KPP (acts as interior wave mixing); numerical diffusion of the advection scheme often larger than 1e-5.
- Fix: diffKrT/S 1e-5 to 1e-4 global (higher res needs less; ~1.5e-5 improved Equatorial thermocline for ECCO-like), viscAr 1e-3-1e-4; Bryan-Lewis (diffKrBL79*) is one profile; for regional local tuning edit vddiff in kpp_calc.F by xC/yC; sea-level/ice melt etc. outside scope.
- Era: 2006-2020.
- Src: mitgcm-support 2017-January 'On the vertical mixing parameterization' ; 2009-March 'upper ocean temperature drift' ; 2006-June 'reducing vertical mixing' ; 2010-April 'KPP scheme and background viscosities'
