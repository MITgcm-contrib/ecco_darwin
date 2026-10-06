# Troubleshooting: pkg/diagnostics, output formats (MDS/MNC/NetCDF), budget closure
Scope: mitgcm-support (2003-2026) + MITgcm GitHub issues/PRs on data.diagnostics, diagnostic names, output timing, MDS/MNC output, rdmds/MITgcmutils, and closing heat/salt/tracer/momentum budgets. Names were checked against origin/master (3 Oct 2026) and the mitgcm-ecco index; "Era" says when a name or bug has since changed. Mailing-list thread URLs: http://mailman.mitgcm.org/pipermail/mitgcm-support/YYYY-Month/thread.html

## Budget closure: tracers (heat / salt / passive)

### Heat/salt budget does not close at k=1 with linear free surface (nonlinFreeSurf=0, no r*)
- Cause: d(eta)/dt is inconsistent with the fixed surface-cell thickness, so a small surface tracer source/sink (WTHMASS/WSLTMASS at k=1) exists. ADVr_TH / ADVr_SLT at k=1 are zero, so they are NOT the right term.
- Fix: (a) linFSConserveTr=.TRUE. in `data` corrects it globally (calc_wsurf_tr.F; uniform TsurfCor/SsurfCor applied at k=1). (b) With the default linFSConserveTr=.FALSE., add the surface term yourself from diagnostics WTHMASS / WSLTMASS level 1: tend_k1 = -( WTHMASS(k=1) - TsurfCor )/(drF(1)*hFacC(k=1)), TsurfCor = SUM(WTHMASS(:,:,1)*RAC)/globalArea. Offline k=1 divergence (Menemenlis/JMC): ADVx_TH(i)-ADVx_TH(i+1) + ADVy_TH(j)-ADVy_TH(j+1) + ADVr_TH(k=2) - WTHMASS(k=1)*RAC. First close k>1 (no surface problem), then revisit k=1.
- Era: 2011-2017, still valid. Ref doc: doc/old_doc/Heat_Salt_Budget_MITgcm.pdf (Chakraborty & Campin) and doc/old_doc/diags_changes.txt.
- Src: mitgcm-support 2011-June 'Heat Budget in tutorial' (JMC/Menemenlis, 2011-May-20 msgs) http://mailman.mitgcm.org/pipermail/mitgcm-support/2011-June/thread.html ; 2014-April 'T, S budget in the surafce layer' http://mailman.mitgcm.org/pipermail/mitgcm-support/2014-April/thread.html ; 2014-November 'Heat budget in MITGCM' http://mailman.mitgcm.org/pipermail/mitgcm-support/2014-November/thread.html

### TOTTTEND / TOTSTEND do not close the budget with nonlinear free surface / r* (ECCO v4 style), or any run where hFac varies
- Cause: TOTTTEND/TOTSTEND = d(theta)/dt (finite difference of theta across the time step, units degC/day, so /86400), not d(h*theta)/dt. With NLFS the volume change term S*d(h)/dt is missing. Time-averaging does not fix it because the weight h is time varying.
- Fix: for NLFS/r* write snapshots (frequency<0, timePhase=0) of THETA/SALT and ETAN (and hFacC/etc. from grid output) and form d(h*S)/dt offline from the snapshots, then compare with the divergence of ADVx/y/r_SLT, DFx/yE_SLT, DFrE/I_SLT (+KPPg_SLT) and surface fluxes (SFLUX, oceFWflx, surForcS) accounting for the time-varying cell volume. Use the Piecuch memo for ECCO v4 (https://dspace.mit.edu/handle/1721.1/111094). Proposed (not merged) new "h*S" tendency diagnostics: issue #1010 and PR #969; Losch suggests simply saving snapshots of h*S.
- Era: 2016-2026. Still open in master (TOTTTEND still filled from diagnostics_fill_state.F as -theta*86400/dTtracerLev, before/after the step; no NLFS correction). Related: PR #720 (layers) notes the same limitation.
- Src: https://github.com/MITgcm/MITgcm/issues/1010 ; mitgcm-support 2016-December 'heat and salt budget with nonlinear free surface' http://mailman.mitgcm.org/pipermail/mitgcm-support/2016-December/thread.html ; 2017-May 'Salinity budget tendencies' (Buckley) http://mailman.mitgcm.org/pipermail/mitgcm-support/2017-May/thread.html ; 2016-May 'Question about closing heat budget' http://mailman.mitgcm.org/pipermail/mitgcm-support/2016-May/thread.html

### Which surface-flux diagnostic closes the surface heat budget (TFLUX vs oceQnet vs surForcT)
- Cause: several overlapping diagnostics; sea ice modifies the ocean fluxes so EXF fields differ from what the ocean feels.
- Fix: TFLUX is the budget-closing one (Menemenlis): if useRealFreshWaterFlux=F or (nonlinFreeSurf=0 & z-coords): TFLUX = oceQnet + TRELAX + oceFreez; if useRealFreshWaterFlux=T and (nonlinFreeSurf>0 or p-coords): TFLUX = oceQnet + TRELAX + PmEpR*temp_EvPrRn*Cp + oceFreez. TFLUX contains the vertically integrated penetrating shortwave; to budget level-by-level remove oceQsw from TFLUX and add back the swfrac-weighted SW at each level. Units: surface flux W/m2 -> degC/s: divide by HeatCapacity_Cp*rUnit2mass (=rhoConst in z coords). To relate to EXF go EXF -> oceQnet/oceFWflx/oceTAUX/oceTAUY -> surForcT/surForcS (model variables surfaceForcingT/S) as an intermediate step.
- Era: 2011-2021, still valid.
- Src: mitgcm-support 2011-June (above); 2021-September 'Understanding gT_Forc and gS_Forc diagnostics' http://mailman.mitgcm.org/pipermail/mitgcm-support/2021-September/thread.html ; 2014-November (Ryan Abernathey msg on TFLUX units)

### Shortwave heating missing in the vertical heat balance in the top ~100 m
- Cause: SW is distributed over depth by swfrac (model/src/swfrac.F; hard-coded Jerlov type IA; only active if CPP SHORTWAVE_HEATING defined in CPP_OPTIONS.h, which is #undef by default).
- Fix: recompute offline from oceQsw: source_k = oceQsw/(rhoConst*Cp) * (swfrac(rF(k)) - swfrac(rF(k+1)))/(drF(k)*hFacC). Dimitris: shortwave is exponentially decaying in the top ~200 m and is modified within the KPP mixing layer. A diagnostic net heat flux divergence = (flux at z=-30 m) already includes the surface cooling; look at divergence, not flux profiles (Losch).
- Era: 2009-2020. Name checks: swfrac.F, SHORTWAVE_HEATING exist today.
- Src: mitgcm-support 2020-August 'Vertical heat balance in the ocean' http://mailman.mitgcm.org/pipermail/mitgcm-support/2020-August/thread.html ; 2009-March 'Absorption of short wave radiation' http://mailman.mitgcm.org/pipermail/mitgcm-support/2009-March/thread.html ; 2014-February 'help with turbidity' http://mailman.mitgcm.org/pipermail/mitgcm-support/2014-February/thread.html

### Subsurface budget residual: ADV*/DF* are fluxes (degC m3/s), plus implicit, KPP and GM terms
- Cause: ADVx/y/r_TH, DFx/y/rE_TH, DFrI_TH are fluxes through cell faces already multiplied by face area; people forget the divergence, cell volume, implicit part or non-local KPP flux. Adams-Bashforth is also in the model tendency.
- Fix: tend(k>1) = -[ (ADVx(i+1)-ADVx(i)) + (ADVy(j+1)-ADVy(j)) + (ADVr(k+1)-ADVr(k)) + same for DFxE, DFyE, DFrE, DFrI, KPPg_TH ] / (rA*drF*hFacC); with Linear FS use fluxes at the same time level and TOTTTEND/86400 for the LHS. KPPg_TH (non-local, ~ KPPghat) is a vertical flux treated like DFrI_TH (zero at k=1; KPPghatK in older code). Do not use centred differences (i+1)-(i-1). With implicitDiffusion=T you must include DFrI_TH. When using GM/Redi see next entry. Always use ADVx_TH (not UVELTH, UTHMASS) for budgets: ADVx_TH includes flux-limiter/numerical diffusion; UVELTH is a plain U*T correlation; UTHMASS is hFac-weighted U*T (Fox-Kemper).
- Era: 2007-2016; names ADVx/y/r_TH, DFxE_TH, DFrE_TH, DFrI_TH, KPPg_TH exist today.
- Src: mitgcm-support 2014-November 'Heat budget in MITGCM' (JMC) http://mailman.mitgcm.org/pipermail/mitgcm-support/2014-November/thread.html ; 2007-February 'advective fluxes and transports' http://mailman.mitgcm.org/pipermail/mitgcm-support/2007-February/thread.html ; 2011-July 'Diagnostics & temperature transport equation' http://mailman.mitgcm.org/pipermail/mitgcm-support/2011-July/thread.html ; 2010-April 'KPP and Heat Budget' http://mailman.mitgcm.org/pipermail/mitgcm-support/2010-April/thread.html

### DFrE_TH non-zero although implicitDiffusion=.TRUE. (or 'DFrE_TH has not been filled')
- Cause: with pkg/gmredi, the Kwx and Kwy (Redi + skew-flux GM) contributions to the vertical flux are always explicit (DFrE_*); only the Kwz part goes into DFrI_* when implicitDiffusion=T. Without gmredi and with implicitDiffusion=T, DFrE_TH is never filled -> "has not been filled (ndiag= 0)" warning and zeros.
- Fix: include both DFrE_* and DFrI_* in the budget when useGMRedi=T. Drop DFrE_* only for runs without GM/Redi.
- Era: 2007-2017 (JMC explanation 2017-Jan), still valid.
- Src: mitgcm-support 2017-January thread of 2007-March 'advective and diffusive diagnostics' http://mailman.mitgcm.org/pipermail/mitgcm-support/2007-March/thread.html ; 2007-July 'diagnostics DFrE_TH' http://mailman.mitgcm.org/pipermail/mitgcm-support/2007-July/thread.html

### WVELTH vs WTHMASS differ when staggerTimeStep=.TRUE.
- Cause: WVELTH is filled before the momentum solve, WTHMASS after it (with updated wVel) when staggerTimeStep=T; identical when staggerTimeStep=F. Neither has hFac or area factor (m degC/s).
- Fix: use WTHMASS/WSLTMASS for tracer budgets (consistent with tracer advection), WVELTH for eddy-flux statistics.
- Era: 2018-Mar; both names exist today (diagnostics_fill_state.F).
- Src: mitgcm-support 2018-March 'WTHMASS vs WVELTH' http://mailman.mitgcm.org/pipermail/mitgcm-support/2018-March/thread.html

### ADVx_TH / ADVy_TH / ADVr_TH off by a factor of the time step with SOM advection (scheme 80/81)
- Cause: gad_som_adv_* returned fluxes multiplied by the tracer time step.
- Fix: update code (JMC changed SOM advection S/R 8 Jan 2008); old code: divide by deltaTtracer. Related: layers budget with scheme 80/81 did not close (issue #60), fixed by PR #988 (LinFSConserveTr=F).
- Era: fixed upstream 2008-Jan; layers fixed 2026.
- Src: mitgcm-support 2008-January 'Strange results from GAD diagnostics using SOM' http://mailman.mitgcm.org/pipermail/mitgcm-support/2008-January/thread.html ; https://github.com/MITgcm/MITgcm/issues/60

### ADVr_TH / DFrE_TH are identically zero at k=1 but WVELTH is not
- Cause: these are fluxes through the upper cell interface ('L' 9th diagCode character); nothing is fluxed through the sea surface. W*THETA at k=1 is just transport by w = d(eta)/dt, not a flux across the surface.
- Fix: no fix needed; see the k=1 entry for the surface correction term.
- Era: 2011-Jun; still valid. 'L'/'U' 9th-character semantics are misleading in 2-D diagnostics (JMC, 2021-Jan: 'L' inherited from p-coordinate atmosphere; k-1/2 is physically upper in ocean).
- Src: mitgcm-support 2011-June 'diagnostics : ADVr_TH,DFrE_TH and WVELTH' http://mailman.mitgcm.org/pipermail/mitgcm-support/2011-June/thread.html ; 2021-January 'Diagnostics package codes' http://mailman.mitgcm.org/pipermail/mitgcm-support/2021-January/thread.html

### Divergence of u not zero / volume budget with nonlinear free surface
- Cause: with finite volumes div(u)=0 exactly only with correct hFac weighting (use UVELMASS/VVELMASS/WVEL, hFacW/S/C, dxG/dyG/rA/drF). With NLFS (no r*) the surface-cell thickness varies in time; with r* div(u) is not zero by design (it balances d(eta)/dt scaling).
- Fix: compute divergence with face areas (flu=u*hFacW*dyG*drF etc.); for NLFS account for time-varying surface thickness; for r* use the r* volume budget (JMC: see Adcroft & Campin 2004 Ocean Modelling 7).
- Era: 2016-Apr, valid.
- Src: mitgcm-support 2016-April 'calculation of divu' http://mailman.mitgcm.org/pipermail/mitgcm-support/2016-April/thread.html

### DIC / passive-tracer budget does not close (ptracers, DIC, darwin)
- Cause: free-surface dilution at k=1 (same issue as for heat/salt); extra DIC terms not diagnosed; OBCS and negative-value resets in DIC/BLING have no diagnostics; calendar-dump output not applied to ptracer dumps.
- Fix: (1) PTRACERS_linFSConserve(n)=.TRUE. (data.ptracers) and/or linFSConserveTr=.TRUE.; darwin has darwin_linFSConserve (data.darwin). Core advice: best conservation with r* (nonlinFreeSurf=4, select_rStar=2 plus CPP NONLIN_FRSURF); exactConserv=.TRUE. may also help. (2) Closing terms: Tp_gTr?? (total transport tendency), ForcTr??, AB_gTr??, ADV*Tr??, DF*Tr??; DIC bio term = Rcp*(-DICBIOA + DICPFLUX + DICRDOP) (Rcp=117) plus DICCARB; air-sea DICTFLX (compare DICCFLX). DICPFLUX is not a diagnostic in the stock code (add one). See Lauderdale 2016 GBC notebook. (3) OBCS contributions to ptracers are not diagnosed; ForcTr includes RBCS.
- Era: 2024-Mar/2025-Jul. Names PTRACERS_linFSConserve, linFSConserveTr, select_rStar, exactConserv exist today.
- Src: mitgcm-support 2024-March 'Closing DIC budget' http://mailman.mitgcm.org/pipermail/mitgcm-support/2024-March/thread.html ; 2025-July 'ptracer/DIC budget closure/calendar pkg' http://mailman.mitgcm.org/pipermail/mitgcm-support/2025-July/thread.html

### Restoring (RBCS) / geothermal / SW forcing not visible in budget
- Cause: RBCS had no diagnostics; surface forcing diagnostics are 2-D.
- Fix: 3-D forcing tendencies gT_Forc and gS_Forc (and ForcTr?? for ptracers) exist since checkpoint65a; they include RBCS relaxation plus at k=1 the surface forcing (subtract surForcT/S to isolate RBCS). Um_Ext/Vm_Ext already include the RBCS momentum relaxation.
- Era: since c65a (Jul 2014); names gT_Forc, surForcT exist today.
- Src: mitgcm-support 2014-August 'rbcs diagnostics' http://mailman.mitgcm.org/pipermail/mitgcm-support/2014-August/thread.html

### ocean-to-sea-ice heat flux: which SI* diagnostic
- Cause: oceQnet (positive down) and SIqnet (positive up) are the same flux (SIqnet = -oceQnet in practice); neither separates the ice part.
- Fix: Losch: heat absorbed by the ice-ocean system = SIatmQnt - SIqnet; SIqneti (atmosphere through ice cover) plus SIaQbOCN (ocean heat that melts ice, 'a_QbyOCN' in seaice_growth.F) are the ice-related terms; for details read seaice_growth.F.
- Era: 2024-Aug; names SIatmQnt, SIaQbOCN exist today.
- Src: mitgcm-support 2024-August 'Diagnostic for Heat flux from ocean to sea ice' http://mailman.mitgcm.org/pipermail/mitgcm-support/2024-August/thread.html

## Budget closure: momentum / energy

### Momentum budget residual (flux-form or vector-invariant)
- Cause: missing terms: surface-pressure gradient, implicit vertical viscosity, AB term, atm/sea-ice loading, NH pressure; Um_Advec already contains Coriolis; units of TOTUTEND are per day.
- Fix (JMC recipe, closes to machine precision over one time step with writeBinaryPrec=64, valid for both formulations): TOTUTEND/86400 = Um_Advec + Um_dPhiX + Um_Diss + Um_ImplD + Um_Ext + AB_gU - g*(ETAN(i)-ETAN(i-1))/dxC [or use PHI_SURF diagnostic] ; also subtract (atmP_load(i)-atmP_load(i-1))/(rhoConst*dxC) and g*(sIceLoad(i)-sIceLoad(i-1))/(rhoConst*dxC) (use rhoConst, not rhoConstFresh; ECCO: SSH=ETAN+sIceLoad/rhoConstFresh in the cited post). Um_Cori only if CD-scheme. Um_ImplD (post-2019) replaces the manual -d(VISrI_Um/(rAw drF hFacW))/dz term. Over long averaging periods AB error is small.
- Era: 2010-2024. Um_dPHdx was renamed Um_dPhiX in PR #219 (Aug 2019; diagnostics_utils.F maps old names). AB_gU/AB_gV were unfilled due to a typo before timestep.F rev 1.55 (Dec 2013). The 5 groups: TOTUTEND/86400 = Um_Advec + Um_dPhiX + Um_Diss + Um_ImplD + Um_Ext + AB_gU (#820 JMC). Finer split: Um_Advec = Um_Cori + Um_AdvZ3 + Um_AdvRe + KE gradient + (Um_Cori3, Um_Metr) (Ferreira 2025; Bernoulli KE term has no diagnostic; compute -d(momKE)/dx offline).
- Src: mitgcm-support 2010-December 'Momentum Budget' http://mailman.mitgcm.org/pipermail/mitgcm-support/2010-December/thread.html ; 2013-December digest threads http://mailman.mitgcm.org/pipermail/mitgcm-support/2013-December/thread.html ; 2014-February 'Momentum Budget' http://mailman.mitgcm.org/pipermail/mitgcm-support/2014-February/thread.html ; https://github.com/MITgcm/MITgcm/issues/820 ; 2025-April 'Diagnostic Um_Advec' http://mailman.mitgcm.org/pipermail/mitgcm-support/2025-April/thread.html

### Momentum budget off by ~10% in non-hydrostatic runs
- Cause: with exactConserv=T, etaN is re-computed after the solve, so g*ETAN != surface pressure used; filters (Shapiro/FFT) after solve_for_pressure; surface and NH pressure are not cleanly separated (zeroPsNH=F hard-coded).
- Fix: use PHI_SURF + PHI_NH diagnostics (not g*ETAN) for the pressure gradient, plus Um_dPhiX; tighten the solver tolerance; turn off GGL90/implicit viscosity warning case when diagnosing. There is no Wm_dPhiZ / TOTWTEND diagnostic for the w equation (Delorme 2020).
- Era: 2017-Aug JMC; 2020-Nov.
- Src: mitgcm-support 2017-August 'non-hydrostatic momentum budget' http://mailman.mitgcm.org/pipermail/mitgcm-support/2017-August/thread.html ; 2020-November 'Non-Hydrostatic Energy Budget' http://mailman.mitgcm.org/pipermail/mitgcm-support/2020-November/thread.html

### Um_Diss misses bottom drag / Smag-3D, UBotDrag diagnostic absent
- Cause: Um_Diss contains explicit dissipation including bottom/side drag only for selectImplicitDrag=0; with selectImplicitDrag>0 drag is in the implicit part (VISrI_Um). Smagorinsky-3D contribution was missing from Um_Diss (JMC agreed to add, 2014). UBotDrag/VBotDrag were replaced (PR #219) by botTauX/botTauY; PR #817/#820 added Um_Diss2/4, UBotDrag, UShIDrag etc. only under CPP ALLOW_MOM_TEND_EXTRA_DIAGS (MOM_COMMON_OPTIONS.h, default #undef).
- Fix: check selectImplicitDrag in STDOUT; define ALLOW_MOM_TEND_EXTRA_DIAGS if you really need UBotDrag (not with selectImplicitDrag=2), else use botTauX/botTauY.
- Era: 2018-2024. botTauX/Y fix for -ur4 in PR #276.
- Src: mitgcm-support 2018-September 'Dose Um_Diss involve UBotDrag term?' http://mailman.mitgcm.org/pipermail/mitgcm-support/2018-September/thread.html ; https://github.com/MITgcm/MITgcm/issues/820 ; mitgcm-support 2014-December 'Smagorinsky 3D tendency diagnostics' http://mailman.mitgcm.org/pipermail/mitgcm-support/2014-December/thread.html

### Kinetic-energy budget / momKE definition
- Cause: momKE (mom_calc_ke.F, KEscheme=2 in flux form) is 0.5*(avg(U^2*hFacW)+avg(V^2*hFacS))/hFacC and excludes w even in NH runs; other KEscheme options exist only with vectorInvariantMomentum. monitor's KE (mon_ke.F) is the better reference. No standard KE budget diagnostics exist; local cancellation of Coriolis terms is not expected on the C-grid.
- Fix: close the momentum budget first, multiply by U/V (Adams-Bashforth matters), or use Klymak's online energy code (github.com/jklymak/MITgcmcode) / bderembl energydiag fork.
- Era: 2014-2019, still true.
- Src: mitgcm-support 2015-January 'KE diagnostics' http://mailman.mitgcm.org/pipermail/mitgcm-support/2015-January/thread.html ; 2015-October 'Kinetics energy budget' http://mailman.mitgcm.org/pipermail/mitgcm-support/2015-October/thread.html ; https://github.com/MITgcm/MITgcm/issues/214

## Output timing: frequency, timePhase, calendarDumps

### Snapshots come out at the middle of the interval (e.g. t=200,600,1000 instead of 0,400,800)
- Cause: default phase for frequency<0 (snapshots) is |frequency|/2 (averages default 0). Snapshot file iteration suffix is myIter-1 (consistent with state variables, not with time-average convention).
- Fix: timePhase(n) = 0. in data.diagnostics (also use for central-day averages: timePhase=-43200 for noon-centred daily averages when startTime is at noon).
- Era: 2005-2023, still valid (documented in verification data.diagnostics).
- Src: mitgcm-support 2005-September 'odd diagnostics package output' http://mailman.mitgcm.org/pipermail/mitgcm-support/2005-September/thread.html ; 2009-February 'snap-shot outputs from diagnostics pkg' http://mailman.mitgcm.org/pipermail/mitgcm-support/2009-February/thread.html ; 2011-March 'what's the difference between U, V, W, salt in data.diagnostic' http://mailman.mitgcm.org/pipermail/mitgcm-support/2011-March/thread.html ; 2017-June 'averaged output frequency' http://mailman.mitgcm.org/pipermail/mitgcm-support/2017-June/thread.html

### Output at unexpected times / extra files with frequency that is not a multiple of deltaT
- Cause: output only when a multiple of freq lies within +/- step/2; if no time step lands on it the nearest step is used, and rounding can switch between .4 and .6 steps. Also, when frequency < deltaT, output is every step and timePhase is ignored.
- Fix: choose frequency as an integer multiple of deltaTclock (and timePhase accordingly); frequency >= deltaT if you want timePhase honored (JMC 2013-Sep).
- Era: 2009-2023, valid.
- Src: mitgcm-support 2009-February 'mnc output problem-multiple time step' http://mailman.mitgcm.org/pipermail/mitgcm-support/2009-February/thread.html ; 2023-May 'Strange output frequency' http://mailman.mitgcm.org/pipermail/mitgcm-support/2023-May/thread.html ; 2013-September 'diagnostics parameter' http://mailman.mitgcm.org/pipermail/mitgcm-support/2013-September/thread.html

### Extra outputs at iter-1 / iter+1 around the right one after long runs (e.g. U.0016783199, ...200, ...201), or 3 files per dump
- Cause: round-off in eesupp DIFFERENT_MULTIPLE (very large myTime vs step); also reported to be caused by sick compute nodes (kipmi0 / power supply) giving odd output times.
- Fix: put eesupp/src/different_multiple.F (and diff_phase_multiple.F) in NOOPTFILES with NOOPTFLAGS=-O0 (or -Kieee / -fp-model precise); check with a different node/platform.
- Era: 2006-2016. DIFFERENT_MULTIPLE exists today. Overflow of integer in diagnostics_out.F timeRec for huge myTime: issue #426 (use startTime-relative timeRec; PR #429).
- Src: mitgcm-support 2016-July 'Uncanny output frequency after some time of simulation' http://mailman.mitgcm.org/pipermail/mitgcm-support/2016-July/thread.html ; 2006-January 'mnc output' http://mailman.mitgcm.org/pipermail/mitgcm-support/2006-January/thread.html ; 2015-November 'Multiple outputfiles' http://mailman.mitgcm.org/pipermail/mitgcm-support/2015-November/thread.html ; https://github.com/MITgcm/MITgcm/issues/426

### Start diagnostics output only after spin-up / for a limited window
- Cause: dumpFreq/taveFreq cannot do it; pkg/diagnostics can via timePhase (write at timePhase + k*|frequency|).
- Fix: timePhase(n)=<seconds> in data.diagnostics (and stat_phase for DIAG_STATIS_PARMS). Caveat (JMC): the first time-average record (frequency>0) is an average from startTime to the first write, not just one period. Different dumpFreq in sequence: use several streams (cannot stop a stream); extra files must be removed. Needs frequency>=deltaT.
- Era: 2013-2024, valid. Does not apply with calendarDumps for timePhase > 1 year until PR #1013 (see below).
- Src: mitgcm-support 2016-September 'timephase() in data.diagnostics' http://mailman.mitgcm.org/pipermail/mitgcm-support/2016-September/thread.html ; 2010-August 'diagnostics - question' http://mailman.mitgcm.org/pipermail/mitgcm-support/2010-August/thread.html ; 2022-August 'different dumpFreqs during a single run' http://mailman.mitgcm.org/pipermail/mitgcm-support/2022-August/thread.html ; 2024-September 'High frequency velocity outputs for a specific period' http://mailman.mitgcm.org/pipermail/mitgcm-support/2024-September/thread.html

### Snapshot of the last time step missing; averaged output at end of run
- Cause: diagnostics are filled at the beginning of the time step, so a snapshot at the last iteration is not available; dumpAtLast (data.diagnostics) only applies to time-averages (frequency>0).
- Fix: run one extra iteration (nTimeSteps+1) and set timePhase to the required multiple. For main-model dumps dumpInitAndLast=.FALSE./dumpFreq=0 turns state files off.
- Era: 2007-2015, valid (DIAGNOSTICS.h: dumpAtLast = 'always write time-ave (freq>0) diagnostics at end of the run').
- Src: mitgcm-support 2007-March 'diagnostics at last timestep' http://mailman.mitgcm.org/pipermail/mitgcm-support/2007-March/thread.html ; 2015-May 'Diagnostics at end of run' http://mailman.mitgcm.org/pipermail/mitgcm-support/2015-May/thread.html ; 2008-August 'No diagnostics output on last timestep' http://mailman.mitgcm.org/pipermail/mitgcm-support/2008-August/thread.html

### Calendar-month / calendar-year averages (calendarDumps) misbehave: half-day offset, snapshot ignored, wrong timeInterval
- Cause: pkg/cal CAL_TIME2DUMP converts frequencies of ~28-31 d / ~360-372 d (posFreq 2592000-2678400 or 31104000-31968000 s) to calendar months/years. Details: (1) only affects dumps (diagnostics freq, chkPt, taveFreq...), not forcing records (monthly forcing needs period=(365*3+366)/4/12 d); (2) month stamps shifted by half a day if startDate_2 is not 000000; (3) snapshots (freq<0) ignored calendarDumps in 2012 (JMC then: 'will try to change'; diagnostics_write.F now calls CAL_TIME2DUMP for both) and default to mid-interval unless timePhase=0; (4) the timeInterval in .meta does not account for calendar rounding (meta file bug, not data); (5) dumpFreq/PTRACERS_dumpFreq and other pkg dumps do not use CAL_TIME2DUMP (only diagnostics, pickups, KPP/seaice/shelfice tave).
- Fix: calendarDumps=.TRUE. in data.cal with a nominal freq (e.g. 2635200 or 2592000), check time stamps not .meta timeInterval; snapshots need timePhase=0; for end-of-month snapshots plus correct tendencies use TOTTTEND averages (Wang 2012).
- Era: 2012-2026. Leap-year bug: with a Gregorian calendar and timePhase > 366 d (e.g. timePhase=946771200 for a 30 yr spin-up, start 19600101), output lands on Jan-2 instead of Jan-1 in some years (issue #992). Fix proposed in PR #1013 (open as of Oct 2026, not in master): ties year/month units to frequency. Workaround until merged: choose timePhase = a date-aligned multiple that does not cross leap years, or change cal_time2dump.F so shTime=myTime when usingGregorianCalendar and phase>31622400.
- Src: https://github.com/MITgcm/MITgcm/issues/992 ; https://github.com/MITgcm/MITgcm/pull/1013 ; mitgcm-support 2014-March 'Using pkg/calendar with diagnostics' http://mailman.mitgcm.org/pipermail/mitgcm-support/2014-March/thread.html ; 2012-June 'caledar package' http://mailman.mitgcm.org/pipermail/mitgcm-support/2012-June/thread.html ; 2018-November 'diagnostics at end of calendar months' http://mailman.mitgcm.org/pipermail/mitgcm-support/2018-November/thread.html ; 2020-September 'Data diagnostic output' http://mailman.mitgcm.org/pipermail/mitgcm-support/2020-September/thread.html ; 2025-July 'ptracer/DIC budget closure/calendar pkg' (above)

### Average over a restart shorter than the averaging frequency gives wrong mean / diagnostics pickup segfault
- Cause: averaging accumulators are not saved across runs unless diagnostics pickups are used (diag_pickup_write/diag_pickup_read, diag_pickup_read_mnc/_write_mnc in data.diagnostics); that facility was never completed or regression-tested.
- Fix: avoid restarts that cut an averaging period (align frequency with run length) or accumulate externally; do not rely on diag_pickup_read (seg fault reported 2025). Leap/phase: first record after restart covers only the current run.
- Era: 2023-2025. diag_pickup_read still exists in diagnostics_readparms.F; JMC: 'never fully completed and tested'.
- Src: mitgcm-support 2025-May 'diagnostics read and write pickups' http://mailman.mitgcm.org/pipermail/mitgcm-support/2025-May/thread.html ; 2023-August 'what happens when run is shorter than diagnostic averaging frequency' http://mailman.mitgcm.org/pipermail/mitgcm-support/2023-August/thread.html

## data.diagnostics: setup errors and warnings

### "WARNING - diag.# NN : NAME has not been filled (ndiag= 0) ... write ZEROS instead" (stats: "write UNDEF instead")
- Cause: the diagnostic is defined but the code that fills it was not executed. Typical: option off (DFrI_TH needs implicitDiffusion=T; momVort3 needs vectorInvariantMomentum=T; DFrE_TH not filled with implicit-only; EXF diagnostics only filled by pkg/exf; UBotDrag only with CPP flag), user diagnostic whose DIAGNOSTICS_FILL is under an #ifdef/branch not run (e.g. dic OMEGAC needs CPP CAR_DISS), packages under adjoint (TAF removed fill call), or iter=0 stats over a region with no data (sea ice at rest; 'CPL_Qic1' warnings in cpl_aim+ocn, benign, issue #900).
- Fix: read available_diagnostics.log (written when debugLevel>=1; lists exactly what exists for this config); STDERR lists which streams are unfilled. For flux-form momentum with no variable viscosity, momVort3 is not computed (comment IF(useVariableVisc) in mom_fluxform.F, or compute offline); for user diagnostics print the array before the FILL call.
- Era: 2005-2025, valid.
- Src: mitgcm-support 2018-July 'Diagnostics have not been filled' http://mailman.mitgcm.org/pipermail/mitgcm-support/2018-July/thread.html ; 2005-October 'diagnostic errors' http://mailman.mitgcm.org/pipermail/mitgcm-support/2005-October/thread.html ; 2015-August 'user-defined diagnostic filled with zeros' http://mailman.mitgcm.org/pipermail/mitgcm-support/2015-August/thread.html ; https://github.com/MITgcm/MITgcm/issues/900

### "*** ERROR *** DIAGNOSTICS_SET_POINTERS: <name> is not a Diagnostic"
- Cause: name not defined for this build: misspelled/not 8 characters (blank name after a missing comma or bad continuation), package not compiled/enabled, custom package's <pkg>_DIAGNOSTICS_INIT never called.
- Fix: names are exactly 8 characters (pad with spaces); put fields on explicit lines, e.g. fields(1:12,1) = 'UVEL    ','VVEL    ',...; check available_diagnostics.log; for a user pkg call <PKG>_DIAGNOSTICS_INIT from <pkg>_init_fixed.F guarded by IF(useDiagnostics) and make sure the package is in packages.conf/data.pkg (also add it to packages_init_fixed.F as other pkgs). Use debugMode=.TRUE. in eedata to see the S/R call sequence.
- Era: 2014-2015, valid.
- Src: mitgcm-support 2014-April 'Diagnostics Package' http://mailman.mitgcm.org/pipermail/mitgcm-support/2014-April/thread.html ; 2015-October 'New added Diagnostics-problem' http://mailman.mitgcm.org/pipermail/mitgcm-support/2015-October/thread.html

### "DIAGNOSTICS_ADDTOLIST: Exceed Max.Number of diagnostics ndiagMax= 500" / numDiags, numperlist, diagSt_size too small / segmentation or memory problems
- Cause: ndiagMax = total available diagnostics for the compiled packages (shown in available_diagnostics.log), not the ones you use; numDiags/diagSt_size are deliberately small default storage for active 2-D/3-D diagnostics and stats.
- Fix: copy pkg/diagnostics/DIAGNOSTICS_SIZE.h to code/ and raise ndiagMax; numDiags >= (sum of 3-D active fields)*Nr (+ 2-D count); diagSt_size for stat diags; numlists/numperlist/numLevels for stream sizes (separate "Exceed Max.Num. of Fields/list numperList=" errors). Large numDiags is double precision storage: for memory-limited runs use more cores (mcmodel=medium related relocation errors, 2015-Oct).
- Era: current defaults ndiagMax=500, numlists=10, numperlist=50, numDiags=1*Nr (spelled numDiags, was numdiags before), nRegions=0, sizRegMsk=1, nStats=4, diagSt_size=10*Nr.
- Src: mitgcm-support 2021-September 'Problem with diagnostics package: Exceed Max.Number of diagnostics ndiagMax= 500' http://mailman.mitgcm.org/pipermail/mitgcm-support/2021-September/thread.html ; 2014-April 'Diagnostics Package' (JMC) http://mailman.mitgcm.org/pipermail/mitgcm-support/2014-April/thread.html ; 2014-August 'Memory issues with regional simulation' http://mailman.mitgcm.org/pipermail/mitgcm-support/2014-August/thread.html ; 2015-October '(no subject)' http://mailman.mitgcm.org/pipermail/mitgcm-support/2015-October/thread.html

### Mixing 2-D (ETAN) and 3-D (UVEL) fields in one output stream gives only level 1 of the 3-D field
- Cause: a stream has one level list; 2-D and 3-D fields cannot share a stream.
- Fix: separate fileName for 2-D and 3-D diagnostics (also for NetCDF output). Use levels(:,n) to select levels; vertical slices or sub-regions are not supported in pkg/diagnostics (see below). To write a field on Nr+1 levels: CALL DIAGNOSTICS_SETKLEV(diagName, Nr+1, myThid) after ADDTOLIST.
- Era: 2007-2019, valid; warning added to tutorial (issue #146).
- Src: https://github.com/MITgcm/MITgcm/issues/146 ; mitgcm-support 2018-August 'Write diagnostics with different vertical levels in same file' http://mailman.mitgcm.org/pipermail/mitgcm-support/2018-August/thread.html ; 2007-January '2d and 3d-diagnostics in one netcdf file?' http://mailman.mitgcm.org/pipermail/mitgcm-support/2007-January/thread.html ; 2015-February 'Diagnostics output on Nr+1 levels' http://mailman.mitgcm.org/pipermail/mitgcm-support/2015-February/thread.html

### How to add your own diagnostic
- Cause: n/a (frequent question).
- Fix: in <pkg>_diagnostics_init.F: diagName (exactly 8 chars), diagTitle, diagUnits, diagCode (e.g. 'SM      MR      '; see diagnostics_main_init.F/diagnostics_init_early.F or doc table: char1 S/U/V/W, char2 U/V/M/Z grid point, char3 'r'/'R' for level-integrated, char9 U/M/L level position, char10 0/1/R/L/M level count), CALL DIAGNOSTICS_ADDTOLIST(diagNum, diagName, diagCode, diagUnits, diagTitle, 0, myThid). Mates (vectors): diagMate = diagNum+2 / diagNum. Fill: IF (useDiagnostics) CALL DIAGNOSTICS_FILL(arr,'SIetop  ',0..Nr,1,2,bi,bj,myThid). Use DIAGNOSTICS_FILL_RS for _RS arrays (RL vs RS mismatch breaks -ur4/real4 RS builds, PRs #276, #838). The diagCode does not change values, only coordinates; set bibjflag/'mate' correctly. No longer edit diagnostics_init_early.F / add2list. User-defined placeholders SDIAG1.., UDIAG1.. exist. To read a diagnostic inside the code: DIAGNOSTICS_GET_POINTERS + DIAGNOSTICS_GET_DIAG (diagnostics_utils.F).
- Era: 2009-2024, valid (names DIAGNOSTICS_ADDTOLIST, DIAGNOSTICS_SETKLEV, DIAGNOSTICS_FILL_RS exist).
- Src: mitgcm-support 2021-July 'Adding new diagnostics to the code' http://mailman.mitgcm.org/pipermail/mitgcm-support/2021-July/thread.html ; 2009-September '(no subject)' http://mailman.mitgcm.org/pipermail/mitgcm-support/2009-September/thread.html ; 2011-August 'accessing diagnostic values' http://mailman.mitgcm.org/pipermail/mitgcm-support/2011-August/thread.html ; https://github.com/MITgcm/MITgcm/pull/276 ; https://github.com/MITgcm/MITgcm/pull/838

### "DIAGNOSTICS_STATUS_ERROR ... DIAGNOSTICS_FILL called from the WRONG place" (after the last DIAGNOSTICS_WRITE in DO_THE_MODEL_IO)
- Cause: DIAGNOSTICS_WRITE (modelEnd=T) was called too early in a modified/coupled driver, so later fills are out of sequence.
- Fix: check DO_THE_MODEL_IO call and modelEnd setting in the coupling driver; test with low optimisation; JMC: not a diagnostics_switch_onoff.F problem.
- Era: 2018-Sep.
- Src: mitgcm-support 2018-September 'diagnostics fill error ...' http://mailman.mitgcm.org/pipermail/mitgcm-support/2018-September/thread.html

### Namelist read failures that look like diagnostics bugs (end-of-file, unit 11, crashes only with useDiagnostics=T)
- Cause: missing namelist terminator in data.diagnostics (no '&'/'/'), a '&' used inside a namelist (terminates it early; also in data.ptracers), missing eedata/mis-sized SIZE.h, or stdout buffering hiding the real failing line.
- Fix: end each namelist with ' &' (or '/'), no stray ampersands; check STDOUT.0000/ STDERR for last line; eedata must exist even if empty (' &EEPARMS' newline ' &'); set debugLevel/-DNML_TERMINATOR as appropriate in optfile; flush buffers (debug run).
- Era: 2009-2018.
- Src: mitgcm-support 2009-July 'output using diagnostics package' http://mailman.mitgcm.org/pipermail/mitgcm-support/2009-July/thread.html ; 2017-October 'MITgcm-support Digest, Vol 172, Issue 9' http://mailman.mitgcm.org/pipermail/mitgcm-support/2017-October/thread.html ; 2018-November 'Error on Cheyenne HPC with diagnostics package' http://mailman.mitgcm.org/pipermail/mitgcm-support/2018-November/thread.html ; 2017-September 'getting started and the eedata file' http://mailman.mitgcm.org/pipermail/mitgcm-support/2017-September/thread.html

### Output precision: diagnostics file is 64-bit when you want 32-bit (or reverse); NetCDF in float
- Cause: precision of pkg/diagnostics output follows writeBinaryPrec (default 32), unless set per file.
- Fix: writeBinaryPrec=32/64 in `data` PARM01 (also governs MNC); per stream in data.diagnostics: fileFlags(n)='D       ' (64-bit) or 'R       ' (32-bit); readBinaryPrec for input. Budget analysis: use 64-bit.
- Era: 2011-2016, valid (fileFlags exists).
- Src: mitgcm-support 2016-October '32-bit diagnostic output' http://mailman.mitgcm.org/pipermail/mitgcm-support/2016-October/thread.html ; 2011-August 'MNC and variable types' http://mailman.mitgcm.org/pipermail/mitgcm-support/2011-August/thread.html

## Statistics diagnostics, regions, sub-domain output

### Regional / single-point statistics (DIAG_STATIS_PARMS, region masks)
- Cause: regions are set via mask files; documentation is thin.
- Fix: compile with #define DIAGSTATS_REGION_MASK (DIAG_OPTIONS.h) and set nRegions and sizRegMsk in DIAGNOSTICS_SIZE.h (missing sizRegMsk gives 'A COMMON block data object must not be an automatic object' compile error, 2010-Sep); in data.diagnostics: diagSt_regMaskFile='regMask.bin' (size Nx*Ny*nSetRegMskFile, real*8 if readBinaryPrec=64), nSetRegMskFile (layers in mask file), set_regMask(i)=layer of region i, val_regMask(i)=mask value of region i, stat_region(:,n)=region indices (integers, no dots), stat_fields/stat_fname/stat_freq/stat_phase. Overlapping regions: use extra mask layers. Single grid point time series every step: region mask of 1 point, stat_freq = deltaT (example aim.5l_cs input.thSI, global_ocean.cs32x15 input.thsice). Request region 0 (global) explicitly: otherwise the first record of NetCDF stats output is junk. Statistics are over ocean points (maskC/W/S), area/volume weighted, and exclude OBCS boundary points; '_vol' includes the number of accumulated steps (use SIheff_ave*SIheff_vol with care).
- Era: 2009-2022, valid. Wrong region count (e.g. 386 regions) and file errors came from not changing stat_region. MITgcmutils.diagnostics.readstats failed with >1 region: fixed in PR #823 (2024).
- Src: mitgcm-support 2022-January 'High-frequency time series output at a single domain point' http://mailman.mitgcm.org/pipermail/mitgcm-support/2022-January/thread.html ; 2016-June 'Defining different regions over the same area' http://mailman.mitgcm.org/pipermail/mitgcm-support/2016-June/thread.html ; 2010-October 'output of regional statistics' http://mailman.mitgcm.org/pipermail/mitgcm-support/2010-October/thread.html ; 2010-September 'compile problems with regional statistics' http://mailman.mitgcm.org/pipermail/mitgcm-support/2010-September/thread.html ; 2009-April 'Regional statistics diag (region-mask)' http://mailman.mitgcm.org/pipermail/mitgcm-support/2009-April/thread.html ; 2010-April 'diagnostics: diagstats' http://mailman.mitgcm.org/pipermail/mitgcm-support/2010-April/thread.html ; https://github.com/MITgcm/MITgcm/issues/818

### Cannot write a sub-region or vertical section with pkg/diagnostics; need only point/section output
- Cause: only full 2-D horizontal levels or 3-D fields; stats machinery gives statistics only, and is slow for many points.
- Fix: write full fields and cut offline; for xz/yz slabs call WRITE_REC_XZ_RL (pkg/rw/write_rec.F; example pkg/obcs/obcs_output.F) from write_state.F; or write 3-D scratch and post-process in chained jobs; for a time series use stats with a 1-point mask.
- Era: 2012-2020, valid.
- Src: mitgcm-support 2018-May 'Output 2D slice from 3D simulation' http://mailman.mitgcm.org/pipermail/mitgcm-support/2018-May/thread.html ; 2013-February 'Advice on writing 2-D vertical sections from 3-D model' http://mailman.mitgcm.org/pipermail/mitgcm-support/2013-February/thread.html ; 2019-December 'How to save output for a small sub-region' http://mailman.mitgcm.org/pipermail/mitgcm-support/2019-December/thread.html

### Snapshot adjoint diagnostics (ADJtheta etc.) have the wrong sign for frequency < 0
- Cause: ADJ snapshot output timing in diagnostics_write_adj.F differed from forward diagnostics (diagnostics_switch_onoff.F fill decision), so the field was filled at the wrong step (sign flip for freq <= -2*deltaT).
- Fix: update to code with PR #783 (autodiff_inadmode_set_ad.F call to DIAGNOSTICS_SWITCH_ONOFF adjusted); workaround: time-average (freq>0) or adjDumpFreq.
- Era: issue #774 (Sep 2023), fixed upstream in PR #783 (Oct 2023).
- Src: https://github.com/MITgcm/MITgcm/issues/774 ; https://github.com/MITgcm/MITgcm/pull/783

## MDS / MNC / NetCDF output

### useMNC=.TRUE. but "pkg/mnc has not been compiled (#undef ALLOW_MNC)", or genmake2 "mnc package was enabled but tests failed to compile NetCDF applications ... DISABLED"
- Cause: genmake2 could not compile/link a NetCDF Fortran test (include path, library, compiler mismatch/name-mangling, only C interface installed) so mnc is silently disabled; also stale build dir/PACKAGES_CONFIG.h.
- Fix: inspect genmake.log (check_netcdf_libs section) in the build dir; set NETCDF_ROOT and export it (genmake2 looks in $NETCDF_ROOT/include,lib) or set INCLUDES/LIBS (-lnetcdff -lnetcdf) in the optfile; use netcdf-fortran built with the same Fortran compiler; add 'mnc' to code/packages.conf (+ useMNC=.TRUE. in data.pkg, data.mnc); then rerun genmake2 and make CLEAN; make depend. Needs HDF5 libs for netCDF4. Checked: if useMNC is TRUE but undef ALLOW_MNC is just a warning, model keeps writing MDS. MNC is disabled when running with >1 thread.
- Era: 2004-2025, valid. Current optfiles take NETCDF_ROOT/ NETCDF_F_ROOT.
- Src: mitgcm-support 2017-May 'Compiling error when using MNC package' http://mailman.mitgcm.org/pipermail/mitgcm-support/2017-May/thread.html ; 2012-July 'mnc package was enabled but tests failed' http://mailman.mitgcm.org/pipermail/mitgcm-support/2012-July/thread.html ; 2014-May 'netcdf on macbook' http://mailman.mitgcm.org/pipermail/mitgcm-support/2014-May/thread.html ; 2011-October 'Run mitgcm with openMP and mnc' http://mailman.mitgcm.org/pipermail/mitgcm-support/2011-October/thread.html ; 2025-December 'NCAR Derecho - Running with MNC on' http://mailman.mitgcm.org/pipermail/mitgcm-support/2025-December/thread.html

### MNC writes one NetCDF file per tile; useSingleCpuIO / globalFiles have no effect on NetCDF
- Cause: pkg/mnc is per tile; singleCpuIO/globalFiles only apply to MDS binary files.
- Fix: glue afterwards with utils/python/MITgcmutils/scripts/gluemncbig (e.g. gluemncbig state.*.t???.nc -o state.all.nc; -2 for 64-bit offset if a file would exceed 2 GB; usage via no args; can subset variables; fast and python3 compatible, PR #90); or, if the run is small, nSx=nSy=1 on one process; or use MDS + xmitgcm. Avoid utils/matlab/mnc_assembly.m (old netcdf toolbox, empty output with new Matlab) and the shell gluemnc/xplodemnc (needs nco; NC_EVARSIZE errors >2 GB: add -4/--64). Global NetCDF patch exists (PR #31, never merged).
- Era: 2005-2020; gluemncbig exists in utils/python/MITgcmutils/scripts.
- Src: mitgcm-support 2018-November 'Matlab scripts for viewing netcdf files' http://mailman.mitgcm.org/pipermail/mitgcm-support/2018-November/thread.html ; 2017-May 'Some problems about SingleCpuIO' http://mailman.mitgcm.org/pipermail/mitgcm-support/2017-May/thread.html ; 2014-June 'gluemnc' http://mailman.mitgcm.org/pipermail/mitgcm-support/2014-June/thread.html ; 2014-June 'Problems' (Losch) http://mailman.mitgcm.org/pipermail/mitgcm-support/2014-June/thread.html ; https://github.com/MITgcm/MITgcm/pull/90 ; https://github.com/MITgcm/MITgcm/issues/364

### MNC restart/pickup errors: "MNC ERROR: variable 'RC' is already defined ... different grid shape", "dimension 'Z'/'Zp1' does not exist", "cannot open either a per-face or a per-tile file"
- Cause: MNC never overwrites existing .nc files; leftover grid.*.nc/state files from a previous run (or output named like a state file, e.g. PHIHYD vs phiHyd, case-insensitive match) collide; MNC pickups need identical decomposition and are incomplete ('incomplete MNC pickup files implementation' warning).
- Fix: rm *.nc before every run, or mnc_use_outdir=.TRUE. with mnc_outdir_str per run; with mnc_use_outdir=.FALSE. remove grid.*.nc or run with debugLevel=-1 (grids only written when nIter0=0 / debugLevel high, see below); use MDS pickups: pickup_write_mnc=.FALSE., pickup_read_mnc=.FALSE. in data.mnc together with useSingleCpuIO=.TRUE. (global pickups allow changing nPx/nSx). Do not name diagnostics streams like state-file names (phiHyd, state); set dumpFreq=0, dumpInitAndLast=.FALSE. to suppress state*.nc. Also increase MNC_MAX_ID if 'MNC_GET_NEXT_EMPTY_IND: array size exceeded'.
- Era: 2005-2019 (still true: pickup_write_mnc warning exists; Losch: 'I would never use mnc for pickups').
- Src: mitgcm-support 2012-August 'Problem using pickup with mnc' http://mailman.mitgcm.org/pipermail/mitgcm-support/2012-August/thread.html ; 2014-February 'problem with MNC pickup/restart' http://mailman.mitgcm.org/pipermail/mitgcm-support/2014-February/thread.html ; 2016-December 'strange mnc behavior' http://mailman.mitgcm.org/pipermail/mitgcm-support/2016-December/thread.html ; 2013-February 'Ice thickness category diagnostics' http://mailman.mitgcm.org/pipermail/mitgcm-support/2013-February/thread.html ; 2019-August 'MNC pickup' http://mailman.mitgcm.org/pipermail/mitgcm-support/2019-August/thread.html

### MNC output slower than MDS (up to 10x) or produces many files
- Cause: pkg/mnc per-tile NetCDF with syncs; filesystem/NetCDF installation dependent.
- Fix: use pkg/diagnostics with MDS (plus useSingleCpuIO=.TRUE. on shared filesystems) and read with xmitgcm/rdmds; fewer fields per .nc file; for coupled runs set cpl_taveFreq (data.cpl) to reduce coupler averages.
- Era: 2014-2020 (jrscott benchmark found no big difference on his system).
- Src: https://github.com/MITgcm/MITgcm/issues/364 ; mitgcm-support 2014-January 'coupler time averages writing frequency' http://mailman.mitgcm.org/pipermail/mitgcm-support/2014-January/thread.html

### MNC output has Xp1/Yp1 one longer than X/Y (halo copy), XG/YG last column duplicates the first
- Cause: by design MNC writes U(1:Nx+1), V(1:Ny+1) with the exchanged neighbour/periodic value; MDS does not. gluemnc (shell) truncates, gluemncbig keeps the extra point.
- Fix: drop the last row/column, e.g. XG(1:end-1,1:end-1); not a bug.
- Era: 2013-2018.
- Src: mitgcm-support 2018-May 'Inconsistent sizes between MNC and binary output' http://mailman.mitgcm.org/pipermail/mitgcm-support/2018-May/thread.html ; 2013-February 'MNC grid boundaries with curvilinear coordinates' http://mailman.mitgcm.org/pipermail/mitgcm-support/2013-February/thread.html ; 2017-March 'Strange U-velocity and V-velocity on North, East' http://mailman.mitgcm.org/pipermail/mitgcm-support/2017-March/thread.html

### MNC missing values (useMissingValue) mask ice-shelf or surface diagnostics wrongly
- Cause: diagnostics_mnc_out.F applies maskC at klev=1 to compressed 2-D fields (shelfice SHIfwFlx, SHIhtFlx, SHIuStar; also p-coordinate and orography cases).
- Fix: only NetCDF diagnostics output; PR/patch by Naughten/Losch in issue #97 (open). Workaround: do not use useMissingValue for shelfice fields, or mask offline with kTopC. Vector components zero (not missing) on land boundary; hFacW cannot be used as missing flag. CPP DIAGNOSTICS_MISSING_VALUE is old name.
- Era: issue #97 (2018), still open at last check.
- Src: https://github.com/MITgcm/MITgcm/issues/97 ; mitgcm-support 2010-January 'mnc missing values' http://mailman.mitgcm.org/pipermail/mitgcm-support/2010-January/thread.html

### MDS: one file per tile (S.0000000000.001.001.data ...), joinmds fails, restart with different processor count
- Cause: default writes per tile; joinmds needs complete .meta; useSingleCpuIO gives one global file per field and global pickups; globalFiles is not robust across file systems.
- Fix: useSingleCpuIO=.TRUE. (best, requires data to fit on one node) or globalFiles=.TRUE.; rdmds/MITgcmutils.rdmds reads tiles transparently (wildcards like '0*/T' work); merge existing tile files by rdmds + wrmds (MITgcmutils mds.py). Restarting on a different decomposition requires global (single-file) pickups. useSingleCpuIO is not implemented for vector/slice output (OBCS control xx_obcs*, use globalFiles=.TRUE. too). Local-disk output (useSingleCpuIO=.FALSE.) can be much faster on big clusters.
- Era: 2007-2025, valid.
- Src: mitgcm-support 2025-July 'Question about generating unified .meta and .data files' http://mailman.mitgcm.org/pipermail/mitgcm-support/2025-July/thread.html ; 2014-July 'Outputting too many files' http://mailman.mitgcm.org/pipermail/mitgcm-support/2014-July/thread.html ; 2015-February 'restart .meta & .data files from a different processor decomposition?' http://mailman.mitgcm.org/pipermail/mitgcm-support/2015-February/thread.html ; 2016-October 'xmitgcm' http://mailman.mitgcm.org/pipermail/mitgcm-support/2016-October/thread.html ; 2008-May 'is this bug in mitgcm code?' http://mailman.mitgcm.org/pipermail/mitgcm-support/2008-May/thread.html

### Grid files (grid.*.nc, XC/YC .data) not written
- Cause: write_grid is only called at the first iteration of a run from time 0 (startTime=baseTime) or when debugLevel >= debLevA/C; on restarts they are skipped.
- Fix: run once from nIter0=0 (or a 0-step run) to produce grids, or raise debugLevel; debugLevel=-1 suppresses grid.* on restarts.
- Era: 2013-2022 (initialise_fixed.F: IF debugLevel>=debLevA .OR. startTime==baseTime).
- Src: mitgcm-support 2022-January 'model_grid_output' http://mailman.mitgcm.org/pipermail/mitgcm-support/2022-January/thread.html ; 2013-February (Losch) as above

### rdmds / MITgcmutils errors (ValueError parsing .meta; NaN import; cannot read missing .meta; dimension reshape)
- Cause: .meta missing for one tile; MITgcmutils old versions failed on numpy 2 (from numpy import NaN) and Python 3.7 re.split; Matlab rdmds.m matched T* to Ttave/TSK (fixed by wildcard fix 2013); Fortran array order vs Matlab; changing Nx instead of sNx.
- Fix: pip install --upgrade MITgcmutils (v0.2+) or update MITgcm/utils/python; each .data needs matching .meta (rdmds('XC')); use rdmds(fname, iter) with explicit iteration; transpose fields in Matlab (contour(xc',yc',eta')); change sNx/nSx/nPx not Nx; rebuild (rm mitgcmuv, make) and delete old .data/.meta after SIZE.h changes; always make CLEAN when headers change.
- Era: 2013-2024; issues #893 (fixed Jun 2024, PyPI 0.2), #139 (PR #140).
- Src: https://github.com/MITgcm/MITgcm/issues/893 ; https://github.com/MITgcm/MITgcm/issues/139 ; mitgcm-support 2024-November 'rdmds & error' http://mailman.mitgcm.org/pipermail/mitgcm-support/2024-November/thread.html ; 2021-April '(no subject)' (changing Nx) http://mailman.mitgcm.org/pipermail/mitgcm-support/2021-April/thread.html ; 2013-May 'Bug in rdmds.m when loading multiple record files' http://mailman.mitgcm.org/pipermail/mitgcm-support/2013-May/thread.html

### Reading MDS output in other tools; offline-run file naming
- Cause: files are big-endian, headerless; names carry the 10-digit iteration; offline pkg reads iterations according to rwSuffixType.
- Fix: np.fromfile(f,'>f4' or '>f8').reshape((ny,nx)); write input with a.astype('>f8').tofile(...) (match readBinaryPrec); reading with Matlab fread(...,'b'); pkg/offline 'MDS_READ_FIELD: Files DO not exist' for u.0000000144: match rwSuffixType (default 0 = 10-digit iteration; ECCO v4 style 'u.0001.data' needs the same rwSuffixType as the run that wrote them) and the period in data.off. FLT output is per-tile without .meta for tiles without floats; read with verification/flt_example/input/read_flt_traj.m and use big-endian convert in Fortran.
- Era: 2015-2025; rwSuffixType exists today.
- Src: mitgcm-support 2025-February 'Inquiry for the offline experiment' http://mailman.mitgcm.org/pipermail/mitgcm-support/2025-February/thread.html ; 2015-April 'open source software for creating binary files from NetCDF' http://mailman.mitgcm.org/pipermail/mitgcm-support/2015-April/thread.html ; 2020-November 'Merging Floaters' Trajectory Output' http://mailman.mitgcm.org/pipermail/mitgcm-support/2020-November/thread.html

## Interpreting specific diagnostics

### PHIHYD / PH / PHIBOT / PHI_NH units and how to get pressure
- Cause: phiHyd is a potential, not pressure; z-coordinate output is m^2/s^2 (p/rho) anomaly relative to g*rhoConst*|z|; default rhoConst=999.8.
- Fix: P = (PHIHYD + g*|RC|)*rhoConst (nothing else needed dynamically); PHIBOT: P_b = PHIBOT*rhoConst + g*rhoConst*H, includes ETAN (and atm pressure) so remove the horizontal mean for bottom-pressure anomalies (REMOVE_MEAN_RL); for NH runs PH + PHI_NH (bottom NH value approximated by phi_nh at last wet level; PHIBOT excludes it). PHIBOT before 2010-01-15 was wrong with nonlinFreeSurf=4 (fixed). PHRefC = g*Z; linear-EOS reference density integral is not stored. In p-coordinates ETAN is bottom pressure anomaly and PHIBOT is g*SSH. rhoNil (EOS scaling, 'rho_c') differs from rhoConst ('rho_0').
- Era: 2004-2022, valid.
- Src: mitgcm-support 2013-November 'What is the unit of PH and PNH in the output files?' http://mailman.mitgcm.org/pipermail/mitgcm-support/2013-November/thread.html ; 2018-November 'Pressure field in MITgcm experiment' http://mailman.mitgcm.org/pipermail/mitgcm-support/2018-November/thread.html ; 2010-February 'bottom pressure diagnostics' http://mailman.mitgcm.org/pipermail/mitgcm-support/2010-February/thread.html ; 2004-August 'phiHydLow' http://mailman.mitgcm.org/pipermail/mitgcm-support/2004-August/thread.html ; 2019-December 'Bottom pressure, non hydrostatic run' http://mailman.mitgcm.org/pipermail/mitgcm-support/2019-December/thread.html ; 2022-October 'How to calculate the phiHyd' http://mailman.mitgcm.org/pipermail/mitgcm-support/2022-October/thread.html

### Buoyancy frequency / stratification from output with a nonlinear EOS
- Cause: in-situ density with a nonlinear EOS needs model pressure at the levels; recomputing offline gives spurious N2<0.
- Fix: use diagnostic DRHODR (model vertical density gradient): N2 = -g/rhoConst*DRHODR; density anomaly RHOAnoma for buoyancy b = g*(RHOAnoma+rhoConst)/rhoConst.
- Era: 2007-2022, valid.
- Src: mitgcm-support 2022-September 'Diagnostic calculation of buoyancy frequency with nonlinear equation of state' http://mailman.mitgcm.org/pipermail/mitgcm-support/2022-September/thread.html ; 2007-April 'Buoyancy anomaly' http://mailman.mitgcm.org/pipermail/mitgcm-support/2007-April/thread.html

### Bolus overturning from GM_PsiX/PsiY (and GM_Kwy) is non-zero at depth
- Cause: GM_PsiY is already vertically integrated from the bottom and masked (no hFacS needed), but zonal integration needs dxG*deepFacF(k); GM_Kwy/Kwx (skew-flux form) live at w-points and must be averaged to V points and masked with maskS before integrating.
- Fix: advective GM form: MOC_bolus(j,k)=SUM_i GM_PsiY(i,j,k)*dxG*deepFacF(k) (matches gmredi_residual_flow.F). Skew form: average Kwy horizontally, mask, then integrate (pkg/layers/layers_fluxcalc.F does it); skew-flux and advective forms legitimately give different eddy MOC (small-slope approximation). 3-D GM/Redi coefficients: only via pkg/ctrl (GM_background_K3dFile, GM_isopycK3dFile) or 1-D*2-D factors K = GM_background_K*GM_bolFac1d(k)*GM_bolFac2d(i,j) from GM_bol1dFile/GM_bol2dFile (GM_iso1dFile/2dFile for Redi).
- Era: 2013-2025.
- Src: mitgcm-support 2025-August 'GM_Psi diagnostic when using HFacC and DeepFacC' http://mailman.mitgcm.org/pipermail/mitgcm-support/2025-August/thread.html ; 2013-June 'GM: Skew flux vs Advective form' http://mailman.mitgcm.org/pipermail/mitgcm-support/2013-June/thread.html ; 2018-October 'GMRedi ptracers Diagnostics' http://mailman.mitgcm.org/pipermail/mitgcm-support/2018-October/thread.html

### EXFuwind/EXFvwind are not the 10-m winds with useRelativeWind=.TRUE.
- Cause: pkg/exf subtracts the ocean surface velocity from uwind/vwind in place; diagnostic names suggest raw winds.
- Fix: add ocean velocity back (uwind + 0.5*(uVel(i)+uVel(i+1)) at k=1) or avoid; also useRelativeWind does not work with constant-in-time wind (uwindperiod=0). Losch/JMC agreed u/vwind should remain unmodified (not yet changed).
- Era: issue #67 closed 2021; behaviour unconfirmed changed.
- Src: https://github.com/MITgcm/MITgcm/issues/67

### EXF diagnostics constant in time / ERA5 forcing sanity (EXFswdn, EXFlwdn, EXFatemp)
- Cause: diagnosing what the model reads is the way to catch bad forcing: ERA5 'net' or accumulated radiation fields (J/m2) instead of downward W/m2 gave lwdown<8 W/m2 and latent flux x3; a data.exf typo made the model read only the first record.
- Fix: need downward (not net) SW and LW in W/m2 (50<lwdown<450); check EXFswdn/EXFlwdn/EXFatemp/EXFhl in diagnostics; rewrite data.exf from a template; use debugLevel>3 to see forcing record reads; calendarDumps does not affect forcing periods.
- Era: 2025-Dec.
- Src: mitgcm-support 2025-December 'unrealistic model results in a regional model' http://mailman.mitgcm.org/pipermail/mitgcm-support/2025-December/thread.html

### pkg/layers diagnostics are all zero or have confusing names/units
- Cause: old namelist syntax (layers_G, NLayers) vs new (layers_name, layers_bounds, layers_bolus); LAYERS_THERMODYNAMICS diagnostics (LTto, LSto, LaTz, LaSz, LaTs, LaSs) were only filled when corresponding non-layer diagnostics were on (PR #720: new layers_useThermo); implicit diffusive fluxes (LaT/Sz1RHO) still miss (2026-Apr); names carry index 1TH/2SLT/3RHO by order in layers_name; MNC not implemented for layers; units were m deg/s (issue #40 fixed).
- Fix: use data.layers with layers_name(1)='TH', layers_bounds(1:N,1)=..., NLayers in LAYERS_SIZE.h; define LAYERS_THERMODYNAMICS in LAYERS_OPTIONS.h (default #undef for memory); output via data.diagnostics e.g. 'LaVH1TH ','LaHs1TH ','LaPs1TH '. Rule of thumb: Nlayers ~ Nr for uniform stratification, > Nr for global; for density, space layers by equal volume; need LAYERS_UFLUX/VFLUX/THICKNESS defined.
- Era: 2014-2026; layers_useThermo exists in master; CPP LAYERS_DIAG_TOTTEND mentioned in PR #720 discussion does not exist in master (final design uses layers_useThermo). Open list: issue #987.
- Src: https://github.com/MITgcm/MITgcm/pull/720 ; https://github.com/MITgcm/MITgcm/issues/987 ; https://github.com/MITgcm/MITgcm/issues/40 ; mitgcm-support 2014-June 'Layers package : nothing but zeros for output?' http://mailman.mitgcm.org/pipermail/mitgcm-support/2014-June/thread.html ; 2016-August 'Layers and NLayers' http://mailman.mitgcm.org/pipermail/mitgcm-support/2016-August/thread.html

### Mixed layer depth diagnostics: MXLDEPTH vs KPPhbl
- Cause: two definitions.
- Fix: MXLDEPTH (model/src/calc_oce_mxlayer.F) uses hMixCriteria: <0 (default -0.8 degC equivalent density, Kara et al. 2000) or >1 (local stratification exceeds mean above by that factor; dRhoSmall); KPPhbl is the KPP boundary layer depth.
- Era: 2015-Jul, valid (hMixCriteria=-0.8 in set_defaults.F).
- Src: mitgcm-support 2015-July 'the calculation of the MXLDEPTH' http://mailman.mitgcm.org/pipermail/mitgcm-support/2015-July/thread.html

### Frictional heating diagnostic (HeatDiss) is all zeros in the ocean
- Cause: addFrictionHeating/ALLOW_FRICTION_HEATING exist but the ocean never computes frictionHeating (only atmosphere code does).
- Fix: not implemented for ocean; add physics yourself in analogy with atm_phys tendency apply.
- Era: 2020-Oct; ALLOW_FRICTION_HEATING still in CPP_OPTIONS.h.
- Src: mitgcm-support 2020-October 'Frictional heating in the ocean model' http://mailman.mitgcm.org/pipermail/mitgcm-support/2020-October/thread.html

## Misc

### Removed: pkg/timeave and CPP EXACT_CONSERV
- Cause: superseded by pkg/diagnostics and always-on exact conservation.
- Fix: replace *tave output (Ttave, uVeltave) by diagnostics; run-time switch exactConserv remains. pkg/timeave averaged state variables with half-weights at first/last step (tave_lastIter), unlike diagnostics (full steps, 'snapshot of t-1').
- Era: removed Nov 2025 (PRs #926, #930, #939; issue #924); Jan 2012 thread documents differences.
- Src: https://github.com/MITgcm/MITgcm/issues/924 ; mitgcm-support 2012-January 'timing of snapshots and time averages in pkg/timeave and pkg/diagnostics' http://mailman.mitgcm.org/pipermail/mitgcm-support/2012-January/thread.html ; 2008-October 'output_mitgcm' http://mailman.mitgcm.org/pipermail/mitgcm-support/2008-October/thread.html

### DIC pkg: dic_pCO2 error; surface-only diagnostics; OMEGAC zero; 3-D silicate
- Cause: dic_pCO2 can only be changed if dic_int1=1 ('cannot change default dic_pCO2 if dic_int1=0'); pH/pCO2/flux are surface diagnostics; calcite saturation at depth needs CPP CAR_DISS (DIC_OPTIONS.h, default #undef, expensive); only surface silicate was read until PR #620 (Dec 2022).
- Fix: data.dic DIC_FORCING: dic_int1=1, dic_pCO2=340.E-6; define CAR_DISS; update to code after PR #620.
- Era: 2015-2022, verified dic_int1/dic_pCO2 in dic_readparms.F.
- Src: mitgcm-support 2015-January 'Assigning atmospheric CO2 levels in DIC module' http://mailman.mitgcm.org/pipermail/mitgcm-support/2015-January/thread.html ; 2021-April 'Regarding surface diagnostic parameters in DIC package' http://mailman.mitgcm.org/pipermail/mitgcm-support/2021-April/thread.html ; https://github.com/MITgcm/MITgcm/issues/619

### Mellor-Yamada (pkg/my82) / MY82 diagnostics location and units
- Cause: pkg/my82 is level 2.0 (MY-2.0, not 2.5); unmaintained; vertical diffusivities at W-points ('L' diagCode).
- Fix: diagnostic names MYHBL, MYVISCAR, MYDIFFKR (upper case after Martin's merge); KPPviscAz/KPPdiffKzT are at W points (interfaces). The tkel unit question: GM/GH are 1/s^2, so MYviscAr is m^2/s (OK).
- Era: 2008-2024.
- Src: mitgcm-support 2024-February 'Mellor-Yamada in MITgcm' http://mailman.mitgcm.org/pipermail/mitgcm-support/2024-February/thread.html ; 2013-March 'Is my82 just the one usually called MY-2.5?' http://mailman.mitgcm.org/pipermail/mitgcm-support/2013-March/thread.html ; 2008-April 'MY82 with netcdf output' http://mailman.mitgcm.org/pipermail/mitgcm-support/2008-April/thread.html

### Huge STDOUT, too many text files, memory and link errors
- Cause: monitor/debug output frequency; grids with static arrays > 2 GB.
- Fix: debugLevel=-1, monitorFreq=21600., useSingleCpuIO=.TRUE. (writes STDOUT only from rank 0); 'relocation truncated to fit: R_X86_64_PC32 against .bss' -> add -mcmodel=medium to FFLAGS/CFLAGS, or use more MPI ranks / smaller DIAGNOSTICS_SIZE arrays; undefine second-order-moment advection in GAD_OPTIONS.h if not used; OpenMP/MPI hybrid reduces per-core memory. Tiles < ~30x30 scale poorly (cg2d).
- Era: 2009-2017.
- Src: mitgcm-support 2009-April 'reduce the size of STDOUT.000X' http://mailman.mitgcm.org/pipermail/mitgcm-support/2009-April/thread.html ; 2017-January 'error' http://mailman.mitgcm.org/pipermail/mitgcm-support/2017-January/thread.html ; 2014-August 'Memory issues with regional simulation' http://mailman.mitgcm.org/pipermail/mitgcm-support/2014-August/thread.html ; 2012-May 'speedup for cs64 on a linux cluster' http://mailman.mitgcm.org/pipermail/mitgcm-support/2012-May/thread.html

### Zero output / zero depth because of wrong input files, grid sign or SIZE.h mismatch
- Cause: bathymetry must be negative (z of the bottom), hFacMin must be <= 1 (hFacMin=2.5 emptied the domain; hFacC all zeros -> monitor stats all zero); all input files big-endian with readBinaryPrec matching; dimension mismatch between SIZE.h and files (sNx*nSx*nPx); Matlab transpose; delX/dxSpacing wrong; U/V initial-condition files are Nx x Ny (model fills Nx+1/Ny+1), never Nx+1 x Ny.
- Fix: check hFacC.data, Depth.data, STDOUT grid summary; set hFacMin ~0.2-0.3, hFacMinDr; use dxSpacing/delX files carefully; recompile after SIZE.h edits (make CLEAN). Stable vertical grid ratio |dz(k+1)/dz(k)| < ~1.4 (Losch); layer 1 must be >> wave height.
- Era: 2004-2017.
- Src: mitgcm-support 2014-June 'Monitor statistics are all zeros except time and vorticity' http://mailman.mitgcm.org/pipermail/mitgcm-support/2014-June/thread.html ; 2011-October 'bathy file with Spherical Polar Grid' http://mailman.mitgcm.org/pipermail/mitgcm-support/2011-October/thread.html ; 2014-July 'Initial Conditions on Staggered Grid' http://mailman.mitgcm.org/pipermail/mitgcm-support/2014-July/thread.html ; 2008-April 'cell size variation' http://mailman.mitgcm.org/pipermail/mitgcm-support/2008-April/thread.html ; 2006-December 'urgent question' http://mailman.mitgcm.org/pipermail/mitgcm-support/2006-December/thread.html

### Anelastic option and RhoRef.data
- Cause: confusion over tRef/sRef/rhoRef/rhoConst.
- Fix: rhoRefFile alone turns on the anelastic formulation (loaded into rho1Ref; at the very top ratio=1, i.e. rhoConst); tRef/sRef initialise T/S (unless 3-D files) and, with nonlinear EOS, define rhoRef (pRef4EOS etc.); RhoRef.data/meta is always written; with linear EOS tRef/sRef have no real effect.
- Era: 2021-Feb (JMC); rhoRefFile in doc.
- Src: mitgcm-support 2021-February 'Confusion about Anelastic Option' http://mailman.mitgcm.org/pipermail/mitgcm-support/2021-February/thread.html
