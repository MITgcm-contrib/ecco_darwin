# Troubleshooting: adjoint / TAF / Tapenade / ctrl / cost / grdchk (MITgcm support history 2003-2026)
Scope: answered mitgcm-support + GitHub threads on AD builds (TAF, Tapenade, OpenAD), tapes/store directives, pkg/ctrl, pkg/ecco/cost, pkg/grdchk, optim, profiles. Verified against origin/master (ae4c03af7, 2026-10-03) and the mitgcm-index; ctrl went through rewrites (old ECCO_CTRL_DEPRECATED code removed PR 631/406 ~2021-22; generic ctrl_map_*genarr/gentim2d; ecco_cost generic gencost). OpenAD is being removed (issue 971).
Mail threads cited as "mitgcm-support YYYY-Month 'Subject'" (index: http://mailman.mitgcm.org/pipermail/mitgcm-support/YYYY-Month/thread.html). JM = Jean-Michel Campin.

## TAF errors / warnings and AD build

### TAF ERROR "keyword REC, KIND, or SHAPE expected" on `CADJ STORE ... kind = isbyte` (TAF -f08 / TAF_FORTRAN_VERS='F08')
- Cause: TAF bug with white space around `kind` in STORE directives, only under f08 mode (Ralf Giering). Not an MITgcm bug.
- Fix: use TAF >= 6.8.10. Stop-gap seen: `kind= isbyte` (works) vs `kind =isbyte` (fails); or comment the STORE.
- Era: 2025-10, TAF 6.8.8; fixed in TAF 6.8.10.
- Src: https://github.com/MITgcm/MITgcm/issues/940

### TAF ERROR "cannot generate correct recomputations for ..." / "unresolvable conflict" then run dies with "!!!!!!! PANIC !!!!!!! in S/R BARRIER myThid = 0"
- Cause: (1) massive recomputation that TAF gives up on (fix store directives first, read taf_ad.log); (2) myThid=0 comes from an argument-list mismatch between TAF-generated and hand-written AD routines, e.g. `#define AUTODIFF_TAMC_COMPATIBILITY` (TAMC-only, never use with TAF) or stale local copies of .F files not matching the code tree.
- Fix: remove AUTODIFF_TAMC_COMPATIBILITY; make local code_ad/*.F compatible with the tree (Mazloff: also check this for segfaults); fix recomputations before running.
- Era: 2010, 2015 (flag still in AUTODIFF_OPTIONS area today); 2018 reminder.
- Src: mitgcm-support 2015-February 'Help with TAF-generated adjoint'; mitgcm-support 2010-January 'Problem building adjoint'; mitgcm-support 2018-June 'segmentation fault' (thread: http://mailman.mitgcm.org/pipermail/mitgcm-support/2018-June/thread.html)

### TAF errors only after upgrading TAF: "TAF ERROR ... store directive for non-existing variable" / "storage directive not initialised"
- Cause: newer TAF checks what older TAF silently ignored: TAF 5.0.9 errors on `CADJ STORE` tapes without `CADJ INIT` (tapelev_ini_bibj, onetape); TAF 5.8.0 errors on STORE for variables that do not exist under the current CPP set (e.g. OmegaC without DIC_CALCITE_SAT).
- Fix: add `CADJ INIT <tape> = ...` ; wrap/remove stale STORE (check your own code_ad/ store files after TAF upgrades).
- Era: TAF 5.0.9 (2020-08, PR 365), TAF 5.8.0 (2023-03, PR 712); both fixed upstream.
- Src: https://github.com/MITgcm/MITgcm/pull/365 ; https://github.com/MITgcm/MITgcm/pull/712

### TAF ERROR "code of function tbar not seen and no flow information set, cannot generate derivative common block variable" (new cost term)
- Cause: custom cost-term edits left a variable without a flow path (partial edit of cost/ctrl files).
- Fix: do not hand-edit cost files; use pkg/ecco generic cost (gencost_barfile = 'm_boxmean_theta' etc., gencost_posproc = 'boxmean') so cost terms switch at run time; start from verification_other/global_oce_cs32 input_ad.sens data.ecco.
- Era: 2023-03 (ecco gencost era).
- Src: mitgcm-support 2023-March 'easier ways to implement cost function' (http://mailman.mitgcm.org/pipermail/mitgcm-support/2023-March/thread.html)

### `make adall`: "staf: Command not found" / ad_taf_output.f not generated
- Cause: TAF is commercial (FastOpt); the `staf` client must be on PATH and licensed. No license -> Tapenade (OpenAD being dropped).
- Fix: put staf in PATH; or build with `genmake2 -tap`/Tapenade instead.
- Era: 2016, 2020; still true.
- Src: https://github.com/MITgcm/MITgcm/issues/335 ; mitgcm-support 2016-May 'Errors occured on tutorial global oce optim' (http://mailman.mitgcm.org/pipermail/mitgcm-support/2016-May/thread.html)

### TAF errors on mixed fixed/free-format source; `-nocat4ad` (-ncad) breaks on *.flowdir / macOS `.for` files
- Cause: since TAF 4.1.1 *.flowdir (non .f/.F suffix) is assumed free-format; with -ncad all individual files go to TAF; macOS case-insensitive FS yields `.for` files read as free form.
- Fix: genmake2 adds `-fixed` to TAF_EXTRA for -ncad (PR 208); mixed .F + .F90 supported for AD (PR 709) and TL (PR 863) with TAF >= 5.8.0; F90 set via `TAF_FORTRAN_VERS` (default in tools/adjoint_options/adjoint_default) / `ALWAYS_USE_F90=1` in genmake_local. Suggested (not done) default `-f90` for TAF: issue 854.
- Era: 2019 (PR 208), 2023-10 (issue 785, gone after PR 709), 2024-08 (PR 863). Fixed upstream.
- Src: https://github.com/MITgcm/MITgcm/pull/208 ; https://github.com/MITgcm/MITgcm/issues/785 ; https://github.com/MITgcm/MITgcm/pull/863 ; https://github.com/MITgcm/MITgcm/issues/854

### `make -j N adtaf` intermittently fails comparing ad_config.template / AD_CONFIG.h (race)
- Cause: two make targets (ad_input_code.f, ad_inpF90_code.f90) used the same temp file ad_config.template.
- Fix: PR 915 uses ad_config.template{0,1,2}. Older trees: use `-ncad` or no -j for adtaf.
- Era: 2025-04; fixed upstream.
- Src: https://github.com/MITgcm/MITgcm/issues/914 ; https://github.com/MITgcm/MITgcm/pull/915

### TAF-generated tangent-linear model fails at run time: "EXF_INTERP_READ: filename ???? File does not exist" (USE_EXF_INTERPOLATION)
- Cause: exf_set_uv_tl.f calls a non-existent exf_interp_uv_tl; missing flow information for exf_getyearlyfieldname.
- Fix (Losch): add `CADJ SUBROUTINE exf_getyearlyfieldname REQUIRED` to exf_ad.flow (exf_getyearlyfieldname exists in pkg/exf today); TLM of global_oce_latlon.w_exf not yet tested.
- Era: 2025-12, issue still open when digested.
- Src: https://github.com/MITgcm/MITgcm/issues/955

### Link error "undefined reference to adexch_uv_3d_rl_ / adexch_xy_rs_" (and mdthe_main_loop_) in the adjoint
- Cause: simple AD set-ups where TAF never generates adjoints of the EXCH routines, but hand-written adjoint-output code (addummy_in_stepping.F, copy_ad_uv_outp.F, monitor_ad.F) calls them.
- Fix: comment out those calls, or `#undef ALLOW_AUTODIFF_MONITOR` in AUTODIFF_OPTIONS.h (then no ADJ* output).
- Era: 2015-11 (c65i); mechanism unchanged.
- Src: mitgcm-support 2015-November 'Issue compiling offline model with TAF on ARCHER' (http://mailman.mitgcm.org/pipermail/mitgcm-support/2015-November/thread.html)

### TAF "*ERROR* : variable precip must be available" in a customised ctrl_map_ini_genarr.F
- Cause: local copy of ctrl file does not include the EXF headers declaring the variable (uwind worked because already in scope).
- Fix: add the exf includes, or better use xx_gentim2d_* for time-varying atmospheric controls (hs94.1x64x5 is the example).
- Era: 2015-10 (pre generic-ctrl rewrite; gentim2d now standard).
- Src: mitgcm-support 2015-October 'Control variables in adjoint sensitivity experiments' (http://mailman.mitgcm.org/pipermail/mitgcm-support/2015-October/thread.html)

### Link error "relocation truncated to fit: R_X86_64_32S" building mitgcmuv_ad (or "relocation truncated" in ifort libs)
- Cause: static executable > 2 GB because of tape arrays (nchklev_1 too big) or large tiles.
- Fix: lower `nchklev_1` in tamc.h (raise nchklev_2/3), or compiler memory model `-mcmodel=medium|large`.
- Era: 2010 (Pleiades), 2018 (ARCHER); unchanged.
- Src: mitgcm-support 2010-January 'Problem building adjoint' ; mitgcm-support 2018-May 'ADWRITE package error "tapeFileCounter > tapeMaxCounter"' (http://mailman.mitgcm.org/pipermail/mitgcm-support/2018-May/thread.html)

### Old checkpoints only: "duplicate symbol _adams_bashforth2_" / "undefined option -nonew_arg" / wrong adjoint after a TAF default change
- Cause: FastOpt changed TAF defaults (v2 naming) or dropped flags; unpinned TAF can silently give wrong code (TAF 1.9.12: log normal, recomputations wrong).
- Fix: pin with `-version X` in AD_TAF_FLAGS or update tools/adjoint_options/*; use gradient/TLM checks after every TAF upgrade.
- Era: 2007-2010 only; historical.
- Src: mitgcm-support 2010-October 'Linking problem (with adjoint, gfortran, OS X)'; mitgcm-support 2009-April 'taf errors'; mitgcm-support 2007-October 'note to TAF users'

## Recomputation, tapes and store directives

### "TAF RECOMPUTATION WARNING" (any): how the core team handles it
- Cause: store directive missing/wrong key, or TAF cannot see that an if-branch is the same; warnings in taf_ad.log. Recomputation can also change results.
- Fix: testreport counts them (e= warnings in the summary line); fix by (a) local tape + store inside the routine, (b) ikey from all loop indices (e.g. `ikey = bi + (bj-1)*nSx + (ikey_dynamics-1)*nSx*nSy`, no per-thread grouping needed), (c) reorder code so scalar temporaries become 1D/2D arrays (k-loop outside i,j). Always compare gradients with and without storing: wrong keys silently corrupt gradients.
- Era: 2020-2023 (Losch/JM clean-ups: PR 395, 496, 628, 731). Keys: issue 608.
- Src: https://github.com/MITgcm/MITgcm/pull/628 ; https://github.com/MITgcm/MITgcm/issues/608 ; https://github.com/MITgcm/MITgcm/issues/763

### Odd result: `if (c) call A; if (c) call B` recomputes but `if (c) {call A; call B}` does not
- Cause: TAF cannot prove the two IF blocks are the same; combine blocks (or call a dummy else-branch).
- Fix: put calls in a single IF; in forward_step the update_cg2d call is always made under ALLOW_AUTODIFF_TAMC for this reason.
- Era: 2020-11.
- Src: https://github.com/MITgcm/MITgcm/issues/391

### Inner (checkpoint level 1) tape computation unstable / wrong with nonlinFreeSurf > 2 (r*), calc_r_star "too SMALL rStarFac"
- Cause: cg2d is declared self-adjoint in cg2d.flow, so TAF dropped update_cg2d from the tape recomputation; for nonlinear free surface the coefficients change each step.
- Fix: PR 392: hand-written fake adjoint (cg2d_sad.F, now pkg/autodiff/cg2d_mad.F, `CADJ SUBROUTINE cg2d ADNAME = cg2d_mad`) + extra storage of the 6 cg2d common-block fields. Quick fix on old trees: `CADJ SUBROUTINE update_cg2d REQUIRED` in cg2d.flow. TAF warning "self adjoint routine has more than one active input" disappears. Monitor each tape level by commenting the monitorFreq reset in turnoff_model_io.F.
- Era: 2020-11 -> fixed upstream (2020-12). All results with nonlinFreeSurf>2 changed.
- Src: https://github.com/MITgcm/MITgcm/issues/391 ; https://github.com/MITgcm/MITgcm/pull/392

### mitgcmuv_ad blows up (CALC_R_STAR too SMALL) 16 days into the REVERSE run, forward run fine; or garbage gradients far from tape start
- Cause: forward recomputation during tape levels is marginal; ECCO-style tapes store single precision.
- Fix: `doSinglePrecTapelev = .FALSE.` in data.ctrl (Ou Wang); then check nonlinFreeSurf/cg2d issue above; turn off problematic pkgs in AD (seaice, ggl90/kpp, gmredi, saltplume as Losch did).
- Era: 2020-08 (doSinglePrecTapelev still read in ctrl_readparms.F).
- Src: mitgcm-support 2020-August 'mitgcmuv_ad explodes during tape computations' (http://mailman.mitgcm.org/pipermail/mitgcm-support/2020-August/thread.html)

### Adjoint run blows up (adxx ~1e16, finite-difference fine, optim "linesearch failed")
- Cause: some physics unstable in reverse (KPP, GMRedi, ptracers, seaice, some advection schemes).
- Fix: turn those off for the adjoint stage only: set useKPP/useGMREDI/usePtracers/useSEAICE=.FALSE. in autodiff_inadmode_set (or data.pkg of the divided-adjoint passes); multiDimAdvection=.FALSE., tempAdvScheme=saltAdvScheme=30; a TLM vs FD test finds store-directive faults; GMREDI_WITH_STABLE_ADJOINT + GM_taper_scheme='stableGmAdjTap' is the ECCO v4r5 recipe; GGL90: mxlMaxFlag=1 (2 is unstable in AD).
- Era: 2009 original; ggl90 advice 2021-2026.
- Src: mitgcm-support 2009-August 'diagnosing problems with the adjoint' ; https://github.com/MITgcm/MITgcm/issues/518

### Wrong-looking gradients with GMREDI_WITH_STABLE_ADJOINT (stableGmAdjTap): tape key in gmredi_slope_limit
- Cause: keys did not include k and the 3 calls from gmredi_calc_tensor, so tape held only slopeX/Y of level 1, last call. ECCO v4r5 uses this option.
- Fix: PR 686 (local tape loctape_gm / proper key with new kPos argument). No change was visible in global_oce_cs32 / llc90 results.
- Era: 2022-11, fixed upstream 2023-02.
- Src: https://github.com/MITgcm/MITgcm/issues/668

### Run dies "ADWRITE: tapeFileCounter > tapeMaxCounter 691 690" (tapelev3_..._gptrnm1)
- Cause: ALLOW_AUTODIFF_WHTAPEIO + ALLOW_WHIO_3D buffer too small for the number of 3D fields stored (BLING tracers).
- Fix: raise `nWh` in pkg/mdsio/MDSIO_BUFF_WH.h (default `nWh=30*Nr`; 60*Nr for BLING; 633 needed with PHYTO_SELF_SHADING).
- Era: 2018-05, 2024-06; code unchanged (nWh still a PARAMETER).
- Src: mitgcm-support 2018-May 'ADWRITE package error "tapeFileCounter > tapeMaxCounter"' (http://mailman.mitgcm.org/pipermail/mitgcm-support/2018-May/thread.html); https://github.com/MITgcm/MITgcm/pull/839

### Wrong adxx (all negative / flipped sign) with divided adjoint + ALLOW_AUTODIFF_WHTAPEIO + ALLOW_ADCTRLBOUND
- Cause: with WHTAPEIO the init tapes (tapelev_init, tapelev_ini_bibj_k) were common blocks, lost between DIVA passes; control bounds then triggered wrongly in adctrl_bound.F.
- Fix: PR 756 (AUTODIFF_WHTAPEIO_SYNC(1,...) around INITIALISE_VARIA in the_main_loop.F; lab_sea.noseaicedyn grdchk on xx_salt, grdchkvarindex=202, tests it). Gradient printed in STDOUT did not reveal it; adxx files did.
- Era: 2023-08, fixed upstream.
- Src: https://github.com/MITgcm/MITgcm/pull/756

### divided adjoint: "end-of-file during read, unit 76, file divided.ctrl" at the end of the run
- Cause: normally benign. divided.ctrl goes `2 1` -> `1 0` -> `0 -1`; `0 -1` means the adjoint cycle is complete.
- Fix: nothing; remove costfinal / divided.ctrl before next run. Genuine EOF in gencost runs comes from time masks that are too short (see gencost mask entry).
- Era: 2023-03/04.
- Src: mitgcm-support 2023-April 'end-of-file during read file divided.ctrl' (http://mailman.mitgcm.org/pipermail/mitgcm-support/2023-April/thread.html)

### ecco_phys.F recomputation after adding BEGIN_MASTER/END_MASTER (multithreaded safety) with ALLOW_PSBAR_STERIC
- Cause: TAF does not understand MITgcm master-thread blocks on non-tiled common-block variables: major recomputation, slow-down.
- Fix: PR 726 reverts them when pkg/autodiff is compiled (keeps local VOLsumGlob_1/RHOsumGlob_1 for thread safety); alternative with store directives judged too intrusive.
- Era: 2023-04; fixed upstream. Adjoint is not multithread-safe in general.
- Src: https://github.com/MITgcm/MITgcm/pull/726

### Major recomputations from pkg/ecco (ecco_phys.F below line ~141) when NONLIN_FRSURF is undefined
- Cause: a store needed outside the NONLIN_FRSURF block.
- Fix: PR 250 (extra `#else` store block in checkpoint_lev1_directives.h: detahdt, gsnm1/gtnm1/gunm1/gvnm1 or gsnm/gtnm..., wvel, etaH). Also keep DISABLE_SIGMA_CODE defined with autodiff.
- Era: 2018-2019; fixed 2019-11.
- Src: https://github.com/MITgcm/MITgcm/issues/68

### pkg/bling AD: wrong store keys / recomputations (bling_light.F, mixedlayer, bio_nitrogen), PHYTO_SELF_SHADING crashes
- Cause: 3D fields stored inside i,j,k loops with tile-only key; k-loop innermost; `chl**(e_rd-1)` with chl=0 gives division by zero in the adjoint.
- Fix: PR 839/628: k-loop outside, 2D stores, global tapes (+1.3% exe size); fix/define PHYTO_SELF_SHADING and ML_MEAN_LIGHT in a local BLING_OPTIONS.h; USE_QSW undef when no exf.
- Era: 2022-04 to 2024-06; fixed upstream.
- Src: https://github.com/MITgcm/MITgcm/issues/763 ; https://github.com/MITgcm/MITgcm/pull/839 ; https://github.com/MITgcm/MITgcm/pull/628

### Sea-ice AD: seaice_dynsolver called from seaice_model_ad; useHB87stressCoupling + SEAICE_deltaTdyn>deltaT; SEAICE_USE_GROWTH_ADX
- Cause: store directives of uice/vice placed before the solver; ALLOW_AUTODIFF_TAMC used where ALLOW_AUTODIFF meant; HB87 coupling with larger dyn step not AD-safe.
- Fix: PR 395 (stores after seaice_dynsolver), PR 496 (fewer stores, seaice_model recompute avoided), PR 731 (stop in seaice_check.F for HB87 case); lab_sea has `#undef SEAICE_USE_GROWTH_ADX` so its gradients stay inexact.
- Era: 2020-12 to 2023-05; fixed upstream.
- Src: https://github.com/MITgcm/MITgcm/pull/395 ; https://github.com/MITgcm/MITgcm/pull/731 ; https://github.com/MITgcm/MITgcm/issues/518

### Convective adjustment (isomip) gradient depends on unrelated code (store of theta/salt at k-1)
- Cause: TAF restores only level k before convectively_mixtracer_ad.
- Fix: PR 457 new store directives in convective_adjustment.F.
- Era: 2021-04; fixed.
- Src: https://github.com/MITgcm/MITgcm/pull/457

### Adjoint 10x slower with an advection or physics option (GAD_SMOLARKIEWICZ_HACK, RBCS, PCELL_MIX_CODE)
- Cause: extensive recomputations (not "adjointed" options); read taf_ad.log.
- Fix: PCELL_MIX_CODE / dwnslope_apply fixed by PR 628 (local tape, reorganised loops); RBCS and Smolarkiewicz never fixed in these threads.
- Era: 2009-2022.
- Src: https://github.com/MITgcm/MITgcm/pull/628 ; mitgcm-support 2009-July 'GAD_SMOLARKIEWICZ_HACK and adjoint code'; mitgcm-support 2010-March 'RBCS adjoint recomputation problem'

### SOLVE_DIAGONAL_LOWMEMORY / SOLVE_DIAGONAL_KINNER with AD (implicit solvers)
- Cause: default undef/undef = AD-safe and vectorises; LOWMEMORY is fast but not AD-suitable; KINNER AD-suitable but slowest.
- Fix: keep both #undef (model/inc/CPP_OPTIONS.h, still present) for adjoint; LOWMEMORY ok for forward-only LLC.
- Era: 2021-01 (JM, Losch).
- Src: mitgcm-support 2021-January 'SOLVE_DIAGONAL_LOWMEMORY' (http://mailman.mitgcm.org/pipermail/mitgcm-support/2021-January/thread.html)

### AD arrays real*4 mismatch: floating overflow in update_cg2d_ad.f with `-ur4` and cg2dFullAdjoint=.TRUE.; adread/adwrite type error with `-warn all`
- Cause: hand-written AD arrays declared _RL while forward arrays are _RS; adread/adwrite only handle real*8 (integer kTopC, _RS fields).
- Fix: PR 527 (cg2d_mad.F types); for integers use a separate user tape (`CADJ INIT tapelev3int = USER,'adint'`); TAF's kind= is ignored by adread/adwrite; ifort -warn all needs adread_adwrite.F in NOOPTFILES.
- Era: 2021-05/08; PR 527 fixed. Integer tape idea only on a branch (issue 474 closed).
- Src: https://github.com/MITgcm/MITgcm/issues/526 ; https://github.com/MITgcm/MITgcm/issues/474

### SIGSEGV in adjoint only (forward fine): ctrl_map_genarr2d_ad / adthe_main_loop / profiles, "address not mapped"
- Cause: usually stack limit or memory violation inside TAF code; here: HPC-dependent (string compare after max_len_fnam change) or stack.
- Fix: `ulimit -s unlimited` / `limit stacksize unlimited` in the batch script (Discover 2024-12, Mac/ifort 2021-02, 2012-02); recompile with -devel; use `genmake2 -ncad` so each *_ad.f file is separate for reading tracebacks; test lower optimisation (-O0) and NOOPTFILES for mom_calc_visc.F.
- Era: 2012-2026.
- Src: mitgcm-support 2024-December 'Running on Discover' ; mitgcm-support 2024-May 'Segmentation Fault in MITgcm Adjoint Simulations' ; mitgcm-support 2018-June 'segmentation fault'

### Segfault in adjoint on 2 nodes with ifort+OpenMPI only (tutorial_tracer_adjsens): fixed with AUTODIFF_USE_OLDSTORE_3D
- Cause: record length of merged 3D tape arrays; option wrote one file per state variable instead.
- Fix: (historical) `#define AUTODIFF_USE_OLDSTORE_3D`. Option no longer in the tree (only tag-index); use ALLOW_AUTODIFF_WHTAPEIO + nWh instead.
- Era: 2012-02; obsolete.
- Src: mitgcm-support 2012-February 'Adjoint tutorial problems' (http://mailman.mitgcm.org/pipermail/mitgcm-support/2012-February/thread.html)

## Performance / memory / disk

### Adjoint runs on Lustre (Pleiades) slow and variable; tapes dominate I/O
- Cause: many small files, tapes and namelists on Lustre.
- Fix: `adTapeDir` (PARM05, fast local/tmp disk; ~8-15%), profilesDir, diagMdsDir, ctrl `tmpDir`; run directory on fast disk if space allows. SINGLE_DISK_IO gave ~0.1%. mdsioLocalDir is ignored by mdsio_write_field/tape when useSingleCpuIO=.TRUE. Dan Kokron's PLEIADES_LUSTRE_OPT1-3 fork (checkpoint67l_fileHash, ~33% on 3-yr adjoint) was never merged (CPP flags absent from master).
- Era: 2021-2022 (open issue).
- Src: https://github.com/MITgcm/MITgcm/issues/535

### mdsioLocalDir + control pack/unpack fails ("attempt to read/write past end of record" in mdsio_gl.f, MDSREADFIELD_3D_GL)
- Cause: ctrl_pack not cleanly parallel: rank 0 reads the global mask; with ALLOW_PACKUNPACK_METHOD2 undefined the read ignores per-rank local dirs.
- Fix: `useSingleCpuIO=.TRUE.` (slow) or `#define EXCLUDE_CTRL_PACK` and pack/unpack in a separate step after linking files into one dir. ctrl_check prints WARNING when useSingleCpuIO is used with mdsio_gl packing.
- Era: 2009, 2014; still open in design (ALLOW_PACKUNPACK_METHOD2 still exists).
- Src: mitgcm-support 2014-March 'issue using mdsioLocalDir with control package' ; mitgcm-support 2009-October 'MDSIO'

### optim.x killed "Exit code -5" on large ecco_cost files
- Cause: serial optim.x needs the whole control vector in memory (several GB for LLC90).
- Fix: run optim on a large-memory node; (Mazloff 2009, Ranger).
- Era: 2009; same mechanism today.
- Src: mitgcm-support 2009-January 'memory on optim.x => Exit code -5'

## ctrl / cost set-up

### Run stops "MDS_READ_FIELD ... Non-existing record number" xx_atemp.effective.NNNN.data when restarting (nIter0 > 0) with time-varying gentim2d controls
- Cause: records in xx_*.effective / adxx_*.effective inconsistent with start/end/diffrec from CTRL_INIT_REC; several generations of bugs.
- Fix: update pkg/ctrl (at least ctrl_get_gen_rec.F) to >= PR 934 (merged 2025-10-22); earlier: PR 380 (loop to diffrec, no meta files in REVERSE_SIMULATION, 2020-21), PR 662/664 (2022: restart from different niter0 than parent xx; adxx length; docs of .effective/.tmp in ocean_state_est chapter). Set xx_gentim2d_startdate1/2 correctly.
- Era: 2020-2025; gentim2d era (ctrl_map_ini_gentim2d.F).
- Src: mitgcm-support 2026-April 'Issue of restarting adjoint run with time-variant controls' (http://mailman.mitgcm.org/pipermail/mitgcm-support/2026-April/thread.html) ; https://github.com/MITgcm/MITgcm/pull/664 ; https://github.com/MITgcm/MITgcm/pull/380

### gentim2d sensitivities appear in adxx_<name>.effective, adxx_<name> is zero
- Cause: reverse mode accumulates in .effective; the adjoint of ctrl_map_ini_gentim2d (weights/smoothing) converts to adxx_<name> only at the very end of a completed run.
- Fix: let the run finish; test with xx_gentim2d_preproc='noscaling' (then adxx and .effective agree); .tmp files show intermediate stages.
- Era: 2020-02 (docs improved 2023 with PR 664).
- Src: https://github.com/MITgcm/MITgcm/issues/330 ; https://github.com/MITgcm/MITgcm/issues/329

### Time-varying ctrl interpolated backwards between records when useCAL=.FALSE.
- Cause: ctrl_get_gen_rec.F passes GET_PERIODIC_INTERVAL weights in swapped order (wght2 used as weight of record count0) for OBCS and gentim2d without pkg/cal.
- Fix: PR 1040 (swap to wght1). Also: gentim2d is intended to run with pkg/cal; ctrl_check.F warns without it (still in master); JM/Losch plan an xx_gentim2d_startTime/forcingCycle (not in code yet; field does not exist).
- Era: 2026-10; checkpoint68u and master up to 2026-09 affected; PR open.
- Src: https://github.com/MITgcm/MITgcm/issues/1039 ; https://github.com/MITgcm/MITgcm/pull/1040

### exf field with fldPeriod = 0 plus xx_gentim2d control: control has no effect or accumulates
- Cause: exf_set_fld.F did not re-initialise the field when fldPeriod==0 while xx_gentim2d is always added.
- Fix: PR 980 (JM): treat like vector-field interpolation, drop the `fldPeriod .NE. 0.` condition so the field is rebuilt every step with ALLOW_GENTIM2D_CONTROL; old hack fldPeriod=repeatPeriod=0. in offline_exf_seaice/input_ad/data.exf.
- Era: 2026-03 to 2026-07; fixed in master.
- Src: https://github.com/MITgcm/MITgcm/pull/980

### Control adjustment non-zero where control weight is zero; mixed scaled/unscaled adxx (genarr3d, doscaling)
- Cause: xx_gen not zeroed where wgenarr3d=0 (introduced by PR 497).
- Fix: PR 615 (set xx_gen to 0 there; log10 controls now `EXP(ln10*x)`, differences at 1e-13 level).
- Era: 2022-03; fixed.
- Src: https://github.com/MITgcm/MITgcm/pull/615

### Compile error ctrl_getobcs[nsew].F when ALLOW_OBCSN_CONTROL defined but ALLOW_OBCS_NORTH is not
- Cause: ctrl options allowed unusable flag combination.
- Fix: PR 889 (compiles; stop if the matching xx_obcs?_file is set). Also PR 815 (non-active reads when pkg/autodiff off).
- Era: 2024-03, 2024-11; fixed upstream.
- Src: https://github.com/MITgcm/MITgcm/issues/888 ; https://github.com/MITgcm/MITgcm/pull/815

### Mixing exf atmospheric ctrl flags (ALLOW_ATEMP_CONTROL, ...) with ALLOW_GENARR3D_CONTROL: exf controls silently not applied / missing xx_uwind.effective
- Cause: with ctrlUseGen=TRUE the old exf_getffields branches are skipped or read .effective files that only generic controls write.
- Fix: use generic controls for everything (atemp/aqh/precip/uwind via xx_gentim2d_file with xx_gentim2d_preproc 'rmcycle','docycle','smooth').
- Era: 2016-08 (before ECCO_CTRL_DEPRECATED code removal, PR 631).
- Src: mitgcm-support 2016-August 'bugs in ctrl pkg?' (http://mailman.mitgcm.org/pipermail/mitgcm-support/2016-August/thread.html)

### Scalar (0-D) control needed
- Cause: genarr/gentim only 2D/3D.
- Fix (Heimbach): declare the variable in CTRL (ctrl.h-like common) and in the AD tool argument list/head, supply your own adjoint output file; do not abuse genarr/gentim (alternative: use first few records of genarr).
- Era: 2023-02.
- Src: mitgcm-support 2023-February 'using (zero-d) scalars as controls with state estimate/optimisation' (http://mailman.mitgcm.org/pipermail/mitgcm-support/2023-February/thread.html)

### CPP options defined in several header files give order-dependent behaviour (ALLOW_PACKUNPACK_METHOD2, ALLOW_GENCOST_CONTRIBUTION)
- Cause: AUTODIFF_OPTIONS.h vs CTRL_OPTIONS.h, COST_OPTIONS.h vs ECCO_OPTIONS.h included in different order per file.
- Fix: PR 689/711 consolidated in reference headers; check your code_ad/*_OPTIONS.h for duplicates (JM has a script).
- Era: 2021-12 to 2023-04; fixed in reference files.
- Src: https://github.com/MITgcm/MITgcm/issues/577

### gencost time mask: "end-of-file during read" at end of forward run; cost term is zero
- Cause: gencost machinery reads the temporal mask (maskT) beyond the forward period; also gencost_datafile read precision mismatch (readBinaryPrec) gives zero cost.
- Fix: make the time-mask vector 2*N long (N records at gencost_avgperiod); put weight 1.0 on the target step, or 1/n on n steps (not auto-normalised); data.profiles empty switches profiles off.
- Era: 2016-2023 (gencost era).
- Src: mitgcm-support 2023-March 'easier ways to implement cost function' ; mitgcm-support 2023-April 'end-of-file during read file divided.ctrl'

### ADJ* files not written; or wrong snapshot sign; adEtaN zero; ADJ values not "per unit"
- Cause/Fix: (a) need `adjDumpFreq` (PARM03, with the decimal point) AND `#define ALLOW_AUTODIFF_MONITOR` in AUTODIFF_OPTIONS.h (default defined in master); with dumpAdByRec check the record counter; (b) snapshot ADJ diagnostics (frequency<0) had the wrong fill step and sign: PR 783 (2023-10); (c) adEtaN needs dummy_for_etan (PR 151, rec counter dumpAdRecEt); (d) ADJtheta includes volume and time-step weighting of the cost (sensitivity of sum, not mean); (e) new cost terms: Forget advises using adxx_* from controls instead of hand-coded ADJ files.
- Era: 2011-2024.
- Src: mitgcm-support 2023-March 'ADJ* missing' ; https://github.com/MITgcm/MITgcm/issues/774 ; https://github.com/MITgcm/MITgcm/pull/151 ; mitgcm-support 2017-May 'Values about adjoint ADJ* results' ; mitgcm-support 2024-December 'Adding Biogeochemistry (BGC) Cost/Objective Function'

### viscFacInAd has no effect
- Cause: store directives skipped mom_calc_visc recomputation (fixed PR 384, 2021); viscFacAdj only multiplies the prescribed 3D field (viscAhDfld etc., needs ALLOW_3D_VISCAH / ALLOW_3D_VISCA4 and a viscA4Dfile); `viscFacInAd` did nothing in mom_vecinv between checkpoint66e and 67v; needs AUTODIFF_ALLOW_VISCFACADJ.
- Fix: define AUTODIFF_ALLOW_VISCFACADJ, supply 3D viscosity file; applying the factor to Alin (JM suggestion) is still open (issue 222).
- Era: 2019-2026.
- Src: https://github.com/MITgcm/MITgcm/issues/222 ; mitgcm-support 2023-April 'viscA4Dfile in adjoint'

### Bottom drag / user control gives zero sensitivity
- Cause: out-of-date local adoptfile (tutorial_global_oce_optim code_ad/adjoint_hfluxm) lacks xx_bottomdrag_dummy.
- Fix: use tools/adjoint_options/adjoint_default and grep taf_ad.log for the dummy; ALLOW_BOTTOMDRAG_COST is not needed for sensitivities. Today use generic controls (xx_genarr2d_file='xx_bottomdrag').
- Era: 2011-07 (pre-generic ctrl).
- Src: mitgcm-support 2011-July 'Adjoint BOTTOMDRAG_CONTROL, adbottomdragfld : can't get non-zero sensitivity' (http://mailman.mitgcm.org/pipermail/mitgcm-support/2011-July/thread.html)

### ecco gencost: m_horflux_vol_N suffix disables the cost; objf_gencost array order; boxmean wrong weights
- Cause: string length typos (ecco_check/readparms, PR 350, 654), objf_gencost/num_gencost index order wrong when nSx*nSy>1 (PR 488), 3D boxmean numerator used current volume but denominator t=0 volume (r*, nonlinear free surface).
- Fix: PR 654, 350, 488; PR 795 (current volume; CPP ECCO_VARIABLE_AREAVOLGLOB selects time-varying areavolGlob); update to those.
- Era: 2020-2023; all fixed upstream (ECCO_VARIABLE_AREAVOLGLOB exists in ECCO_OPTIONS.h).
- Src: https://github.com/MITgcm/MITgcm/pull/795 ; https://github.com/MITgcm/MITgcm/pull/488 ; https://github.com/MITgcm/MITgcm/pull/654

### Observation records misaligned: cost_gencal assumes obs start with run start; single multi-year gencost_datafile treated as climatology
- Cause: first-year record from days since run start, not Jan 1; no flag for multi-year files.
- Fix: use yearly obs files starting 1 Jan (Mazloff patch) / gencost_startdate; issues stayed open.
- Era: 2018-12 (open in the digest).
- Src: https://github.com/MITgcm/MITgcm/issues/185 ; https://github.com/MITgcm/MITgcm/issues/183

### Cost = time-average of a field (SSH, box) in the pre-gencost code
- Fix: extend cost_accumulate_mean.F (pre-computed means) as cost_atlantic_heat does; today use gencost (avgperiod + posproc).
- Era: 2014-10.
- Src: mitgcm-support 2014-October 'How to write a cost function with time-averaged SSH?'

## Gradient check, TLM, profiles

### grdchk summary RMS of ratios bad (1e-3 ... 1e-1) although AD is right
- Cause: grdchk_eps too large/small for the experiment; inexact AD (lab_sea seaice growth adx; ggl90 mxlMaxFlag).
- Fix: tune `grdchk_eps` (e.g. 1e-2 -> 1e-4 gave RMS 5.9e-2 -> 2.8e-6 in global_ocean.90x40x15); use mxlMaxFlag=1/3 with ggl90; check with useCentralDiff.
- Era: 2021-08 (issue 518, closed 2026-01).
- Src: https://github.com/MITgcm/MITgcm/issues/518

### Gradient check fails with r* (time-varying hFacC) on the 2nd perturbation / central difference
- Cause: forward perturbation runs not re-initialised: recip_hFac wrong after INITIALISE_VARIA.
- Fix: PR 438 (fix in ini_nlfs_vars.F; do not just add INI_MASKS_ETC in grdchk_main.F, which breaks OBCS masks).
- Era: 2021-03; fixed.
- Src: https://github.com/MITgcm/MITgcm/pull/438

### grdchk_main: recglo ignored; grdchk output missing with useSingleCpuIO; gradient index depends on data.ctrl order
- Cause/Fix: recglo fixed (PR 528, 2021-09; static controls with recglo>1 give zero); AD/TL MDS output not written with useSingleCpuIO (issue 874, 2024-10: condition removed for MPI-safe case); grdchk works on one tile only; select the control by name (grdchkvarname in GRDCHK.h) rather than by index (index depends on order in data.ctrl).
- Era: 2021-2024.
- Src: https://github.com/MITgcm/MITgcm/pull/528 ; https://github.com/MITgcm/MITgcm/issues/874 ; https://github.com/MITgcm/MITgcm/issues/786

### TLM monitor `g__dynstat_*` all zeros; WARNING CTRL_INIT_CTRLVAR "could not find g_xx_theta ... initialize this file with all zeros"
- Cause: TLM is forced only by g_xx_* perturbation files; in verification only grdchk_main sets individual points to 1; if the grdchk control is a passive tracer (xx_ptr2) dynamics stay zero.
- Fix: perturb g_xx_* (or xx_theta in data.grdchk) to see dynamics; TL output g_fc (g_objc_state_final with adjoint_state_final adoptfile).
- Era: 2018-02, 2024-10.
- Src: https://github.com/MITgcm/MITgcm/issues/874 ; mitgcm-support 2018-February 'MITgcm tangent linear' (http://mailman.mitgcm.org/pipermail/mitgcm-support/2018-February/thread.html)

### profiles: gradient check wrong with NetCDF model equivalents (profilesDoNcOutput); segfault in profiles_interp_ad
- Cause: varid(10) too small in profiles_init_ncfiles.F and netcdf file not re-opened for read (issue 884, PR 887); segfault when a profile has more levels than NLEVELMAX.
- Fix: PR 887; profiles_init_fixed.F now stops when levels exceed NLEVELMAX (set NLEVELMAX in PROFILES_SIZE.h); TLM profiles init PR 873.
- Era: 2024-11, 2026-04 (cp69m).
- Src: https://github.com/MITgcm/MITgcm/issues/884 ; https://github.com/MITgcm/MITgcm/issues/985 ; https://github.com/MITgcm/MITgcm/issues/717

### Salt (or other) sensitivities all zero when nSx > 1 (e.g. nSx=6,nPx=4) in global_oce_cs32 ad.sens, fine with nSx=1,nPx=6
- Cause: not diagnosed in thread (adjoint with several tiles per process; note TAF key-index/tile issues were fixed later, e.g. PR 488 objf_gencost order).
- Fix: use nSx=1 with more processes to confirm; re-test on current master.
- Era: 2018-06 (checkpoint67b); unresolved.
- Src: mitgcm-support 2018-June "Can't replicate verification/global_oce_cs32/ad.sens" (http://mailman.mitgcm.org/pipermail/mitgcm-support/2018-June/thread.html)

## Optimisation (optim / lsopt / optim_m1qn3)

### optim_m1qn3: omode=5/4/6, terminates early or "nsim/niter changed"; line search "mlis3" stalls
- Cause: omode=1 means epsg reached (default epsg 1e-6), omode=6 = no further improvement i.e. gradient not accurate enough for epsg, omode=7 non-descent. nsim/niter are overwritten on exit and saved to OPWARM, so a script that blindly warm-restarts stops at once.
- Fix: stop the driver loop when m1qn3_output.txt says "output mode"; cold restart by removing OPWARM* or `&M1QN3 coldstart=.TRUE.` (forgets Hessian approx.; fmin must be below the last cost or leave unset and use fminFrac); check the gradient with grdchk first; bad mlis3 usually means wrong gradient.
- Era: 2014-2020.
- Src: mitgcm-support 2020-October 'optim_m1qn3: maximum iterations reached?' ; mitgcm-support 2019-February 'optim_m1qn3 line search getting stuck' ; mitgcm-support 2014-October 'coldstart option for m1qn3'

### optim.x EOF / NaN on Cray (ARCHER): "End of file" in optim_readdata.f
- Cause: ecco_cost/ecco_ctrl files written with pkg/ctrl, not mdsio, so -D_BYTESWAPIO is not in effect; endianness differs from optim.x; initialisation bug gave NaNs.
- Fix: same byte-order compiler flag (e.g. `-h byteswapio`) for mitgcmuv_ad and optim.x; use current optim_m1qn3 (NaN init bug fixed 2019-04); myX/YglobalLo=1 for column models (PR 338).
- Era: 2019-04.
- Src: mitgcm-support 2019-April 'optim_m1qn3 on ARCHER supercomputer' (http://mailman.mitgcm.org/pipermail/mitgcm-support/2019-April/thread.html); https://github.com/MITgcm/MITgcm/pull/338

### Bundled lsopt/optim: build errors, "ifail = 4 the search direction is not a descent one"
- Cause: optim/lsopt are unmaintained (-lblas1, makedepend, missing example); ifail=4 on the first line search is normal output in the tutorial, but cost may not drop if data.optim optimcycle was not incremented.
- Fix: `./testreport -t tutorial_global_oce_optim -adm -ncad -devel; cd ../lsopt && make; cd ../optim && make` (PR 206); Losch recommends optim_m1qn3 (github.com/mjlosch/optim_m1qn3); namelist &M1QN3 needed in data.optim for m1qn3.
- Era: 2018-2021 (README updated 2021-12).
- Src: mitgcm-support 2018-April 'tutorial_global_oce_optim optimisation failed' ; https://github.com/MITgcm/MITgcm/pull/206

## Tapenade (AD tool for pkg/tapenade, genmake2 -tap)

### genmake2: Tapenade run fails "flow_tap" not found when the build dir is outside the MITgcm tree
- Cause: PR 949 hard-wired a relative path `../../../tools/TAP_support/flow_tap` in tools/adjoint_options/adjoint_tap.
- Fix: PR 973: reference it as `$(TOOLSDIR)` (escaped so Makefile, not bash, expands it; same for diffsizes.F90). Updated genmake2 and adjoint_tap needed.
- Era: 2026-02; fixed upstream.
- Src: https://github.com/MITgcm/MITgcm/issues/972 ; https://github.com/MITgcm/MITgcm/pull/973

### Tapenade Makefile rebuilds everything / race / unusable on macOS (case-insensitive FS)
- Cause: touch hacks on *_b.f, adStack.c/adBinomial.c recompiled each time.
- Fix: PR 935: generated files concatenated into adj/tlm/_tap_all.f (default) or `-ncad` appends forward source to each Tapenade file via append_fwd_src; symlinks to ADFirstAidKit files in pkg/tapenade; case-insensitive fix; docker recipe + optfile darwin_arm64_gfortran for macOS.
- Era: 2025-09 (checkpoint after that).
- Src: https://github.com/MITgcm/MITgcm/pull/935 ; https://github.com/MITgcm/MITgcm/issues/735

### Tapenade: "System: Not a tree operator" for every file / FortranParser crash
- Cause: fortranParser crashes at start on a system with outdated glibc.
- Fix: build FortranParser from the new tarball for that system (Hascoet); test `tapenade program.f` on hello-world first.
- Era: 2024-10.
- Src: mitgcm-support 2024-October 'Compiling Problems with Tapenade' (http://mailman.mitgcm.org/pipermail/mitgcm-support/2024-October/thread.html)

### Tapenade streamice test/adjoint 50x slower or wrong; fixed-point; PETSC compile errors
- Cause: old Tapenade binary (2023-05) vs 2025-07-29+ (144 min vs 3 min); STREAMICE needs short-circuited CG solve (flow_tap + stubs_tap_adj.F) and FP-LOOP fixed-point; single-solve cost on surface velocity needs one extra streamice_vel_phistage call.
- Fix: update Tapenade; PR 927 (fixed point), 974 (ALLOW_PETSC compile), 1004 (extra phistage call).
- Era: 2025-07 to 2026-07; fixed upstream.
- Src: https://github.com/MITgcm/MITgcm/pull/927 ; https://github.com/MITgcm/MITgcm/pull/974 ; https://github.com/MITgcm/MITgcm/pull/1004

### Tapenade/TAF compile: CTRL_SIZE.h / CTRL_DUMMY.h included twice (ALLOW_GENTIM2D_CONTROL), wrong _RL vs _RS in exch_tap_b.F
- Cause: header include left over from PR 740; EXCH1_RS_B arguments declared _RL (breaks `-ur4`).
- Fix: PR 843, PR 737. Also TLM with gentim2d and profiles missing init: issues 717/873.
- Era: 2023-05 to 2024-06; fixed.
- Src: https://github.com/MITgcm/MITgcm/pull/843 ; https://github.com/MITgcm/MITgcm/pull/737

### OpenAD is being removed (issue 971); build questions
- Cause: maintenance burden; Tapenade covers the open-source route.
- Fix: C.I. step done (PR 1000, 2026-07); a tagged checkpoint before the removal keeps OpenAD; docs mention it. Do not start new work on OpenAD.
- Era: 2026-01 to 2026-07.
- Src: https://github.com/MITgcm/MITgcm/issues/971

## ECCO-specific notes

### ECCO v4 (llc90) "rStarFac too SMALL" after altering wind forcing
- Cause: sea ice pushed into corners (r* with real freshwater flux).
- Fix: `SEAICEpressReplFac = 0.` and no sea-ice diffusion (ECCO v4r2 data.seaice has SEAICEdiffKh*=400); tested config: verification_other/global_oce_llc90 (input.ecmwf data.seaice update). Did not cure all cases.
- Era: 2024-02 (v4r2 namelists).
- Src: mitgcm-support 2024-February 'ECCOv4 Configurations: R-Star Coordinate Error with Modified Wind Forcing' (http://mailman.mitgcm.org/pipermail/mitgcm-support/2024-February/thread.html)

### ecco adjoint verification set-ups (what to run)
- Fix: verification_other/global_oce_cs32 (input_ad, input_ad.sens) and global_oce_llc90 (.ecco_v4, .ecmwf, .core2) are the daily ECCO-like AD tests; `../verification/testreport -t global_oce_cs32 -adm`; a run that says "adm_boxmean_theta ... File does not exist" was prepared without autodiff (wrong prepare_run/pkg).
- Era: 2020-09; set-ups renamed over time.
- Src: mitgcm-support 2020-September 'verification_other case global_oce_cs32 fails at runtime in adjoint mode' (http://mailman.mitgcm.org/pipermail/mitgcm-support/2020-September/thread.html)

### Monthly/averaged adjoint output with ECCO v4
- Fix: adxx_ are control-period snapshots (xx_gentim2d_period, 2 weeks in v4); change the period or add zero-filled controls; ADJ* diagnostics via pkg/diagnostics for state variables (snapshots: frequency<0, see ADJ entry).
- Era: 2019-02.
- Src: mitgcm-support 2019-February 'Monthly averaged output of adjoint sensitivity with ECCOv4' (http://mailman.mitgcm.org/pipermail/mitgcm-support/2019-February/thread.html)
