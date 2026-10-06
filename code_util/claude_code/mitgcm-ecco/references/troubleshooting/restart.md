# Troubleshooting: pickups, restarts, nIter0, checkpoints
Distilled from answered mitgcm-support threads (2003-2026) and MITgcm GitHub issues/PRs on restart/pickup problems; names verified against origin/master (Oct 2026) unless Era says otherwise.
Thread URLs are month indexes: http://mailman.mitgcm.org/pipermail/mitgcm-support/YYYY-Month/thread.html

## Choosing / naming the pickup (nIter0, pickupSuff, ckptA/B)
### Restart fails: "MDS_READ_FIELD: File does not exist" / reads wrong pickup / only pickup.ckptA exists
- Cause: pickup file name must carry the iteration: pickup.<nIter0 as 10 digits>.<tile>.data+.meta (+ pickup_<pkg>.* for each package). Rolling pickups (chkptFreq>0) are pickup.ckptA / ckptB, written alternately starting with ckptA after every (re)start, so ckptA is NOT always the latest; each is overwritten. The other half of the "file does not exist" message is in STDOUT.0000, not STDERR.
- Fix: read the iteration in pickup.ckptX.meta (`timeStepNumber`), then either rename ALL pickup*.ckptX.* files (data, meta, every pkg and every tile) to that iteration and set `nIter0` (PARM03), or keep names and set `pickupSuff='ckptA'` (remove it again for the next restart). Better: `chkptFreq=0.` and `pChkptFreq` = multiple of run length (permanent pickups named by iteration), and `writePickupAtEnd=.TRUE.`. STDOUT prints a grep-able line every time a pickup is written. Specify EITHER `nIter0` OR `startTime` (startTime = nIter0*deltaTClock) not both inconsistently.
- Era: 2004-2022 (all current).
- Src: mitgcm-support 2011-January 'how to restart simulation with only pickup.ckptA file' http://mailman.mitgcm.org/pipermail/mitgcm-support/2011-January/thread.html ; 2004-August 'Restarting the MITgcm' http://mailman.mitgcm.org/pipermail/mitgcm-support/2004-August/thread.html ; 2022-August 'rolling pickups: nchecklev value on restarting model' http://mailman.mitgcm.org/pipermail/mitgcm-support/2022-August/thread.html ; 2010-June 'A very strange problem'/'Restart MITGCM model' http://mailman.mitgcm.org/pipermail/mitgcm-support/2010-June/thread.html ; 2012-December 'restart MITGCM' (pickupSuff='ckptA') http://mailman.mitgcm.org/pipermail/mitgcm-support/2012-December/thread.html

### Restarted run deviates from the continuous run right after restart (forcing/phase wrong, SST jump)
- Cause: inconsistent time specification: `baseTime`/`startTime` set (e.g. =0) while nIter0>0; or deltaT changed between runs (model time = baseTime + iter*deltaTClock); or wrong nIter0 (1569600 x 3600 s = 179 y, not 50 y).
- Fix: when deltaT is unchanged, change ONLY nIter0 (remove baseTime/startTime from the namelist). Output times, forcing records and calendar all derive from model time.
- Era: 2018-2020.
- Src: mitgcm-support 2018-December 'restart simulations using pickup files' http://mailman.mitgcm.org/pipermail/mitgcm-support/2018-December/thread.html ; 2020-February 'SST problem when restarting a run' http://mailman.mitgcm.org/pipermail/mitgcm-support/2020-February/thread.html

### Restart with a different deltaT (pickup named for old iteration, time axis wrong, EXF reads wrong forcing)
- Cause: model time = baseTime + nIter0*deltaTClock; by default baseTime=0 so a changed deltaT shifts the model clock (netcdf time axis = iter*new dt). Default not changed because it would break existing set-ups; deltaTClock is not stored in the pickup.
- Fix: (a) rename pickup to nIter0' = nIter0_old*dt_old/dt_new and use nIter0' (works, EXF follows time), or (b) keep nIter0 and set `baseTime = nIter0*(dt_old - dt_new)`, or (c) set both nIter0 and startTime. Use synchronous deltaT only (not deltaTClock != deltaT). Rolling/changed pChkptFreq phase: output times are multiples of freq on the model clock.
- Era: 2005-2013 (baseTime added 2005; behaviour unchanged).
- Src: mitgcm-support 2013-January 'restart from pickup files but change deltaT problem with EXF' http://mailman.mitgcm.org/pipermail/mitgcm-support/2013-January/thread.html ; 2012-November 'Time variable and changing deltaT' http://mailman.mitgcm.org/pipermail/mitgcm-support/2012-November/thread.html ; 2005-May 'Restart after changing deltatT' http://mailman.mitgcm.org/pipermail/mitgcm-support/2005-May/thread.html ; 2009-April 'Restarting with multiple time steps' http://mailman.mitgcm.org/pipermail/mitgcm-support/2009-April/thread.html

### Pickup iteration after switching calendar (360-day "model" vs 365/Gregorian) or calendarDumps
- Cause: iteration number is tied to model time; "TheCalendar='model'" is 360-day; with `calendarDumps=.TRUE.` pickups/dumps fall at calendar month ends and the file name uses the actual iteration (e.g. 26784 for 31 d with dt=100), not freq/deltaT. Changing calendar mid-run needs pickup edits.
- Fix: rename pickup to the iteration matching the new calendar (e.g. 39500 y x 360 d, dt=14400 -> nIter0=85320000), or keep old file with `pickupSuff='0086505000'` plus the new nIter0 and then remove pickupSuff. In data.exf keep start dates fixed across restarts; use `repeatPeriod` (data.exf) for periodic repeating forcing; one-file-per-year naming only works up to year 9999.
- Era: 2021-2024.
- Src: mitgcm-support 2024-January '{Disarmed} Re: Re: MITgcm-support Digest, Vol 247, Issue 3' http://mailman.mitgcm.org/pipermail/mitgcm-support/2024-January/thread.html ; 2024-January '... Issue 1' ; 2021-May 'calendarDumps and pickup file names' http://mailman.mitgcm.org/pipermail/mitgcm-support/2021-May/thread.html

## Exact (bit-for-bit) restart; results differ after restart
### Restarted run differs from non-restart run (T,S,U diverge after the first step)
- Cause: MITgcm is designed to restart exactly (nightly "2+2=4" tests: tst_2+2 in testreport; pickup files compared), but not with every compiler/option: aggressive optimisation or higher internal precision than the 64-bit pickup (x87/vector), untested option combinations, or NetCDF (MNC) pickup. Multi-processor runs are not bit-reproducible run-to-run, but a restart from the same pickup with the same layout should match for at least a few steps.
- Fix: build with `genmake2 -ieee` / exact IEEE flags or -O0 and repeat 2+2 test (`testreport` option for restart experiments; scripts fixed for POSIX sed in PR #371 so they work on macOS); use plain MDS pickups (`pickup_write_mnc=.FALSE.`, `pickup_read_mnc=.FALSE.`); check pkgs used are tested (`GLOBAL_SUM_TILE` vs GLOBAL_SUM etc.). Restart from last pickup before a blow-up often runs further (round-off chaos), not a bug by itself.
- Era: 2014-2021 (PR #371 Sep 2020).
- Src: mitgcm-support 2021-October 'Deviated results after restart' http://mailman.mitgcm.org/pipermail/mitgcm-support/2021-October/thread.html ; 2021-April 'pickup question' http://mailman.mitgcm.org/pipermail/mitgcm-support/2021-April/thread.html ; 2014-December 'Different output shortly after restart' http://mailman.mitgcm.org/pipermail/mitgcm-support/2014-December/thread.html ; 2018-June 'Reproducibility of blowup' http://mailman.mitgcm.org/pipermail/mitgcm-support/2018-June/thread.html (in instability.md sources) ; https://github.com/MITgcm/MITgcm/pull/371

### Results change when the number/size of tiles changes between pickups
- Cause: global sums (cg2d, monitors) change order; at the truncation level, growing with time/resolution. No 'bug'.
- Fix: keep nSx*nPx (and tile size) unchanged when you only change processors; compile with `#define GLOBAL_SUM_SEND_RECV` (CPP_EEOPTIONS.h) and -O0 to get identical results for the same tile count; `#define CG2D_SINGLECPU_SUM` (CPP_EEOPTIONS.h; not NH) when tile size changes. Editing SIZE.h requires recompiling; "No. of processes not equal to nPx*nPy" means your mpirun count mismatches SIZE.h.
- Era: 2008-2010.
- Src: mitgcm-support 2008-March 'restart from pickup files (cube-sphere)' http://mailman.mitgcm.org/pipermail/mitgcm-support/2008-March/thread.html ; 2009-December 'pickup after changing number of processors' http://mailman.mitgcm.org/pipermail/mitgcm-support/2009-December/thread.html

## NetCDF (MNC) pickups
### WARNING MNC_READPARMS: incomplete MNC pickup files implementation / MNC_CW_RL_R: cannot open either a per-face or a per-tile file pickup...t001.nc / variable 'RC' already defined in file 'grid.t001.nc'
- Cause: MNC pickups are poorly maintained: one file per MPI process (t001...), cannot change domain decomposition, cannot read a global NetCDF, cannot be inspected easily, restart not exact in some set-ups, old versions had a name mismatch (pickup.0000nnnnnn.0000.000001.nc vs read name) and "grid.t001.nc already defined" errors arise because MNC never overwrites existing .nc files (also stale grid.*.nc/state.*.nc in the run dir). Warning is still printed by mnc_readparms.F.
- Fix: use MDS pickups (`pickup_write_mnc=.FALSE.`, `pickup_read_mnc=.FALSE.` in data.mnc, like all verification experiments), with `useSingleCpuIO=.TRUE.` for one global file per field usable with any tiling; MNC can still be used for other output. If you must: copy pickup.*.t???.nc into the run dir, `mnc_use_indir=.FALSE.`, remove existing *.nc (or `mnc_use_outdir=.TRUE.` with new directory per run), keep the same number of processes. MNC pickup also fails when a field (e.g. phi_nh) is missing.
- Era: 2005-2019 (warning current).
- Src: mitgcm-support 2019-August 'MNC pickup' http://mailman.mitgcm.org/pipermail/mitgcm-support/2019-August/thread.html ; 2017-July 'using netcdf pickup files' http://mailman.mitgcm.org/pipermail/mitgcm-support/2017-July/thread.html ; 2014-February 'problem with MNC pickup/restart' http://mailman.mitgcm.org/pipermail/mitgcm-support/2014-February/thread.html ; 2012-August 'Problem using pickup with mnc' http://mailman.mitgcm.org/pipermail/mitgcm-support/2012-August/thread.html ; 2011-July 'gluemnc and nc file restart' http://mailman.mitgcm.org/pipermail/mitgcm-support/2011-July/thread.html ; 2015-January 'Pickup Output Times' http://mailman.mitgcm.org/pipermail/mitgcm-support/2015-January/thread.html

## Changing decomposition, resolution or initial fields
### Restart on a different number of MPI processes (or different tile size)
- Cause: tiled pickup files (one per tile) are tied to the tiling.
- Fix: write global pickups (`useSingleCpuIO=.TRUE.` preferred; `globalFiles=.TRUE.` is "not safe in MPI runs" and failed on Cray (lib-5058 'read system call read less data', cured by useSingleCPUIO)); a single pickup.<iter>.data (no .001.001) is read by any decomposition (mdsreadfield tries FILENAME, FILENAME.data, FILENAME.00?.00?.data). Combine tiled files with rdmds+fwrite, or restart from init-files. NOT possible with MNC.
- Era: 2009-2017 (current).
- Src: mitgcm-support 2014-June 'using low resolution spin up as input to a high resolution model' http://mailman.mitgcm.org/pipermail/mitgcm-support/2014-June/thread.html ; 2017-April 'error while writing pickup files with Cray compilers' http://mailman.mitgcm.org/pipermail/mitgcm-support/2017-April/thread.html ; 2004-August 'Restarting the MITgcm' (Heimbach: mdsreadfield search order) http://mailman.mitgcm.org/pipermail/mitgcm-support/2004-August/thread.html ; 2009-December 'pickup after changing number of processors'

### Start a new/finer run from a coarser run's state, or replace theta/salt/eta/u/v at restart
- Cause: interpolated pickups carry inconsistent tendencies (guNm1, gtNm1...), unbalanced fields; often unstable (and pickups fail if netcdf).
- Fix: do NOT interpolate the pickup. Cold-start from initial-condition files in PARM05: `hydrogThetaFile`, `hydrogSaltFile`, `uVelInitFile`, `vVelInitFile`, `pSurfInitFile` (readBinaryPrec, ieee-be; w from continuity; tendency terms recomputed), plus seaice `HeffFile`, `AreaFile`, `HsnowFile`, `HsaltFile`, `UiceFile`, `ViceFile` in data.seaice; or use time-averaged output (smoother -> more stable); nearest-neighbour interpolation; short deltaT initially, higher viscosity. Equivalent: edit theta/salt/eta inside the pickup (tendencies stay old). Extract IC from MNC pickups with rdmnc.m. Do not expect a continuously running model to re-initialise mid-run: restart.
- Era: 2009-2024.
- Src: mitgcm-support 2014-June 'using low resolution spin up as input to a high resolution model' (Menemenlis list of PARM05 + seaice files) http://mailman.mitgcm.org/pipermail/mitgcm-support/2014-June/thread.html ; 2024-November 'Restart a simulation using finer spatial simulation (interpolation)' http://mailman.mitgcm.org/pipermail/mitgcm-support/2024-November/thread.html ; 2023-November 'How to use pickups' (Menemenlis) http://mailman.mitgcm.org/pipermail/mitgcm-support/2023-November/thread.html ; 2020-April 'How to start from given initial conditions without restarting the simulation' http://mailman.mitgcm.org/pipermail/mitgcm-support/2020-April/thread.html ; 2009-August 'Read theta(salt)&eta from outer files when the model restarted' http://mailman.mitgcm.org/pipermail/mitgcm-support/2009-August/thread.html

### Restart with changed physics: seaice initial values ignored; nonhydrostatic from hydrostatic pickup; EOS change
- Cause: on a restart, state variables come from the pickup (SEAICE_initialHEFF, HeffFile unused); NH needs `phi_nh`/gW, MDJWF and JMD95P need hydrostatic geopotential phiHyd (pickup field), so a pickup written with another EOS/mode lacks them.
- Fix: error "field phi-NHyd is missing": add a zero field to the pickup and list it in the .meta fldList (or `pickupStrictlyMatch=.FALSE.`), restarting NH from hydrostatic pickup is OK; for EOS change append phiHyd record to pickup.data and edit pickup.meta fldList (use order found in a pickup that has it). NetCDF pickups cannot do this.
- Era: 2009-2013.
- Src: mitgcm-support 2009-March 'restart from a pickup file' http://mailman.mitgcm.org/pipermail/mitgcm-support/2009-March/thread.html ; 2013-January "restart from pickup when using 'MDJWF' eos" http://mailman.mitgcm.org/pipermail/mitgcm-support/2013-January/thread.html ; 2009-June '(no subject)' (SEAICE_initialHEFF) http://mailman.mitgcm.org/pipermail/mitgcm-support/2009-June/thread.html

## Old/mismatched pickups, .meta handling
### Reading an old pickup with new code (garbled fields, CALC_R_STAR at first step, salt looks like temperature)
- Cause: pickup contents/order changed over time (c54 2004: new order; c76-78 'cube' era; seaice format changes); if the .meta is missing the code assumes the NEW order (nbFields=-1) -> wrong fields; MDS_READ_META warns "no field-list found ... try to read pickup as currently written". Corrupted/foreign-endian pickups also give absurd EtaH (checked via rdmds EtaH).
- Fix: keep the .meta; `pickupStrictlyMatch=.FALSE.` (PARM03) to continue despite missing/extra fields; for pickups older than checkpoint54 (2004-07-02) add `usePickupBeforeC54=.TRUE.` (PARM01), turn it off for new pickups; support lower bound checkpoint35 (usePickupBeforeC35 is retired and now traps with an error). Without a meta use a minimal .meta from a verification experiment. Check pickup with `rdmds('pickup',iter)`, EtaH is last record; try 'ieee-be'/'ieee-le'. Menemenlis' mk_pickup78.m documents pickup layouts and mapping from older formats. Reading old netcdf checkpoint formats is not supported.
- Era: 2004-2015; flags exist in ini_parms.F (usePickupBeforeC54 present, C35 retired).
- Src: mitgcm-support 2008-September 'clever pickup too clever for me' http://mailman.mitgcm.org/pipermail/mitgcm-support/2008-September/thread.html ; 2015-January 'Reading an old pickup file' http://mailman.mitgcm.org/pipermail/mitgcm-support/2015-January/thread.html ; 2004-July 'pickup c52 to c54' http://mailman.mitgcm.org/pipermail/mitgcm-support/2004-July/thread.html ; 2014-November 'old pickups, new model version -- fails on first timestep with STOP in CALC_R_STAR' http://mailman.mitgcm.org/pipermail/mitgcm-support/2014-November/thread.html ; 2003-November 'Re: c51 (fwd)' (new->old not supported; pickup_cd changed) http://mailman.mitgcm.org/pipermail/mitgcm-support/2003-November/thread.html

### Pickup read/write fixed upstream: missing-field list index typo, thsice MNC pickup, cheapaml passive tracer
- Cause: typos in the missFldList index in 7 pkg *_read_pickup.F, thsice MNC pickup and cheapaml `useCheapTracer` pickup (gCheaptracerm) partially broken.
- Fix: PR #1022 (closed Aug 2026; check tag-index for merge); update, or copy the corrected files.
- Era: through checkpoint69x; fixed in master 2026-08.
- Src: https://github.com/MITgcm/MITgcm/pull/1022

### New field added to the pickup (e.g. QH stagger tendency, NH) breaks old restart
- Cause: new code reads a field absent from old pickup.
- Fix: `pickupStrictlyMatch=.FALSE.` allows restart (not exact); QH+staggerTimeStep field under `ALLOW_QHYD_STAGGER_TS` (PR #433, Feb 2021).
- Era: 2021.
- Src: https://github.com/MITgcm/MITgcm/pull/433

## Reading, editing and writing pickups by hand
### Edit a pickup in Matlab/Python (which record is which field; "floating invalid" after writing)
- Cause: pickup.data is a stack of 2-D (nx,ny) records, listed in pickup.meta `fldList` (3-D fields take Nr records each); always float64; `readBinaryPrec` does NOT apply to pickups (Matlab needs 'real*8', otherwise NaN).
- Fix: `p = rdmds('pickup',iter)` gives (nx,ny,nrec); slice via fldList and Nr (e.g. last record = EtaH); write back with fwrite(fid,p,'real*8') big-endian ('b') or gcmfaces writebin; keep the meta consistent; pkg/flt pickups (pickup_flt.*) are per-tile (no singleCpuIO) with first 9-element record special. For the FLT `tend` (max integration time) rewrite field 9 of arr; the pickup of Prather tracers includes pickup_somTRAC01 (1st/2nd moments) which can be left unchanged when editing means (switching off scheme 80 drops the file). Invalid values written into hand-made pickups produce "floating invalid" at the first step (traceback with -g, debugLevel=4).
- Era: 2014-2020.
- Src: mitgcm-support 2016-September 'Variables in pickup files' (etaN vs etaH, gXNm1 definitions) http://mailman.mitgcm.org/pipermail/mitgcm-support/2016-September/thread.html ; 2016-January 'R/W pickup' http://mailman.mitgcm.org/pipermail/mitgcm-support/2016-January/thread.html ; 2014-September 'matlab script for reading pickup files' http://mailman.mitgcm.org/pipermail/mitgcm-support/2014-September/thread.html ; 2015-August 'MITgcm reading binary restarts' http://mailman.mitgcm.org/pipermail/mitgcm-support/2015-August/thread.html ; 2016-August 'floating invalid with new pickup' http://mailman.mitgcm.org/pipermail/mitgcm-support/2016-August/thread.html ; 2017-May 'pickup file for ptracer' http://mailman.mitgcm.org/pipermail/mitgcm-support/2017-May/thread.html ; 2020-May 'changing one entry in pickup files' http://mailman.mitgcm.org/pipermail/mitgcm-support/2020-May/thread.html

### Nonlinear free surface restart: do not shift eta; Eta fields consistent
- Cause: restart needs etaN, etaH and dEtaHdt consistently; subtracting the mean SSH from etaN only (not etaH) gave divergence; only perfect with consistent triple.
- Fix: if you must shift SSH subtract the same constant from etaN and etaH; otherwise leave as is (volume conserved; OBCS may drift).
- Era: 2005.
- Src: mitgcm-support 2005-April 'Re: pickup problem' http://mailman.mitgcm.org/pipermail/mitgcm-support/2005-April/thread.html

## Forcing / boundary conditions at restart
### Restart and forcing files: do I need to change OBCS/forcing file names or start records?
- Cause/Answer: the model derives the forcing/OBCS record from model time since 0 (not nIter0), interpolating linearly between records; with pkg/exf the records are located from `*startdate1/2` + `*period` relative to data.cal startDate. Start dates must NOT be changed between restarts. startDate_1 in data.cal corresponds to model time 0 (iteration 0), not to nIter0.
- Fix: keep data.cal, data.exf and forcing files identical across restarts; only nIter0 changes. Errors "Non-existing record number" after a restart mean the forcing/xx files are shorter than model time needs (N+1 rule), not a pickup problem; changing `externForcingPeriod` w/o Cycle also does it. Link shifted year files (e.g. EOG_rain_2009 -> 2010) when a 3-hr offset needs the previous year file.
- Era: 2015-2024.
- Src: mitgcm-support 2015-December 'How to set the OB files when restart the simulation' http://mailman.mitgcm.org/pipermail/mitgcm-support/2015-December/thread.html ; 2015-March 'Pickup model in the middle of forcing file' http://mailman.mitgcm.org/pipermail/mitgcm-support/2015-March/thread.html ; 2016-September 'About pickup' http://mailman.mitgcm.org/pipermail/mitgcm-support/2016-September/thread.html ; 2016-December 'build errors for llc_1080' http://mailman.mitgcm.org/pipermail/mitgcm-support/2016-December/thread.html ; 2024-January 'Re: Re: MITgcm-support Digest, Vol 247, Issue 1' http://mailman.mitgcm.org/pipermail/mitgcm-support/2024-January/thread.html

### Starting tracers / packages from scratch on top of a physics pickup (ptracers, DIC, floats)
- Cause: package pickups (pickup_ptracers, pickup_dic_co2atm, pickup_flt, pickup_seaice, pickup_cd, ...) are expected for nIter0>0.
- Fix: ptracers: `PTRACERS_Iter0 = nIter0` (data.ptracers) with `PTRACERS_initialFile` (or PTRACERS_ref -> zero fields) to start tracers while dynamics come from the pickup; `PTRACERS_Iter0 > nIter0` delays tracers; `PTRACERS_Iter0 < nIter0` (cfc_example style) requires pickup_ptracers.<Iter0> to exist. `PTRACERS_EvPrRn`, DIC_bounds optional; atmospheric-CO2 pickup (dic_int1=3) not in any verification experiment (#191 open). Floats: no FLT_Iter0 (idea 2012): float initial positions are only read at nIter0=0; options: nIter0=0 with `pickupSuff`, build your own pickup_flt, or modify flt_init_varia.F. diagnostics pickup (`diag_pickup_read`) never completed: segfault.
- Era: 2007-2025.
- Src: mitgcm-support 2020-November 'Restart simulation with biogeochemical packages' http://mailman.mitgcm.org/pipermail/mitgcm-support/2020-November/thread.html ; 2007-April 'offline simulation' http://mailman.mitgcm.org/pipermail/mitgcm-support/2007-April/thread.html ; 2012-April 'MITgcm floats.................' http://mailman.mitgcm.org/pipermail/mitgcm-support/2012-April/thread.html ; 2025-May 'diagnostics read and write pickups' http://mailman.mitgcm.org/pipermail/mitgcm-support/2025-May/thread.html ; https://github.com/MITgcm/MITgcm/issues/191

### Orlanski open boundary: W or boundary phase speed zero after restart
- Cause: (2008) Orlanski state was not checkpointed.
- Fix: now implemented: pkg/obcs obcs_read_pickup.F states only Orlanski (and Stevens) need pickup files; older checkpoints restart with zero phase speed.
- Era: 2008 (fixed upstream since; verify in your checkpoint).
- Src: mitgcm-support 2008-February 'Non-hydrostatic Orlanski pickup and W velocity' http://mailman.mitgcm.org/pipermail/mitgcm-support/2008-February/thread.html

## Adjoint / control restarts
### Adjoint (or forward) restart with time-varying gentim2d controls: "Non-existing record number" in xx_*.effective.*.data / "attempting to access non-existing record in xx_precip"
- Cause: ECCO/gentim2d control records indexed assuming start at record 1 (nIter0 = parent start): startrec/endrec/diffrec mismatch, adxx files of wrong length; fixed in steps: #380 (2021), #662/#664 (Oct 2022, forward reads right xx records; adxx length), #934 (merged 2025-10-22: further start/restart fixes, ctrl_get_gen_rec.F).
- Fix: use current code (replace ctrl_get_gen_rec.F etc. with master); xx_gentim2d_startdate1/2 and period must match parent; old workaround was to rename pickups to iter 1 (ECCOv4 does start at nIter0=1). Breaking: results change vs. code between #380 and #664 if xx nonzero.
- Era: 2020-2026.
- Src: mitgcm-support 2026-April 'Issue of restarting adjoint run with time-variant controls' http://mailman.mitgcm.org/pipermail/mitgcm-support/2026-April/thread.html ; https://github.com/MITgcm/MITgcm/issues/377 ; https://github.com/MITgcm/MITgcm/pull/664 ; https://github.com/MITgcm/MITgcm/pull/934

### Divided adjoint (restart of adjoint) and pkg/profiles segfault
- Cause/Answer: restarting the adjoint in pieces is the "divided adjoint" (DIVA, manual adjoint chapter); pkg/profiles adjoint (checkpoint 69m) can segfault when profiles have more levels than NLEVELMAX.
- Fix: DIVA recipe in docs; for profiles add a check in profiles_init_fixed.F or raise NLEVELMAX (PR by averdy, Apr 2026).
- Era: 2009, 2026.
- Src: mitgcm-support 2009-October 'Adjoint Restart' http://mailman.mitgcm.org/pipermail/mitgcm-support/2009-October/thread.html ; https://github.com/MITgcm/MITgcm/issues/985

## Platform / scripting
### Segmentation fault when writing pickup (Darwin/DIC, big runs)
- Cause: stack size limit; or bad I/O option (globalFiles, Cray).
- Fix: `ulimit -s unlimited` (bash) / `limit stacksize unlimited` before mpirun; compile with -g + traceback; use useSingleCpuIO.
- Era: 2017.
- Src: mitgcm-support 2017-December 'Failure to create pickup file for DIC' http://mailman.mitgcm.org/pipermail/mitgcm-support/2017-December/thread.html ; 2017-April 'error while writing pickup files with Cray compilers' http://mailman.mitgcm.org/pipermail/mitgcm-support/2017-April/thread.html

### Segfault on restart because of NaN in bathymetry/other input (not the pickup)
- Cause: NaN in bathymetry/topography file.
- Fix: clean inputs; bathyFile negative, no NaN.
- Era: 2014.
- Src: mitgcm-support 2014-February 'Segmentation Fault' http://mailman.mitgcm.org/pipermail/mitgcm-support/2014-February/thread.html

### Automating restarts under PBS/SLURM (find most recent pickup, update nIter0, resubmit)
- Cause/Answer: community scripts: update `nIter0` from the latest complete pickup (pChkpt or ckpt) and resubmit; ECCO-Darwin `modpickup` picks latest available complete pickup set; model must finish (STDOUT 'Execution ended Normally') before the resubmit, otherwise the next job starts from a missing/partial pickup.
- Fix: Abernathey's gist `most_recent_pickup.sh`; modpickup scripts at MITgcm_contrib/high_res_cube/input/modpickup and https://github.com/MITgcm-contrib/ecco_darwin/blob/master/v03/cs510_latest/input/modpickup; use `writePickupAtEnd=.TRUE.`.
- Era: 2017-2022.
- Src: mitgcm-support 2017-October 'Bash scripts to restart' http://mailman.mitgcm.org/pipermail/mitgcm-support/2017-October/thread.html ; 2022-August 'rolling pickups: nchecklev value on restarting model' http://mailman.mitgcm.org/pipermail/mitgcm-support/2022-August/thread.html

### Vector-computer / aggressive optimisation: crash in seaice_growth after unrelated change (stale restart behaviour)
- Cause: aggressive optimisation on SX-ACE changed behaviour after unrelated diagnostics change (PR #348).
- Fix: add the file to NOOPTFILES in the optfile (tools/build_options/SUPER-UX_SX-ACE_sxf90_awi).
- Era: 2020 (checkpoint67s).
- Src: mitgcm-support 2020-October 'Bug in seaice code since checkpoint67s?' http://mailman.mitgcm.org/pipermail/mitgcm-support/2020-October/thread.html

### Misc setup errors seen alongside restarts
- "useOASIS set but pkg/oasis not compiled" / strange parameter reads after code update: local copies of headers (PARAMS.h, EEPARAMS.h, ini_parms.F) in code/ are stale after updating the checkout (common-block length mismatch); diff them against the new originals. (2005-December viscC4smag; 2011-January oasis.)
- Src: mitgcm-support 2011-January 'oasis ????' http://mailman.mitgcm.org/pipermail/mitgcm-support/2011-January/thread.html ; 2005-December 'bug with viscC4smag ?' http://mailman.mitgcm.org/pipermail/mitgcm-support/2005-December/thread.html
