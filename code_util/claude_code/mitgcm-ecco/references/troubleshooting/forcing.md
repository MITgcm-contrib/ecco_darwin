# Troubleshooting: surface forcing (pkg/exf, cal, restoring, runoff)
Distilled from mitgcm-support (2003-2026) and MITgcm GitHub issues/PRs, core-developer answers only; names checked against upstream master (Oct 2026) and the source index. "Martin" = Martin Losch, "J-M" = Jean-Michel Campin. Month-thread URLs are http://mailman.mitgcm.org/pipermail/mitgcm-support/YYYY-Month/thread.html.

## EXF set-up: CPP combinations, namelist errors, unrealistic heat flux

### Temperature goes wild / hflux far larger than input, e.g. exf_hflux_max 1200 when the file max is 679 (EXF_OPTIONS.h flags inconsistent with forcing files)
- Cause: `ALLOW_ATM_TEMP` / `ALLOW_DOWNWARD_RADIATION` / `ALLOW_BULKFORMULAE` defined, but atemp/aqh/lwdown are not supplied, so exf uses atemp = 0 K, aqh = 0 and computes enormous sensible/latent/longwave fluxes. `useExfCheckRange=.FALSE.` hides it.
- Fix: pick one of the combinations in the table at the top of pkg/exf/EXF_OPTIONS.h. Read-in hflux+swflux+sflux: `#undef ALLOW_ATM_TEMP`, `#undef ALLOW_DOWNWARD_RADIATION`, `#define ALLOW_ATM_WIND` (or stress), `#define ALLOW_RUNOFF` as needed. Generic hindcast (also required with pkg/seaice): `ALLOW_ATM_TEMP`, `ALLOW_DOWNWARD_RADIATION`, `ALLOW_BULKFORMULAE`. Always dump EXFhs, EXFhl, EXFlwnet, EXFswnet, EXFqnet, EXFtaux/y, EXFlwdn, EXFswdn with pkg/diagnostics and compare to hflux typical range -250..600, swflux -350..0 W/m^2.
- Era: 2009-2025, one of the most frequent EXF questions.
- Src: mitgcm-support 2013-January 'strange hflux value in MON exf_hflux'; 2011-December 'Query regarding "exf package".'; 2009-September 'Shortwave radiation'; 2014-April 'something wrong with exf'

### Ocean cools indefinitely (-20 C coastal, lake, "simple EXF setup") with downward radiation missing or wrong
- Cause: upward longwave = emissivity*sigma*T^4 (~330 W/m^2 at 5 C) is always computed; with lwdown/swdown absent, too small (ERA5 accumulated J/m^2 not converted to W/m^2 rate), or wrong sign, the ocean freezes/cools. exf wants DOWNWARD sw/lw, not net, and rates not accumulations.
- Fix: check lwdown 50-450 and swdown 0-450+ W/m^2; convert accumulations by dividing by the accumulation period (value belongs to mid-interval). Use EXFlwdn/EXFswdn diagnostics, `debugLevel>3`, then 'print' debugging. Reduces to: sanity-check input files first.
- Era: 2019-2026.
- Src: mitgcm-support 2025-December 'unrealistic model results in a regional model'; 2019-November 'Weird cooling with simple EXF setup'; 2022-February 'interpolate atmospheric forcing data'; https://github.com/MITgcm/MITgcm/issues/978

### STOP "EXF_CHECK: use u,v_wind components but not wind-stress" (or the reverse)
- Cause: `useAtmWind=.TRUE.` (set by `ALLOW_ATM_WIND`) while you give ustress/vstress files.
- Fix: `useAtmWind=.FALSE.` in data.exf (exf then derives wind speed from stress, see exf_wind.F); the older CPP `EXF_ATM_WIND` is largely replaced by the run-time flag. The 2008 typo (ustressfile tested twice) is fixed.
- Era: 2008-2018, message still in exf_check.F.
- Src: mitgcm-support 2014-January 'error with exf "use u, v_wind components but not wind-stress..."'; 2017-December 'how to use windstress to calculate wind speed'

### Namelist error "variable not in namelist EXF_NML_0x" / "read unexpected character" / "Fortran runtime error: syntax error in NAMELIST"
- Cause: no comma after an entry, '#' not in column 1, blank line, a name that does not exist in that namelist (e.g. `*_lon0` without `USE_EXF_INTERPOLATION`), `_lat_inc` given as one value.
- Fix: add trailing commas; first column only ' ', '&' or '#'; `#define USE_EXF_INTERPOLATION` in EXF_OPTIONS.h for lon0/lat0 etc. (EXF_NML_04 must still exist, empty, if you do not interpolate). Also check which namelist was last opened in STDOUT (scratch copy name).
- Era: 2008-2018 repeatedly.
- Src: mitgcm-support 2008-August 'input variables not recognised'; 2013-December 'error with exf forcing "variable not in namelist"'; 2018-April 'exf package'; 2008-September 'parameter set'

### "EXF_INTERP_UV: input grid must encompass output grid" / striped or distorted forcing / "Non-existing record number" with interpolation
- Cause: `*_lat_inc` must have nlat-1 entries (vector, e.g. `720*0.25` -> `719*0.25`), `*_lon_inc` is a single scalar; input grid must span the model domain (no extrapolation). Latitude 90N edge or too-small `exf_max_nLon/Lat` (default 520/260) break the read. A file with fewer records than needed gives "Non-existing record number".
- Fix: `precip_lat_inc = 719*0.25D0,` (nlat-1 values; irregular Gaussian grids: list real spacings); `precip_lon_inc = 0.25D0,`; lon0/lat0 are the centre of the first input cell; define `EXF_INTERP_USE_DYNALLOC` or raise `exf_max_nLon/nLat` in EXF_INTERP_SIZE.h (different compilers behave differently with the dynamic version). Different grid per field: set `fld_nlon/nlat/lon0/lat0` for each. Rotated or non-lat/lon grids: use `uvInterp_*` and see EXF_NML_04 docs; Cartesian grids cannot use exf interpolation (stop added in 2006).
- Era: 2006-2023; still valid.
- Src: mitgcm-support 2023-February 'Problems with the EXF package and binary input files'; 2023-September 'EXF configuration'; 2016-January 'EXF Interpolation problem.'; 2022-April 'monthly forcing data'; 2016-June '(no subject)' (different sizes); 2012-May 'Definition of spatial boundary in EXF package'

### Bicubic vs bilinear in exf_interp: spurious values, rectangular patterns in wind-stress curl
- Cause/Fix: defaults are bicubic for vectors (`*_interpMethod = 12/22`, avoids derivative artefacts) and bilinear for scalars (=1, avoids negative rain/humidity). Piecewise-linear winds give piecewise-constant curl, so keep bicubic for winds. Interpolation near the poles with a data last-latitude very close to 90N can blow up (LAGRAN denominator); south-pole row treatment was fixed in 2016.
- Era: 2009-2017.
- Src: mitgcm-support 2009-October 'exf interp'; 2017-May 'numerical issues with exf_interp_uv?'; 2016-September 'Bug in exf_interp routine?'

### EXF values wrong with compiler optimisation on (atemp of constant 300 K shows max 1200)
- Fix: compile exf_interpolation.f with lower optimisation: `NOOPTFILES='exf_interpolation.f'`, `NOOPTFLAGS='-O0 -g'` (or a lower -O level) in the optfile; then `make makefile && make depend`.
- Era: 2017.
- Src: mitgcm-support 2017-December 'EXF'

### Constant-in-time forcing with exf; stress read, wrong grid
- Fix: `*period = 0.` (default) in data.exf for all constant fields (one record); `readStressOnCgrid=.TRUE.` (or `readStressOnAgrid`) must match where the stress is defined; without exf set `zonalWindFile` (C-grid) and no `periodicExternalForcing`. Quirk: constant wind-stress with spatial interpolation has special code in exf_init_varia.F. A 2015 STOP for period 0 with interpolation in exf_set_uv.F was removed.
- Era: 2011-2024.
- Src: mitgcm-support 2024-July 'Setting up external wind stress in ocean model'; 2012-January (thread 'Query regarding "exf package".'); 2015-October 'EXF period with interpolation'

### Wrong values in generic (non-exf) forcing: taux0/tauy0 empty, fu/fv constant
- Cause: without `periodicExternalForcing` stress is read directly into `fu`,`fv`; with it into `taux0/tauy0` then interpolated (issue #702).
- Fix: use `fu`, `fv` (surface flux of momentum) for constant forcing; in `data` PARM03 `periodicExternalForcing=.TRUE., externForcingPeriod=..., externForcingCycle=...` also for non-periodic series (just keep integration time below the cycle). Those params are ignored when exf is used.
- Era: 2012-2023.
- Src: https://github.com/MITgcm/MITgcm/issues/702 ; mitgcm-support 2012-November 'add the time dependent wind by zonalWindFile'; 2005-August 'Re: single-layer configuration problem update'

### Forcing in non-exf mode is linearly interpolated, "suddenly changing" wind not possible
- Fix: edit weights in model/src/external_fields_load.F (set aWght/bWght so only one record is used). Records are assumed mid-interval for non-exf fields (use start date mid-month in exf to match). Period-0 weights in get_periodic_interval.F are 0.5 at t=0.
- Era: 2013.
- Src: mitgcm-support 2013-February 'Wind stress without interpolation'

### exf relative wind: is ocean current subtracted for latent/sensible heat too?
- Fix/answer: `useRelativeWind=.TRUE.` (data.exf) subtracts the surface ocean velocity in exf_getffields.F, so wind stress AND ustar/latent/sensible computations use uwind-uVel; over ice, seaice_get_dynforcing.F uses uwind-uIce.
- Era: 2016.
- Src: mitgcm-support 2016-May 'if the ocean surface current is subtracted...'; 2016-June 'how to calculate the relative velocity...'

### exf requires pkg cal (build fails with missing cal headers; STOP CAL_FULLDATE called too early)
- Cause: exf depends on cal (pkg_depend); in packages added after Apr 2012 pkg/cal setup moved to cal_init_fixed.F, so cal_fulldate called from your own *_readparms is too early (`cal_setStatus=0`).
- Fix: add `cal` to packages.conf (pkg_depend handles it today); move early cal calls to your package's init_fixed. pkg/exf can run without cal (as in offline_exf_seaice) which also allows a time step that is not a multiple of 1 s.
- Era: 2006-2019.
- Src: mitgcm-support 2006-September 'Exf depends on cal?'; 2014-September 'Problem with cal_fulldate'; 2019-October 'Sub-1s Timestep with EXF'

## Calendar, start dates, periods, restarts

### Model reads the wrong record / crashes with "Fortran runtime error: Non-existing record number" (forcing file too short)
- Cause: the model linearly interpolates between records and needs one record past the end. Daily forcing for 30 days starting at 00:00 needs 31 records; hourly data needs `hfluxperiod = 3600.` (not deltaT). Period must equal the real spacing between records.
- Fix: add records; `*period` in seconds between records; `*startdate1/2` = time of record 1; data.cal `startDate_1/2` = time of nIter0=0. Mid-interval start (e.g. 120000) needs an extra record before. Turn on `debugLevel=3` to see records read.
- Era: 2007-2023, perennial.
- Src: mitgcm-support 2020-December 'Getting Backtrace error'; 2018-January 'wind stress file'; 2010-January 'Exf package and interpolation in time'; 2016-June '(no subject)'

### Monthly/climatological forcing: period -12, mid-month start, leap years, "how long is a year"
- Fix: `*period = -12.` = 12 monthly records that repeat (even in Gregorian calendar); start date mid-month (`00000115`); with `useExfYearlyFields=.TRUE.` only fields with period -12 are climatologies, others look for `file_YYYY`. Regular monthly period in a Gregorian run: 2629800 s = (365*3+366)/4*86400/12 (calendarDumps does NOT apply to forcing). Gregorian has 366-day years; a leap day needs an extra Feb 28 in a file, or `cal_isleap.F` returns false / `nDaysLeap` 365 hack. `startdate` year is replaced by the file year with yearly fields, so start on Jan 01 (extra file for the previous year), otherwise yearly OBCS may look constant.
- Era: 2011-2022.
- Src: mitgcm-support 2020-October 'using EXF with monthly data and pkg cal gregorian'; 2014-April 'how long is a year? & monthly mean forcing'; 2013-August 'problem with data.exf and dates'; 2011-April 'useOBCSYearlyFields'; 2021-January 'issues with SST/SSS restoring and OBCS in EXF package'

### Starting year/phase issues: restart does not read forcing from the right time; repeatPeriod, yearly fields, very long runs
- Fix: put all years in one file per field if possible and keep `startdate*` fixed during restarts; the model chooses the record from the model date. `repeatPeriod` > 0 with `useExfYearlyFields=.TRUE.` is not implemented. With model calendar and a 360-day year, switching calendars at restart needs a renamed pickup or `pickupSuff`. For periodic non-exf forcing restarts the first record depends on nIter0 (e.g. DIC perturbation file mismatched by 50 days): set `DIC_forcingCycle` explicitly. Record selection rounding bug in external_fields_load (2011) fixed.
- Era: 2009-2024.
- Src: mitgcm-support 2024-January 'Time-varying forcing and EXF package' / 'Re: Re: Digest, Vol 247, Issue 1'; 2009-July 'EXF package: useExfYearlyFields=.FALSE.'; 2015-May 'Model/DIC_forcing time mismatch'

### calendarDumps / month-end output and pickups
- Fix: `calendarDumps=.TRUE.` (data.cal) converts approximate months/years in `chkPtFreq/pChkPtFreq/taveFreq/freq` of diagnostics (not forcing) to exact calendar ones; works with Gregorian only meaningfully. dumpFreq/`$PKG_dumpFreq` output does not use calendar dumps (only diagnostics, KPP_taveFreq, seaice_output, pickups); pickup name = actual iteration at month end. To run exactly one month use days*86400 (+1 step) endTime, `rn_pickup` script; `dumpInitAndLast` controls first/last writes. `timePhase > 366 days` with `calendarDumps` and leap years gives wrong output times (issue #992, open, proposed fix not merged as of Oct 2026). Time-interval stamps in .meta do not account for calendar rounding.
- Era: 2014-2026.
- Src: https://github.com/MITgcm/MITgcm/issues/992 ; mitgcm-support 2017-September 'how to run for exactly 1 month with calendarDumps'; 2018-November 'diagnostics at end of calendar months'; 2021-May 'calendarDumps and pickup file names'; 2014-March 'Using pkg/calendar with diagnostics'; 2025-July 'ptracer/DIC budget closure/calendar pkg'

### EXF fldPeriod=0 with ctrl xx_gentim2d; ctrl time-varying weights backwards with useCAL=F
- Cause: exf_set_fld.F does not re-initialise when fldPeriod==0; and `CTRL_GET_GEN_REC` returns the second-record weight when `useCAL=.FALSE.` (issue #1039: time-varying controls interpolated backwards within each interval).
- Fix: PR #980 (merged 2026) removes the `fldPeriod .NE. 0.` condition; PR #1040/issue #1039 open as of Oct 2026; workaround: use pkg/cal. gentim2d controls without cal are periodic or constant (use externForcingCycle logic).
- Era: 2023-2026.
- Src: https://github.com/MITgcm/MITgcm/pull/980 ; https://github.com/MITgcm/MITgcm/issues/1039 ; mitgcm-support 2023-February 'nonperiodic gentim2d controls w/out calendar'

### Interannual monthly forcing broke OBCS/EXF in checkpoint68d
- Cause: PR #545 computed `${fld}startTime` even when `${fld}period` is zero.
- Fix: PR #562 (Nov 2021) skips it. Use master after Nov 2021 or revert.
- Era: checkpoint68d, fixed upstream.
- Src: https://github.com/MITgcm/MITgcm/pull/562 ; mitgcm-support 2021-November 'OBCS/EXF in Checkpoint68d'

### TAF TL code with USE_EXF_INTERPOLATION: "EXF_INTERP_READ: File does not exist"
- Fix: add `CADJ SUBROUTINE exf_getyearlyfieldname REQUIRED` to exf_ad.flow (issue #955, open as of Dec 2025). Also adjoint ADJ* files need `ALLOW_AUTODIFF_MONITOR` in AUTODIFF_OPTIONS.h (on by default now) and adjDumpFreq.
- Era: 2023-2025.
- Src: https://github.com/MITgcm/MITgcm/issues/955 ; mitgcm-support 2023-March 'ADJ* missing'

## Heat and freshwater fluxes, runoff, free surface

### Negative salinity / crazy salinity (-200..90) at coastal cells or big river runoff
- Cause: virtual salt flux with linear free surface; isolated surface cells (single connection) cannot advect away the forcing; large runoff (1e-4 m/s) into thin cell.
- Fix: `useRealFreshWaterFlux=.TRUE.` with `exactConserv=.TRUE.` (on old checkpoints also CPP `EXACT_CONSERV`; that CPP option was removed in PR #930), better still `nonlinFreeSurf=4` (+`select_rStar=2`); alternatives `convertFW2Salt=-1.` (local salinity) or remove/open the isolated points. nonlinFreeSurf requires exactConserv (config_check). EmPmRfile of cold freshwater input: `temp_EvPrRn` / NLFS needed.
- Era: 2004-2021.
- Src: mitgcm-support 2016-June 'Large runoff causes negative salinity'; 2019-November 'Bulk forcing and partial cells'; 2021-January 'Forcing a surface fresh water flux with EmPmRFile'; 2007-April 'EmPmRfile, EmPmR'

### Runoff: how to put river water in; at depth; heat content of runoff
- Fix: default is surface level, subtracted from EmP into EmPmR (only wet cells count): `runoffperiod=-12`, `runoffFile`; interpolation off with `runoff_interpMethod=0` and a model-grid file; 360x180x12 file needs `runoff_nlat-1` values of `runoff_lat_inc`. Depth: `ALLOW_ADDFLUID` + `selectAddFluid=1` + `addMassFile` (static only), or pkg/obcs. Runoff temperature: exf has no entry; heat content depends on deltaT (issue #719, open); for coupled aim/ocn set `temp_EvPrRn=0.` with NLFS; `runOffMapFile` (cpl_aim+ocn) is runOffMapSize float64 triples [atm index i+(j-1)Nx, ocean index, catchment area m^2].
- Era: 2013-2025.
- Src: mitgcm-support 2016-November 'river runoff'; 2016-February 'Error with Arctic configuration and runoff file'; 2025-January 'runoff file of cpl_aim+ocn'; https://github.com/MITgcm/MITgcm/issues/719

### Prescribing heat/salt flux (hflux, sflux) with seaice on
- Cause: the sea-ice model recomputes fluxes where there is ice; code stops/ignores.
- Fix: not supported. Over ice, `Qnet = Qnet_ocean*(1-AREA) + Qice*AREA` after seaice and lwdown is unused; to add a heat anomaly dQ, add `dQ*(1-AREA)` after seaice (do_oceanic_phys). Alternative: use diagnosed oceQnet/oceFWflx (opposite sign to hflux/sflux) from a spun-up run. With seaice use EXF option (3) bulk. `EXF_READ_EVAP` with thsice stops in thsice_get_exf.F.
- Era: 2015-2022.
- Src: mitgcm-support 2022-April 'forcing model with hflux and sflux'; 2015-July 'EXF/SEAICE Surface Heat Fluxes'; 2017-October 'Prescribed atmospheric temperatures...'

### Net heat flux sign and components (Qnet, Qsw, shortwave penetration)
- Fix: `Qnet` (hflux) is UPWARD positive = latent+sensible+net lw+net sw (`>0` cools). Qsw is subtracted from Qnet at surface and re-added through swfrac.F (Jerlov water type hard-coded jwtype=2); to get penetration provide `qswfile` as well as `surfQnetFile`. With `ALLOW_BULKFORMULAE` cloud effects must already be in swdown (no cloud file in exf). exf_albedo is open-ocean water albedo (~0.1); ice albedos in data.seaice. Apply `balanceEmPmR`/`balanceQnet` (CPP `ALLOW_BALANCE_FLUXES`) to remove global mean drift in long runs.
- Era: 2005-2022; config_summary printing of `balanceEmPmR` fixed by PR #229 (2019).
- Src: mitgcm-support 2013-June 'Qnet and Qsw heat fluxes after ceckpointc54'; 2017-May 'The Net Heat Flux'; 2010-April 'paleo runs and model drift'; 2010-January 'Understanding the shortwave radiation'; https://github.com/MITgcm/MITgcm/issues/226

### Freshwater flux, sea-ice melt vs exf budget (docs gap)
- Fix: no clean manual section for combos of `useRealFreshWaterFlux`, `convertFW2Salt`, `temp_EvPrRn`, linear/NLFS (issues #184, #700). Martin's recommended sets: rigid lid; linear free surface with `implicitFreeSurface=.TRUE.`, `exactConserv`; NLFS 4 only with `exactConserv=.TRUE., useRealFreshWaterFlux=.TRUE.`. FW flux with salinity restoring in adjoint = mixed boundary conditions, avoid.
- Era: 2018-2022.
- Src: https://github.com/MITgcm/MITgcm/issues/184 ; mitgcm-support 2022-September 'why fresh water flux not recommened together with salinity restoring for adjoint'

### Evaporation/precip input precision: salt drifts although sflux integrates to zero
- Cause: online fields carry 64-bit precision; fields written to disk are rounded (real*4).
- Fix: compare using the on-disk (read back) fields; use `exf_iprec=64`, `readBinaryPrec=64`.
- Era: 2018.
- Src: mitgcm-support 2018-January 'precision difference between freshwater flux bin file...'

### Wrong friction velocity ustar in exf_wind.F / wrong stable-case term psimh in exf_bulkformulae.F (old code)
- Fix: `exf_wind.F` ustar: `ustar = SQRT(wStress)*recip_sqrtRhoA` now correct (2018 discussion). `exf_bulkformulae.F` psimh bracket bug (missing parenthesis, -2*ATAN(x)+pi/2 added to both stable and unstable cases) fixed 2005; exf_check.F `ustressfile` typo fixed 2008.
- Era: fixed upstream.
- Src: mitgcm-support 2018-April 'Friction velocity ustar in S/R exf_wind.F'; 2005-June 'bug in exf_bulkformulae.F ?'; 2008-June 'typo in exf_check.F?'

### Frictional heating (ALLOW_FRICTION_HEATING) gives zeros in ocean
- Cause: `frictionHeating` is never computed for the ocean, only atmosphere physics.
- Fix: implement it yourself (in analogy to atm_phys_tendency_apply).
- Era: 2020; may still be true (check source).
- Src: mitgcm-support 2020-October 'Frictional heating in the ocean model'

## Restoring, RBCS, sponge, relaxation

### Surface SST/SSS restoring: set-up, time scale, regional restoring, restarts
- Fix: with exf `climsstTauRelax`, `climsssTauRelax`, `climsst*period`, without exf `tauThetaClimRelax`, `thetaClimFile`. Spatially varying strength: `lambdaThetaFile/lambdaSaltFile` (PARM05) but you still need a positive tau (doThetaClimRelax prints in STDOUT); `latBandClimRelax` limits latitude. Tau must be >= deltaT; stability: fully explicit, `deltaT/tau < 1` (<1/2 with AB-2 for non-oscillatory); tau<deltaT gives noise. Piston velocity logic: scale tau with surface-layer thickness (tau ~ dz/piston). Crash on adding SSS relaxation: check the relaxation field for NaN/unrealistic values over ocean, and look at SRELAX diagnostic. Restoring to monthly climatology with long tau adds a phase shift. Restore only some part of exf: put all forcing specs in either data or data.exf, not both.
- Era: 2009-2026.
- Src: mitgcm-support 2026-March 'Problem regarding how to add SSS relaxation'; 2014-October 'SSS restoring in selected regions'; 2014-February 'Adding constant submerged and surface ice.'; 2012-June 'tauSaltClimRelax for global ocean simulation'; 2014-May '"how to proper use both net heat flux and sstrelaxation??"'; 2010-February 'Bulk Force and SSS restoring'

### pkg/rbcs: restoring velocities, tracers, ice; maskLEN; time-scale; restart offset
- Fix: 3D restoring of T, S, U, V and ptracers with 3D masks (also regional/sponge); `#undef DISABLE_RBCS_MOM` for velocities; masks for u/v are staggered, a T-mask can mask velocity at boundary; `maskLEN` in RBCS_SIZE.h and irbc = 2+iTr; timing params `rbcsForcingPeriod`, `rbcsForcingCycle`, `rbcsForcingOffset` (seconds; replaced old `rbcsIniter`; default 0 relaxes from restart start); do not use relaxation time < deltaT (T(n+1)=T_r at tau=deltaT); no support for sea ice (extend yourself; relaxing ice velocity is pointless). With `debugLevel>=3` RBCS writes masks; Um_Ext/Vm_Ext diagnostics contain RBCS tendency. pkg/ctrl cannot restore. A non-rbcs 'sponge' of OBCS is not in gT_Forc diagnostics, rbcs is.
- Era: 2009-2026.
- Src: mitgcm-support 2014-January 'RBCS relax u velocities problem'; 2021-March 'Change relaxation forcing after the spin up'; 2023-December 'Trying to relax ptracer using RBCS'; 2018-March 'RBCS relaxation time scale limits'; 2026-April 'How to restore the surface currents...'; 2024-September 'Lateral heat flux due to temperature restoring'

### Restoring tracer RBCS: asymmetric solution under symmetric forcing
- Cause: C-grid staggering when a user criterion uses theta at velocity points.
- Fix: average theta to velocity points before testing: `0.5*(theta(i,j)+theta(i,j-1))`.
- Era: 2023.
- Src: mitgcm-support 2023-October 'Asymmetric surface temperature under symmetric forcing when rbcs package is used'

### 3D forcing/tendency in momentum or tracer equations (body force, drag, tidal body force, bottom heat flux)
- Fix: `model/src/apply_forcing.F` (APPLY_FORCING_U/V/T) is the active place; `external_forcing.F` only with `USE_OLD_EXTERNAL_FORCING`. Or pkg/mypackage (`MYPACKAGE_TENDENCY`, `myPa_applyTendU/V`, tendency apply T for bottom heat flux with kLowC). In apply_forcing_u/v `theta` is available. gchem forcing of theta: do it in gchem_calc_tendencies, not gchem_forcing_sep (called at end of step). Use XC/YC (not indices) or global indices `ig=myXGlobalLo-1+(bi-1)*sNx+i` to localise forcing under MPI. Non-hydrostatic is not needed for tidal body force.
- Era: 2013-2023.
- Src: mitgcm-support 2023-June 'How to add a linear drag for U and V when temperature is lower than a value'; 2020-February 'apply_forcing.F vs external_forcing.F'; 2022-February '3D forcing for the momentum equations'; 2015-April 'Forcing theta within gchem'; 2016-November 'Process tracking'

### Tidal forcing options
- Fix: surface tidal potential: `tidePot` in pkg/exf (J-M added, prescribe separately from atmospheric pressure) or JPL SPICE-based tide potential (Oliver Jahn); regional domains need both OBCS tidal velocities (`useOBCStides=.TRUE.`, amplitude and phase as in verification/seaice_obcs/input.tides, phase in seconds relative to startdate) and potential. Do not double count tides with `useOBCSprescribe` fields and `useOBCSbalance` (balance suppresses all non-tidal transport). Old: `ATMOSPHERIC_LOADING` with pLoadFile/apressurefile. OBC tidal implementation improvements: issue #617.
- Era: 2008-2025.
- Src: mitgcm-support 2024-December 'How to include tides by the atmospheric pressure like the tidal forcing in LLC4320'; 2025-July 'Question About Sea Surface Height (Eta) Behavior in MITgcm'; 2017-August 'tidal period'; https://github.com/MITgcm/MITgcm/issues/617

## OBCS interactions (forcing-related)

### OBCS + EXF: STOP EXF_CHECK_RANGE after some time, or noise at open boundary
- Cause: unbalanced boundary data produce strange SST/hflux; large convergence/divergence at the OB; bulk fluxes respond.
- Fix: `useExfCheckRange=.FALSE.` to look at the fields; inspect snapshots (not 3 h means); Stevens BC reduces vertical motion at the OB; viscosity; test exf interpolation offline. Stevens BC crash (EXTREME Pot.Temp): normal velocities are vertically averaged, try without tangential velocities (`OBWvFile=' '`).
- Era: 2019-2025.
- Src: mitgcm-support 2025-June 'Questions about EXF and OBCS'; 2019-July 'Grid-scale noise in a nested simulation with EXF'; 2021-May 'Stevens BC'

### Prescribed OBCS time control (with and without exf); yearly OBCS fields
- Fix: without exf, OB files follow externForcingPeriod/Cycle (constant if not periodic); with exf, time specs in EXF_NML_OBCS (`obcsNperiod=0.` gives constant in time; verification/obcs_ctrl/input_ad tests it); OB files shape (N or S: Nx*Nr*nt; E or W: Ny*Nr*nt); `OB_Iwest=ny*1`, `OB_Ieast=ny*-1` integer arrays; Eta OBC is not read (zero). For constant OBC but pickup restart: `pickupSuff`. Monthly OB data with `obcsperiod=-12`. wVel OBC reading added 2009.
- Era: 2005-2021.
- Src: mitgcm-support 2014-October 'OBCS_PRESCRIBE_READ'; 2011-April 'useOBCSYearlyFields'; 2021-January 'issues with SST/SSS restoring and OBCS in EXF package'; 2007-March 'Fw: OB*file interpolation'

### OBCS sponge not tile-proof; domain decomposition changes results
- Cause: only tiles that contain an OB are processed, so a sponge wider than sNx/sNy is truncated (tile with sNy=12). Reduction order in cg2d (and cg3d) changes with decomposition (also exch2 blank tiles).
- Fix: sponge thickness < tile size, sNx/sNy >= ~30; `#define CG2D_SINGLECPU_SUM` in CPP_EEOPTIONS.h (slower) for bit-reproducibility; cg3d not covered; chaotic runs differ anyway. Prefer pkg/rbcs over OBCS sponge (3D mask, no overlap issues at OB corners).
- Era: 2024.
- Src: mitgcm-support 2024-May 'Domain decompositions affecting simulations outcome'; 2014-July 'Logarithmic decrease in Eta without OBCS/EXF'

### Periodicity assumptions with walls and OBCS
- Fix: model is always doubly periodic in exchanges; put a wall (zero depth) at one end or OBCS to break it; OBCS overwrites values (u duplicated at ib and ib+1, flat topography across the boundary, `maskInC` avoids cg2d seeing OB tracer cells). No `notUsingXPeriodicity` option in current code. Tiny sNx with one point in a direction acts as 2D (no gradients).
- Era: 2005-2021.
- Src: https://github.com/MITgcm/MITgcm/issues/362 ; mitgcm-support 2019-November 'notUsingX/YPeriodicity'; 2020-August 'basic questions about periodicity in MITgcm domains'

## Input files, binary format, tools

### Input file format: always unblocked IEEE big-endian; precision; orientation
- Fix: write ieee-be, real*4 for `readBinaryPrec=32`/`exf_iprec=32`, real*8 for 64; fastest index x first (Fortran order; numpy `[t,y,x]`, matlab `(x,y,t)`); (1,1) must be the SW corner; Python: `data.astype('>f4').tofile(fname)` (do not mix byteswap and astype); meshgrid `indexing='ij'` pitfalls (Y,X order). No netcdf forcing input (initial fields only). Compilers without a big-endian flag: `-D_BYTESWAPIO` in the optfile (see eesupp DEF_IN_MAKEFILE.h). Single-CPU I/O with mnc not supported, use gluemncbig. Namelist terminator `/` vs `&`: `-DNML_TERMINATOR` (offline example 2018). Cartesian grids: rotate vectors to grid i/j yourself.
- Era: 2005-2021.
- Src: mitgcm-support 2019-January 'file runoff-360x180x12.bin'; 2021-April 'Barotropic ocean gyre tutorial'; 2018-December 'Using NetCDF input files'; 2010-March 'Error in reading the input files'; 2018-December 'Offline model help'

### 'make' after header changes / "relocation truncated to fit" / makedepend out of space
- Fix: `make Clean` (then `make makefile; make depend; make`) whenever a source file is moved/added in a -mods dir, a package is added, or a *_OPTIONS.h changes; "relocation truncated to fit" means >2 GB static memory: ifort `-mcmodel=medium` or more MPI ranks; makedepend MAXFILES: use tools/cyrus-imapd-makedepend and `genmake2 -makedepend`.
- Era: 2007-2014.
- Src: mitgcm-support 2014-August 'EXF and Calendar Package'; 2013-February 'Error Message'; 2007-February 'makedepend: error'

### Forcing is read but nothing happens / ocean at rest / no response
- Fix: check that the model uses the files: `zonalWindFile` etc. must be named in `data`; `periodicExternalForcing`; stress diagnostics oceTAUX/oceTAUY vs EXFtaux. Without bottom drag/friction a wind-forced jet accelerates until NaN. Stratified salinity with diffusion on sloped bathymetry or unbalanced ice-shelf `SHELFICEloadAnomalyFile` creates spurious currents.
- Era: 2019.
- Src: mitgcm-support 2019-January 'Prescribing winds using stress versus bulk formulae'; 2019-November 'Stratification creating current without forcing'

## pkg/cheapaml, aim, thsice forcing (few entries)

### CheapAML on llc90 / with seaice
- Fix: defaults fit lat-lon basin scale; specify all cheapaml input files (verification_other/offline_cheapaml/input/data.cheapaml); `useDLongWave=.TRUE.` on Cartesian grids; `useFreshWaterFlux=.TRUE.`; cheapaml advects over land, so llc corners limit its time step; seaice+thsice+cheapaml example offline_cheapaml `input.dyn`.
- Era: 2014-2018.
- Src: mitgcm-support 2016-August 'cheapAML + global_oce_llc90'; 2014-December 'THSIce noisy/unphysical output'; 2018-June 'cheapaml with sea ice'

### aim_v23 with thsice: SST from file, orbital parameters, insolation
- Fix: `stepFwd_oceMxL=.TRUE.` with `tauRelax_MxL=5184000.` (verification/aim.5l_cs/input.thSI) relaxes the mixed-layer SST to the aim file; `aim_useFMsurfBC` vs `aim_useMMsurfFc`; obliquity via `ALLOW_INSOLATION` and `OBLIQ` in data.aimphys (circular orbit only).
- Era: 2025.
- Src: mitgcm-support 2025-October 'Digest Vol 267, Issue 7'; 2025-May 'Digest Vol 263, Issue 1'

## Misc forcing physics Q&A with clear answers

### Adding mass/evaporation without salinity change (addmass)
- Fix: use real freshwater mode: prescribe E and P (useRealFreshWaterFlux, NLFS); to suppress salt effect set salinity to zero; addMass is for 3D sources.
- Era: 2025.
- Src: mitgcm-support 2025-April 'About "addmass" in MITgcm'

### Atmospheric pressure: sign/offset in EOS; setting p_top > 0
- Fix: EOS expects pressure relative to surface (sea pressure; TEOS-10 subtract 10.1325e4 Pa); issue #199 corrected pressure_for_EOS. Atmosphere top pressure >0: adjust Ro_SeaLevel/delR (p* code assumes pTop=0, breaks if pTop ~ delR(Nr)).
- Era: 2015-2020.
- Src: https://github.com/MITgcm/MITgcm/issues/199 ; mitgcm-support 2015-October 'p /= 0 at top of atmosphere model'

