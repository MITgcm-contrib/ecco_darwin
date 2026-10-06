# Configuration sanity: code_* + run_* pairs

MITgcm validates surprisingly little across files. Mismatches show up as crashes, or
worse as silently wrong indexing. When asked to check, change or port a configuration,
walk this list and report each item as OK / mismatch / can't tell.

## Compile-time ↔ run-time couplings

| Quantity | Set in | Must agree with |
|---|---|---|
| `sNx,sNy,OLx,OLy,nSx,nSy,nPx,nPy,Nr` | `SIZE.h` | grid dims in `data`/`data.exch2`; `mpirun -np`; `delZ`/`tRef`/`sRef` length = `Nr` |
| Overlap `OLx,OLy` | `SIZE.h` | ≥ what advection schemes/packages need (e.g. pkg/bbl needs ≥2; high-order advection/`useCDscheme` more) |
| `PTRACERS_num` | `PTRACERS_SIZE.h` | `PTRACERS_numInUse` and the name/units/initialFile lists in `data.ptracers` |
| Darwin `nplank`,`nGroup`,`nPhoto`,`nopt` | `DARWIN_SIZE.h` | `grp_nplank` sum, `grp_photo` count, `grp_names` length in `data.darwin`; one Chl tracer per photo type if `DARWIN_ALLOW_CHLQUOTA` |
| `nlam` | `RADTRANS_SIZE.h` | counts of `RT_EdFile`/`RT_EsFile`/`RT_wbRefWLs`, `darwin_scatSlope*`; `RT_wbEdges` = `nlam+1`; `RT_kmax ≤ Nr` |
| Diagnostics counts | `DIAGNOSTICS_SIZE.h` | `numlists`, `numperlist`, `numDiags` vs `data.diagnostics` |
| MNC / mdsio list sizes | `MNC_SIZE.h`, `RW_MFLDS.h` | many tracers × streams can overflow the defaults |
| `ALLOW_<PKG>` (packages.conf) | genmake2 | `use<PKG>` in `data.pkg`; presence of `data.<pkg>` |
| CPP options | `CPP_OPTIONS.h`, `<PKG>_OPTIONS.h` | namelist params that only exist when a flag is defined (e.g. `USE_EXF_INTERPOLATION` for `*_lon0` blocks) |

For darwin3, generate the ptracers block instead of hand-typing it:
`<darwin3>/tools/darwin/mkdarwintracers -r -f -u` run in the build dir (needs current
`DARWIN_SIZE.h` and `DARWIN_INDICES.h`). Reordering tracers breaks `PTRACERS_initialFile`
mapping and pickups silently.

## Time and calendar

- Check `deltaT × nTimeSteps` (or `nEndIter − nIter0`) against the intended period, and
  every `*Freq`/`*Period` against the calendar: 360-day years (31104000 s, 30-day months
  2592000 s) vs Gregorian with pkg/cal (`data.cal` `startDate_1`, `calendarDumps`).
- pkg/offline keeps its own clock: `deltaToffline`, `offlineIter0`, `offlineForcingPeriod`,
  `offlineForcingCycle` (`0.` = non-repeating archive). `offlineIter0 × deltaToffline` must
  equal `nIter0 × deltaT`; moving `nIter0` mid-cycle requires moving `offlineIter0` too.
- Record timing: offline record k is centred at day k−0.5 (one-day-early pitfall); EXF
  fields load at the start of a step while OBCS values apply at the end (one-step lag);
  EXF reads ahead, so pad forcing files with an extra record. Each source's
  `*startdate1` must be that file's first real record date.
- `data_org`-style copies of namelists are not includes — they are stale backups.
- **Every forcing start date vs the model start.** For each EXF/darwin/radtrans field (ice, wind,
  iron, pCO2, OASIM, runoff, ...), compute `*startdate1/2` against `data.cal` `startDate_1/2`. With
  `repeatPeriod`/`*RepCycle = 0` (the EXF default), a first record *later* than the model start
  makes `EXF_GetFFieldRec` STOP at t=0 ("myTime ... earlier than 1rst record"). Also check the last
  record covers the end time (EXF needs the record after it too). Real incident: re-centring
  `icestartdate` to 19920102 00Z with the model starting 19920101 12Z.
- **Native-grid forcing read with interpolation on.** With `USE_EXF_INTERPOLATION` defined, every
  field defaults to `*_interpMethod=1` (lat-lon). Fields already on the model grid (e.g. LLC90
  compact ice or wind files) need `*_interpMethod=0`, or they are silently garbled.

## Namelist traps

- Inline comments use `!`. A trailing `#` kills the read. Keep lines < ~200 chars.
- Name families can't be mixed: `viscAz` vs `viscAr` → "Cannot mix z, p and r".
- Packages read their own namelist even if empty; a compiled-in package with a missing
  `data.<pkg>` (ggl90, exf, …) stops the run.
- `diag_mnc` defaults to `useMNC` (diagnostics_readparms.F), so `useMNC=.TRUE.` with mnc compiled
  and `diag_mnc` unset writes diagnostics as NetCDF, not MDS. Pin `diag_mnc=.FALSE.` for `.data/.meta`.
- `offlineForcingCycle` equal to the full archive length hits a boundary case in
  GET_PERIODIC_INTERVAL; use `0.` for a non-repeating archive.
- Under pkg/offline: `momStepping`/`tempStepping` flags in `data` don't make dynamics
  prognostic; GGL90 complains it "needs implicitViscosity" — a custom `ggl90_check.F` is
  carried in `ECCO/offline/code_offline_ggl90`.
- `GAD_SMOLARKIEWICZ_HACK` (GAD_OPTIONS.h): "Do not use with Adams-Bashforth (for ptracers)! Do not use with OBCS!"
- `EXF_INTERP_READ` wants one combined file; `MDS_READ_FIELD` needs a real `.data/.meta` pair.
- Duplicate inputs (e.g. sea-ice area in both `data.darwin` `icefile` and `data.radtrans`
  `RT_icefile`) must change together.

## Offline runs must mirror their forward run

When a pkg/offline run consumes a forward run's archive, compare the two configurations directly
(e.g. `run_v4r6_forward/data*` vs the offline `data*`). Anything the forward run used that the
offline run doesn't reproduce is a silent error:
- Grid: `hFacMin`, `hFacMinDr`, `bathyFile`, and ice-shelf cavities (`useShelfIce` forward but not
  offline turns cavities into open surface columns).
- Eddy transport: forward `GM_AdvForm=.TRUE.` means archived `Kwx/Kwy` are Redi-only and bolus
  velocity isn't in UVEL/VVEL, so offline loses GM bolus transport unless residual velocities or
  the full tensor are archived.
- Mixing: 3-D `diffKrFile` (`ALLOW_3D_DIFFKR`) and 3-D Redi (`xx_kapredi`) fields.
- Free surface (z*/NLFS vs linear), real freshwater vs virtual salt flux (which flux drives
  DIC/ALK dilution).
- Archive timing: which iteration each record's suffix holds vs where pkg/offline centres it, and
  whether the archive covers the final step.

## Paths and inputs

- Grep every `data*` file for absolute paths (`grep -n "/" data* | grep -v "^.*!"`) and
  list which ones don't resolve on the target machine.
- Check input binaries' size = nx·ny·(nz)·nrec·(4 or 8 bytes) for `readBinaryPrec`, and
  big-endian.
- Diagnostics: check each requested field name exists (`available_diagnostics.log` from a
  short run); bad names appear only in `STDERR.0000`. `fileName` dirs must exist.
- Per-timestep stats (`stat_freq` = `deltaT`) produce huge output on long runs.

## Before a long run

1. Short run (temporary `nTimeSteps`, then restore it) → `NORMAL END`, sane `%MON`
   stats, expected record counts.
2. Restart test if pickups changed: N+N steps from pickup equals 2N steps straight.
3. Re-sync local edits to the HPC copy and diff.
