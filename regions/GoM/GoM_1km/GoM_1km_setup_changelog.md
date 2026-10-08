# GoM_1km setup changelog

2026-10-08

The ERA5 coastal fix closed the temperature runaway (`theta_max` drift 0.527 to ~0.00 degC/day);
the OBCS period and record count were both wrong and are now repaired; river runoff still
dilutes one cell toward zero salinity; and at 4.5 model days per wall day the 26-year run
cannot finish as configured.

---

## Surface forcing — ERA5

The ERA5 cell covering the Belize lagoon is land, so EXF's bilinear stencil pulled a land wind
of 2.76 m/s onto the ocean cell instead of 6.46. That halved evaporative cooling and drove an
equilibrium SST of 32.47 degC. Flood-filling the land cells brought it to 29.06 degC, against
an offshore control of 29.28.

| Change | Before to after | State | Note |
| --- | --- | --- | --- |
| Land/sea mask derived from the mean diurnal `atemp` range | none available to derived | Done | `fix_era5_coastal_wind.py mask`; exact against synthetic truth |
| `u10m`, `v10m` flood-filled, `--erode 2` | 2.76 to 6.46 m/s at the Belize cell | Done | worth only ~10 W/m2 at model level, not the 85 predicted from file values |
| `aqh` flood-filled | — | Done | ~21 W/m2, the larger of the two levers |
| `atemp` flood-filled | — | Done | same fill, same erode |
| Filename convention | `ERA5_u10m_filled_1992` | Done | suffix before the year, for `useExfYearlyFields` |
| Years filled | 1992, 1993 | Partial | 1994-2019 outstanding |

`theta_max` drift fell from +0.527 degC/day to +0.0098, and in the current run it is flat:
28.218 to 28.196 over 2.25 days.

`ERA5/` in the run directory is a symlink into `/nobackupp17/hzhang1/forcing/era5`, another
user's read-only space, so filled files are written to
`/nobackupp27/rsavelli/LOAC/GoM_1km/forcings/era5_filled` instead.

---

## River forcing — GloFAS runoff and temperature

Runoff was cut by a factor of 4.2 and the salinity collapse survived it. `exf_runoff_max` is
now 4.134e-5 m/s, down from 1.729e-4, which bought an e-folding time of 17.1 h at the affected
cell instead of roughly 3 h: the same outcome, slower.

| Change | Before to after | State | Note |
| --- | --- | --- | --- |
| Regridding routine | `crop_glofas_for_exf.py` to `glofas_to_mitgcm.py` | Done | grid-agnostic: NetCDF grid reading, spherical kd-tree, great-circle Gaussian footprint, participation-ratio area |
| Runoff magnitude | 1.729e-4 to 4.134e-5 m/s | Partial | `salt_min` still decays to zero |
| River temperature | added, max 32.435 degC | Partial | `make_glofas_river_temp_forcing.py`; years past 1997 still needed, to 2026 |
| Conservation check in the regridder | kept-sources-only to offered vs delivered, plus a halo report | Done | the old form accumulated only kept sources, so a dropped river could never show |
| `EXFroff` / `EXFroft` diagnostics | on to off | Open | re-enable before production |

### Why reducing runoff did not fix it

The budget closes exactly. `empmr_min` is -0.0400 kg/m2/s = 4.004e-5 m/s = 3.46 m/day, and it
matches `exf_runoff_max` almost exactly, so it is runoff and not precipitation (~1e-6) or
evaporation (~2e-7). Solving `dS/dt = -S*F/h` for the observed 17.1 h e-folding gives h = 2.47 m,
recovering the top-cell thickness from the decay rate alone.

The cell is net-**divergent**: volume leaves, so no salty water flows in to replace it, and the
only thing restoring salt is lateral diffusion, of which there is almost none at 1 km. So the
equilibrium is near zero rather than near 5 or 10, and reducing R alone only slows the approach.

Three levers actually change it: enlarge the footprint (area x ~14, so sigma x ~3.7), distribute
the runoff vertically over ~10 m instead of 2.5 m (x4 on its own), or raise the shallow-cell
cutoff so the collapsing cell stops existing — it survived the `maskshallow` pass at an effective
depth of 2.47 m.

A suspended-sediment / TSS tracer from GlobalNEWS was evaluated for its effect on light
attenuation and dropped before implementation.

---

## Bathymetry and vertical grid

The bathymetry was edited three ways, all recorded in the file name
`GoM_1km_bathymetry_8_maskshallow_lagoons_deepenNE.bin`:

- the NE side of the domain was deepened
- some lagoons were filled and blocked
- very shallow cells were removed

`diffKr` was changed as well: `diffKrFile = 'GoM_1km_it42_diffkr_r8.data'` supplies a 3-D field
over a background of `diffKrT = diffKrS = 1.e-5`. What was separately established, by reading
`ggl90_calc_diff.F`, is that `diffKr` acts as a working floor rather than a dead parameter under
GGL90 — that is why editing it has an effect, not a reason it was left alone.

The vertical grid and partial-cell settings were examined and left as they stand:

- `delR(1) = 1.00` m, then 1.14, 1.30, 1.49, 1.70. The top cell is 1 m, and river water mixes
  over the top two to three levels; the salinity decay fit puts the effective depth at 2.47 m.
- `hFacMin = 0.3`, `hFacMinDr = 1.`, `hFacInf = 0.1`, `hFacSup = 5.`. rStar can expand a column
  fivefold, so the river cell is not hitting the clip, and `eta` confirms it: domain max 0.55 m
  and stable, rather than ballooning by the metres of water being added.

One connection worth making: filling and blocking lagoons and removing very shallow cells
targeted exactly the class of cell that is now collapsing in salinity, one too enclosed to flush.
So either the river cell survived that pass or it is a cell the pass did not cover. `empmr_min`
has held essentially constant at -0.040 kg/m2/s through the current run, so whichever cell it is,
it is still there and still being fed.

Outstanding here: 303 wet cells carry `diffKr = 0` and still need re-patching.

---

## Open boundary conditions

Two independent faults, found in this order: the period in `data.exf` was monthly when the files
are daily, and the generator had been silently dropping 32 of 34 years.

### The merge bug

```python
if itr == itrs[0]: dv_diag_all = dv_diag
else:              dv_diag = concat([dv_diag_all, dv_diag])
```

`dv_diag_all` is assigned once and never updated, so the loop ends holding `concat(first, last)`.
Running that pattern over 34 segments returns exactly segments 0 and 33, which is 730 records:
`ETAN_south.bin` and `UVEL_south.bin` held 1992 glued to roughly 2025, with 32 years cut out of
the middle and nothing in a flat binary to show it. The same function also globbed without
sorting, so even a corrected merge could assemble the time axis out of order.

### The cadence, established twice

The parent's dv segments are spaced 26,280 iterations apart. At LLC270's 1200 s timestep that is
72 iterations per day and exactly 365.000 days per segment, so 730 records is 2 x 365 and the
cadence is daily.

Separately, from the data: the autocorrelation of `ETAN_south` gives lag 1 = +0.983 and
lag 182 = -0.788. At 14-day cadence lag 182 would be 6.98 years and come out strongly positive,
so bimensual is ruled out; monthly would run 730 records to 2053. One corroboration:
`_merge_dv_chunks` dates record 0 of segment 0 at (2-1) x 1200 s = 20 minutes, and
`obcsSstartdate2 = 002000` is exactly that.

### Changes

| Change | Before to after | State | Note |
| --- | --- | --- | --- |
| `obcsN/S/Eperiod` | 2629800.0 to 86400.0 | Done | verified in the STDOUT |
| Physics OBCS regenerated | 730 to **12,321** records | Done | all 114 files agree on the count, covering 1992-12-31 to 2026-09-24 |
| `ETAN_east.bin` | truncated to whole | Done | now 12,321 records like every other boundary file; the 904-byte fragment is gone |
| `read_dv_diags` | unfixed to fixed | Done | the regenerated record count is the proof it took: 12,321 is a full 34-segment merge, not first-plus-last |
| `UTILS_DIR` in `gen_obcs_fast.py` | hardcoded Point Dume path | Open | still imports Point Dume's `gen_obcs.py`, so Point Dume's own OBCS may carry the same hole |

### Two notes

**The record count is resolved.** The regenerated files hold 12,321 daily records, all 114 in
agreement, running 1992-12-31 to 2026-09-24 — that is 2,822 records, or 7.7 years, past the end
of the run. The 95-record shortfall against the 12,416 estimate is explained exactly by a final
partial dv segment of 270 days: the parent run itself only reaches 2026-09-24. OBCS is no longer
a ceiling on this configuration.

**Segment 22 is not a defect.** `578161 -> 604873` spans 26,712 iterations = 371 days, and it
holds 371 records, so the axis is continuous across it. It looked like a silent 6-day defect from
the spacing alone, before checking the record count against it.

Tools written for this: `obcs_cadence.py` (recovers the cadence by counting annual cycles,
validated against synthetic daily / 14-day / monthly files of identical size),
`check_obcs_records.py`, `gen_obcs_patch.py`.

---

## Initial conditions and pickups

The cold start from iteration 0 segfaulted on negative tracer concentrations, and that is fixed.
The initial state comes from the parent's pickup files at iteration 26352 — `hydrogThetaFile`,
`hydrogSaltFile`, `uVelInitFile`, `vVelInitFile`, `pSurfInitFile` and `GGL90TKEFile` — with
`nIter0 = 0`. One temperature artifact in that state is still unexplained.

| Change | Before to after | State | Note |
| --- | --- | --- | --- |
| Negative PTRACERS clipped | `ptracer19` (O2) min -8.858; 17 of 31 tracers negative | Done | `PTRACERS_initialFile` points pTr08-11 and pTr19-31 at `forcings/ptr_clipped/`, exactly 17 files. The other 14 still use the originals, which is what `clip_ptracers.py --tracers auto` does |
| Write target for the clipped files | in place to a separate directory | Done | `clip_ptracers.py` refuses to write into `--dir`: the inputs are symlinks into a shared pickups directory |
| `theta_min` = 0.7685 degC | unchanged | Open | too cold for this domain. Deep water in the Gulf and Caribbean is about 4 degC; 0.77 degC is Antarctic Bottom Water. Set at iteration 0 and barely moving, so an initial-condition artifact rather than a model instability |
| `Zinterp` bottom-gap fill | `tmp1[:,idB[0]]` to `tmp1[:,idB]` | Partial | fixed in the fast drivers via `_zinterp_fixed`, but the iteration-26352 fields may well predate that fix |

The original form filled only the first bottom-gap level, leaving the rest NaN, which the
unconditional `tmp1[np.isnan(tmp1)] = 0` then zeroed: spurious zero layers mid-water-column.
That is the most likely family of cause for the cold `theta_min`, though 0.7685 is not itself a
zero, so the connection is a hypothesis rather than a finding.

---

## Input file manifest

251 file references across the namelists. `check_run_inputs.py` parses them, expands the EXF
names over the years the run needs, and reports the earliest date at which forcing runs out.
Run it in the run directory before submitting.

| Group | Paths | Where they live |
| --- | --- | --- |
| Grid and bathymetry | `GoM_1km_bathymetry_8_maskshallow_lagoons_deepenNE.bin`, `GoM_1km_it42_diffkr_r8.data`, `delYFile` | run dir |
| Initial state (iter 26352) | `pickup_THETA`, `pickup_SALT`, `pickup_U`, `pickup_V`, `pickup_ETAN`, `pickup_ggl90` | symlinks to `GoM_highres/grid2/forcings/pickups/` |
| PTRACERS, 14 of 31 | `pickup_pTr{01-07,12-18}.0000026352.data` | symlinks to `GoM_highres/grid2/forcings/pickups/` |
| PTRACERS, 17 of 31 | `pickup_pTr{08-11,19-31}.0000026352.data` | `/nobackup/rsavelli/LOAC/GoM_1km/forcings/ptr_clipped/` |
| EXF, filled | `ERA5_filled/ERA5_u10m_filled`, `_v10m_filled`, `_spfh2m_filled`, `_tmp2m_degC_filled` | `LOAC/GoM_1km/forcings/era5_filled`, **1992-93 only** |
| EXF, unfilled | `ERA5/ERA5_pres`, `_dlw`, `_dsw`, `_rain` | symlink to `/nobackupp17/hzhang1/forcing/era5` (read-only) |
| EXF, rivers | `GloFas_GoM`, `GloFas_temp_GoM` | run dir |
| EXF, other | `iron_monthly_clim_Hamilton_kgFem2s_GoM_1km` (climatology, no `_YYYY`), `apCO2` | run dir / `GoM_highres/grid/forcings/` |
| OBCS | 114 files: `{SALT,THETA,UVEL,VVEL,ETAN}_{north,south,east}.bin` + `TRAC01-31_{north,south,east}.bin` | run dir |
| Diagnostics output | 10 `diags/<stream>/` directories | must exist; MITgcm will not create them |

`useExfYearlyFields = .TRUE.`, so every EXF name above gets `_YYYY` appended, except a
climatology read on a `*RepCycle`.

Every large input is currently a symlink into another tree. 37 files at 2.82 GB each (diffKr,
the five 3-D initial fields, all 31 tracers) plus two 31 MB 2-D files is **~104 GB** if they are
to be physically copied rather than linked. Provenance is currently split across three trees:
`GoM_1km/darwin3/run`, `GoM_highres/grid2/forcings/pickups`, and `LOAC/GoM_1km/forcings`.

### Where this run actually stops

`startDate_1 = 19930101`, so the run begins 1993-01-01. **One ceiling remains, and it is the
forcing, not the OBCS.**

**1994-01-01, model day 365.** The filled ERA5 exists for 1992 and 1993 only, so
`ERA5_filled/ERA5_u10m_filled_1994` does not exist, nor the `v10m`, `aqh` and `atemp`
equivalents, nor `GloFas_GoM_1994` and `GloFas_temp_GoM_1994`. At 4.5 model days per wall day
that arrives after about 81 wall days.

OBCS used to be the second ceiling at 1994-12-30. It is not any more: 12,321 records reach
2026-09-24, 7.7 years past the run's end.

### What finishing the forcing costs

ERA5 at N320 is 1280 x 640 hourly, so one field-year in real*4 is **28.7 GB**:

| Approach | Size for 4 fields x 26 yr |
| --- | --- |
| Fill the global fields, as done for 1992-93 | **2.99 TB** |
| Crop to the GoM domain first, 3-cell margin (74 x 59 of 1280 x 640) | **18 GB** |
| Crop with a 5-cell margin (78 x 69) | **20 GB** |

A 150-fold saving, so cropping is almost certainly the route. Two constraints on the crop: the
margin must exceed `--erode` plus EXF's own bilinear stencil, so at least 3 cells, and `data.exf`
needs its `*_lon0`, `*_lat0`, `*_nlon`, `*_nlat` updated for the four cropped fields, since they
will no longer be on the global grid the unfilled four still use.

One observation worth a decision: the four fields that were flood-filled are the four whose land
contamination comes from surface properties (wind, 2 m temperature, 2 m humidity). `swdown`,
`lwdown`, `precip` and `apressure` still come from the unfilled shared `ERA5/`, and the case for
filling them is weaker, since downward radiation depends on the atmosphere above rather than the
surface below. If the coastal cloud banding turns out to matter, `swdown` and `lwdown` are where
to look.

---

## data file settings

Three things changed in the course of this work: `bathyFile` and `diffKrFile` (both covered under
Bathymetry above) and `obcsN/S/Eperiod`, from 2629800.0 to 86400.0. Everything below is as read
from the 2026-10-08 STDOUT.

| Parameter | Value |
| --- | --- |
| `bathyFile` | `GoM_1km_bathymetry_8_maskshallow_lagoons_deepenNE.bin` |
| `diffKrFile` | `GoM_1km_it42_diffkr_r8.data` |
| `diffKrT` / `diffKrS` | 1.e-5 / 1.e-5 (background under the 3-D file) |
| `startDate_1` / `startDate_2` | 19930101 / 000000 |
| `nIter0` | 0 |
| `deltaT` | 30. s |
| `endTime` | 820540800. s (26.0 yr) |
| `nonlinFreeSurf` / `select_rStar` | 4 / 2 |
| `implicitFreeSurface` / `implicitDiffusion` / `implicitViscosity` | .TRUE. |
| `eosType` | JMD95Z |
| `rhoConstFresh` | 999.8 |
| `tempAdvScheme` / `saltAdvScheme` | 7 / 7 (OS7MP, monotonicity-preserving) |
| `PTRACERS_advScheme` | 31*33 (DST3 with flux limiter) |
| `multiDimAdvection` / `staggerTimeStep` / `vectorInvariantMomentum` | .TRUE. |
| `viscAr` | 5.6614e-04 |
| `viscC4Leith` / `viscC4Leithd` / `viscA4GridMax` | 2.15 / 2.15 / 0.8 |
| `useGGL90` | .TRUE. |
| `useRealFreshWaterFlux` / `convertFW2Salt` | .TRUE. / -1. |
| `hFacMin` / `hFacMinDr` / `hFacInf` / `hFacSup` | 0.3 / 1. / 0.1 / 5. |
| `delR(1..5)` | 1.00, 1.14, 1.30, 1.49, 1.70 m |
| `useOBCS` / `useOBCSprescribe` / `useOBCSbalance` | .TRUE. |
| `obcsN,S,Eperiod` | 86400.0 s (was 2629800.0) |
| `obcs*startdate1` / `startdate2` | 19921231 / 002000 |
| `useExfYearlyFields` / `useRelativeWind` | .TRUE. / .TRUE. |
| `exf_iprec` | 32 |
| `readBinaryPrec` / `writeBinaryPrec` | 64 / 32 |
| `globalFiles` | .TRUE. |

### Known-bad, not yet fixed

- `data.diagnostics` has 10 active streams, 1 to 10, with blocks 11 and 12 (`EXFroff`, `EXFroft`)
  commented out. That is consistent, so the `DIAGNOSTICS_SET_POINTERS: 0 is not a Diagnostic`
  fault is **not** present in this configuration. It returns the moment 11 and 12 are re-enabled
  without raising `numlists` and rebuilding `DIAGNOSTIC_SIZE.h` — those two are coupled.
- The old THETA/SALT-into-`vel_3D_mon_mean` mismatch is gone. Block 10 is now `diags/exf_flux`
  carrying EXFwspee, EXFhl, EXFhs, EXFlwnet, EXFswnet, EXFatemp, EXFaqh and EXFqnet, and every
  block's fields match its file name.
- `data.darwin`: `darwin_strict_check = .FALSE.`, which turns a tracer-range message into a
  segfault. Still `.FALSE.` in the current run, and it is why the negative-oxygen cold start
  crashed instead of reporting.
- `deltaT = 30.` against a maximum advective CFL of 0.119 leaves a factor of four to five unused.

---

## Watch list for the job now running

`sss_del2` is the one to watch: it rose three orders of magnitude in the first six hours and has
since levelled at about 2.7e-5, which says the fresh anomaly is confined to one cell. If it
climbs again, or `uvel_max` jumps past ~2 m/s (the density artifact from near-zero salinity), the
bad cell has started contaminating its neighbours and the run is no longer telling you anything
useful.

| Watch | Healthy | Tripwire |
| --- | --- | --- |
| `theta_max` | flat near 28.2 | any resumed upward trend |
| `uvel_max` / `vvel_max` | <= ~1.2 m/s | a jump past ~2 m/s |
| `sss_del2` | plateaued ~2.7e-5 | growing again |
| `salt_mean` | 35.163 | any drift |
| `advcfl_W_hf` | ~0.12 | > 0.3 |
