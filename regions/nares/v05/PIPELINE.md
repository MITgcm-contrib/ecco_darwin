# How the Nares inputs were generated

The repository holds configuration text only. Everything the model reads as
data is produced by the `ecco_darwin/regions/downscaling` toolkit plus the
Nares scripts in `scripts/pipeline/`. This file records what was run for Nares
and the choices that are specific to it. For the generic description of each
step see `regions/downscaling/STEP{1,2,3,4}.md`.

All paths below are relative to the Nares configuration directory on Pleiades
(`<config>` = `/nobackup/dcarrol2/downscaling/nares/`), which is separate from
this repository. Dates refer to the decision log.

**DRAFT: pending the frozen physics (run_arch_D10).** The inputs themselves do
not depend on the sea-ice tuning, but the final OBCS and IC sets for a
2016-2025 production run do not exist yet.

---

## STEP 1 — regional grid, bathymetry, shelf, tiling

Produces `<config>/grid/`.

| script | output |
|---|---|
| `grid/gen_nares_mitgrid.py` | `Nares.mitgrid` (16 fields, 1121 x 361 real*8, ONE facet, also linked as `tile001.mitgrid`) |
| `grid/gen_bathy_nares.py` | `Nares_bathymetry.bin`, `Nares_shelficeTopo.bin`, `Nares_bathymetry_stats.json` |
| `grid/wall_fill_bridges_nares.py` | `Nares_bathymetry_walled.bin` (the bathymetry the runs use) |
| `gen_ncgrid/` (code, input, 2 PBS jobs) | `Nares_ncgrid.nc` (1-step dump on all 252 tiles, no blankList) |
| not sourced (DECISIONS 2026-10-15) | `Nares_grid18.face001.bin` = the 16 mitgrid fields + AngleCS/AngleSN from `Nares_ncgrid.nc` |
| `grid/gen_exch2_blanklist.py` | `input/data.exch2` (28 x 9 = 252 tiles of 40 x 40, 60 blank, x-fastest ids) |
| `grid/gen_ob_index.py` | the `OB_I{east,west}` / `OB_J{north,south}` arrays in `input/data.obcs` |
| `utils/gen_dvmasks.py` (toolkit; no Nares wrapper found) | `parent/inputs/{east,west,north,south}_BC_mask.bin` |
| `grid/gen_humboldt_plumes.py` | iceplume mask + length files (baseline, sheet, no-discharge variants) and `Humboldt_Qsg_daily.npz` |

Grid: a rotated lat-lon grid whose rotated equator is the great circle through
A = (83.6 N, 55 W, Lincoln Sea) and B = (76.0 N, 73.5 W, North Water), 180 km
each side, extended 60 km past A and 120 km past B (to reach Devon Island and
cover the whole North Water), padded to nx = 1120; ny = 360. All metrics are
analytic on the sphere (dx 999.6 to 1000.0 m). i runs S to N along the strait.
`usingCurvilinearGrid = .TRUE.`; `delX`/`delY` are not used.

Bathymetry: GEBCO 2026 exact cell-mean binning in rotated coordinates, wet
where the mean elevation < 0, 10 m minimum depth (deepen only), largest
4-connected component, side-edge wet runs shorter than 30 cells (fjord
mouths) walled off. GEBCO's ice-surface grid reads floating tongues as land,
so the Petermann cavity comes from BedMachine Greenland v6: floating cells
with a column >= 20 m become shelfice cells (1200 cells, draft 9 to 586 m)
in `Nares_shelficeTopo.bin`. Other floating cells become land. No wetting and
drying.

Walling (2026-10-15): the 10 m rule made 571 near-dry GEBCO cells wet; they
formed thin "fill bridges" that short-circuited two parts of the domain and
carried 2 to 5 m/s spurious currents with tides. They are set to land in
`Nares_bathymetry_walled.bin`. No new wet cells, so the pickups generated on
the unwalled bathymetry are used unchanged. A 20 m minimum depth was tested
(`deepen_min_depth_nares.py`) and not adopted.

Tiling: a tile is active if it holds any wet cell (bathymetry < 0 or shelfice
cavity). 60 of 252 tiles are blank, so 192 MPI ranks. `SIZE.h` (nPx = 192),
`data.exch2` (blankList) and the bathymetry must change together; the build
scripts assert this.

Grid angles: the run reads `Nares_grid18.face001.bin` (`horizGridFile =
'Nares_grid18'`). Without it the model computes the angles itself and gets the
east and north edge rows wrong by up to 167 degrees (found 2026-10-15; the
first 8 days of tide tests ran with wrong edge winds).

---

## STEP 2 — parent-side extraction (shared with west_AO)

Produces `/nobackupp19/dcarrol2/downscaling/parent_V5r1_westAO_nares/run/dv/`
(4272 files on 2026-11-16, growing).

A `pkg/diagnostics_vec`-enabled LLC270 V5r1 Darwin run (767 ranks, 28 bro_ele),
restarted from the 2015-01-01 pickup (`nIter0 = 604872`) and chained in 48 h
segments (`chain_next.sh`) to 2025-12-31, writing hourly boundary values at the
dv mask points only. `scripts/pipeline/parent_run/` holds its `code_darwin/`,
`input/`, `job_build.sh`, `job_ECCO_darwin` (csh), `chain_next.sh` and
`check_dv.py`.

* `DIAGNOSTICS_VEC_SIZE.h` and `data.diagnostics_vec` are a **matched pair**:
  8 masks (4 west_AO + 4 Nares) x 5 groups = 40 `nml_vecFiles` entries.
* Field names are diagnostics_vec's internal names (`UICE`, `AREA`,
  `PTRACEnn`), not `SIuice`/`SIarea`.
* Output is one file per field per mask, e.g. `nares_east_BC_mask_THETA.<iter>.bin`,
  big-endian **float64**. A restart writes a separate file set tagged
  `<nIter0+1>` (e.g. `0000604873`).
* The parent has 50 levels; the Nares east mask has 144 parent points, the
  west 57, north 60, south 41 (counts from the first snapshot, 2026-10-07).

---

## STEP 3 — regional boundary conditions and initial conditions

### Physics OBCS — `obcs_ic/job_obcs_2015JanAug.sh` -> `<config>/forcings/OBCS2_rot_tf_JanAug/` (40 files)

Chain, as run (2026-10-23, same chain as `job_obcs_v2.sh` of 2026-10-09):

1. `snapshot_dv.py`: copy the first 5856 hourly records (2015-01-01 00:00 to
   2015-09-01) of each field into a fixed-length float32 file set, because
   `gen_obcs_fast.py` infers record counts from file size and the extraction
   is still growing.
2. `gen_obcs_fast.py -seaice -disc off`: fast replacement for the toolkit's
   `gen_obcs.py` (2-D `griddata` horizontally, 1-D linear vertically onto the
   regional `drF`, nearest-neighbor fallback). Fields: THETA, SALT, UVEL, VVEL,
   ETAN (not wired in `data.obcs`), SIarea, SIheff, SIhsnow, SIuice, SIvice,
   for east, west, north, south. The 2026-10-09 fix of its box-path masks
   matters: the first set had about 80 % of the north-edge V and 31 % of the
   Baffin-edge U zeroed.
3. `rotate_obcs_velocity.py`: the toolkit rotates parent LLC (u, v) to
   geographic east/north only. The Nares grid is rotated 43 to 83 degrees, so
   the components are rotated again into the grid axes (`uE = CS*u - SN*v`
   inverted), for ocean and sea-ice velocities. (`-> OBCS_rot`)
4. `obcs_boundary_budget.py` + `nares_throughflow_offset.py`: net transport by
   strait end. The parent carries 0.519 Sv southward through the boundaries
   (its own strait section: 0.39 Sv) against about 0.85 Sv observed. A uniform
   barotropic offset Q = 0.331 Sv southward is added over the wet east
   (Lincoln) and west (Baffin) faces, net zero, so `useOBCSbalance` leaves it
   alone. (`-> OBCS_rot_tf`; the "2" in `OBCS2_` marks the fixed masks.)
5. Jan to Jun of the Jan to Sep set is asserted identical to the earlier
   `OBCS2_rot_tf` (4344 records).

### Tides — `tides/` -> `<config>/forcings/OBCS_tides_2015Jan_c2/` (16 files)

`nares_ob_points.py` writes the open-boundary point list for the rotated
grid; `gen_obcs_tides.py` (shared tool, `downscaling/python/tides`) writes
`OB{W,E,S,N}_{u,v}Tid{Am,Ph}.bin` for M2 S2 N2 K2 K1 O1 P1 Q1 from FES2014
currents, phase zero at 2015-01-01 00:00, nodal corrections at 2015-01-16,
real*8; `check_tide_files.py` checks the net boundary inflow against the
FES2022 prism budget. Variant C2 (2026-10-15): amplitudes scaled per
constituent to the FES2022 prism (M2 0.917, S2 1.02, N2 0.99, K2 0.943,
K1 0.909, O1 0.99, P1 0.893, Q1 1), tangential amplitudes 0, normal amplitude
tapered 0, 1/3, 2/3 over the 3 cells next to each coast. Regenerate for any
other start year.

### Initial conditions — `obcs_ic/job_nares_arch_inputs.sh` -> `<config>/forcings/pickups/` (22 files)

* `gen_pickups.py -oh raise` (toolkit): parent state at iteration 604872 onto
  the regional grid: `pickup_{THETA,SALT,U,V,ETAN}`, `pickup_ggl90`,
  `pickup_seaice` (real*8).
* `nares_seaice_ic.py`: splits `pickup_seaice` into the single-field files
  `data.seaice` reads (AreaFile ... vIceFile). Without them the model
  cold-starts with no ice.
* `make_ic_thin_mouth.py` (2026-10-20, "ICT"): `pickup_seaice_siHEFF_ICT`, the
  HEFF field with the strait-mouth ice thinned to 2.0 m floe thickness, so the
  mouth is mobile on Jan 1 as in MODIS. Adopted into the baseline.
* `nares_tref_sref.py`: horizontal wet-cell means of the IC per level ->
  `tRef`/`sRef` in `input/data`.
* `nares_shelfice_load.py` (+ `jmd95.py`): `Nares_shelficeLoadAnomaly.bin`,
  the ISOMIP `gendata.m` port (JMD95Z, rhoConst 1027.5, partial-cell drloc).

### iceplume forcing — `iceplume/` (package: https://github.com/mhwood/iceplume at e4d4ac9 + local fixes) -> `<config>/forcings/iceplume/runoffQsg_YYYY` (12 files)

`gen_runoffQsg_files.py` (PBS `job_runoffQsg.sh`) expands the Mankoff v8
Humboldt daily discharge (`Humboldt_Qsg_daily.npz`, 59 marine outlets lumped
onto the 5 largest plume cells) into yearly 1120 x 360 real*4 files, m3/s at
the plume cells, read through EXF (`runoffQsgperiod = 86400`). Humboldt is
NOT also in EXF runoff. `useICEPLUME = .FALSE.` in every run so far.

### BGC — `bgc/job_bgc_inputs.sh` -> `/nobackupp19/dcarrol2/downscaling/nares_bgc_inputs/forcings/` (2026-10-29/30)

Same steps and rules as west_AO prod3/prod4:

1. `gen_pickups.py -bgc` from parent pickup 604872 -> `pickup_pTr01..31`
   (plus `pickup_darwin` for pH).
2. `gen_derived_tracer_pickup.py`: rDOC = 0.50 x DOC, CDOM = 0.02 x DOC
   (carbon CDOM, `CDOMcoeff 0.20`), CDOM clipped at 0.
3. `clip_bgc_negatives_nares.py`: inherited negative DOC, rDOC, DON, plankton
   and Chl clipped at 0, mass taken from the inorganic pools of the same cell.
4. `seed_plankton_nares.py`: where c_k < 1e-5 in the upper 100 m, set c_k to
   Xmin (c01-05 0.0125, c06-07 0.001), Chl by the median Chl:C; IC and ALL
   four edges (the parent has no plankton under the Lincoln Sea pack).
5. `snapshot_dv.py` for PTRACE01..31 (50 levels), `gen_obcs_fast.py -bgc`,
   `gen_derived_tracer_obcs.py`, clip and seed in OB mode -> 152 OBCS files
   (`TRACnn_{east,west,north,south}.bin`, rDOC, CDOM; about 206 GB).
6. `nutfix_nares.py`: the pan_AO nutrient-fix rule v3 (parent NO3/PO4/SiO2
   low in Nares vs WOA18 and GLODAP: SiO2 about 2x, PO4 about 1.5x). IC factor
   = min(FOB, 3), OB factor = min(FOB, WOA23 winter / Jan-Apr max) >= 1.

The BGC run reads the physics pickups and OBCS from STEP 3 above unchanged.

---

## STEP 4 — regional run

See `readme_pleiades.txt`. Physics built from `code`, run with
`scripts/job_run_arch_long.sh` (27.4 h for Jan 1 to Sep 1 2015 on 192 ranks).
BGC built from `code_darwin` + `code` (DRAFT, not yet built in that form) and
run as 3 segments with `scripts/stage_bgc_prod.sh` (about 227 h in all).

---

## Regeneration order

The dependencies are strictly sequential, and **a change upstream invalidates
everything downstream silently** (the arrays keep their shapes and the model
still starts):

```
mitgrid -> bathymetry + shelf -> ncgrid -> grid18 angles, exch2 blankList, OB index, dv masks
        -> parent extraction -> snapshot -> OBCS (interpolate, rotate, offset) + tides + pickups
        -> BGC pickups + OBCS -> run
```

If the regional grid or the bathymetry's wet mask changes, every later product
must be regenerated, including `data.exch2` and the `OB_I*`/`OB_J*` arrays in
`data.obcs`. The walling of 2026-10-15 was safe only because it removed wet
cells and added none.
