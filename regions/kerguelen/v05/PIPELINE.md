# How the Kerguelen inputs were generated

The repository holds configuration text only. Everything the model reads as
data is produced by the `ecco_darwin/regions/downscaling` toolkit, in four
steps. This file records what was run for Kerguelen and the choices that are
specific to it. For the generic description of each step see
`regions/downscaling/STEP{1,2,3,4}.md`.

All paths below are relative to the Kerguelen configuration directory
(referred to as `<config>`), which is separate from this repository.

---

## STEP 1 — regional grid and bathymetry

Produces `<config>/grid/`.

| script | output |
|---|---|
| `run_step1_mitgrid.sh` | `kerguelen.mitgrid` |
| `job_step1_bathy.sh` | `kerguelen_bathymetry.bin` |
| `mitgrid2tiles` | `grid/tiles/` |
| `job_step1_stitch.sh` | `kerguelen_ncgrid.nc` (~7 GB) |
| `job_step1_dvmasks.sh` | `parent/inputs/{east,west,north,south}_BC_mask.bin` |

The grid is a plain uniform lat/lon box (`delX`/`delY` are both constant), so
unlike GoM_1km no `delYFile` is needed.

The dv masks select the parent (LLC270) cells that feed each regional open
boundary. They are the link between STEP1 and STEP2: they must be regenerated
whenever the regional grid changes, and the parent extraction must then be
re-run, because the extracted files are indexed by mask point and carry no
record of the geometry they came from.

---

## STEP 2 — parent-side extraction

Produces `<config>/parent/outputs/OBCS/` (~60 GB).

A `pkg/diagnostics_vec`-enabled LLC270 Darwin run, restarted from the
2020-01-01 pickup (`nIter0 = 736344`) and integrated one year, writing hourly
boundary values at the dv mask points only.

* `parent_run/code_darwin/DIAGNOSTICS_VEC_SIZE.h` and
  `parent_run/input/data.diagnostics_vec` are a **matched pair**: `nVEC_mask`
  must equal the number of `nml_vecFiles` entries (4 boundaries x 5
  mask-entries = 20). Changing one without the other is a build/run mismatch.
* `VEC_points` bounds the **per-process, per-mask** point count, not the
  global count. Sizing it against the global total is harmless but is the
  wrong quantity.
* Output is one file per field per boundary, named from the mask file's own
  basename, e.g. `east_BC_mask_THETA.<iter>.bin`.
* A restart writes a **separate** file set tagged `<nIter0+1>` whose record
  numbering restarts at 1. Segments must be stitched, not concatenated, and a
  segment that overran its successor's pickup leaves overlap records at its
  tail to discard.

---

## STEP 3 — regional boundary conditions and initial conditions

### OBCS files — `job_gen_obcs_fast.sh` -> `<config>/forcings/OBCS/` (164 files, ~100 GB)

`gen_obcs_fast.py` interpolates each parent boundary point onto the regional
boundary: 2-D `griddata` horizontally, 1-D linear vertically onto the regional
`drF`, with a nearest-neighbour fallback outside the parent convex hull.
Output is `{VAR}_{bnd}.bin`, e.g. `THETA_east.bin`.

### Barotropic transport correction — `job_gen_obcs_transport_correction.sh`

Rewrites the `UVEL_*`/`VVEL_*` files in place.

Nothing in the interpolation enforces that the volume transport crossing the
regional boundary matches the transport that crossed the same physical line in
the parent. Across a ~23x jump in along-boundary point density the residual
accumulates, and under `nonlinFreeSurf`/`exactConserv` it appears as secular
ETAN drift over a 3-year integration. The correction adds a uniform
(barotropic) velocity increment at each output time so the net transport
matches, preserving the interpolated baroclinic structure.

One Kerguelen-specific subtlety, worth knowing before reusing this: each
parent sample point's share of the boundary length is derived from its
position along the **true boundary coordinate** (latitude for east/west,
longitude for north/south), normalised to the regional edge length — *not*
from the parent's own native DYG/DXG. On this rotated LLC face the native
metrics overcount total boundary length by 2.16x–3.11x, differently per
boundary, which produced a net-transport bias that did not cancel between
opposite boundaries.

### Initial conditions — `job_gen_pickups.sh` -> `<config>/forcings/pickups/` (78 files, ~50 GB)

`gen_pickups.py` interpolates the parent state at iteration 736344 onto the
regional grid, producing `pickup_{THETA,SALT,U,V,ETAN,...}.0000736344.data`
plus the ptracer and sea-ice fields. `input/data` consumes these through
`hydrogThetaFile`/`hydrogSaltFile` etc. as a cold start at `nIter0 = 0`.

Regional wet cells with no parent ocean within the search radius are filled by
the `-oh nearest` fallback, with `-mwo` asserting the exact expected count so
the job fails rather than back-filling somewhere new if that count ever moves.

---

## STEP 4 — regional run

See `readme_pleiades.txt`. Built from `code_darwin` + `code`; run with
`scripts/job_production_chain.sh`, which chains ~6–7 segments under the
`long` queue's 120 h ceiling to cover the ~28 days of compute.

---

## Regeneration order

The dependencies are strictly sequential, and the important property is that
**a change upstream invalidates everything downstream silently** — the arrays
keep their shapes and the model still starts:

```
grid  ->  dv masks  ->  parent extraction  ->  OBCS + pickups  ->  run
```

If the regional grid is rebuilt, every one of the later products must be
regenerated. Dimension and file-size checks will not detect the staleness,
because only the coordinates change, not the shapes.
