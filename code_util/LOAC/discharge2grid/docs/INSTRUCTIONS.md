# Instructions: a new grid, a new dataset

Steps 1–3 take GloFAS onto a grid. Step 4 adds another discharge product. Step 5 checks a
multi-year run. Step 6 collects the ECCO/LLC specifics.

## 0. Install and test

```bash
pip install -e .
python tests/test_mapping.py && python tests/test_glofas.py && python tests/test_readers_checks.py
python examples/demo_synthetic.py
```

Use a compute node for real data. A global GloFAS year is about 2.4 GB of compressed netCDF
(31 GB once decompressed), so on an HPC system run the steps below inside a batch job, not
on a login node.

## 1. The grid

```python
from discharge2grid import ModelGrid
g = ModelGrid.from_netcdf("grid.nc")             # guesses XC/YC/rA/hFacC-style names
g = ModelGrid.from_netcdf("grid.nc", lon="XC", lat="YC", area="rA", mask="maskC", level=0)
g = ModelGrid.from_mds("run/")                   # MITgcm grid dump (XC, YC, RAC, hFacC, DRF)
print(g.summary())                               # shape, dims, wet cells, cell size
```

- The arrays can have any rank. `(j, i)`, `(tile, j, i)` and compact `(1170, 90)` all work,
  because only great-circle distances are used.
- `mask` may be `hFacC`/`maskC` (wet = > 0; with a vertical dimension, `level` picks the
  layer) or a depth/bathymetry (wet = > 0, or < 0 if the file is negative-down).
- **Use the mask the model really treats as wet at the surface.** Our WestAO grid has
  wetting-drying, and its `model_mask` includes floodable land. With that mask, rivers
  landed on mudflats, so we use `mask`, which excludes them.
- Grids with blank (all-land, unwritten) tiles whose XC/YC are 0 or missing make land under
  them look "outside the domain". Check the "dropped as outside the domain" list.

## 2. GloFAS outlets (once per dataset version)

**Preferred: the GloFAS drainage map.** GloFAS publishes its static maps alongside the
model. The local drain direction map (`ldd`, PCRaster codes 1–9, 5 = pit) is in the same
collection as the pixel-area map (`GLOFAS_pixarea_Global_03min.nc` on our system); the
file name depends on where you get it. The variable name is detected (`ldd`, `Band1`) or
can be passed with `--ldd-var`.

```bash
python -m discharge2grid outlets --ldd glofas_ldd.nc --lat-range -60 90 \
    --raw '/data/glofas_original_{year}.nc' --years 2015 --out outlets.csv
```

Every river ends in exactly one pit, so this gives one outlet per river: no duplicates and
no dead-end triage. `kind` in the output:
- `mouth`: the outlet touches the sea.
- `near_sea`: a pit within 250 km of the sea. This is a river that stops short of the coast
  (the Ob); it is delivered at the pit and snapped to water.
- `inland`: an endorheic pit; `use=0`.

**This path is tested on synthetic maps only.** Before relying on it, check on real data:
- The number of mouths is plausible.
- The Ob, Yukon, Mackenzie and Copper come out right.
- Domain totals agree with the coast rule below to within ~10 %.

**Without the drainage map: the coast rule.**

```bash
python -m discharge2grid outlets --raw '/data/glofas_original_{year}.nc' \
    --years 1990,2000,2010,2020 --lat-range -60 90 \
    --out outlets.csv --dead-ends dead_end_candidates.csv
```

Every land cell that touches the sea becomes an outlet, small coastal catchments included.
The sea is where GloFAS is NaN on every scanned day. The cells are found once, from the
years you scan, and that list is then read for every year (a fixed list; METHODS §2.3).
- Scan several years spread over the period. A river whose last cell moves between years
  then has all its cells in the list. In any given year only one of them carries the river,
  so nothing is counted twice.
- **Along-coast duplicates are removed by default** (`dup_of` ≥ 0 → `use=0`). Where a
  channel runs along the coast, two or three coastal cells carry the same river. Keeping
  them counted the Mackenzie twice, and WestAO came out 43 % above our production forcing
  instead of 10.6 % (METHODS V4). Use `--keep-duplicates` only to inspect them.
- The coast rule keeps every small coastal cell. Whether neighbouring small cells are
  separate catchments or drain into each other cannot be told without the drainage map.
  That is about 10 % of the water at 40–90 N, concentrated on fjord coasts (METHODS §2.2).
  Report it as an uncertainty.
- `--min-q` drops tiny cells. The command prints how much of the water each size class
  carries; decide from that.

**Dead ends (needs your eyes).** These are rivers whose GloFAS channel stops inland, so no
coastal cell carries them; the Ob, for example, ends 9 km from the sea. On 2015 at
40–90 N the list had 136 candidates ≥ 100 m³ s⁻¹, including all 27 rivers we had found by
hand. It also contains false alarms:
- Deltas and estuaries where discharge dips slightly downstream. The Yukon at
  -159.975/62.275 is **not** a dead end.
- Mid-river points of big rivers. These show up as "sea within" 100–250 km.
- Repeats of one river.

**Triage the whole list, not just our 27.** We built our list against our production
input, which already contained some rivers that the raw coast rule treats as dead ends:
the Indigirka (1,900 m³ s⁻¹), a James Bay river (1,244), the Alsek (827) and the Kemi (574).
Each stops 7–11 km short of the raw coastline and is in the candidate list (METHODS §2.2).
Triage each row against a map and against GloFAS itself (for example: is there a coastal
outlet with the same discharge nearby?). Set `use=1` on real dead ends, and pass the file
to `run` with `--extra`. Set `deliver_lon/deliver_lat` where the water should enter the sea
if that is not the nearest water; we deliver the Thelon at Chesterfield Inlet.
`examples/glofas_extra_rivers_arctic_alaska.csv` is our triaged list for the Arctic and
Alaska/BC (27 rivers).

## 3. Run

```bash
python -m discharge2grid run --grid grid.nc --grid-mask maskC --source glofas-raw \
    --raw '/data/glofas_original_{year}.nc' --outlets outlets.csv --extra my_dead_ends.csv \
    --years 1992-2025 --sigma-km 18 --max-snap-km 95 --exclude-outside-domain \
    --tag MYGRID --mitgcm-dir out_bin --sparse-dir out_nc --tables-dir out_tables
```

Choosing the parameters (METHODS §3–4):

| parameter | what it does | our value | guidance |
|---|---|---|---|
| `--sigma-km` | Gaussian footprint scale | 18 km (1–4 km grids) | at least the cell size; size it on the peak-flow residence time |
| `--max-snap-km` | farthest a source may be moved to water | 95 km | large enough for estuaries and fjords the grid does not resolve |
| `--exclude-outside-domain` | drop sources beyond the grid's edge | on | needed for regional grids once max-snap is large |
| `--domain-tol-km` | how far off the grid a source may sit | 1.5 × mean cell | raise it if a real river on the edge is excluded |
| `--cutoff-sigma` | footprint truncation | 3.5 | leave it |

Read the report for each year:
- `conservation: delivered / kept` must print 1.000000000000.
- **snapped > 25 km**: check that each of the largest lands in the right basin.
- **dropped as too far**: raise `--max-snap-km` if these rivers belong to the domain.
- **dropped as outside the domain**: raise `--domain-tol-km` if one of them drains into the
  domain.
- **footprint area** far below 4πσ² means the footprint is clipped by the coast. Volume is
  still conserved, but it is concentrated into fewer cells.

## 3b. River temperature (optional)

GloFAS has no water temperature, so it comes from a gridded river-temperature product.
Build a daily calendar climatology from that product once, then sample it at the sources in
every run (METHODS §2.9 and §3.5):

```bash
# binary daily files; the default layout is the JRA55-do-grid wattmp files we used
# (little-endian, north first, Kelvin, 1440 x 720 at 0.25 deg); --grid-spec, --dtype,
# --south-first and --celsius describe other layouts
python -m discharge2grid temp-clim --bin '/data/wattmp*.bin' --lat-range 35 90 --out wattmp_clim.npz
# or netCDF files with dated records (units read from the file, K or deg C)
python -m discharge2grid temp-clim --nc '/data/rivertemp_*.nc' --var tw --out clim.npz

python -m discharge2grid run ...same as step 3... --temp-clim wattmp_clim.npz
```

- `run` then writes `runoff_temp_<tag>_<year>` (°C, 0 where no water; prefix with
  `--temp-prefix`) next to the runoff file, and `temperature(time, point)` in the sparse
  netCDF.
- `--lat-range` keeps the climatology small. Pad it a few degrees beyond the grid so sources
  near the edge still find their nearest usable cell.
- Sources without a temperature (beyond `--temp-max-km`, or rivers given with `--extra` in
  the legacy path) take the temperature of the largest known source within
  `--temp-fill-km` (150 km).
- **In the report check that `heat out/in` is 1.000000000000.** If not, some water still
  has no temperature; the report says how much.
- **A new temperature product.** Check byte order, row order, longitude origin and
  Kelvin/Celsius by mapping its valid cells and comparing known rivers. Each of these was
  wrong once with our files and still gave plausible numbers.
- For the table source, give temperatures with `--temp-table` (same layout as `--series`).
- For the legacy map, give the matching temperature map with `--temp-map`.

## 4. A new dataset

Write a function that returns a `SourceSet`:

```python
import numpy as np
from discharge2grid import SourceSet

def read_my_product(path, year):
    ...                                  # open the file
    return SourceSet(
        lon=mouth_lon, lat=mouth_lat,    # (n,) where the water ENTERS the sea
        q=discharge_m3s,                 # (ntime, n) m3/s; NaN/negative -> 0
        time=times,                      # (ntime,) numpy datetime64 (daily, monthly, ...)
        name=river_names, kind=["my_product"] * n,
        attrs=dict(dataset="my product v1", file=path),
        temp=None)                       # optional (ntime, n) deg C, NaN = unknown
```

Checklist for a new product:
1. **Units → m³ s⁻¹.** Flux products (kg m⁻² s⁻¹, m s⁻¹, mm day⁻¹) need the cell area:
   `readers.gridded.from_field` does it. Check the file's units attribute against the
   numbers. On our system the JRA55-do river files were in m s⁻¹, not the kg m⁻² s⁻¹ we
   expected.
2. **Locations are mouths.** Gauges are inland. Move them to the mouth (and scale for the
   ungauged area) or the snap will be long.
3. **One source per river.** Routing products carry the river at every cell along its
   course. Take outlets only (as `readers.glofas` does), never the whole field.
4. **A fixed source list for the whole run.** Zeros are fine; appearing and disappearing
   sources are not (METHODS §2.3).
5. **Time.** Records can be daily or monthly. The writers write one record per time step
   and the model namelist sets the period. Leap days come from the product.
6. **Temperature (optional).** `temp` (ntime, n) in °C, NaN where unknown. Or leave it out
   and add it from a gridded product with `readers.temperature.add_temperature`.

Existing readers: `readers.points.from_csv` / `from_netcdf` (river tables, Dai &
Trenberth-style station files), `readers.gridded.from_field` / `from_netcdf` (coastal runoff
fields such as JRA55-do `friver`), and `readers.glofas` (raw outlets, and the legacy runoff
map).

## 5. Checks over a multi-year run

```bash
python -m discharge2grid dropouts --tables 'out_tables/sources_MYGRID_*.csv'
```

This lists rivers ≥ 20 m³ s⁻¹ whose yearly mean falls below 30 % of their median in some
years.
- `lost` ≈ 1: the water is really missing (the river vanished from the product).
- `lost` ≈ 0: it only moved to a nearby cell (no action).

With the fixed-outlet GloFAS path this list should be empty or nearly so; with other
products it is how you find problems. Also compare domain totals with published values
(e.g. Arctic runoff ≈ 4,400 km³ yr⁻¹).

## 6. ECCO / LLC notes

- **Grid.** ECCO grid netCDF files (`XC, YC, rA` as `(tile, j, i)`, `maskC`/`hFacC` as
  `(k, tile, j, i)`) read directly: `ModelGrid.from_netcdf(f, lon="XC", lat="YC", area="rA",
  mask="hFacC")`. Dimension names are carried to the sparse output as `tile, j, i` indices.
- **Sparse netCDF → model.** MITgcm EXF reads dense binaries, not the sparse file. Expand a
  record with `read_sparse_netcdf(p).to_dense(k)`, or write the binary directly. For an LLC
  binary, build the grid from the model's own MDS output (`--grid-mds`) so the records come
  out in the model's compact layout. The tile order of the netCDF grid file is not that
  layout.
- **Coarse grids (LLC90 ≈ 1°, LLC270 ≈ ⅓°).** A sigma well below the cell size puts each
  river into its nearest wet cell. A sigma of about one cell spreads big rivers over a few
  cells, which avoids near-zero salinity at the mouth. Check the report's footprint area and
  the residence-time argument in METHODS §3.3.
- **Salinity.** If mouth salinity goes negative, the cause is the salt-flux namelist, not
  the forcing. Use `useRealFreshWaterFlux = .TRUE.` with a nonlinear free surface, or
  `convertFW2Salt = -1` (METHODS §3.3).
- **Speed.** The operator is built once per grid and source list, so later years reuse it.
  Reading GloFAS dominates the run time: one hyperslab per storage chunk
  (`readers.glofas.chunk_points`), never point-by-point indexing.
