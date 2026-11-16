# discharge2grid

Put river discharge from any dataset onto any ocean-model grid. It works on regional MITgcm
grids, ECCO LLC grids and any curvilinear grid, and it conserves volume exactly. It can also
put river water temperature on the same cells (flux-weighted, conserving heat). It handles
freshwater and temperature only, not biogeochemistry.

This is the regridding we use to make the river forcing for two MITgcm/DARWIN set-ups
downscaled from ECCO-Darwin: pan_AO (pan-Arctic) and WestAO (Bering, Gulf of Alaska, Beaufort).
It is cleaned up so it can be shared, and split into two parts:

```
dataset reader  ->  SourceSet  ->  Mapper (snap + Gaussian footprint)  ->  writers
(GloFAS, tables,    n rivers,       grid-agnostic, volume-conserving       MITgcm binary
 gridded fields)    Q(t) in m3/s                                           sparse netCDF
```

Readers know about one dataset and nothing about the grid. The mapper knows about the grid
and nothing about the dataset. To support a new discharge product, write one function that
returns a `SourceSet` (see [docs/INSTRUCTIONS.md](docs/INSTRUCTIONS.md), step 4).

## Install and test

```bash
pip install -e .                 # numpy, scipy, netCDF4 (MITgcmutils only for MDS grids)
# or no install at all: export PYTHONPATH=/path/to/discharge2grid
python tests/test_mapping.py     # or: pytest tests/
python tests/test_glofas.py
python tests/test_readers_checks.py
python tests/test_ldd.py
python tests/test_temperature.py
python examples/demo_synthetic.py   # 10-second end-to-end demo on a toy grid
```

## Quick start

### Python

```python
from discharge2grid import ModelGrid, Mapper, write_mitgcm_runoff, write_sparse_netcdf
from discharge2grid.readers import glofas, points

grid = ModelGrid.from_netcdf("WestAO_grid.nc", mask="mask")       # any shape, incl. (tile, j, i)
src  = points.from_csv("rivers.csv", "discharge.csv")              # or a GloFAS / gridded reader
d    = Mapper(grid, sigma_km=18, max_snap_km=95, exclude_outside=True).map_sources(src)
d.report()                                   # what was kept, dropped, how far it was moved
write_mitgcm_runoff("runoff_WestAO_2015", d)       # dense >f4 records, m s-1
write_sparse_netcdf("WestAO_runoff_2015.nc", d)    # only the cells that receive water

# river temperature (optional): sample a gridded product at the sources, before mapping
from discharge2grid.readers import temperature
clim = temperature.TemperatureClimatology.load("wattmp_clim.npz")   # built once, see INSTRUCTIONS 3b
src  = temperature.add_temperature(src, clim)
# then map as above; d.temperature_points(k) / write_mitgcm_temperature(path, d)
```

### Command line (GloFAS)

```bash
# 1. find outlets once. Preferred: from the GloFAS drainage map (one exact outlet per river)
python -m discharge2grid outlets --ldd glofas_ldd.nc --lat-range 40 90 --out outlets.csv
#    without the drainage map: coast rule from several years of raw GloFAS (along-coast
#    duplicates are removed by default), plus dead-end candidates to check by hand
python -m discharge2grid outlets --raw 'glofas_original_{year}.nc' --years 1990,2005,2015 \
    --lat-range 40 90 --out outlets.csv --dead-ends dead_end_candidates.csv

# 2. map them, year by year (--extra: the dead ends you triaged as real, METHODS §2.5;
#    examples/glofas_extra_rivers_arctic_alaska.csv is our list, not complete for this path)
python -m discharge2grid run --grid grid.nc --source glofas-raw \
    --raw 'glofas_original_{year}.nc' --outlets outlets.csv \
    --extra my_dead_ends.csv --years 1992-2025 \
    --sigma-km 18 --max-snap-km 95 --exclude-outside-domain \
    --tag MYGRID --mitgcm-dir out_bin --sparse-dir out_nc --tables-dir out_tables

#    + river temperature: build a climatology once, then add --temp-clim to `run`
python -m discharge2grid temp-clim --bin 'wattmp*.bin' --lat-range 35 90 --out wattmp_clim.npz

# 3. check that no river drops out in some years
python -m discharge2grid dropouts --tables 'out_tables/sources_MYGRID_*.csv'
```

## Outputs

- **MITgcm binary**: `<prefix><tag>_<year>`, big-endian float32 records of the grid's shape,
  runoff in m s⁻¹. For EXF, set `runoff_interpMethod = 0` and use yearly files. For an LLC
  model, read the grid from the model's own MDS output (`--grid-mds`). Records are written in
  the shape the grid was read in, so a binary built from a tiled `(tile, j, i)` netCDF is in
  tile order, which is not MITgcm's compact LLC layout.
- **Sparse netCDF**: `runoff(time, point)` (m s⁻¹) on the cells that ever receive water, plus
  `point_index` (flat index into the grid) and one 0-based index per grid dimension (for
  example `tile, j, i`). It also stores the sources and the sparse operator
  `runoff = W · discharge` (`link_*`), so any value can be traced back to its rivers.
  On WestAO 2015 it holds 14 % of the cells: 137 MB with sources and operator, against
  854 MB dense (METHODS V4).
- **Temperature** (when the sources carry it): `runoff_temp_<tag>_<year>` in the same binary
  layout (°C, 0 where no water), and `temperature(time, point)` in the sparse file (NaN where
  no water). It is flux-weighted per cell, so the heat the rivers carry is conserved.

## Documentation

- [docs/METHODS.md](docs/METHODS.md): what the code does and why, GloFAS pitfalls (dead-end
  rivers, outlets that drop out in some years, duplicates), parameter choices, and validation.
- [docs/INSTRUCTIONS.md](docs/INSTRUCTIONS.md): step-by-step for a new grid and a new
  dataset, with ECCO/LLC notes.

## Status

Version 0.2.0 (November 2026): river temperature added to 0.1.0.

- The mapping core reproduces our production runoff files bit for bit, on synthetic input
  and on the real 2015 GloFAS input for WestAO (METHODS, V1).
- The GloFAS raw-file outlet path is new in this package. Real-data checks against our
  production input are in METHODS (V2–V6).
  - The plain coast rule counts rivers that run along the coast two or three times.
    Duplicates are removed by default.
  - The 27 rivers that end inland in GloFAS, which we found by hand, are all among the
    automatic candidates.
  - The fixed outlet list keeps rivers our production input dropped in some years.
  - Without the drainage map, small coastal cells leave about 10 % uncertainty in the
    water at 40–90 N (METHODS §2.2).
  - The drainage-map (LDD) outlet finder is tested on synthetic maps only; we did not have
    the GloFAS LDD file.
  - Logs, tables and scripts of the real-data checks are in `validation/`.
- River temperature: the legacy path reproduces our production temperature files bit for
  bit (METHODS V7). The climatology builder repeats our production recipe.
- Biogeochemistry is part of our production pipeline but not of this package.
