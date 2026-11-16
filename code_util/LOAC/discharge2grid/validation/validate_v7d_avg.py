#!/usr/bin/env python3
"""V7d (2026-11-16): is the V7b difference (package climatology vs production temperature map) only the
production map's 3x3 averaging? Emulate it: each day every runoff cell gives its nearest-cell climatology
value to itself and its 8 neighbours, each cell takes the mean of what it receives (as in
make_glofas_river_temp_forcing.py), then compare with the production map at the runoff cells."""
import os
import sys
import time

import numpy as np

PKG = "$DATA/LOAC/discharge2grid_pkg"
sys.path.insert(0, PKG)
from discharge2grid import ModelGrid, Mapper  # noqa
from discharge2grid.readers import glofas, temperature as rt  # noqa

OUT = os.environ.get("OUT", "$DATA/LOAC/d2g_validation")
MAP = "$DATA/LOAC/glofas/river_runoff/GloFas_2015"
TMAP = "$DATA/LOAC/glofas/river_temp/GloFas_temp_2015"
GRID = "$DATA/LOAC/westAO/WestAO_grid.nc"
T0 = time.time()

grid = ModelGrid.from_netcdf(GRID, mask="mask", verbose=False)
mapper = Mapper(grid, sigma_km=18, max_snap_km=95, exclude_outside=True)
src = glofas.runoff_map_sources(MAP, mapper, year=2015, temp_path=TMAP, verbose=False)
clim = rt.TemperatureClimatology.load(os.path.join(OUT, "v7_wattmp_clim_35N.npz"))
ts, _ = clim.sample(src.lon, src.lat, src.time)
# every GloFAS cell that is non-zero in the map on some day in the band (production spreads from ALL of
# them, not only the ones near this grid): read the map band directly
mm = np.memmap(MAP, dtype=">f4", mode="r").reshape(-1, 3000, 7200)
j = np.rint((src.lat - glofas.GLO_LAT0) / glofas.GLO_INC).astype(int)          # south-to-north index
i = np.rint((src.lon - glofas.GLO_LON0) / glofas.GLO_INC).astype(int)
row = 3000 - 1 - j                                                             # map row (north first)
r0, r1 = row.min() - 2, row.max() + 3
eq = tot = 0
worst = 0.0
for k in range(src.nt):
    band = np.asarray(mm[k, r0:r1, :], dtype=np.float32)
    ar, ac = np.nonzero(band != 0)                                             # production "active" cells
    # their climatology value: nearest usable cell, this day's slot
    lat_a = glofas.GLO_LAT0 + glofas.GLO_INC * (3000 - 1 - (ar + r0))
    lon_a = glofas.GLO_LON0 + glofas.GLO_INC * ac
    va, _ = clim.sample(lon_a, lat_a, src.time[k:k + 1])
    va = va[0]
    s = np.zeros(band.shape)
    c = np.zeros(band.shape, dtype=np.int32)
    for dy in (-1, 0, 1):
        for dx in (-1, 0, 1):
            rr, cc = ar + dy, (ac + dx) % 7200
            m = (rr >= 0) & (rr < band.shape[0])
            np.add.at(s, (rr[m], cc[m]), va[m])
            np.add.at(c, (rr[m], cc[m]), 1)
    emul = np.divide(s, c, out=np.zeros_like(s), where=c > 0).astype(np.float32)
    act = src.q[k] > 0
    e = emul[row[act] - r0, i[act]]
    p = src.temp[k, act].astype(np.float32)
    eq += int((e == p).sum())
    tot += int(act.sum())
    worst = max(worst, float(np.abs(e - p).max()) if act.any() else 0.0)
print(f"V7d RESULT: production temperature map reproduced by package climatology + 3x3 averaging on "
      f"{eq:,} of {tot:,} source-days ({100 * eq / max(tot, 1):.3f} %); max |diff| {worst:.5f} C "
      f"({time.time() - T0:.0f} s)", flush=True)
