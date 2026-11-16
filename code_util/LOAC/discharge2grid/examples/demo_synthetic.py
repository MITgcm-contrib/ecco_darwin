#!/usr/bin/env python3
"""End-to-end demo on a toy grid -- no data needed. Run:  python examples/demo_synthetic.py

A 0.05-degree box off a straight coast, three rivers given as a station table (one far
inland and therefore dropped), 31 daily records, with river temperatures. Writes dense
MITgcm runoff and temperature files and a sparse netCDF to ./demo_out and checks that both
hold the same water and the same heat.
"""
import os
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from discharge2grid import (ModelGrid, Mapper, SourceSet, daily_axis, write_mitgcm_runoff,  # noqa
                            write_mitgcm_temperature, write_sparse_netcdf, read_sparse_netcdf)
from discharge2grid.geometry import regular_cell_area  # noqa

out = "demo_out"
os.makedirs(out, exist_ok=True)

# --- a grid: land west of 2 E, sea east of it (any ModelGrid works the same way)
lon = np.arange(0.025, 6.0, 0.05)
lat = np.arange(60.025, 64.0, 0.05)
LON, LAT = np.meshgrid(lon, lat)
grid = ModelGrid(LON, LAT, regular_cell_area(lon, lat), LON > 2.0, dims=("j", "i"), name="toy")
print(grid.summary())

# --- sources: what any reader returns (here written by hand)
t = daily_axis(2015, 31)
season = 1 + 0.5 * np.sin(np.arange(31) / 31 * np.pi)
q = np.stack([2000 * season, 300 * season, 50 * season], axis=1)          # m3/s
temp = np.stack([4 + 6 * season, 9 + 0 * season, 3 + 0 * season], axis=1)  # deg C
src = SourceSet(lon=[1.98, 1.92, 0.5], lat=[61.0, 62.5, 63.0], q=q, time=t,
                name=["Big River", "Small River", "Inland Creek"], temp=temp)

# --- map: 10 km footprint, rivers up to 30 km from water, nothing from outside the box
d = Mapper(grid, sigma_km=10.0, max_snap_km=30.0, exclude_outside=True).map_sources(src)
d.report()

# --- write and check
nb = write_mitgcm_runoff(os.path.join(out, "runoff_toy_2015"), d)
write_mitgcm_temperature(os.path.join(out, "runoff_temp_toy_2015"), d)
p = write_sparse_netcdf(os.path.join(out, "toy_runoff_2015.nc"), d, title="demo")
r = read_sparse_netcdf(p)
dense = np.fromfile(os.path.join(out, "runoff_toy_2015"), ">f4").reshape((nb,) + grid.shape)
dtemp = np.fromfile(os.path.join(out, "runoff_temp_toy_2015"), ">f4").reshape((nb,) + grid.shape)
for k in (0, 15, 30):
    tot_dense = float((dense[k].astype(float) * grid.area).sum())
    heat_in = float((q[k, :2] * temp[k, :2]).sum())
    heat_file = float((dense[k].astype(float) * dtemp[k] * grid.area).sum())
    print(f"day {k:2d}: rivers kept {q[k, :2].sum():8.1f} m3/s | dense file {tot_dense:8.1f} | "
          f"sparse file {r.total_m3s(k):8.1f} || heat (m3/s x C) rivers {heat_in:9.1f} | "
          f"files {heat_file:9.1f}")
print(f"sparse file: {d.n_point:,} points instead of {grid.size:,} cells "
      f"({os.path.getsize(p) / 1e3:.0f} kB vs {os.path.getsize(os.path.join(out, 'runoff_toy_2015')) / 1e3:.0f} kB)")
