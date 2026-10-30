#!/usr/bin/env python3
"""Crop global FES tide files (FES2022b elevations, FES2014 elevations or
currents) to the LR17 region and delete the global file.

Usage:  python3 fes_crop.py FILE [FILE ...]
Writes FILE with '_lr17' before '.nc' (lat 45.4 .. 46.8, lon 358.0 .. 359.5,
i.e. 2.0 .. 0.5 W), keeping all variables, then removes FILE.
"""
import os
import sys

import xarray as xr

BOX = dict(lat=slice(45.4, 46.8), lon=slice(358.0, 359.5))

for f in sys.argv[1:]:
    out = f.replace('.nc', '_lr17.nc')
    d = xr.open_dataset(f)
    sub = d.sel(**{k: v for k, v in BOX.items() if k in d.dims})
    sub.load().to_netcdf(out)
    d.close()
    os.remove(f)
    print(f'{os.path.basename(f)} -> {os.path.basename(out)} {dict(sub.sizes)}')
