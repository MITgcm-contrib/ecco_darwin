#!/usr/bin/env python3
"""Correct the socal 2-km nest's upper-ocean initial temperature with MUR SST (2026-10-10).

The parent (V5r1 LLC270) 2025-01-01 state interpolated by gen_pickups.py is too warm near the
surface: vs 38 CalCOFI 2501RL CTD casts (3-16 Jan 2025) +1.2-1.6 C at 0-60 m on the shelf/slope
(< 500 m) and +0.7-1.1 C offshore, but within 0.1 C at 60-200 m; MUR differs from the same casts by
only +0.14-0.19 C at the surface. Correction: dT = MUR(2025-01-01) - THETA(level 1), MUR bilinear to
the cell centers, coastal cells without MUR filled from the nearest cell with MUR; applied with
weight 1 down to 30 m, linearly to 0 at 60 m, wet cells only. Salinity is left as it is.

Usage: socal_ic_murfix.py <pickup_THETA.0000867888.data> <mur nc> <out file>
"""
import os
import sys

import numpy as np
import xarray as xr
from scipy import ndimage
from scipy.interpolate import RegularGridInterpolator

HERE = os.path.dirname(os.path.abspath(__file__))
g = xr.open_dataset(os.path.join(HERE, '..', 'work', 'socal', 'socal_ncgrid.nc'))
xc, yc, hf, dz = g.XC.values, g.YC.values, g.HFacC.values, g.drF.values
zc = np.cumsum(dz) - dz / 2
nr, ny, nx = hf.shape
src, murf, out = sys.argv[1:4]
T = np.fromfile(src, '>f8').reshape(nr, ny, nx)
m = xr.open_dataset(murf).sel(time='2025-01-01')
mur = RegularGridInterpolator((m.latitude.values, m.longitude.values), m.analysed_sst.values[0],
                              bounds_error=False, fill_value=np.nan)(np.c_[yc.ravel(), xc.ravel()]).reshape(ny, nx)
wet = hf[0] > 0
dT = mur - T[0]
ok = wet & np.isfinite(dT)
_, (jj, ii) = ndimage.distance_transform_edt(~ok, return_indices=True)
dT = dT[jj, ii]
w = np.clip((60.0 - zc) / 30.0, 0.0, 1.0)
T2 = np.where(hf > 0, T + w[:, None, None] * dT[None], T)
T2.astype('>f8').tofile(out)
shelf = wet & (g.Depth.values < 100)
print(f'dT surface: mean {dT[wet].mean():+.2f}, shelf {dT[shelf].mean():+.2f}, offshore '
      f'{dT[wet & (g.Depth.values > 500)].mean():+.2f}, p5/p95 {np.percentile(dT[wet], 5):+.2f}/{np.percentile(dT[wet], 95):+.2f} C; '
      f'filled {(wet & ~ok).sum()} coastal cells; weight 1 to 30 m, 0 at 60 m; wrote {out}')
