#!/usr/bin/env python3
"""
Add the fields the downscaling STEP3 tools need to work/socal/socal_ncgrid.nc.

gen_obcs.py / gen_pickups.py / gen_obcs_transport_correction.py (ecco_darwin
regions/downscaling) read XC, YC, AngleCS, AngleSN, HFacC, HFacS, HFacW, drF,
dxG, dyG from <region>_ncgrid.nc. make_socal_grid.py wrote only XC, YC, XG, YG,
Depth, so this script computes the rest for the 1/48 deg spherical-polar grid:

- AngleCS = 1, AngleSN = 0 (grid lines follow lon/lat).
- hFacC follows MITgcm ini_masks_etc.F (no shelfice): hFacMin = 0.3 and
  hFacMinDr = 1 m, the Malibu nest values (input_pace2024_*/data), so the
  regional `data` must use the same pair.
- HFacW (Nr, ny, nx+1) = min of the two neighboring cells; HFacS (Nr, ny+1, nx)
  likewise. Outer edge faces take the adjacent cell's hFacC (no periodic wrap).
- dxG (ny+1, nx) and dyG (ny, nx+1) on a sphere of radius 6370 km (MITgcm
  rSphere default).

Writes work/socal/socal_ncgrid.nc in place (keeps the original variables).
"""
import json
import os

import numpy as np
import xarray as xr

HERE = os.path.dirname(os.path.abspath(__file__))
D = os.path.join(HERE, '..', 'work', 'socal')
HFACMIN, HFACMINDR, RSPHERE = 0.3, 1.0, 6370.0e3

g = json.load(open(os.path.join(D, 'grid_info.json')))
ds = xr.load_dataset(os.path.join(D, 'socal_ncgrid.nc'))
drF = np.loadtxt(os.path.join(D, 'delR.txt'))
depth = ds['Depth'].values                      # positive down, 0 on land
ny, nx, nr = depth.shape[0], depth.shape[1], drF.size
assert (ny, nx, nr) == (g['ny'], g['nx'], g['nr'])

rF = -np.concatenate([[0.0], np.cumsum(drF)])   # interface depths, negative down
hFacC = np.zeros((nr, ny, nx))
for k in range(nr):
    mnsz = max(HFACMIN, min(HFACMINDR / drF[k], 1.0))
    h = np.clip((rF[k] + depth) / drF[k], 0.0, 1.0)
    small = h < mnsz
    h[small] = np.where(h[small] < 0.5 * mnsz, 0.0, mnsz)
    hFacC[k] = h

hFacW = np.empty((nr, ny, nx + 1))
hFacW[:, :, 1:-1] = np.minimum(hFacC[:, :, :-1], hFacC[:, :, 1:])
hFacW[:, :, 0], hFacW[:, :, -1] = hFacC[:, :, 0], hFacC[:, :, -1]
hFacS = np.empty((nr, ny + 1, nx))
hFacS[:, 1:-1, :] = np.minimum(hFacC[:, :-1, :], hFacC[:, 1:, :])
hFacS[:, 0, :], hFacS[:, -1, :] = hFacC[:, 0, :], hFacC[:, -1, :]

yg = g['yg0'] + g['dlat'] * np.arange(ny + 1)
yc = np.asarray(g['yc'])
dxG = np.repeat((RSPHERE * np.cos(np.deg2rad(yg)) * np.deg2rad(g['dlon']))[:, None], nx, 1)
dyG = np.full((ny, nx + 1), RSPHERE * np.deg2rad(g['dlat']))

ds['AngleCS'] = (('j', 'i'), np.ones((ny, nx)))
ds['AngleSN'] = (('j', 'i'), np.zeros((ny, nx)))
ds['HFacC'] = (('k', 'j', 'i'), hFacC)
ds['HFacW'] = (('k', 'j', 'ig'), hFacW)
ds['HFacS'] = (('k', 'jg', 'i'), hFacS)
ds['drF'] = (('k',), drF)
ds['dxG'] = (('jg', 'i'), dxG)
ds['dyG'] = (('j', 'ig'), dyG)
ds.attrs['hfac'] = 'make_socal_ncgrid_full.py: hFacMin %.1f, hFacMinDr %.1f m' % (HFACMIN, HFACMINDR)
ds.to_netcdf(os.path.join(D, 'socal_ncgrid.nc'))

wet = hFacC[0] > 0
print('wet surface cells', wet.sum(), 'wet 3D cells', (hFacC > 0).sum())
print('boundary wet points W E S N:', wet[:, 0].sum(), wet[:, -1].sum(), wet[0, :].sum(), wet[-1, :].sum())
print('model depth vs Depth, max |diff| m:', np.abs((hFacC * drF[:, None, None]).sum(0) - depth).max())
