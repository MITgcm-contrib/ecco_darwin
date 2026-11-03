#!/usr/bin/env python3
"""Intermediate ~2 km Southern California Bight nest ("socal") between the
ECCO-Darwin v05 LLC270 V5r1 parent (~30 km) and the Malibu nest (~460 m).

Horizontal: uniform 1/48 deg in lon and lat (dx ~1.93 km, dy ~2.32 km at
33.5 N), cell edges lon -121.0 .. -117.5 (168 cells), lat 32.5 .. 32.5 + 101/48
(101 cells), so the Malibu nest (lon -119.1 .. -118.46, lat 33.8 .. 34.1) sits
>= 1.4 deg inside the west, east and south edges and 0.5 deg inside the north.
Bathymetry: NOAA Coastal Relief Model vol. 6 (3 arc-sec; the Malibu nest's
source), read every 4th point (12 arc-sec) and block-averaged per cell
(mean elevation of all samples; land where the mean is >= 0 m), cells
shallower than 4 m deepened to 4 m (as Malibu), closed basins not connected
to the open boundary removed (scipy.ndimage.label).
Vertical: 60 levels; twenty 2-m levels to 40 m (identical to the Malibu nest's
top 20 levels), then 40 levels growing geometrically to cover the deepest cell.
Writes ../work/socal/: BATHY_168x101_socal (>f4, negative down), delR.txt,
grid_info.json and socal_ncgrid.nc (XC, YC, XG, YG, Depth; what
regions/downscaling/utils/gen_dvmasks.py reads).

Usage: python3 make_socal_grid.py [--crm-cache ../work/socal/crm6_socal.npz]
"""
import argparse
import json
import os

import netCDF4
import numpy as np
from scipy import ndimage
from scipy.optimize import brentq

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, '..', 'work', 'socal')
CRM_URL = 'https://www.ngdc.noaa.gov/thredds/dodsC/crm/crm_vol6.nc'
D = 1 / 48
LON0, LAT0, NX, NY = -121.0, 32.5, 168, 101
NR, NTOP, DZTOP = 60, 20, 2.0
HMIN = 4.0


def fetch(cache, stride=4):
    if cache and os.path.exists(cache):
        z = np.load(cache)
        return z['x'], z['y'], z['z']
    nc = netCDF4.Dataset(CRM_URL)
    x, y = nc.variables['x'][:], nc.variables['y'][:]
    i = np.where((x >= LON0 - 0.01) & (x <= LON0 + NX * D + 0.01))[0]
    j = np.where((y >= LAT0 - 0.01) & (y <= LAT0 + NY * D + 0.01))[0]
    i, j = i[::stride], j[::stride]
    z = nc.variables['z'][j[0]:j[-1] + 1:stride, i[0]:i[-1] + 1:stride]
    x, y, z = x[i], y[j], np.asarray(z, float)
    if cache:
        np.savez_compressed(cache, x=x, y=y, z=z)
    return x, y, z


def levels(maxdepth):
    top = NTOP * DZTOP
    n = NR - NTOP
    f = lambda r: DZTOP * r * (r ** n - 1) / (r - 1) - (maxdepth - top)
    r = brentq(f, 1.0001, 2.0)
    dz = np.r_[np.full(NTOP, DZTOP), DZTOP * r ** np.arange(1, n + 1)]
    return dz, r


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--crm-cache', default=os.path.join(OUT, 'crm6_socal.npz'))
    a = ap.parse_args()
    os.makedirs(OUT, exist_ok=True)
    xg = LON0 + D * np.arange(NX + 1)
    yg = LAT0 + D * np.arange(NY + 1)
    xc, yc = 0.5 * (xg[1:] + xg[:-1]), 0.5 * (yg[1:] + yg[:-1])
    x, y, z = fetch(a.crm_cache)
    ii = np.clip(np.searchsorted(xg, x) - 1, -1, NX)
    jj = np.clip(np.searchsorted(yg, y) - 1, -1, NY)
    S = np.zeros((NY, NX)); N = np.zeros((NY, NX))
    J, I = np.meshgrid(jj, ii, indexing='ij')
    ok = (I >= 0) & (I < NX) & (J >= 0) & (J < NY)
    np.add.at(S, (J[ok], I[ok]), z[ok]); np.add.at(N, (J[ok], I[ok]), 1)
    elev = S / np.maximum(N, 1)
    depth = np.where(elev < 0, -elev, 0.0)
    depth = np.where((depth > 0) & (depth < HMIN), HMIN, depth)
    # keep only water connected to an open (domain-edge) wet cell
    lab, n = ndimage.label(depth > 0)
    edge = np.unique(np.r_[lab[0], lab[-1], lab[:, 0], lab[:, -1]])
    keep = np.isin(lab, edge[edge > 0])
    removed = int(((depth > 0) & ~keep).sum())
    depth = np.where(keep, depth, 0.0)
    dz, r = levels(depth.max() + 1.0)
    bathy = -depth
    tag = f'{NX}x{NY}'
    bathy.astype('>f4').tofile(os.path.join(OUT, f'BATHY_{tag}_socal'))
    np.savetxt(os.path.join(OUT, 'delR.txt'), dz, fmt='%.3f')
    XC, YC = np.meshgrid(xc, yc)
    with netCDF4.Dataset(os.path.join(OUT, 'socal_ncgrid.nc'), 'w') as nc:
        nc.createDimension('j', NY); nc.createDimension('i', NX)
        nc.createDimension('jg', NY + 1); nc.createDimension('ig', NX + 1)
        for nm, v, dims in (('XC', XC, ('j', 'i')), ('YC', YC, ('j', 'i')), ('Depth', depth, ('j', 'i'))):
            var = nc.createVariable(nm, 'f8', dims); var[:] = v
        XG, YG = np.meshgrid(xg, yg)
        for nm, v in (('XG', XG), ('YG', YG)):
            var = nc.createVariable(nm, 'f8', ('jg', 'ig')); var[:] = v
        nc.source = 'make_socal_grid.py: 1/48 deg, NOAA CRM vol. 6 block mean'
    wet = depth > 0
    ob = dict(west=int(wet[:, 0].sum()), east=int(wet[:, -1].sum()), south=int(wet[0].sum()), north=int(wet[-1].sum()))
    info = dict(nx=NX, ny=NY, nr=NR, dlon=D, dlat=D, xg0=LON0, yg0=LAT0, xc=xc.tolist(), yc=yc.tolist(),
                delr=dz.round(3).tolist(), stretch=float(r), maxdepth=float(depth.max()),
                bathy=f'BATHY_{tag}_socal', wet_cells=int(wet.sum()), closed_basin_cells_removed=removed,
                open_boundary_wet_points=ob, crm_stride_arcsec=12)
    json.dump(info, open(os.path.join(OUT, 'grid_info.json'), 'w'), indent=1)
    print(json.dumps({k: v for k, v in info.items() if k not in ('xc', 'yc', 'delr')}))
    print('delR top/bottom', dz[:3], dz[-3:], 'sum', dz.sum())


if __name__ == '__main__':
    main()
