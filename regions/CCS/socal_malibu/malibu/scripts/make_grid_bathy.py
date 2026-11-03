#!/usr/bin/env python3
"""Grid, vertical levels and bathymetry for the Malibu / Pt. Dume kelp nest.

Horizontal: spherical-polar, dlon = 0.005 deg (~460 m at 34N), delY scaled by
cos(lat) so cells are ~square (same delYfile approach as CCS v05).
Vertical: 20 x 2 m levels over the kelp canopy depth range (0-40 m), then
geometric stretching to the deepest point of Santa Monica Basin in the domain.
Bathymetry: NOAA Coastal Relief Model vol. 6 (3 arc-sec), via NCEI OPeNDAP,
cell-averaged onto the model grid.

Writes (big-endian float32, MITgcm convention) into ../work/:
  BATHY_<nx>x<ny>_malibu   (negative = ocean depth)
  delYFile_<ny>_malibu
  grid_info.json           (read by the other generators)
and prints the &PARM04 / OBCS snippets for input/data and input/data.obcs.

Usage: python3 -I make_grid_bathy.py [--crm-cache ../work/crm6_subset.npz]
"""
import argparse
import json
import os

import numpy as np
from scipy import ndimage

HERE = os.path.dirname(os.path.abspath(__file__))
WORK = os.path.join(HERE, '..', 'work')

# ---- domain --------------------------------------------------------------
LON0, DLON, NX = -119.10, 0.005, 128      # -119.10 .. -118.46
LAT0, NY = 33.80, 72                      # 33.80 .. ~34.10
NR, NTOP, DZTOP = 50, 20, 2.0             # 20 x 2 m, then stretched
MINDEPTH = 4.0                            # shallowest wet column [m]
CRM_URL = 'https://www.ngdc.noaa.gov/thredds/dodsC/crm/crm_vol6.nc'


def lat_edges():
    """Isotropic-ish rows: dy_j = dlon*cos(lat_j), marched from LAT0."""
    yg = [LAT0]
    for _ in range(NY):
        y = yg[-1]
        dy = DLON * np.cos(np.deg2rad(y + 0.5 * DLON * np.cos(np.deg2rad(y))))
        yg.append(y + dy)
    return np.array(yg)


def fetch_crm(lon_rng, lat_rng, cache):
    if cache and os.path.exists(cache):
        d = np.load(cache)
        return d['x'], d['y'], d['z']
    import netCDF4
    nc = netCDF4.Dataset(CRM_URL)
    x = nc['x'][:]
    y = nc['y'][:]
    i = np.where((x >= lon_rng[0]) & (x <= lon_rng[1]))[0]
    j = np.where((y >= lat_rng[0]) & (y <= lat_rng[1]))[0]
    z = np.asarray(nc['z'][j[0]:j[-1] + 1, i[0]:i[-1] + 1], dtype='f4')
    x, y = np.asarray(x[i]), np.asarray(y[j])
    nc.close()
    if cache:
        np.savez_compressed(cache, x=x, y=y, z=z)
    return x, y, z


def vertical_levels(maxdepth):
    """NTOP levels of DZTOP, then geometric growth reaching maxdepth."""
    nlow = NR - NTOP
    rest = maxdepth - NTOP * DZTOP
    lo, hi = 1.0001, 2.0
    for _ in range(200):            # bisection on growth factor r
        r = 0.5 * (lo + hi)
        tot = DZTOP * sum(r ** k for k in range(1, nlow + 1))
        lo, hi = (r, hi) if tot < rest else (lo, r)
    dz = [DZTOP] * NTOP + [DZTOP * r ** k for k in range(1, nlow + 1)]
    return np.round(np.array(dz), 2)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--crm-cache', default=os.path.join(WORK, 'crm6_subset.npz'))
    a = ap.parse_args()
    os.makedirs(WORK, exist_ok=True)

    xg = LON0 + DLON * np.arange(NX + 1)
    yg = lat_edges()
    xc, yc = 0.5 * (xg[1:] + xg[:-1]), 0.5 * (yg[1:] + yg[:-1])
    dely = np.diff(yg)

    x, y, z = fetch_crm((xg[0] - .01, xg[-1] + .01), (yg[0] - .01, yg[-1] + .01),
                        a.crm_cache)
    # cell average of depth (land counted as 0) and wet fraction
    ix = np.searchsorted(xg, x) - 1
    jy = np.searchsorted(yg, y) - 1
    ok_x, ok_y = (ix >= 0) & (ix < NX), (jy >= 0) & (jy < NY)
    zz = z[np.ix_(ok_y, ok_x)]
    II, JJ = np.meshgrid(ix[ok_x], jy[ok_y])
    wet = np.isfinite(zz) & (zz < 0)
    dep = np.where(wet, -zz, 0.0)
    flat = JJ.ravel() * NX + II.ravel()
    n = np.bincount(flat, minlength=NX * NY).reshape(NY, NX)
    sw = np.bincount(flat, wet.ravel(), NX * NY).reshape(NY, NX)
    sd = np.bincount(flat, dep.ravel(), NX * NY).reshape(NY, NX)
    frac = sw / np.maximum(n, 1)
    depth = np.where(frac >= 0.5, sd / np.maximum(sw, 1), 0.0)
    depth = np.where(depth > 0, np.maximum(depth, MINDEPTH), 0.0)

    # close the north edge (only a few wet cells off Pt. Mugu reach it);
    # open boundaries are then W, S and E
    depth[-1, :] = 0.0
    # keep only ocean connected to the open boundary (drop lagoons / ponds)
    lab, _ = ndimage.label(depth > 0)
    edge = np.unique(np.r_[lab[0], lab[-1], lab[:, 0], lab[:, -1]])
    depth[~np.isin(lab, edge[edge > 0])] = 0.0
    # fill 1-cell-wide inlets (3 land neighbours) to keep the solver happy
    for _ in range(5):
        w = depth > 0
        nb = (np.roll(w, 1, 0).astype(int) + np.roll(w, -1, 0)
              + np.roll(w, 1, 1) + np.roll(w, -1, 1))
        nb[0], nb[-1], nb[:, 0], nb[:, -1] = 4, 4, 4, 4
        kill = w & (nb <= 1)
        if not kill.any():
            break
        depth[kill] = 0.0

    maxdepth = float(np.ceil(depth.max() / 10.0) * 10.0)
    dz = vertical_levels(maxdepth)

    bathy = -depth  # (NY, NX) -> MITgcm reads x fastest
    fb = os.path.join(WORK, f'BATHY_{NX}x{NY}_malibu')
    fy = os.path.join(WORK, f'delYFile_{NY}_malibu')
    bathy.astype('>f4').tofile(fb)
    dely.astype('>f4').tofile(fy)

    # open boundaries: any wet cell on a domain edge
    ob = {'north': bool((depth[-1] > 0).any()), 'south': bool((depth[0] > 0).any()),
          'west': bool((depth[:, 0] > 0).any()), 'east': bool((depth[:, -1] > 0).any())}
    info = dict(nx=NX, ny=NY, nr=NR, lon0=LON0, dlon=DLON, lat0=LAT0,
                xc=xc.tolist(), yc=yc.tolist(), xg=xg.tolist(), yg=yg.tolist(),
                dely=dely.tolist(), delr=dz.tolist(), maxdepth=float(depth.max()),
                bathy=os.path.basename(fb), delyfile=os.path.basename(fy), open=ob)
    with open(os.path.join(WORK, 'grid_info.json'), 'w') as f:
        json.dump(info, f)

    km = 111.2 * np.cos(np.deg2rad(34.0)) * DLON
    print(f'grid {NX}x{NY}x{NR}, dx ~{km*1e3:.0f} m, lat {yg[0]:.4f}..{yg[-1]:.4f}')
    print(f'wet columns {int((depth>0).sum())}/{NX*NY}, max depth {depth.max():.0f} m,'
          f' sum(delR) {dz.sum():.1f} m, deepest dz {dz[-1]:.1f} m')
    print('open boundaries:', ob)
    print('\n # &PARM04 snippet')
    print(f' xgOrigin = {LON0},\n ygOrigin = {LAT0},\n delX = {NX}*{DLON},')
    print(f" delYfile = '{os.path.basename(fy)}',\n delR =")
    for k in range(0, NR, 8):
        print('   ' + ', '.join(f'{v:.2f}' for v in dz[k:k + 8]) + ',')


if __name__ == '__main__':
    main()
