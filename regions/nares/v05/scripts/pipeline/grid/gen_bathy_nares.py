#!/usr/bin/env python3
"""gen_bathy_nares.py -- model bathymetry for the Nares Strait downscale.

Same policy as the WestAO no-WAD bathymetry (west_AO/grid/gen_bathy_westAO.py
--nowad, itself following pan_AO): GEBCO 2026 cell-mean, wet where the mean
elevation < 0, 10 m minimum depth (deepen-only; pan_AO's CALC_R_STAR crash),
largest 4-connected wet component only, depth clipped to the 90-level column.

Binning is exact: the grid is regular in rotated coordinates, so every GEBCO
pixel is mapped to (lam', phi') with gen_nares_mitgrid.to_rot and floored into
its cell -- no interpolation, no distance metric.

Ice-shelf cavities (2026-10-06, user: Petermann as a pkg/shelfice cavity): with a
BedMachine Greenland file as 3rd argument, BedMachine pixels flagged floating ice
(mask = 3) are binned into the same cells; a cell >= FLOAT_FRAC floating becomes a
shelfice cell with bed = mean BedMachine bed and ice draft = mean (surface - thickness)
over its floating pixels, kept wet only if the cavity column is >= MIN_CAVITY. The
draft is written to Nares_shelficeTopo.bin (SHELFICEtopoFile, m, negative, 0 elsewhere).
BedMachine heights are relative to the EIGEN-6C4 geoid (~MSL), as GEBCO's are.

Usage: python3 gen_bathy_nares.py <grid dir> <gebco_crop.npz> [BedMachineGreenland.nc]
"""
import os
import sys
import json
import numpy as np
from scipy import ndimage

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import gen_nares_mitgrid as G

MIN_DEPTH = 10.0
MIN_OB_RUN = 30     # side-edge wet runs shorter than this are fjord mouths: wall them off
ZMAX = 7357.0     # sum of the 90-level delR used for WestAO/pan_AO
FLOAT_FRAC = 0.5  # fraction of a cell's BedMachine pixels that must be floating ice
MIN_CAVITY = 20.0 # m: thinner cavities are grounded (land)
# cavities kept as shelfice (user, 2026-10-06: Petermann only; Ryder and Steensby tongues -> land)
KEEP_CAVITIES = {'Petermann': (80.83, -60.97)}
KEEP_RADIUS_KM = 60.0


def bedmachine_float(path, nx, ny, lamG, phiG, d, e1, e2, n):
    """Per-cell floating-ice fraction, mean bed and mean draft from BedMachine."""
    import netCDF4 as nc
    from pyproj import Transformer
    tr_inv = Transformer.from_crs(4326, 3413, always_xy=True)
    tr = Transformer.from_crs(3413, 4326, always_xy=True)
    ds = nc.Dataset(path)
    bx, by = ds['x'][:], ds['y'][:]
    # bounding box of the model grid in EPSG:3413 (sample the grid outline densely)
    L = np.r_[lamG, np.full(phiG.size, lamG[-1]), lamG[::-1], np.full(phiG.size, lamG[0])]
    P = np.r_[np.full(lamG.size, phiG[0]), phiG, np.full(lamG.size, phiG[-1]), phiG[::-1]]
    lo, la = G.to_geo(L, P, e1, e2, n)
    X, Y = tr_inv.transform(lo, la)
    i0, i1 = np.searchsorted(bx, X.min()), np.searchsorted(bx, X.max())
    jy = np.where((by >= Y.min()) & (by <= Y.max()))[0]
    nf = np.zeros(nx * ny); ntot = np.zeros(nx * ny); sbed = np.zeros(nx * ny); sdr = np.zeros(nx * ny)
    xs = bx[i0:i1]
    for j0 in range(jy.min(), jy.max() + 1, 500):
        j1 = min(j0 + 500, jy.max() + 1)
        m = ds['mask'][j0:j1, i0:i1]
        XX, YY = np.meshgrid(xs, by[j0:j1])
        lon, lat = tr.transform(XX, YY)
        lam, phi = G.to_rot(lon, lat, e1, e2, n)
        ci = np.floor((lam - lamG[0]) / d).astype(int); ri = np.floor((phi - phiG[0]) / d).astype(int)
        ok = (ci >= 0) & (ci < nx) & (ri >= 0) & (ri < ny)
        idx = ri * nx + ci
        ntot += np.bincount(idx[ok], minlength=nx * ny)
        fl = ok & (np.ma.getdata(m) == 3)
        if fl.any():
            bed = np.ma.getdata(ds['bed'][j0:j1, i0:i1]).astype(float)
            dr = (np.ma.getdata(ds['surface'][j0:j1, i0:i1]) - np.ma.getdata(ds['thickness'][j0:j1, i0:i1])).astype(float)
            nf += np.bincount(idx[fl], minlength=nx * ny)
            sbed += np.bincount(idx[fl], bed[fl], nx * ny); sdr += np.bincount(idx[fl], dr[fl], nx * ny)
    with np.errstate(invalid='ignore', divide='ignore'):
        frac = np.where(ntot > 0, nf / ntot, 0.0).reshape(ny, nx)
        mbed = (sbed / nf).reshape(ny, nx); mdr = (sdr / nf).reshape(ny, nx)
    return frac, mbed, mdr


def main(gdir, crop, bedmachine=None):
    meta = json.load(open(os.path.join(gdir, 'Nares_grid_meta.json')))
    nx, ny = meta['nx'], meta['ny']
    g = np.fromfile(os.path.join(gdir, 'Nares.mitgrid'), '>f8').reshape(16, ny + 1, nx + 1)
    e1, e2, n = (np.array(meta[k]) for k in ('e1', 'e2', 'pole'))
    lamG, _ = G.to_rot(g[5][0, :], g[6][0, :], e1, e2, n)
    _, phiG = G.to_rot(g[5][:, 0], g[6][:, 0], e1, e2, n)
    d = np.mean(np.diff(lamG))

    D = np.load(crop)
    elev, glon, glat = D['elev'], D['lon'], D['lat']
    ssum = np.zeros(nx * ny); cnt = np.zeros(nx * ny); wsum = np.zeros(nx * ny)
    for j0 in range(0, glat.size, 400):                   # row blocks keep memory small
        LON, LAT = np.meshgrid(glon, glat[j0:j0 + 400])
        lam, phi = G.to_rot(LON, LAT, e1, e2, n)
        ci = np.floor((lam - lamG[0]) / d).astype(int)
        ri = np.floor((phi - phiG[0]) / d).astype(int)
        ok = (ci >= 0) & (ci < nx) & (ri >= 0) & (ri < ny)
        idx = (ri * nx + ci)[ok]
        e = elev[j0:j0 + 400][ok].astype(float)
        ssum += np.bincount(idx, e, nx * ny); cnt += np.bincount(idx, minlength=nx * ny)
        wsum += np.bincount(idx, e < 0, nx * ny)
    cnt = cnt.reshape(ny, nx)
    assert cnt.min() > 0, 'cell with no GEBCO pixels -- crop too small'
    mean_elev = ssum.reshape(ny, nx) / cnt

    bathy = np.where(mean_elev < 0, np.maximum(mean_elev, -ZMAX), 0.0)
    icetopo = np.zeros((ny, nx))
    dropped_cav = 0
    cav = np.zeros((ny, nx), bool)
    if bedmachine:
        frac, mbed, mdr = bedmachine_float(bedmachine, nx, ny, lamG, phiG, d, e1, e2, n)
        flt = frac >= FLOAT_FRAC
        col = mdr - mbed
        cav = flt & (col >= MIN_CAVITY) & (mdr < 0)
        g8 = np.fromfile(os.path.join(gdir, 'Nares.mitgrid'), '>f8').reshape(16, ny + 1, nx + 1)
        xc, yc = g8[0][:-1, :-1], g8[1][:-1, :-1]
        clab, ncl = ndimage.label(cav, structure=np.ones((3, 3)))
        kept = np.zeros_like(cav)
        for kk in range(1, ncl + 1):
            sel = clab == kk
            clat, clon = yc[sel].mean(), xc[sel].mean()
            for name, (plat, plon) in KEEP_CAVITIES.items():
                dkm = np.hypot((clat - plat) * 111.2, (clon - plon) * 111.2 * np.cos(np.deg2rad(plat)))
                if dkm < KEEP_RADIUS_KM:
                    kept |= sel
        dropped_cav = int((cav & ~kept).sum())
        cav = kept
        bathy[flt] = 0.0                       # floating but not kept (or too thin): land
        bathy[cav] = np.maximum(mbed[cav], -ZMAX)
        icetopo[cav] = mdr[cav]
    lab, nlab = ndimage.label(bathy < 0)
    sizes = np.bincount(lab.ravel())[1:]
    main_lab = 1 + int(np.argmax(sizes))
    removed = int(((lab != main_lab) & (bathy < 0)).sum())
    bathy[lab != main_lab] = 0.0
    # Side edges (j = 0, ny-1) cross the open sea only at the Baffin Bay corner and
    # the Lincoln Sea; short wet runs there are fjord mouths (N Greenland, one on
    # Ellesmere). An open boundary across a fjord mouth has no parent data worth
    # having, so wall the edge cell and keep the fjord as a closed inlet.
    walled = 0
    for row in (0, ny - 1):
        e = bathy[row] < 0
        i = 0
        while i < nx:
            if e[i]:
                k = i
                while k < nx and e[k]:
                    k += 1
                if k - i < MIN_OB_RUN:
                    bathy[row, i:k] = 0.0
                    walled += k - i
                i = k
            else:
                i += 1
    lab, nlab2 = ndimage.label(bathy < 0)
    sizes = np.bincount(lab.ravel())[1:]
    keep = 1 + int(np.argmax(sizes))
    removed += int(((lab != keep) & (bathy < 0)).sum())
    bathy[lab != keep] = 0.0
    icetopo[bathy >= 0] = 0.0
    cav = cav & (bathy < 0)
    shallow = (bathy < 0) & (bathy > -MIN_DEPTH) & ~cav
    bathy[shallow] = -MIN_DEPTH
    wet = bathy < 0
    stats = dict(nx=nx, ny=ny, gebco_px_per_cell=[int(cnt.min()), int(cnt.max())],
                 wet=int(wet.sum()), wet_fraction=float(wet.mean()), components=int(nlab),
                 removed_cells=removed, walled_fjord_mouth_cells=walled, deepened_to_min_depth=int(shallow.sum()),
                 max_depth_m=float(-bathy.min()),
                 shelfice_cells=int(cav.sum()), floating_cells_made_land=dropped_cav if bedmachine else 0,
                 shelfice_draft_m=[float(-icetopo[cav].max()), float(-icetopo[cav].min())] if cav.any() else None,
                 shelfice_column_m=[float((icetopo-bathy)[cav].min()), float((icetopo-bathy)[cav].max())] if cav.any() else None,
                 ob_wet_cells=dict(south_x0=int(wet[:, 0].sum()), north_xmax=int(wet[:, -1].sum()),
                                   side_y0=int(wet[0, :].sum()), side_ymax=int(wet[-1, :].sum())))
    bathy.astype('>f8').tofile(os.path.join(gdir, 'Nares_bathymetry.bin'))
    icetopo.astype('>f8').tofile(os.path.join(gdir, 'Nares_shelficeTopo.bin'))
    json.dump(stats, open(os.path.join(gdir, 'Nares_bathymetry_stats.json'), 'w'), indent=1)
    print(json.dumps(stats, indent=1))


if __name__ == '__main__':
    main(*sys.argv[1:4])
