#!/usr/bin/env python3
"""Wave forcing for a run window: ERA5 (fetch_era5_waves.py), or WAVEWATCH III
offshore (fetch_ww3_wcoast.py) blended with CDIP MOP nearshore (fetch_mop.py).

--source era5 (default): as below.
--source mop_ww3: NOAA GFS-Wave WAVEWATCH III wcoast 1/6-deg fields (nearest
   ocean value, as for ERA5) offshore; CDIP MOP v1.1 10-m-depth points
   (O'Reilly et al. 2016, Coast. Eng. 116:118) in the nearshore: weight
   w = clip((MOP_D1 - depth)/(MOP_D1 - MOP_D0), 0, 1) (1 shallower than 20 m,
   0 deeper than 40 m; the 20/40-m blend depths are ASSUMED), from the nearest
   MOP point within MOP_RMAX km; direction blended as unit vectors. MOP is the
   10-m-depth sea state; it is applied at the local depth with linear theory.
   EXF WvHeight/WvPeriod are then written on a regular 0.005 x 0.004-deg grid
   covering the nest (so the kelp wave mortality sees the MOP nearshore field),
   prefix 'mopww3'.
--source waverys (2026-10-26): Copernicus Marine WAVERYS (GLOBAL_MULTIYEAR_WAV_001_032, dataset
   cmems_mod_glo_wav_my_0.2deg_PT3H-i, 0.2 deg, 3-hourly; VHM0, VTPK, VMDR) from
   ../work/waves/waverys/waverys_*.nc (copernicusmarine subset), as Raphael Savelli uses in PointDume_20m.
   Each nest cell takes the nearest WAVERYS ocean value (no nearshore transformation; WAVERYS has few wet
   points along the Malibu coast); the 3-hourly values are held for each hour (searchsorted, as ERA5).
   EXF files on the same regular 0.005 x 0.004-deg nest grid as mop_ww3, prefix 'waverys', so data.exf
   only needs the file prefix changed.

1. pkg/sediment sedWaveMode = 3 (nest grid, hourly records from the run start):
   ../work/<case>/waves/sed_{Hs,Tp,Dir}_<nx>x<ny>  [m, s, deg toward, ccw from +x]
   ERA5 has no values in coastal/land cells of its 0.5-deg grid: each nest wet
   cell takes the nearest ERA5 ocean value (no nearshore transformation;
   PLACEHOLDER for a nearshore product such as CDIP MOP).
2. EXF WvHeight / WvPeriod for the kelp wave-mortality term (darwin_macroalgae,
   v05_kelp), as yearly files on the ERA5 grid, hourly from Jan 1:
   ../work/<case>/waves/era5_WvHeight_<YYYY>, era5_WvPeriod_<YYYY>
   Hours outside the downloaded window are filled with the nearest available
   record (only the run window is used).
Writes waves.json with the ERA5 grid (lon0, dlon, lat0, dlat, nx, ny) for data.exf.

Usage: python3 make_waves.py <run start YYYY-MM-DD> <days> <out dir> [--source era5|mop_ww3|waverys]
"""
import calendar
import datetime as dt
import glob
import json
import os
import sys

import numpy as np
import xarray as xr
from scipy import ndimage

HERE = os.path.dirname(os.path.abspath(__file__))
WORK = os.path.join(HERE, '..', 'work')
MOP_D0, MOP_D1, MOP_RMAX = 20.0, 40.0, 3.0


def load(source='era5'):
    if source == 'era5':
        ds = xr.open_mfdataset(sorted(glob.glob(os.path.join(WORK, 'waves', 'era5_waves_*.nc'))),
                               combine='by_coords')
        names = ('swh', 'pp1d', 'mwd')
    elif source == 'waverys':
        ds = xr.open_mfdataset(sorted(glob.glob(os.path.join(WORK, 'waves', 'waverys', 'waverys_*.nc'))),
                               combine='by_coords')
        names = ('VHM0', 'VTPK', 'VMDR')
    else:
        ds = xr.open_mfdataset(sorted(glob.glob(os.path.join(WORK, 'waves', 'ww3', 'ww3_wcoast_*.nc'))),
                               combine='by_coords')
        names = ('swh', 'perpw', 'dirpw')
    tname = 'valid_time' if 'valid_time' in ds.dims else 'time'
    ds = ds.sortby(tname).sortby('latitude')
    t = ds[tname].values.astype('datetime64[s]')
    lat, lon = ds['latitude'].values, ds['longitude'].values
    lon = np.where(lon > 180, lon - 360, lon)
    return (t, lat, lon) + tuple(ds[n].values.astype(float) for n in names)


def load_mop():
    ds = xr.open_mfdataset(sorted(glob.glob(os.path.join(WORK, 'waves', 'mop', 'mop_*.nc'))),
                           combine='nested', concat_dim='waveTime').sortby('waveTime')
    ds = ds.isel(waveTime=np.unique(ds.waveTime.values, return_index=True)[1])
    return (ds.waveTime.values.astype('datetime64[s]'), ds.metaLatitude.values, ds.metaLongitude.values,
            ds.waveHs.values.astype(float), ds.waveTp.values.astype(float), ds.waveDp.values.astype(float))


def toward_ccw_east(frm):
    """'coming from, cw from north' (deg) -> unit vector of the direction of travel."""
    th = np.deg2rad(frm)
    return -np.sin(th), -np.cos(th)


def fill_nearest(a):
    """Fill NaNs (land/coast) with the nearest valid value, per record."""
    bad = np.isnan(a[0])
    _, (jj, ii) = ndimage.distance_transform_edt(bad, return_indices=True)
    return a[:, jj, ii]


def main():
    source = 'era5'
    if '--source' in sys.argv:
        k = sys.argv.index('--source')
        source = sys.argv[k + 1]
        del sys.argv[k:k + 2]
    start = np.datetime64(sys.argv[1], 's')
    days = float(sys.argv[2])
    out = sys.argv[3]
    os.makedirs(out, exist_ok=True)
    t, lat, lon, hs, tp, md = load(source)
    hs, tp = fill_nearest(hs), fill_nearest(tp)
    # direction: 'coming from, cw from north' -> 'toward, ccw from east'
    cx, cy = toward_ccw_east(md)
    cx, cy = fill_nearest(cx), fill_nearest(cy)
    dtoward = np.rad2deg(np.arctan2(cy, cx)) % 360.0

    # 1. nest grid, hourly from the run start
    g = json.load(open(os.path.join(WORK, 'grid_info.json')))
    xc, yc = np.array(g['xc']), np.array(g['yc'])
    nrec = int(days * 24) + 25
    tt = start + np.arange(nrec) * np.timedelta64(3600, 's')
    it = np.clip(np.searchsorted(t, tt), 0, len(t) - 1)
    jj = np.clip(np.searchsorted(lat, yc) - 1, 0, lat.size - 1)
    ii = np.abs(lon[None, :] - xc[:, None]).argmin(1)
    jn = np.abs(lat[None, :] - yc[:, None]).argmin(1)
    nest = {name: f[it][:, jn][:, :, ii] for name, f in [('Hs', hs), ('Tp', tp), ('cx', cx), ('cy', cy)]}
    info = {}
    if source == 'mop_ww3':
        depth = -np.fromfile(os.path.join(WORK, g['bathy']), '>f4').reshape(len(yc), len(xc))
        tm, plat, plon, mhs, mtp, mdp = load_mop()
        X, Y = np.meshgrid(xc, yc)
        kx = 111.32 * np.cos(np.deg2rad(34.0))
        dist = np.hypot((X[..., None] - plon) * kx, (Y[..., None] - plat) * 111.32)
        near = dist.argmin(-1)
        w = np.clip((MOP_D1 - depth) / (MOP_D1 - MOP_D0), 0, 1) * (dist.min(-1) < MOP_RMAX) * (depth > 0)
        im = np.clip(np.searchsorted(tm, tt), 0, len(tm) - 1)
        mx, my = toward_ccw_east(mdp)
        for name, m in [('Hs', mhs), ('Tp', mtp), ('cx', mx), ('cy', my)]:
            mv = m[near][:, :, im].transpose(2, 0, 1)          # (nrec, ny, nx) from the nearest MOP point
            ww = np.where(np.isfinite(mv), w, 0.0)
            nest[name] = ww * np.nan_to_num(mv) + (1 - ww) * nest[name]
        info = dict(mop_points=int(plat.size), mop_cells=int((w > 0).sum()), mop_blend_depths=[MOP_D0, MOP_D1],
                    mop_rmax_km=MOP_RMAX, t_mop=[str(tm[0]), str(tm[-1])])
    nest['Dir'] = np.rad2deg(np.arctan2(nest['cy'], nest['cx'])) % 360.0
    for name in ('Hs', 'Tp', 'Dir'):
        nest[name].astype('>f4').tofile(os.path.join(out, f"sed_{name}_{g['nx']}x{g['ny']}"))
    # 2. EXF yearly hourly files: ERA5 grid (era5) or a regular grid over the nest (mop_ww3)
    years = sorted({int(str(x)[:4]) for x in tt})
    if source == 'era5':
        prefix, glon, glat = 'era5', lon, lat
        fields, tsrc = {'WvHeight': hs, 'WvPeriod': tp}, t
    else:
        prefix = 'mopww3' if source == 'mop_ww3' else source
        glon = np.round(xc[0] + 0.005 * np.arange(len(xc)), 6)
        glat = np.round(np.array(g['yg'])[0] + 0.004 * np.arange(int(np.ceil((g['yg'][-1] - g['yg'][0]) / 0.004)) + 1), 6)
        gj = np.abs(yc[None, :] - glat[:, None]).argmin(1)
        fields = {'WvHeight': nest['Hs'][:, gj], 'WvPeriod': nest['Tp'][:, gj]}
        tsrc = tt
    for y in years:
        n = (366 if calendar.isleap(y) else 365) * 24
        ty = np.datetime64(f'{y}-01-01T00:00:00') + np.arange(n) * np.timedelta64(3600, 's')
        iy = np.clip(np.searchsorted(tsrc, ty), 0, len(tsrc) - 1)
        for f, v in fields.items():
            v[iy].astype('>f4').tofile(os.path.join(out, f'{prefix}_{f}_{y}'))
    json.dump(dict(source=source, exf_prefix=prefix, lon0=float(glon[0]), dlon=float(np.diff(glon).mean()),
                   lat=glat.tolist(), nx=int(glon.size), ny=int(glat.size), nrec_sed=nrec, sedWaveDt=3600.0,
                   start=str(start), t_grid=[str(t[0]), str(t[-1])], **info),
              open(os.path.join(out, 'waves.json'), 'w'), indent=1)
    print(source, f'nest Hs {nest["Hs"].min():.2f}-{nest["Hs"].max():.2f} m, median {np.median(nest["Hs"]):.2f}; grid Hs {np.nanmin(hs):.2f}-{np.nanmax(hs):.2f} m, Tp {np.nanmin(tp):.1f}-{np.nanmax(tp):.1f} s; '
          f'{nrec} hourly sediment records from {start}; EXF years {years}; source grid {lon.size}x{lat.size}')


if __name__ == '__main__':
    main()
