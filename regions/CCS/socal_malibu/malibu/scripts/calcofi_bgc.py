#!/usr/bin/env python3
"""Observation-based biogeochemistry for the PACE-year ICs and OBCS.

The ECCO-Darwin v05 V5r1 parent has surface Chl ~0.07-0.09 mg/m3 and NO3 ~0
around Santa Monica Bay in Nov 2024 - Jan 2025, while CalCOFI bottles give
surface Chl ~0.6-5 mg/m3 and NO3 median ~2 uM (raw_data/obs_2024-25/README.txt).
This replaces the upper-ocean NO3, PO4, SiO2 and plankton of the V5r1-derived
ICs/OBCS (make_ic_obcs.py --mode v5r1) with CalCOFI profiles; physics, DIC,
ALK, O2, Fe, DOM, POM, NO2, NH4 and CDOM stay as they are.

Data: CalCOFI PRELIMINARY bottle data (non-publication use; "oxygen, nitrate and
chlorophyll may change significantly after post-cruise calibrations"), cruises
2411SR (15-24 Nov 2024; Chl only, nutrients not yet imported) and 2501RL +
2502RL (3-19 Jan 2025; Chl, NO3, PO4, SIL), casts in 33.0-34.3 N, 120.0-117.5 W.
QC: 2501_037u 3-m bottle NO3/PO4/SIL (deep-water values in 14.3 C surface water)
and 2501_039u SIL at 3 and 12 m excluded (as in the README).

Method:
- each cast interpolated (linear) to the nest level centres between its
  shallowest and deepest bottle, held constant above the shallowest bottle;
  shelf (water depth < 500 m and < 20 km from land, station table 'zone') and
  offshore casts give median profiles (>= 3 casts per level, gaps filled from
  the level above);
- nutrients: Jan 2025 profiles for all times; Chl: Nov 2024 profile at the 2411
  median cast time, Jan 2025 at the 2501 median, linear in between, constant
  outside;
- columns: shelf weight w = clip((500 - depth)/400, 0, 1);
- depth: observations over 0-200 m, linear taper to the parent at 300 m;
- plankton: total Chl from CalCOFI split over Chl01-05 with the parent's
  domain-mean class fractions per level, carbon c01-c05 = Chl / (parent Chl:C
  per class and level); grazers c06-c07 kept at the parent's values (default).
  --scale-grazers scales them by the same factor as total phytoplankton carbon
  (keeps the parent's grazer:phyto ratio; ASSUMED). In the 2024-25 runs this
  overloaded the grazers: run_pace2024_obs (parent C:Chl ~268) ended with grazer
  C 6.7 vs phytoplankton C 0.14 mmol/m3 and surface Chl ~0.01 mg/m3.

Usage: python3 calcofi_bgc.py <in dir (work/pace2024)> <out dir> <start YYYY-MM-DD> [C:Chl g/g] [--scale-grazers]
       [--with-2504] [--rec-days D] [--rec-offset-h H]
--with-2504 (2026-10-10): adds CalCOFI 2504SH (spring 2025, 20 Mar - early Apr, Bell M. Shimada;
  calcofi_prelim_bottle_SCB_2504SH.csv). Nutrients: Jan 2025 profiles up to the 2501 median cast time,
  linear to the 2504 profiles at the 2504 median time, constant after; Chl: Nov -> Jan -> Apr piecewise
  linear. 2504 casts have no station-table zone: shelf = water depth < 500 m and < 20 km from land,
  both from the socal 2-km grid (work/socal/socal_ncgrid.nc), the rule of the derived station table.
--rec-days / --rec-offset-h: OBCS record r is at start + offset + r * rec_days (default 30 d, 0 h;
  the daily socal-mode OBCS of make_ic_obcs.py use 1 d and -12 h).
With C:Chl given, phytoplankton carbon = Chl x C:Chl / 12.011 (mmol C/m3) for
every class instead of the parent's Chl:C, but not below each class's darwin
C:Chl floor 12.011/chl2cmax (data.traits CHL2CMAX: 0.5, 0.5, 0.2, 0.2, 0.2 mg Chl
(mmol C)-1 -> floors 24, 24, 60, 60, 60 g/g). Literature: southern CCE inshore
winter mixed-layer C:Chl 18-63 g/g, mean ~38 (Taylor, Landry, Selph & Wokuluk
2015, Deep-Sea Res. II 112:117, Table 1, Line 90 inshore stations); ~40 with
shallow nitraclines (same, Sect. 3.2); inshore 2004-13 baseline 48 +/- 4 (Landry
et al. 2024, L&O 69:1683, Table 1). Recommended 40 (range 20-60).
"""
import json
import os
import shutil
import sys

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
WORK = os.path.join(HERE, '..', 'work')
OBS = os.path.join(HERE, '..', '..', '..', 'raw_data', 'obs_2024-25')
BOX = dict(lat=(33.0, 34.3), lon=(-120.0, -117.5))
Z_OBS, Z_PARENT = 200.0, 300.0
VARS = {'Chl-a': 'chl', 'NO3': 'NO3', 'PO4': 'PO4', 'SIL': 'SiO2'}


def zone_2504(b):
    """shelf/offshore for casts without a station-table zone: depth < 500 m and < 20 km from land."""
    import xarray as xr
    from scipy import ndimage
    g = xr.open_dataset(os.path.join(WORK, 'socal', 'socal_ncgrid.nc'))
    xc, yc, dep = g.XC.values, g.YC.values, g.Depth.values
    dkm = ndimage.distance_transform_edt(dep > 0) * 0.0208333 * 111.2 * np.cos(np.deg2rad(33.5))
    z = {}
    for cid, c in b.groupby('Cast_ID'):
        j = int(np.clip((c.Lat_Dec.iloc[0] - yc[0, 0]) / 0.0208333 + 0.5, 0, yc.shape[0] - 1))
        i = int(np.clip((c.Lon_Dec.iloc[0] - xc[0, 0]) / 0.0208333 + 0.5, 0, xc.shape[1] - 1))
        z[cid] = 'shelf' if (dep[j, i] < 500 and dkm[j, i] < 20) else 'offshore'
    return z


def casts(with_2504=False):
    b = pd.read_csv(os.path.join(OBS, 'calcofi_prelim_bottle_SCB_2411SR_2501RL_2502RL.csv'), low_memory=False)
    st = pd.read_csv(os.path.join(OBS, 'derived_calcofi_station_table_Nov2024-Jan2025.csv'))
    if with_2504:
        b4 = pd.read_csv(os.path.join(OBS, 'calcofi_prelim_bottle_SCB_2504SH.csv'), low_memory=False)
        b = pd.concat([b, b4], ignore_index=True)
    b = b[(b.Lat_Dec.between(*BOX['lat'])) & (b.Lon_Dec.between(*BOX['lon']))].copy()
    for c in VARS:
        b[c] = pd.to_numeric(b[c], errors='coerce')
    b['BTL_Depth'] = pd.to_numeric(b['BTL_Depth'], errors='coerce')
    # QC exclusions (README)
    bad = (b.Cast_ID == '2501_037u') & (b.BTL_Depth < 5)
    b.loc[bad, ['NO3', 'PO4', 'SIL']] = np.nan
    b.loc[(b.Cast_ID == '2501_039u') & (b.BTL_Depth < 15), 'SIL'] = np.nan
    zone = dict(zip(st.cast, st.zone))
    if with_2504:
        zone.update(zone_2504(b[b.Cast_ID.str.startswith('2504')]))
    b['zone'] = b.Cast_ID.map(zone).fillna('offshore').map(lambda z: 'shelf' if 'shelf' in z else 'offshore')
    b['time'] = pd.to_datetime(b.Date_Time_UTC, format='%d-%b-%Y %H:%M:%S')
    b['cruise'] = b.Cast_ID.str[:4]
    return b


def profiles(b, zc, cruises, var):
    out = {}
    for zone in ('shelf', 'offshore'):
        P = []
        for cid, c in b[(b.zone == zone) & (b.cruise.isin(cruises))].groupby('Cast_ID'):
            c = c[np.isfinite(c[var]) & np.isfinite(c.BTL_Depth)].sort_values('BTL_Depth')
            if len(c) < 2:
                continue
            z, v = c.BTL_Depth.values, c[var].values
            p = np.interp(zc, z, v, left=v[0], right=np.nan)
            p[zc > z[-1]] = np.nan
            P.append(p)
        P = np.array(P)
        n = np.isfinite(P).sum(0)
        med = np.where(n >= 3, np.nanmedian(np.where(np.isfinite(P), P, np.nan), 0), np.nan)
        for k in range(1, med.size):                     # gaps from the level above
            if not np.isfinite(med[k]):
                med[k] = med[k - 1]
        out[zone] = (med, len(P))
    return out


def main():
    scale_grazers = '--scale-grazers' in sys.argv
    with_2504 = '--with-2504' in sys.argv
    argv = sys.argv[1:]
    opt = {}
    for k in ('--rec-days', '--rec-offset-h'):
        if k in argv:
            i = argv.index(k)
            opt[k] = float(argv[i + 1])
            del argv[i:i + 2]
    rec_days, rec_off = opt.get('--rec-days', 30.0), opt.get('--rec-offset-h', 0.0)
    args = [a for a in argv if a not in ('--scale-grazers', '--with-2504')]
    src, dst, start = args[:3]
    c2chl = float(args[3]) if len(args) > 3 else None
    g = json.load(open(os.path.join(WORK, 'grid_info.json')))
    nx, ny, nr = g['nx'], g['ny'], g['nr']
    dz = np.array(g['delr'])
    zc = np.cumsum(dz) - 0.5 * dz
    depth = -np.fromfile(os.path.join(WORK, g['bathy']), '>f4').reshape(ny, nx)
    if os.path.exists(dst):
        sys.exit(f'{dst} exists: remove it first')
    shutil.copytree(src, dst)
    ob = json.load(open(os.path.join(dst, 'obcs_files.json')))
    tag = ob['tag']
    b = casts(with_2504)
    tN = b[b.cruise == '2411'].time.median()
    tJ = b[b.cruise.isin(['2501', '2502'])].time.median()
    prof = {'chl_N': profiles(b, zc, ['2411'], 'Chl-a'), 'chl_J': profiles(b, zc, ['2501', '2502'], 'Chl-a')}
    for v in ('NO3', 'PO4', 'SIL'):
        prof[v] = profiles(b, zc, ['2501', '2502'], v)
    tA = None
    if with_2504:
        tA = b[b.cruise == '2504'].time.median()
        prof['chl_A'] = profiles(b, zc, ['2504'], 'Chl-a')
        for v in ('NO3', 'PO4', 'SIL'):
            prof[v + '_A'] = profiles(b, zc, ['2504'], v)
    info = {k: {z: dict(n_casts=n, surface=float(p[0]), z50=float(np.interp(50, zc, p)))
                for z, (p, n) in v.items()} for k, v in prof.items()}
    info.update(t_nov=str(tN), t_jan=str(tJ), t_apr=str(tA), rec_days=rec_days, rec_offset_h=rec_off)
    # weights
    wz = np.clip((Z_PARENT - zc) / (Z_PARENT - Z_OBS), 0, 1)        # 1 above 200 m, 0 below 300 m
    shelf_w = lambda d: np.clip((500.0 - d) / 400.0, 0, 1)

    def obs_field(key, d2, wchl=None, wapr=0.0):
        """(nr, ...) observation profile field for depths d2 (any shape). wapr: weight of the 2504 profiles."""
        sh = shelf_w(d2)
        if key == 'chl':
            if wapr > 0:
                ps = (1 - wapr) * prof['chl_J']['shelf'][0] + wapr * prof['chl_A']['shelf'][0]
                po = (1 - wapr) * prof['chl_J']['offshore'][0] + wapr * prof['chl_A']['offshore'][0]
            else:
                ps = (1 - wchl) * prof['chl_N']['shelf'][0] + wchl * prof['chl_J']['shelf'][0]
                po = (1 - wchl) * prof['chl_N']['offshore'][0] + wchl * prof['chl_J']['offshore'][0]
        else:
            ps = (1 - wapr) * prof[key]['shelf'][0] + (wapr * prof[key + '_A']['shelf'][0] if wapr > 0 else 0)
            po = (1 - wapr) * prof[key]['offshore'][0] + (wapr * prof[key + '_A']['offshore'][0] if wapr > 0 else 0)
        return ps.reshape((-1,) + (1,) * d2.ndim) * sh + po.reshape((-1,) + (1,) * d2.ndim) * (1 - sh)

    def time_w(t):
        return float(np.clip((t - tN) / (tJ - tN), 0, 1))

    def apr_w(t):
        return 0.0 if tA is None else float(np.clip((t - tJ) / (tA - tJ), 0, 1))

    # parent class fractions and Chl:C per level (domain means of the parent-derived IC)
    rd = lambda f, shp: np.fromfile(f, '>f4').reshape(shp).astype(float)
    icf = lambda v: os.path.join(dst, 'init', f'{v}_{tag}')
    chl_names = [f'Chl{i}' for i in range(1, 6)]
    c_names = [f'c{i}' for i in range(1, 8)]
    CH = np.stack([rd(icf(v), (nr, ny, nx)) for v in chl_names])
    CC = np.stack([rd(icf(v), (nr, ny, nx)) for v in c_names])
    wet = depth[None] > (np.cumsum(dz) - dz)[:, None, None]
    s = lambda X: np.array([[X[i][k][wet[k]].sum() for k in range(nr)] for i in range(X.shape[0])])
    chs, ccs = s(CH), s(CC[:5])
    frac = chs / np.maximum(chs.sum(0), 1e-30)                    # (5, nr)
    theta = chs / np.maximum(ccs, 1e-30)                           # mg Chl / mmol C, (5, nr)
    if c2chl:     # C:Chl [g/g], but not below each class's darwin floor 12.011/chl2cmax
        import re
        tr = open(os.path.join(HERE, '..', 'input', 'data.traits')).read()
        cm = []
        for tok in re.search(r'CHL2CMAX\s*=\s*([^\n]*)', tr, re.I).group(1).split(','):
            tok = tok.strip()
            if tok:
                n, v = (tok.split('*') + [None])[:2] if '*' in tok else (1, tok)
                cm += [float(v)] * int(n)
        chl2cmax = np.array(cm[:5])                                # mg Chl / mmol C, phyto 1-5
        theta = np.repeat(np.minimum(12.011 / c2chl, chl2cmax)[:, None], nr, axis=1)
        info['c2chl_by_class_g_per_g'] = (12.011 / theta[:, 0]).round(1).tolist()
    info['c2chl_g_per_g'] = c2chl if c2chl else 'parent'
    info['parent_c2chl_surface_g_per_g'] = float(12.011 * ccs[:, 0].sum() / max(chs[:, 0].sum(), 1e-30))
    info['parent_class_fraction_surface'] = frac[:, 0].round(3).tolist()

    def apply(fields, d2, t):
        """fields: dict name -> array (nr, ...) of parent values; returns modified."""
        out = dict(fields)
        w = wz.reshape((-1,) + (1,) * d2.ndim)
        for key, name in (('NO3', 'NO3'), ('PO4', 'PO4'), ('SIL', 'SiO2')):
            out[name] = w * obs_field(key, d2, wapr=apr_w(t)) + (1 - w) * fields[name]
        chl_obs = obs_field('chl', d2, time_w(t), apr_w(t))
        phy_c_old = sum(fields[f'c{i}'] for i in range(1, 6))
        new_c = []
        for i in range(5):
            fr = frac[i].reshape((-1,) + (1,) * d2.ndim)
            th = theta[i].reshape((-1,) + (1,) * d2.ndim)
            chl_i = w * chl_obs * fr + (1 - w) * fields[f'Chl{i + 1}']
            c_i = w * (chl_obs * fr / np.maximum(th, 1e-12)) + (1 - w) * fields[f'c{i + 1}']
            out[f'Chl{i + 1}'], out[f'c{i + 1}'] = chl_i, c_i
            new_c.append(c_i)
        scale = sum(new_c) / np.maximum(phy_c_old, 1e-12)
        scale = np.where(phy_c_old > 1e-9, scale, 1.0)
        if scale_grazers:
            for i in (6, 7):
                out[f'c{i}'] = fields[f'c{i}'] * scale
        return out

    names = ['NO3', 'PO4', 'SiO2'] + chl_names + c_names
    t0 = pd.Timestamp(start)
    # IC
    F = {v: rd(icf(v), (nr, ny, nx)) for v in names}
    F = apply(F, depth, t0)
    for v in names:
        F[v].astype('>f4').tofile(icf(v))
    # OBCS (records every 30 d from the start)
    for side in ob['sides']:
        d_line = {'west': depth[:, 0], 'east': depth[:, -1], 'south': depth[0], 'north': depth[-1]}[side]
        fn = lambda v: os.path.join(dst, 'obcs', ob['files'][v][side])
        nrec = os.path.getsize(fn('NO3')) // (4 * nr * d_line.size)
        Fo = {v: rd(fn(v), (nrec, nr, d_line.size)) for v in names}
        for r in range(nrec):
            t = t0 + pd.Timedelta(hours=rec_off) + pd.Timedelta(days=rec_days * r)
            m = apply({v: Fo[v][r] for v in names}, d_line, t)
            for v in names:
                Fo[v][r] = m[v]
        for v in names:
            Fo[v].astype('>f4').tofile(fn(v))
    info['note'] = ('CalCOFI PRELIMINARY bottles; obs 0-200 m, taper to parent at 300 m; shelf weight '
                    'clip((500-depth)/400); grazers ' + ('scaled with phyto C (assumed)' if scale_grazers else 'kept at parent'))
    json.dump(info, open(os.path.join(dst, 'calcofi_bgc.json'), 'w'), indent=1, default=str)
    for k in [k for k in ('chl_N', 'chl_J', 'chl_A', 'NO3', 'PO4', 'SIL', 'NO3_A', 'PO4_A', 'SIL_A') if k in info]:
        print(k, {z: (v['n_casts'], round(v['surface'], 2), round(v['z50'], 2)) for z, v in info[k].items()})
    print('parent surface class fractions', info['parent_class_fraction_surface'])
    sfc = rd(icf('Chl1'), (nr, ny, nx))[0]
    print('new IC surface Chl01 median', float(np.median(sfc[depth > 0])), '| wrote', dst)


if __name__ == '__main__':
    main()
