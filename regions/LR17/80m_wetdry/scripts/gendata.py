#!/usr/bin/env python3
"""Inputs for the LR17 80 m wetting-drying prototype (Pertuis Charentais).

Usage:  python3 gendata.py START DAYS [OUTDIR] [--zb=FILE] [--sed] [--transport] [--stamp-times]
                          [--tide=fes2022|fes2014] [--fes-boundary[=eta] [--band-download=START_END]]
  LR17_ROOT (environment, optional): the wetting_drying checkout that holds
         data/ and MITgcm/ (default: the checkout of this script), so a git
         worktree without them can still run gendata.py
  START  first time, YYYY-MM-DD (00 UTC; model time 0), within the IBI and
         ERA5 downloads (data/LR17/ibi, data/LR17/era5: 2022-02-24 .. 03-18)
  DAYS   run length in days
  OUTDIR default: run_START next to this script
Framework: Adeola's LR17 80 m Darwin P1 (grid, EXF interpolation, prescribed
IBI open boundaries), physics only here (Darwin off), with pkg/wetdry, KPP and
Manning drag. Writes bathy.bin, eta.bin, T.bin, S.bin, the open-boundary
files (obcs/), the ERA5 forcing (exf/) and the namelists.

Grid: 936 x 875 cells of 0.0008 deg (62 x 89 m), lon -1.7001 + ..., lat
45.7499 + ... (Adeola's LR17). Bed: data/LR17/work/LR17_zb.bin
(scripts/LR17_bathy.py: tidal flats and marsh seaward of the sea defences,
m MSL, land = 99). Model rest level DATUM = 5.5 m above MSL (floodable bed
up to +5 m), 14 levels (8 x 2.5, 2 x 5, 3 x 10, 8.5 m; 68.5 m in all, the
bed reaches -60.9 m). Edge cells of the western and southern open
boundaries shallower than 2 m are closed (a prescribed inflow on a cell
that dries does not work).
Mean sea level: IBI zos is relative to its own reference; here MSL is the
time mean of zos over the whole download at each IBI point (removes the
mean dynamic topography of the 3-km grid too; the tide and surge stay).
Initial state (START): IBI daily T, S interpolated to the model grid and
levels (nearest IBI sea point where IBI is land), sea level from hourly zos,
at rest; dry cells carry a 5 cm film.
Open boundaries (west column, south row; hourly, obcs/): sea level zos
(tide + surge), normal and tangential velocity = daily 3-D IBI velocity
minus its depth mean plus the hourly depth-mean ubar, vbar (the barotropic
tide), T and S daily, interpolated to the hour.
Atmosphere (exf/): ERA5 hourly on its 0.25 deg grid, interpolated by EXF:
10 m wind, 2 m air temperature, specific humidity (from dewpoint and sea-
level pressure), downward short- and long-wave radiation, precipitation.
Rivers (data.rivers, code/lr17_rivers.F): daily Hub'Eau discharge
(scripts/LR17_rivers.py) as point volume sources (addMass, selectAddFluid=1)
at the upstream tip of each river in the model: the model cell of the
carved river (OSM centreline, scripts/LR17_bathy.py) farthest from the open
boundaries along wet cells; the Charente (+ Boutonne, which joins it
upstream of the domain) at the head of its estuary channel, where the
bathymetry ends. Fresh (S = 0) at RIVER_T. Records at noon of each day,
linear in time.
IBI time stamps: the hourly fields (zos, ubar, vbar) are hourly means
stamped at the start of their hour (like the daily means, stamped at
00 UTC and centred at 12 UTC), so a record stamped HH:00 is the value at
HH:30. The value at the full hour T (initial sea level, open-boundary
records) is the mean of the records stamped T-1h and T. Taking the stamps
as instants made IBI ~20 min and LR17 ~30 min early against the SHOM
gauges (2026-10-03); --stamp-times restores that (old) behaviour.
--tide=fes2022 (or fes2014): the tide on the open boundaries and in the
initial sea level comes from the FES atlas instead of IBI: sea level from
FES2022b (or FES2014b) and depth-mean current from FES2014, predicted from
20 constituents (M2 S2 N2 K2 2N2 MU2 NU2 L2 T2 K1 O1 P1 Q1 M3 M4 MS4 MN4 N4 S4
M6) with nodal corrections (scripts/fes_tide.py, UTide); IBI keeps the
non-tidal part (surge, mean flow), its hourly sea level and depth-mean
current low-pass filtered with the Godin filter (24, 24, 25 h running means).
Reason: IBI has a third of the observed MS4 at the boundary and ~6 % less M2
than FES2022b (LR17_fes_compare.py, 2026-10-05).
--fes-boundary (LR17 v3 baseline, 2026-10-21): on the open boundaries the IBI
semidiurnal and quarter-diurnal bands (1.80-2.06 and 3.6-4.2 cycles per day,
sea level and depth-mean current) are replaced by the FES tide of M2 S2 N2
K2 2N2 MU2 NU2 L2 T2 M4 MS4 (sea level FES2022b, current FES2014, nodal
corrections; scripts/LR17_bc_overtides.py swap_terms). The IBI bands are
band-passed (FFT) over the 50-day IBI download ibi_*_2022-02-20_2022-04-10,
so the run must lie inside it, at least 8 days from both ends (checked);
--band-download=START_END names another ssh+cur download for the band-pass,
e.g. 2022-01-25_2022-05-05 for runs that reach into April. Everything else
(diurnal tide, surge, mean flow, M6, T, S, initial state) stays IBI. Reason:
IBI M2 is ~4 % below FES2022b at the boundary and IBI MS4 a third of the
observed; v3 (run_fesd2b_16) cut the gauge M2 error by a third and the low-water
bias from +0.22 to +0.14 m against v2 (docs/LR17/LR17_tidal_analysis.pptx).
--fes-boundary=eta swaps the bands in the boundary sea level only and keeps
IBI's depth-mean current (test of the v3 current bias at the AWAC, 2026-10-24).
Not with --tide or --stamp-times.
--transport: the hourly IBI depth-mean current on the open boundaries is
scaled at each boundary cell by (IBI depth + eta)/(model depth + eta),
limited to 0.5 .. 2, so that each cell passes IBI's volume transport
rather than IBI's velocity times the model's (different) depth (IBI
deptho, ibi_static.nc). Without it the model's mean sea level drifted by
+-0.8 m with the spring-neap cycle (2026-10-02).
--sed: rain and evaporation on (useRealFreshWaterFlux, ERA5); mud and fine
sand (pTracers 1, 2) with pkg/sediment, layered bed
and IBI waves (sedWaveMode=3): hourly Hs (VHM0), mean period (VTM02) and
direction (VMDR) interpolated from IBI sea points only onto the model grid
(sheltered bays keep their own IBI values), Hs damped as exp(-d/WAVE_SHELTER)
with the distance d from the IBI sea and as exp(-d/WAVE_FLAT) with the
distance from water deeper than FLAT_DEPTH (swell dissipates across the
flats, which the 2.5 km IBI grid does not resolve); pkg/sediment caps Hs
at 0.5 h. Initial bed: the EMODnet Geology seabed substrate (Folk 7,
1:250 000) on the model grid, data/LR17/work/LR17_bed_{mud,sand}.bin
(scripts/LR17_bedmap.py); without it mud on the flats (bed above -5 m MSL)
and sand below -2 m MSL. Each river carries its own mud concentration
(RIVER_MUD_EACH, map sedRiverCFile), open-boundary water is clean.
"""
import datetime as dt
import json
from collections import deque
import os
import re
import sys

import numpy as np
import xarray as xr
from scipy import ndimage
from scipy.interpolate import RegularGridInterpolator

HERE = os.path.dirname(os.path.abspath(__file__))
TOP = os.path.dirname(HERE)                   # 80m_wetdry/ (scripts/, input/, data/)
ROOT = os.environ.get('LR17_ROOT', TOP)       # holds data/LR17/
D = os.path.join(ROOT, 'data', 'LR17')
NX, NY, DG = 936, 875, 0.0008
LON0, LAT0 = -1.7001, 45.7499
DATUM = 5.5
DELR = [2.5] * 8 + [5.0] * 2 + [10.0] * 3 + [8.5]
NR = len(DELR)
DT = 10.0
ZOB = -2.0                     # open-boundary cells must be deeper than this
TAG = '2022-02-24_2022-03-18'
RIVER_FILE = os.path.join(D, 'hubeau', 'LR17_rivers_2022-02-01_2022-03-31.json')
WAVE_SHELTER = 2000.0           # m, e-folding of Hs away from the IBI sea
WAVE_FLAT, FLAT_DEPTH = 1000.0, 2.0   # m: e-folding of Hs over the flats, from water deeper than 2 m
MUD0, SAND0 = 20.0, 50.0        # kg/m2, initial bed without the EMODnet map (illustrative)
RIVER_MUD = 0.03                # kg/m3 mud in river water without RIVER_MUD_EACH
# suspended solids (MES, Hub'Eau qualite_rivieres, parameter 1305) in the
# rivers, February-March 2022, at the most seaward non-tidal station:
# Charente at Taillebourg (5-17 mg/l), Sevre Niortaise at Marans (9-10),
# Lay at La Claye (8-22), Seudre at Saujon (4-11); kg/m3 of mud
RIVER_MUD_EACH = {'Charente': 0.010, 'Sevre': 0.010, 'Lay': 0.015, 'Seudre': 0.008}
GAUGES = {'La Rochelle-La Pallice': (46.1585, -1.2206), "Ile d'Aix": (46.0074, -1.1743)}
RIVER_T = 9.0                    # river temperature, deg C (Feb-Mar)
# model river: (OSM name of the carved river or None, gauged rivers it carries)
RIVERS = {'Charente': (None, ('Charente', 'Boutonne')),
          'Sevre': ('La Sèvre Niortaise', ('Sevre',)),
          'Lay': ('Le Lay', ('Lay',)),
          'Seudre': ('La Seudre', ('Seudre',))}
lon = LON0 + (np.arange(NX) + 0.5) * DG
lat = LAT0 + (np.arange(NY) + 0.5) * DG
zc = np.cumsum(DELR) - 0.5 * np.array(DELR)            # level centres below rest
zmsl = np.maximum(zc - DATUM, 0.5)                     # depth below MSL (IBI)


def ibi(key):
    return xr.open_dataset(os.path.join(D, 'ibi', f'ibi_{key}_{TAG}.nc'))


def fill_nearest(a):
    """fill NaNs (IBI land) with the nearest valid value, horizontally"""
    bad = ~np.isfinite(a)
    if bad.all() or not bad.any():
        return a
    idx = ndimage.distance_transform_edt(bad, return_distances=False, return_indices=True)
    return a[tuple(idx)]


def to_grid(field2d, glon, glat, plon, plat):
    """bilinear from the IBI grid to points (plon, plat), after filling land"""
    f = RegularGridInterpolator((glat, glon), fill_nearest(field2d), bounds_error=False,
                                fill_value=None)
    return f(np.c_[plat.ravel(), plon.ravel()]).reshape(plat.shape)


def sea_grid(field2d, glon, glat, plon, plat, model):
    """IBI sea values only (normalised bilinear interpolation, so a
    sheltered bay keeps its own values instead of the nearest strait's);
    model cells beyond the IBI sea take the nearest interpolated value."""
    ok = np.isfinite(field2d)
    num = RegularGridInterpolator((glat, glon), np.where(ok, field2d, 0.0), bounds_error=False,
                                  fill_value=None)(np.c_[plat.ravel(), plon.ravel()])
    den = RegularGridInterpolator((glat, glon), ok.astype(float), bounds_error=False,
                                  fill_value=None)(np.c_[plat.ravel(), plon.ravel()])
    out = np.where(den > 0.05, num / np.maximum(den, 0.05), np.nan).reshape(plat.shape)
    return fill_nearest(out)


def vert(profile_stack, depth):
    """IBI levels -> model levels: (ndepth, ...) at zmsl, nearest above the bed"""
    out = np.empty((NR,) + profile_stack.shape[1:])
    for k, z in enumerate(zmsl):
        j = np.searchsorted(depth, z)
        j = min(max(j, 1), len(depth) - 1)
        w = (z - depth[j - 1]) / (depth[j] - depth[j - 1])
        lo, hi = profile_stack[j - 1], profile_stack[j]
        v = (1 - w) * lo + w * hi
        # below the IBI bed: carry the deepest valid value down
        out[k] = np.where(np.isfinite(v), v, np.nan)
    for k in range(1, NR):
        out[k] = np.where(np.isfinite(out[k]), out[k], out[k - 1])
    return out


def hourly_from_daily(daily, dtimes, htimes):
    """linear in time, daily means centred at 12 UTC"""
    td = (dtimes - dtimes[0]) / np.timedelta64(1, 'h') + 12.0
    th = (htimes - dtimes[0]) / np.timedelta64(1, 'h')
    i1 = np.clip(np.searchsorted(td, th), 1, len(td) - 1)
    w = np.clip((th - td[i1 - 1]) / (td[i1] - td[i1 - 1]), 0, 1)
    return (1 - w)[:, None, None] * daily[i1 - 1] + w[:, None, None] * daily[i1]


def wet_distance(model):
    """steps from the open boundaries (west column, south row) through model cells"""
    dist = np.full(model.shape, -1)
    q = deque()
    for j, i in list(zip(*np.where(model[:, :1]))) + list(zip(*np.where(model[:1, :]))):
        if dist[j, i] < 0:
            dist[j, i] = 0
            q.append((j, i))
    ny, nx = model.shape
    while q:
        j, i = q.popleft()
        for jj, ii in ((j + 1, i), (j - 1, i), (j, i + 1), (j, i - 1)):
            if 0 <= jj < ny and 0 <= ii < nx and model[jj, ii] and dist[jj, ii] < 0:
                dist[jj, ii] = dist[j, i] + 1
                q.append((jj, ii))
    return dist


def river_tips(model, lon, lat):
    """source cell (j, i) of each model river"""
    sys.path.insert(0, HERE)
    argv, sys.argv = sys.argv, sys.argv[:1]        # LR17_bathy reads argv on import
    from LR17_bathy import rasterize
    sys.argv = argv
    dist = wet_distance(model)
    els = json.load(open(os.path.join(D, 'osm', 'osm_mainrivers.json')))['elements']
    tips = {}
    for name, (osm, _) in RIVERS.items():
        if osm is None:          # Charente: OSM way stops upstream of the model
            area = (lon[None, :] > -1.1) & (lat[:, None] > 45.88) & (lat[:, None] < 46.0)
        else:
            area = ndimage.binary_dilation(rasterize([e for e in els if e['tags'].get('name') == osm]),
                                           np.ones((3, 3), bool))
        tips[name] = np.unravel_index(np.argmax(np.where(area & model, dist, -1)), model.shape)
    return tips


def w32(path, a):
    np.asarray(a, dtype='>f4').tofile(path)


def at_hours(a, times, targets, stamp_times=False):
    """IBI hourly means (stamped at the start of the hour) at the full
    hours `targets`: mean of the records stamped T-1h and T; with
    stamp_times the record stamped T (stamps taken as instants)."""
    k = np.searchsorted(times, targets)
    if not np.all(times[k] == targets) or (not stamp_times and k.min() < 1):
        sys.exit('IBI hourly records do not cover the requested hours')
    return a[k] if stamp_times else 0.5 * (a[k - 1] + a[k])


def godin(a):
    """Godin low-pass (24, 24, 25 h running means) along axis 0 of an hourly
    series: removes the tides, keeps surge and mean flow"""
    for w in (24, 24, 25):
        a = ndimage.uniform_filter1d(a, w, axis=0, mode='nearest')
    return a


def atlas_tide(tide, kind, times, plat, plon):
    """tide (ntimes, npt) at points from the FES atlas: kind 'eta' (m) or
    'u', 'v' (m/s)"""
    sys.path.insert(0, HERE)
    import fes_tide
    fes_tide.D = D
    src = 'FES2022b' if (tide == 'fes2022' and kind == 'eta') else 'FES2014'
    A = fes_tide.atlas(src, kind)
    if len(A) < len(fes_tide.CONS):
        sys.exit(f'{src} {kind}: only {len(A)} of {len(fes_tide.CONS)} constituents found')
    amp, pha, names = fes_tide.constants(A, plat, plon)
    return fes_tide.predict_fast(times, amp, pha, names, float(np.mean(plat))).T


def main(start, days, out, sed=False, transport=False, stamp_times=False, tide=None,
         fes_boundary=False, band_tag=None, zb_file=None):
    t0 = np.datetime64(start + 'T00:00')
    t1 = t0 + np.timedelta64(int(round(days * 24)) + 1, 'h')
    os.makedirs(os.path.join(out, 'obcs'), exist_ok=True)
    os.makedirs(os.path.join(out, 'exf'), exist_ok=True)

    # --- bathymetry
    zb_file = zb_file or os.path.join(D, 'work', 'LR17_zb.bin')
    zb = np.fromfile(zb_file, '>f4').reshape(NY, NX).astype(float)
    print(f'bed: {zb_file}')
    model = zb < 90
    model[:, 0] &= zb[:, 0] < ZOB
    model[0, :] &= zb[0, :] < ZOB
    r_low = np.where(model, zb - DATUM, 0.0)
    # no sliver bottom cells: a bottom cell thinner than 0.125 m (layers 2-4)
    # or 0.25 m (deeper) is removed (bed raised to the interface above).
    # Slivers (hFac ~ 0.01) squeezed by r* at low water blew up the first
    # test. Layer 1 is left exact.
    rf = -np.concatenate([[0.0], np.cumsum(DELR)])
    k = np.clip(np.searchsorted(-rf, -r_low, side='left') - 1, 0, NR - 1)
    dz_bot = rf[k] - r_low
    thin = model & (k >= 1) & (dz_bot < np.where(k <= 3, 0.125, 0.25))
    r_low = np.where(thin, rf[k], r_low)
    print(f'sliver bottom cells raised to the layer above: {thin.sum()}')
    w32(os.path.join(out, 'bathy.bin'), r_low)
    obw = model[:, 0]
    obs = model[0, :]

    # --- sea level: MSL = time mean of zos over the download at each point
    ssh = ibi('ssh')
    glon, glat = ssh.longitude.values, ssh.latitude.values
    zos = ssh.zos.values
    msl = np.nanmean(zos, axis=0)
    zos = zos - msl
    if tide:                                   # IBI keeps only the non-tidal part
        zos = np.where(np.isfinite(zos), godin(np.where(np.isfinite(zos), zos, 0.0)), np.nan)
    stime = ssh.time.values
    sel = (stime >= t0) & (stime <= t1)
    hts = stime[sel]
    LO, LA = np.meshgrid(lon, lat)
    eta0 = to_grid(at_hours(zos, stime, np.array([t0]), stamp_times)[0], glon, glat, LO, LA)
    if tide:
        eta0 = eta0 + atlas_tide(tide, 'eta', np.array([t0]), LA.ravel(), LO.ravel())[0] \
            .reshape(LA.shape)
        print(f'initial sea level: IBI non-tidal + {tide} tide')
    # wet only where connected to the open boundaries: low pockets behind the
    # flats (below the start level but cut off) start dry, else they drain in
    # a burst at t = 0
    wet = model & (eta0 - DATUM > r_low + 0.05)
    lab, _ = ndimage.label(wet)
    keep = np.setdiff1d(np.union1d(lab[:, 0], lab[0, :]), [0])
    sea = np.isin(lab, keep)
    print(f'initial wet cells cut off from the sea, started dry: {(wet & ~sea).sum()}')
    eta = np.where(sea, eta0 - DATUM, np.where(model, r_low + 0.05, 0.0))
    w32(os.path.join(out, 'eta.bin'), eta)

    # --- initial T, S from the daily fields of the start day
    tem, sal = ibi('temp3d'), ibi('sal')
    dtimes = tem.time.values
    kd = int(np.argmin(np.abs(dtimes - (t0 + np.timedelta64(12, 'h')))))
    dep = tem.depth.values
    for name, da in (('T', tem.thetao), ('S', sal.so)):
        prof = vert(da.values[kd], dep)                      # (NR, lat, lon) on IBI grid
        f3 = np.array([to_grid(prof[k], glon, glat, LO, LA) for k in range(NR)])
        w32(os.path.join(out, f'{name}.bin'), np.where(model[None], f3, 0.0))

    # --- open boundaries, hourly
    cur, c3 = ibi('cur'), ibi('cur3d')
    htime = cur.time.values
    hsel = (htime >= t0) & (htime <= t1)
    if tide:
        ub0, vb0 = [np.where(np.isfinite(a), godin(np.where(np.isfinite(a), a, 0.0)), np.nan)
                    for a in (cur.ubar.values, cur.vbar.values)]
    else:
        ub0, vb0 = cur.ubar.values, cur.vbar.values
    ubar = at_hours(ub0, htime, htime[hsel], stamp_times)
    vbar = at_hours(vb0, htime, htime[hsel], stamp_times)
    edges = {'West': (np.full(NY, lon[0] - DG / 2), lat, obw),
             'South': (lon, np.full(NX, lat[0] - DG / 2), obs)}
    zcell = {'West': r_low[:, 0] + DATUM, 'South': r_low[0, :] + DATUM}   # bed, m MSL
    if transport:
        st = xr.open_dataset(os.path.join(D, 'ibi', 'ibi_static.nc'))
        hibi = (st.latitude.values, st.longitude.values, fill_nearest(st.deptho.values))
    fields = {'T': tem.thetao.values, 'S': sal.so.values, 'U': c3.uo.values, 'V': c3.vo.values}
    for side, (plon, plat, open_) in edges.items():
        pts = np.c_[plat, plon]
        n = len(plon)

        def along(a2d):
            return RegularGridInterpolator((glat, glon), fill_nearest(a2d), bounds_error=False,
                                           fill_value=None)(pts)
        eta_b = np.array([along(z) for z in at_hours(zos, stime, hts, stamp_times)]) - DATUM
        ub = np.array([along(u) for u in ubar])
        vb = np.array([along(v) for v in vbar])
        if tide:
            eta_b = eta_b + atlas_tide(tide, 'eta', hts, plat, plon)
            ub = ub + atlas_tide(tide, 'u', htime[hsel], plat, plon)
            vb = vb + atlas_tide(tide, 'v', htime[hsel], plat, plon)
            o = open_
            print(f'{side}: {tide} tide added; eta std {eta_b[:, o].std():.2f} m, '
                  f'u std {ub[:, o].std():.3f}, v std {vb[:, o].std():.3f} m/s')
        if fes_boundary:
            sys.path.insert(0, HERE)
            import LR17_bc_overtides as bco
            T = bco.swap_terms(D, pts, hts, bco.V3_CONS, bco.V3_BANDS, tag_band=band_tag)
            eta_b = eta_b - T['b_eta'] + T['f_eta']
            if fes_boundary != 'eta':
                ub = ub - T['b_u'] + T['f_u']
                vb = vb - T['b_v'] + T['f_v']
            o = open_
            print(f'{side}: FES boundary ({fes_boundary}): IBI band rms {np.sqrt((T["b_eta"][:, o] ** 2).mean()):.3f} m '
                  f'replaced by FES {np.sqrt((T["f_eta"][:, o] ** 2).mean()):.3f} m')
        if transport:
            hi = RegularGridInterpolator(hibi[:2], hibi[2], bounds_error=False,
                                         fill_value=None)(pts)
            em = eta_b + DATUM                                  # m MSL, hourly
            fac = np.clip((hi[None] + em) / np.maximum(em - zcell[side][None], 0.5), 0.5, 2.0)
            ub, vb = ub * fac, vb * fac
            o = np.broadcast_to(open_[None], fac.shape)
            print(f'{side}: transport factor median {np.median(fac[o]):.3f}, '
                  f'5-95% {np.percentile(fac[o], 5):.2f}..{np.percentile(fac[o], 95):.2f}, '
                  f'clipped {(fac[o] >= 2.0).mean() * 100:.1f}% at 2')
        prof = {}
        for key, arr in fields.items():
            dly = np.array([[along(p) for p in vert(arr[d], dep)] for d in range(len(dtimes))])
            prof[key] = hourly_from_daily(dly.reshape(len(dtimes), NR, n), dtimes, hts)
        for key in ('U', 'V'):                                # baroclinic part only
            prof[key] = prof[key] - prof[key].mean(axis=1, keepdims=True)
        prof['U'] = prof['U'] + ub[:, None, :]
        prof['V'] = prof['V'] + vb[:, None, :]
        m = open_[None, None, :]
        for key in ('T', 'S', 'U', 'V'):
            w32(os.path.join(out, 'obcs', f'OB{side[0]}{key.lower()}.bin'),
                np.where(m, prof[key], 0.0))
        w32(os.path.join(out, 'obcs', f'OB{side[0]}eta.bin'), np.where(open_[None], eta_b, 0.0))
    nrec = len(hts)

    # --- ERA5 on its grid (EXF interpolation), hourly from START
    ei = xr.open_dataset(os.path.join(D, 'era5', f'era5_{TAG}_instant.nc'))
    ea = xr.open_dataset(os.path.join(D, 'era5', f'era5_{TAG}_accum.nc'))
    et = ei.valid_time.values
    es = (et >= t0) & (et <= t1)
    elat = ei.latitude.values[::-1]                       # EXF: south to north
    elon = ei.longitude.values
    flip = lambda a: a[:, ::-1, :]
    td = flip(ei.d2m.values[es]) - 273.15
    pa = flip(ei.msl.values[es])
    e = 6.112 * np.exp(17.67 * td / (td + 243.5)) * 100.0      # vapour pressure, Pa
    q = 0.622 * e / (pa - 0.378 * e)
    exf = {'uwind': flip(ei.u10.values[es]), 'vwind': flip(ei.v10.values[es]),
           'atemp': flip(ei.t2m.values[es]), 'aqh': q,
           'swdown': flip(ea.ssrd.values[es]) / 3600.0, 'lwdown': flip(ea.strd.values[es]) / 3600.0,
           'precip': flip(ea.tp.values[es]) / 3600.0}
    for k, a in exf.items():
        w32(os.path.join(out, 'exf', f'ERA5_{k}.bin'), a)

    # --- rivers: daily discharge at noon, from the day before START
    riv = json.load(open(RIVER_FILE))
    tips = river_tips(model, lon, lat)
    days_r = [dt.date.fromisoformat(d[:10]) for d in riv['Charente']['dates']]
    d0 = dt.date.fromisoformat(start) - dt.timedelta(days=1)
    k0 = days_r.index(d0)
    nday = int(np.ceil(days)) + 3
    Q = {n: sum(np.array(riv[g]['q'][k0:k0 + nday]) for g in gauged) for n, (_, gauged) in RIVERS.items()}
    rows = [f' rivN={len(RIVERS)}, rivNrec={nday}, rivT0={-12 * 3600.:.1f}, rivDt=86400.,',
            ' rivI=' + ', '.join(str(tips[n][1] + 1) for n in RIVERS) + ',',
            ' rivJ=' + ', '.join(str(tips[n][0] + 1) for n in RIVERS) + ',']
    for r, n in enumerate(RIVERS, 1):
        rows.append(f' rivQ({r},1:{nday})=' + ', '.join(f'{v:.2f}' for v in Q[n]) + ',')
    open(os.path.join(out, 'data.rivers'), 'w').write(
        '# LR17 rivers (' + ', '.join(RIVERS) + '), Hub\'Eau daily discharge, m3/s; '
        'record 1 = noon of the day before START\n &RIVERS_PARM\n' + '\n'.join(rows) + '\n &\n')
    for n in RIVERS:
        j, i = tips[n]
        print(f'river {n:9s} cell i={i + 1} j={j + 1} ({lon[i]:.3f}E {lat[j]:.3f}N, bed '
              f'{zb[j, i]:+.2f} m MSL): {Q[n].min():.1f}..{Q[n].max():.1f} m3/s')

    # tide-gauge cells (regions 1, 2 of the 10-min statistics 'gauges'):
    # nearest model cell deeper than 2 m below MSL to each SHOM gauge
    gmask = np.zeros(model.shape)
    for k, (gla, glo) in enumerate(GAUGES.values(), 1):
        d2 = np.where(model & (zb < -2.0), (lat[:, None] - gla) ** 2
                      + ((lon[None, :] - glo) * np.cos(np.radians(gla))) ** 2, np.inf)
        j, i = np.unravel_index(np.argmin(d2), d2.shape)
        gmask[j, i] = k
        print(f'gauge {k} {list(GAUGES)[k - 1]}: cell i={i + 1} j={j + 1} ({lon[i]:.4f}E {lat[j]:.4f}N, '
              f'bed {zb[j, i]:+.1f} m MSL)')
    w32(os.path.join(out, 'gauge_mask.bin'), gmask)
    write_namelists(out, start, days, nrec, obw, obs, elon, elat)
    if fes_boundary:
        fn = os.path.join(out, 'data')
        body = open(fn).read()
        if fes_boundary == 'eta':
            head = ('# v3 eta-only: FES boundary tide in sea level only (--fes-boundary=eta): IBI '
                    '1.80-2.06 and 3.6-4.2 cpd sea-level bands replaced by FES2022b M2 S2 N2 K2 2N2 MU2 '
                    'NU2 L2 T2 M4 MS4, IBI depth-mean current kept\n')
        else:
            head = ('# v3: FES boundary tide (--fes-boundary): IBI 1.80-2.06 and 3.6-4.2 cpd '
                    'bands replaced by FES2022b/FES2014 M2 S2 N2 K2 2N2 MU2 NU2 L2 T2 M4 MS4 '
                    '(LR17_bc_overtides.swap_terms, 50-day IBI band-pass)\n')
        open(fn, 'w').write(head + body)
    if zb_file != os.path.join(D, 'work', 'LR17_zb.bin'):
        fn = os.path.join(out, 'data')
        body = open(fn).read()
        open(fn, 'w').write(f'# bed: {os.path.basename(zb_file)} (--zb; LR17_bathy.py thalweg options, '
                            'v4 = --zmouth=-4 --ztip=-3)\n' + body)
    if sed:
        write_sed(out, t0, t1, zb, model, lon, lat)
        dpath = os.path.join(out, 'data')
        txt = open(dpath).read()
        assert 'useRealFreshWaterFlux=.FALSE.' in txt
        open(dpath, 'w').write(txt.replace('useRealFreshWaterFlux=.FALSE.',
                                           'useRealFreshWaterFlux=.TRUE.'))
    print(f'{out}: {model.sum()} model cells; OB west {obw.sum()}, south {obs.sum()} cells; '
          f'{nrec} hourly OB records; ERA5 {es.sum()} records {len(elon)} x {len(elat)}; '
          f'eta0 {np.nanmin(eta0[model]):+.2f} .. {np.nanmax(eta0[model]):+.2f} m MSL')


def write_namelists(out, start, days, nrec, obw, obs, elon, elat):
    d0 = start.replace('-', '')
    nsteps = int(round(days * 86400 / DT))
    dlist = ', '.join(f'{z:g}' for z in DELR)
    open(os.path.join(out, 'data'), 'w').write(f"""# LR17 80 m wetting-drying prototype: {start} + {days} d (model time 0 = {start} 00 UTC; rest level {DATUM} m MSL)
 &PARM01
 tRef={NR}*10., sRef={NR}*34.,
 no_slip_sides=.FALSE., no_slip_bottom=.FALSE.,
 viscAr=1.E-4, viscAh=2., viscAhGrid=1.E-2,
 viscC4Leith=1.5, viscC4Leithd=1.5, viscA4GridMax=0.5, useAreaViscLength=.TRUE.,
 diffKhT=1., diffKhS=1., diffKrT=1.E-5, diffKrS=1.E-5,
 bottomDragQuadratic=2.5E-3, selectImplicitDrag=2, implicitViscosity=.TRUE.,
 implicitDiffusion=.TRUE.,
 rhoConst=1025., gravity=9.81, eosType='JMD95Z',
 rigidLid=.FALSE., implicitFreeSurface=.TRUE., exactConserv=.TRUE.,
 nonlinFreeSurf=4, select_rStar=2,
 hFacMin=0.01, hFacMinDr=0.01, hFacInf=0.001, hFacSup=2.0,
 momStepping=.TRUE., vectorInvariantMomentum=.TRUE., upwindShear=.TRUE.,
 tempStepping=.TRUE., saltStepping=.TRUE., tempAdvScheme=33, saltAdvScheme=33,
 tempVertAdvScheme=3, saltVertAdvScheme=33,
 momImplVertAdv=.TRUE., tempImplVertAdv=.TRUE., saltImplVertAdv=.FALSE.,
 staggerTimeStep=.TRUE., useRealFreshWaterFlux=.FALSE.,
 selectAddFluid=1, temp_addMass={RIVER_T}, salt_addMass=0.,
 readBinaryPrec=32, writeBinaryPrec=32, useSingleCpuIO=.TRUE.,
 &
 &PARM02
 cg2dMaxIters=500, cg2dTargetResidual=1.E-12,
 &
 &PARM03
 nIter0=0, nTimeSteps={nsteps}, deltaT={DT}, abEps=0.1,
 dumpFreq=0., monitorSelect=2, monitorFreq=3600.,
 pChkptFreq={days * 86400:.1f}, chkptFreq=86400.,
 &
 &PARM04
 usingSphericalPolarGrid=.TRUE., xgOrigin={LON0}, ygOrigin={LAT0},
 delX={NX}*{DG}, delY={NY}*{DG},
 delR={dlist},
 &
 &PARM05
 bathyFile='bathy.bin', pSurfInitFile='eta.bin',
 hydrogThetaFile='T.bin', hydrogSaltFile='S.bin',
 &
""")
    open(os.path.join(out, 'data.pkg'), 'w').write(
        ' &PACKAGES\n useCAL=.TRUE.,\n useEXF=.TRUE.,\n useOBCS=.TRUE.,\n useKPP=.TRUE.,\n'
        ' useWETDRY=.TRUE.,\n useDiagnostics=.TRUE.,\n &\n')
    open(os.path.join(out, 'data.cal'), 'w').write(
        f" &CAL_NML\n TheCalendar='gregorian',\n startDate_1={d0},\n startDate_2=000000,\n &\n")
    # open boundaries: every model cell of the west column and south row
    def rle(mask):
        vals = [1 if m else 0 for m in mask]
        out_, prev, n = [], None, 0
        for v in vals + [None]:
            if v == prev:
                n += 1
                continue
            if prev is not None:
                out_.append(f'{n}*{prev}')
            prev, n = v, 1
        return ','.join(out_)
    open(os.path.join(out, 'data.obcs'), 'w').write(f""" &OBCS_PARM01
 OB_Iwest={rle(obw)},
 OB_Jsouth={rle(obs)},
 useOBCSprescribe=.TRUE.,
 OBWuFile='obcs/OBWu.bin', OBWvFile='obcs/OBWv.bin', OBWtFile='obcs/OBWt.bin',
 OBWsFile='obcs/OBWs.bin', OBWetaFile='obcs/OBWeta.bin',
 OBSuFile='obcs/OBSu.bin', OBSvFile='obcs/OBSv.bin', OBStFile='obcs/OBSt.bin',
 OBSsFile='obcs/OBSs.bin', OBSetaFile='obcs/OBSeta.bin',
 &
 &OBCS_PARM05
 &
""")
    grid = lambda k: (f" {k}_lon0={elon[0]}, {k}_lon_inc={elon[1] - elon[0]},\n"
                      f" {k}_lat0={elat[0]}, {k}_lat_inc={len(elat)}*{elat[1] - elat[0]},\n"
                      f" {k}_nlon={len(elon)}, {k}_nlat={len(elat)},\n")
    fl = ('uwind', 'vwind', 'atemp', 'aqh', 'swdown', 'lwdown', 'precip')
    nml2 = ''.join(f" {k}file='exf/ERA5_{k}.bin', {k}startdate1={d0}, {k}startdate2=000000, "
                   f"{k}period=3600.0,\n" for k in fl)
    nml4 = ''.join(grid(k) for k in fl)
    open(os.path.join(out, 'data.exf'), 'w').write(f""" &EXF_NML_01
 useAtmWind=.TRUE., exf_iprec=32, useExfCheckRange=.FALSE.,
 exf_albedo=0.1, ocean_emissivity=0.97,
 &
 &EXF_NML_02
{nml2} &
 &EXF_NML_03
 &
 &EXF_NML_04
{nml4} &
 &EXF_NML_OBCS
 obcsWstartdate1={d0}, obcsWstartdate2=000000, obcsWperiod=3600.0,
 obcsSstartdate1={d0}, obcsSstartdate2=000000, obcsSperiod=3600.0,
 &
""")
    open(os.path.join(out, 'data.kpp'), 'w').write(open(os.path.join(TOP, 'input', 'data.kpp')).read())
    open(os.path.join(out, 'data.wetdry'), 'w').write(
        ' &WETDRY_PARM01\n wdMinDepth=0.05,\n wdCritDepth=0.10,\n wdMonFreq=3600.,\n'
        ' wdManningN=0.02,\n wdMaxFroude=1.0,\n &\n')
    # maps every 15 min (tests), hourly (up to 2 days), 2-hourly (longer runs)
    fr = 900.0 if days < 1 else (3600.0 if days <= 2 else 7200.0)
    open(os.path.join(out, 'data.diagnostics'), 'w').write(f""" &DIAGNOSTICS_LIST
 fields(1:3,1)='ETAN    ','WDdepth','WDdryC ',
 fileName(1)='sfc', frequency(1)=-{fr}, timePhase(1)=0.,
 fields(1:4,2)='UVEL    ','VVEL    ','SALT    ','THETA   ',
 levels(1,2)=1., fileName(2)='top', frequency(2)=-{fr}, timePhase(2)=0.,
 &
 &DIAG_STATIS_PARMS
 diagSt_regMaskFile='gauge_mask.bin',
 nSetRegMskFile=1,
 set_regMask(1:2)=1,1,
 val_regMask(1)=1., val_regMask(2)=2.,
 stat_fields(1,1)='ETAN    ',
 stat_fName(1)='gauges', stat_freq(1)=600., stat_phase(1)=0.,
 stat_region(1:2,1)=1,2,
 &
""")
    open(os.path.join(out, 'eedata'), 'w').write(' &EEPARMS\n nTx=1,\n nTy=1,\n &\n')


def write_sed(out, t0, t1, zb, model, lon, lat):
    """IBI waves, initial bed, pTracers and pkg/sediment namelists (--sed)."""
    wv = ibi('wav')
    wt = wv.time.values
    ws_ = (wt >= t0) & (wt <= t1)
    glon, glat = wv.longitude.values, wv.latitude.values
    LO, LA = np.meshgrid(lon, lat)
    # IBI sea on the model grid (nearest IBI point is sea) and distance from it
    ii = np.clip(np.round((lon - glon[0]) / np.diff(glon).mean()).astype(int), 0, len(glon) - 1)
    jj = np.clip(np.round((lat - glat[0]) / np.diff(glat).mean()).astype(int), 0, len(glat) - 1)
    ibisea = np.isfinite(wv.VHM0.values[np.argmax(ws_)])[jj][:, ii]
    dx = DG * 111.2e3 * np.cos(np.radians(lat.mean()))
    shelter = np.exp(-ndimage.distance_transform_edt(~ibisea) * dx / WAVE_SHELTER)
    # swell loses its energy crossing the flats (friction, breaking), which
    # the 2.5 km IBI grid does not resolve: Hs damped with the distance from
    # water deeper than FLAT_DEPTH below MSL
    deep = model & (zb < -FLAT_DEPTH)
    shelter = shelter * np.exp(-ndimage.distance_transform_edt(~deep) * dx / WAVE_FLAT)
    fh = {k: open(os.path.join(out, f'wav_{k}.bin'), 'wb') for k in ('hs', 'per', 'dir')}
    for k in np.where(ws_)[0]:
        hs = sea_grid(wv.VHM0.values[k], glon, glat, LO, LA, model) * shelter
        # mean period (Tm02) for the orbital velocity, not the peak period of
        # the swell
        per = sea_grid(wv.VTM02.values[k], glon, glat, LO, LA, model)
        # VMDR: coming from, deg clockwise from north -> toward, ccw from +x (east),
        # interpolated as a unit vector
        th = np.radians(wv.VMDR.values[k] + 180.0)
        ue = to_grid(np.sin(th), glon, glat, LO, LA)
        un = to_grid(np.cos(th), glon, glat, LO, LA)
        dr = np.mod(np.degrees(np.arctan2(un, ue)), 360.0)
        for n, a in (('hs', hs), ('per', per), ('dir', dr)):
            np.where(model, np.maximum(a, 0.0), 0.0).astype('>f4').tofile(fh[n])
    for f in fh.values():
        f.close()
    nw = int(ws_.sum())
    fm = os.path.join(D, 'work', 'LR17_bed_mud.bin')
    if os.path.exists(fm):                       # EMODnet substrate (LR17_bedmap.py)
        mud = np.where(model, np.fromfile(fm, '>f4').reshape(NY, NX), 0.0)
        sand = np.where(model, np.fromfile(fm.replace('_mud', '_sand'), '>f4')
                        .reshape(NY, NX), 0.0)
        print(f'bed: EMODnet map, mud {mud[model].mean():.1f}, sand {sand[model].mean():.1f} kg/m2 mean')
    else:
        mud = np.where(model & (zb > -5.0), MUD0, 0.0)
        sand = np.where(model & (zb < -2.0), SAND0, 0.0)
    # mud concentration of the river water, at each river's source cell
    rv = open(os.path.join(out, 'data.rivers')).read()
    ri = [int(v) for v in re.search(r'rivI=([\d, ]+),', rv).group(1).split(',')]
    rj = [int(v) for v in re.search(r'rivJ=([\d, ]+),', rv).group(1).split(',')]
    rc = np.zeros((NY, NX))
    for name, i, j in zip(RIVERS, ri, rj):
        rc[j - 1, i - 1] = RIVER_MUD_EACH[name]
    w32(os.path.join(out, 'river_mud.bin'), rc)
    w32(os.path.join(out, 'bed_mud.bin'), mud)
    w32(os.path.join(out, 'bed_sand.bin'), sand)
    open(os.path.join(out, 'data.ptracers'), 'w').write(""" &PTRACERS_PARM01
 PTRACERS_numInUse=2,
 PTRACERS_Iter0=0,
 PTRACERS_advScheme(1:2)=2*77,
 PTRACERS_diffKh(1:2)=2*0.5,
 PTRACERS_diffKr(1:2)=2*1.E-5,
 PTRACERS_useKPP(1:2)=2*.TRUE.,
 PTRACERS_names(1)='mud',
 PTRACERS_names(2)='sand',
 PTRACERS_units(1)='kg/m^3',
 PTRACERS_units(2)='kg/m^3',
 &
""")
    open(os.path.join(out, 'data.sediment'), 'w').write(f"""# LR17: mud and fine sand, layered bed (as the Colville river-tuned set-up),
# IBI hourly waves (gendata.py --sed). Initial bed from EMODnet (LR17_bedmap.py),
# mud per river (Hub'Eau suspended solids, Feb-Mar 2022), rain on.
 &SEDIMENT_PARM01
 sedNcl=2,
 sedPtrIdx=1,2,
 sedWs=1.E-4,1.E-2,
 sedTauCE=0.05,0.15,
 sedTauCD=0.2,0.,
 sedEroM=5.E-4,5.E-3,
 sedRiverC={RIVER_MUD:g},0.,
 sedRiverCFile='river_mud.bin',' ',
 sedBulkRho=500.,1600.,
 sedBedFile='bed_mud.bin','bed_sand.bin',
 &
 &SEDIMENT_PARM02
 sedNlay=8,
 sedCohesive=.TRUE.,.FALSE.,
 sedPorosity=0.9,
 sedFcLo=0.03,
 sedFcHi=0.2,
 sedTcbNew=0.05,
 sedTcbOffset=4.0236,
 sedTcbSlope=2.5,
 sedTcbMin=0.05,
 sedTcbMax=2.,
 sedTcons=2592000.,
 sedTswell=864000.,
 sedZaCoef=0.007,
 sedZaMin=0.5E-3,
 sedLayMass=1.,
 &
 &SEDIMENT_PARM03
 sedWaveMode=3,
 sedWaveHsFile='wav_hs.bin',
 sedWavePerFile='wav_per.bin',
 sedWaveDirFile='wav_dir.bin',
 sedWaveNrec={nw},
 sedWaveDt=3600.,
 sedWaveGamma=0.5,
 sedZ0=2.E-4,
 &
""")
    pk = open(os.path.join(out, 'data.pkg')).read()
    open(os.path.join(out, 'data.pkg'), 'w').write(
        pk.replace(' &\n', ' usePTRACERS=.TRUE.,\n useSEDIMENT=.TRUE.,\n &\n'))
    dg = open(os.path.join(out, 'data.diagnostics')).read()
    fr = re.search(r'frequency\(1\)=(-[\d.]+)', dg).group(1)
    dg = dg.replace("fields(1:4,2)='UVEL    ','VVEL    ','SALT    ','THETA   ',",
                    "fields(1:6,2)='UVEL    ','VVEL    ','SALT    ','THETA   ','TRAC01  ','TRAC02  ',")
    dg = dg.replace(" &\n &DIAG_STATIS_PARMS", f""" fields(1:6,3)='SEDBED01','SEDBED02','SEDTAUB ','SEDTAUW ','SEDUBW  ','SEDHSW  ',
 fileName(3)='sed', frequency(3)={fr}, timePhase(3)=0.,
 &
 &DIAG_STATIS_PARMS""")
    open(os.path.join(out, 'data.diagnostics'), 'w').write(dg)
    print(f'sediment: {nw} hourly IBI wave records; initial bed mud {mud[model].mean():.1f}, '
          f'sand {sand[model].mean():.1f} kg/m2 (model mean); river mud ' + ', '.join(
              f'{k} {v * 1000:g}' for k, v in RIVER_MUD_EACH.items()) + ' mg/l')


if __name__ == '__main__':
    sed = '--sed' in sys.argv
    transport = '--transport' in sys.argv
    stamp_times = '--stamp-times' in sys.argv
    tide = next((a.split('=', 1)[1] for a in sys.argv if a.startswith('--tide=')), None)
    if tide not in (None, 'fes2022', 'fes2014'):
        sys.exit('--tide=fes2022 or --tide=fes2014')
    fes_boundary = next(('eta' if a == '--fes-boundary=eta' else 'all' for a in sys.argv
                         if a in ('--fes-boundary', '--fes-boundary=eta')), False)
    band_tag = next((a.split('=', 1)[1] for a in sys.argv if a.startswith('--band-download=')), None)
    if band_tag and not fes_boundary:
        sys.exit('--band-download needs --fes-boundary')
    if fes_boundary and (tide or stamp_times):
        sys.exit('--fes-boundary does not combine with --tide or --stamp-times')
    args = [a for a in sys.argv[1:] if a not in ('--sed', '--transport', '--stamp-times',
                                                 '--fes-boundary', '--fes-boundary=eta')
            and not a.startswith(('--tide=', '--band-download=', '--zb='))]
    zb_file = next((os.path.abspath(a.split('=', 1)[1]) for a in sys.argv if a.startswith('--zb=')), None)
    start, days = args[0], float(args[1])
    out = args[2] if len(args) > 2 else os.path.join(HERE, f'run_{start}' + ('_sed' if sed else ''))
    main(start, days, out, sed, transport, stamp_times, tide, fes_boundary, band_tag, zb_file)
