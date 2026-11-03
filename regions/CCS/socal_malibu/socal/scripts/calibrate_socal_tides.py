#!/usr/bin/env python3
"""Calibrate the socal boundary tides against FES2022b elevations inside the domain (2026-10-11).

The 15-day tide test (run_v3_ramp15d, FES2014 currents on the normal velocity only) gave tides ~15-20%
too large and ~10-12 deg late at the NOAA gauges. This finds one complex correction per constituent:
- model: hourly ETAN (hours after the 2-day ramp) at interior wet points (depth > 50 m, > 10 cells from
  the open boundaries, every STRIDE-th point);
- reference: FES2022b elevation (nearest FES ocean cell) predicted for the same hours, all 8
  constituents with nodal corrections (UTide FUV, as gen_obcs_tides.py);
- both series fitted with UTide (M2, S2, K1, O1: the pairs a 13-day record separates), so the same
  aliasing (N2 into M2, K2 into S2, P1 into K1, Q1 into O1) acts on both;
- correction = median over points of A_fes/A_model and of g_model - g_fes; applied as am *= ratio,
  ph -= dg/360 * T (ph is the phase lag in seconds, obcs_add_tides: am cos(2 pi (t - ph)/T)).
  N2 takes M2's correction, K2 S2's, P1 K1's, Q1 O1's.
The NOAA gauges are left out of the calibration and stay an independent check.

Usage: calibrate_socal_tides.py <run dir> <fes2022 dir> <tide dir (OB?am/ph_tides.bin)> <out dir> [--skip-hours 48]
"""
import datetime as dt
import glob
import json
import os
import shutil
import sys

import numpy as np
import utide
import xarray as xr
from utide import harmonics, ut_constants

HERE = os.path.dirname(os.path.abspath(__file__))
GRID = os.path.join(HERE, '..', 'work', 'socal', 'socal_ncgrid.nc')
NX, NY, DT, NCOMP = 168, 101, 120.0, 10
T0 = dt.datetime(2025, 1, 1)
CONS8 = ['M2', 'S2', 'N2', 'K2', 'K1', 'O1', 'P1', 'Q1']
FIT = ['M2', 'S2', 'K1', 'O1']
MAP = {'M2': 'M2', 'N2': 'M2', 'S2': 'S2', 'K2': 'S2', 'K1': 'K1', 'P1': 'K1', 'O1': 'O1', 'Q1': 'O1'}
STRIDE = 7


def fes_pred(fdir, pts, t, lat_ref):
    names = [n.strip() for n in ut_constants.const.name]
    lind = np.array([names.index(c) for c in CONS8])
    from utide._time_conversion import _normalize_time          # UTide's internal time (as gen_obcs_tides.py)
    tu = _normalize_time(np.array([(x - T0).total_seconds() / 86400.0 for x in t]), '2025-01-01')
    F, U, V = harmonics.FUV(tu, tu[0], lind, lat_ref, [0, 0, 0, 0])   # exact nodal f, u and astronomical V
    h = np.zeros((len(pts), len(t)))
    for k, c in enumerate(CONS8):
        d = xr.open_dataset(os.path.join(fdir, f'{c.lower()}_fes2022.nc'))
        lat, lon = d.lat.values, d.lon.values % 360
        LA, LO = np.meshgrid(lat, lon, indexing='ij')
        A, G = d.amplitude.values / 100.0, d.phase.values             # cm -> m, deg
        ok = np.isfinite(A)
        for n, (la, lo) in enumerate(pts):
            dd = np.where(ok, (LA - la) ** 2 + ((LO - lo % 360) * np.cos(np.radians(la))) ** 2, np.inf)
            j, i = np.unravel_index(np.argmin(dd), dd.shape)
            h[n] += F[:, k] * A[j, i] * np.cos(2 * np.pi * (V[:, k] + U[:, k]) - np.radians(G[j, i]))
    return h


def main():
    run, fdir, tdir, out = sys.argv[1:5]
    skip = int(sys.argv[sys.argv.index('--skip-hours') + 1]) if '--skip-hours' in sys.argv else 48
    g = xr.open_dataset(GRID)
    xc, yc, dep, wet = g.XC.values, g.YC.values, g.Depth.values, g.HFacC.values[0] > 0
    jj, ii = np.meshgrid(np.arange(NY), np.arange(NX), indexing='ij')
    sel = wet & (dep > 50) & (ii >= 10) & (ii < NX - 10) & (jj >= 10) & (jj < NY - 10)
    pj, pi = np.where(sel)
    pj, pi = pj[::STRIDE], pi[::STRIDE]
    fs = sorted(glob.glob(os.path.join(run, 'diags', 'hourly', 'EtaN.*.data')))[skip:]
    t = [T0 + dt.timedelta(seconds=int(f.split('.')[-2]) * DT - 1800) for f in fs]
    eta = np.stack([np.fromfile(f, '>f4').reshape(NY, NX)[pj, pi] for f in fs], 1)
    pts = list(zip(yc[pj, pi], xc[pj, pi]))
    lat_ref = float(np.mean(yc[pj, pi]))
    ref = fes_pred(fdir, pts, t, lat_ref)
    tt = np.array(t, dtype='datetime64[s]')
    kw = dict(constit=FIT, method='ols', conf_int='none', verbose=False, nodal=True, trend=False)
    R = {c: [] for c in FIT}
    for n in range(len(pts)):
        cm = utide.solve(tt, eta[n] - eta[n].mean(), lat=pts[n][0], **kw)
        cf = utide.solve(tt, ref[n] - ref[n].mean(), lat=pts[n][0], **kw)
        for c in FIT:
            km, kf = list(cm.name).index(c), list(cf.name).index(c)
            R[c].append((cf.A[kf] / cm.A[km], (cm.g[km] - cf.g[kf] + 180) % 360 - 180))
    corr = {}
    for c in FIT:
        a = np.array(R[c])
        corr[c] = dict(amp_ratio=float(np.median(a[:, 0])), amp_ratio_iqr=np.percentile(a[:, 0], [25, 75]).round(3).tolist(),
                       dphase_deg=float(np.median(a[:, 1])), dphase_iqr=np.percentile(a[:, 1], [25, 75]).round(1).tolist())
    print(f'{len(pts)} points, {len(t)} hours ({t[0]} .. {t[-1]})')
    for c, v in corr.items():
        print(f'  {c}: FES/model amplitude {v["amp_ratio"]:.3f} (IQR {v["amp_ratio_iqr"]}), model - FES phase '
              f'{v["dphase_deg"]:+.1f} deg (IQR {v["dphase_iqr"]})')
    # apply to the boundary files
    os.makedirs(out, exist_ok=True)
    snip = open(os.path.join(tdir, 'data.obcs_tides.snippet')).read()
    T = {c: float(l.split('=')[1].strip().rstrip(','))
         for c, l in zip(CONS8, [l for l in snip.splitlines() if l.strip().startswith('tidalPeriod(')])}
    for X, n in (('W', NY), ('E', NY), ('S', NX), ('N', NX)):
        am = np.fromfile(os.path.join(tdir, f'OB{X}am_tides.bin'), '>f8').reshape(NCOMP, n)
        ph = np.fromfile(os.path.join(tdir, f'OB{X}ph_tides.bin'), '>f8').reshape(NCOMP, n)
        for k, c in enumerate(CONS8):
            v = corr[MAP[c]]
            am[k] *= v['amp_ratio']
            ph[k] = np.where(am[k] != 0, np.mod(ph[k] - v['dphase_deg'] / 360.0 * T[c], T[c]), 0.0)
        am.astype('>f8').tofile(os.path.join(out, f'OB{X}am_tides.bin'))
        ph.astype('>f8').tofile(os.path.join(out, f'OB{X}ph_tides.bin'))
    shutil.copy(os.path.join(tdir, 'data.obcs_tides.snippet'), out)
    json.dump(dict(points=len(pts), hours=len(t), first=str(t[0]), last=str(t[-1]), run=run,
                   correction=corr, mapping=MAP), open(os.path.join(out, 'calibration.json'), 'w'), indent=1)
    print('wrote', out)


if __name__ == '__main__':
    main()
