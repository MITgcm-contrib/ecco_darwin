#!/usr/bin/env python3
"""Tide predictions from the FES atlases for LR17 (module).

atlas(source, kind) reads the cropped harmonic constants (scripts/fes_crop.py)
  source 'FES2022b' (elevations, 1/30 deg) or 'FES2014' (currents, 1/16 deg;
  elevations once data/LR17/fes2014_elevations is unpacked);
  kind 'eta' (m) or 'u', 'v' (eastward, northward depth-mean current, m/s).
constants(A, lat, lon) interpolates amplitude and phase to points (bilinear on
  the complex amplitude; atlas land cells are filled from the nearest sea cell).
predict(times, amp, pha, names, lat) predicts the tide at the points with
  UTide's reconstruction (nodal corrections, Greenwich phase lags), so atlas
  and harmonic analyses (LR17_harmonics.py) share one convention.

Constituent names follow UTide (M2, 2N2, MU2, ...); FES file names are lower
case (m2, 2n2, mu2, ...).
"""
import glob
import os

import numpy as np
import xarray as xr
import utide
from scipy import ndimage
from scipy.interpolate import RegularGridInterpolator

HERE = os.path.dirname(os.path.abspath(__file__))
D = os.path.join(os.path.dirname(HERE), 'data', 'LR17')
CONS = ['M2', 'S2', 'N2', 'K2', '2N2', 'MU2', 'NU2', 'L2', 'T2', 'K1', 'O1', 'P1', 'Q1',
        'M3', 'M4', 'MS4', 'MN4', 'N4', 'S4', 'M6']


def _file(source, kind, c):
    c = c.lower()
    if source == 'FES2022b':
        return os.path.join(D, 'fes2022b', f'{c}_fes2022_lr17.nc')
    sub = {'u': 'fes2014_currents/eastward_velocity', 'v': 'fes2014_currents/northward_velocity',
           'eta': 'fes2014_elevations'}[kind]
    f = glob.glob(os.path.join(D, sub, '**', f'{c}_lr17.nc'), recursive=True)
    return f[0] if f else os.path.join(D, sub, f'{c}_lr17.nc')


def atlas(source, kind, names=CONS):
    """{name: (lat, lon_east, complex amplitude array)} for the constituents
    present; amplitude in m (elevation) or m/s (current)."""
    out = {}
    for c in names:
        f = _file(source, kind, c)
        if not os.path.exists(f):
            continue
        d = xr.open_dataset(f)
        if kind == 'eta':
            a, p = d['amplitude'], d['phase']
        else:
            a = d['Ua' if kind == 'u' else 'Va']
            p = d['Ug' if kind == 'u' else 'Vg']
        z = a.values * 0.01 * np.exp(1j * np.radians(p.values))      # cm -> m
        bad = ~np.isfinite(z)
        if bad.any():                                   # fill land from the nearest sea cell
            idx = ndimage.distance_transform_edt(bad, return_distances=False, return_indices=True)
            z = z[tuple(idx)]
        out[c] = (d.lat.values, d.lon.values, z)
    return out


def constants(A, lat, lon):
    """amplitude (npt, ncons), phase (deg), names at points lat, lon (deg E,
    negative west)"""
    lat, lon = np.atleast_1d(lat), np.atleast_1d(lon) % 360
    names = list(A)
    amp = np.zeros((lat.size, len(names)))
    pha = np.zeros_like(amp)
    for k, c in enumerate(names):
        la, lo, z = A[c]
        f = lambda w: RegularGridInterpolator((la, lo), w, bounds_error=False,
                                              fill_value=None)(np.c_[lat, lon])
        zi = f(z.real) + 1j * f(z.imag)
        amp[:, k] = np.abs(zi)
        pha[:, k] = np.degrees(np.angle(zi)) % 360
    return amp, pha, names


_TEMPLATE = {}


def _template(t, names, lat):
    key = (len(t), float(t[0]), float(t[-1]), tuple(names), round(float(lat), 1))
    if key not in _TEMPLATE:
        s = sum(np.cos(2 * np.pi * t * k / (0.5 + k)) for k in range(1, 4))   # any series
        _TEMPLATE[key] = utide.solve(t, s, lat=lat, constit=names, method='ols',
                                     conf_int='none', nodal=True, trend=False, verbose=False,
                                     epoch='2022-01-01')
    return _TEMPLATE[key]


def predict(times, amp, pha, names, lat):
    """tide (npt, ntimes) at datetime64 times from amplitude/phase per point"""
    t = (np.asarray(times) - np.datetime64('2022-01-01')) / np.timedelta64(1, 'D')
    coef = _template(t, names, lat)
    order = [names.index(n) for n in coef.name]
    out = np.zeros((amp.shape[0], t.size))
    for p in range(amp.shape[0]):
        coef['A'] = amp[p, order]
        coef['g'] = pha[p, order]
        coef['mean'] = 0.0
        coef['slope'] = 0.0
        out[p] = utide.reconstruct(t, coef, verbose=False, min_SNR=0, min_PE=0,
                                   epoch='2022-01-01').h
    return out


def predict_fast(times, amp, pha, names, lat):
    """as predict, for many points at once: the nodal-corrected basis E
    (ntimes, ncons) of UTide's reconstruction is built once and applied to
    the complex amplitudes of all points; returns (npt, ntimes)"""
    from utide._reconstruct import ut_E
    t = (np.atleast_1d(np.asarray(times)) - np.datetime64('2022-01-01')) / np.timedelta64(1, 'D')
    tt = t if t.size > 1 else np.r_[t, t + 1.0]
    coef = _template(np.arange(tt[0] - 8, tt[0] + 8, 1 / 24.), names, lat)
    aux = coef['aux']
    opt = aux['opt']
    order = [names.index(n) for n in coef['name']]
    from utide._time_conversion import _normalize_time
    E = ut_E(_normalize_time(t, '2022-01-01'), aux['reftime'], aux['frq'], aux['lind'], aux['lat'],
             [opt['nodsatlint'], opt['nodsatnone'], opt['gwchlint'], opt['gwchnone']],
             opt['prefilt'])
    ap = 0.5 * amp[:, order].T * np.exp(-1j * np.radians(pha[:, order].T))      # (ncons, npt)
    return (2.0 * np.real(E @ ap)).T
