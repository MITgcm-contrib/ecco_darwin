#!/usr/bin/env python3
"""gen_nares_mitgrid.py -- rotated lat-lon grid for the Nares Strait downscale.

The rotated equator is the great circle through two points on the strait axis
(Lincoln Sea -> North Water), so the strait runs along the model x-axis and the
grid is near-uniform (cos of the rotated latitude stays > 0.998 for +-200 km).
All metrics are analytic on the sphere in rotated coordinates:
    dx = R cos(phi') dlam',  dy = R dphi',  rA = R^2 dlam' (sin phi'_N - sin phi'_S)
Written as ONE facet .mitgrid (16 fields, (ny+1, nx+1), real*8 big-endian), the
layout exch2 reads as tile001.mitgrid for a single facet (pan_AO, 2026-10-02),
so there are no per-tile files and no tile-numbering trap.

Usage: python3 gen_nares_mitgrid.py <outdir> [dx_km]
"""
import os
import sys
import json
import numpy as np

R = 6371.0e3
# strait axis end points (lat, lon) -- north end in the Lincoln Sea beyond the
# Robeson Channel arch, south end in the North Water beyond the Smith Sound arch
A = (83.6, -55.0)
B = (76.0, -73.5)
EXT_N_KM, EXT_S_KM = 60.0, 120.0   # extra along-axis length past A and B (S: +120 km to Devon I. for the whole North Water, Clement 2026-10-06)
HALF_W_KM = 180.0                   # cross-axis half width (ny = 360 at 1 km)
NX_MULT = 40                        # round nx up to a multiple of this (tiling), padding both ends


def unit(lat, lon):
    la, lo = np.deg2rad(lat), np.deg2rad(lon)
    return np.array([np.cos(la) * np.cos(lo), np.cos(la) * np.sin(lo), np.sin(la)])


def frame():
    pa, pb = unit(*A), unit(*B)
    n = np.cross(pb, pa); n /= np.linalg.norm(n)     # rotated pole (x runs B -> A)
    e1 = pa + pb; e1 /= np.linalg.norm(e1)            # rotated (0, 0) at the midpoint
    e2 = np.cross(n, e1)
    return e1, e2, n


def to_geo(lam, phi, e1, e2, n):
    p = (np.cos(phi)[..., None] * (np.cos(lam)[..., None] * e1 + np.sin(lam)[..., None] * e2)
         + np.sin(phi)[..., None] * n)
    lat = np.rad2deg(np.arcsin(np.clip(p[..., 2], -1, 1)))
    lon = np.rad2deg(np.arctan2(p[..., 1], p[..., 0]))
    return lon, lat


def to_rot(lon, lat, e1, e2, n):
    p = unit(lat, lon)
    p = np.moveaxis(p, 0, -1)
    return np.arctan2(p @ e2, p @ e1), np.arcsin(np.clip(p @ n, -1, 1))


def main(outdir, dx_km=1.0):
    e1, e2, n = frame()
    la, _ = to_rot(np.array(A[1]), np.array(A[0]), e1, e2, n)
    lb, _ = to_rot(np.array(B[1]), np.array(B[0]), e1, e2, n)
    lam0 = float(lb) - EXT_S_KM * 1e3 / R
    lam1 = float(la) + EXT_N_KM * 1e3 / R
    d = dx_km * 1e3 / R
    nx = int(np.ceil((lam1 - lam0) / d / NX_MULT) * NX_MULT)
    pad_l = 0.5 * (nx * d - (lam1 - lam0))
    lam0 -= pad_l
    ny = int(round(2 * HALF_W_KM / dx_km))
    lamG = lam0 + d * np.arange(nx + 1)
    phiG = -0.5 * ny * d + d * np.arange(ny + 1)
    lamC = 0.5 * (lamG[:-1] + lamG[1:]); phiC = 0.5 * (phiG[:-1] + phiG[1:])

    def pad(a):  # (ny, nx) -> (ny+1, nx+1) with zero last row/col, mitgrid layout
        o = np.zeros((ny + 1, nx + 1)); o[:ny, :nx] = a; return o

    LC, PC = np.meshgrid(lamC, phiC)
    LG, PG = np.meshgrid(lamG, phiG)
    xC, yC = to_geo(LC, PC, e1, e2, n)
    xG, yG = to_geo(LG, PG, e1, e2, n)
    dxF = R * np.cos(PC) * d
    dyF = np.full_like(PC, R * d)
    rA = R ** 2 * d * (np.sin(PC + d / 2) - np.sin(PC - d / 2))
    # u-points (west faces, lamG, phiC) and v-points (south faces, lamC, phiG)
    LU, PU = np.meshgrid(lamG, phiC)
    LV, PV = np.meshgrid(lamC, phiG)
    dxC = R * np.cos(PU) * d               # between centres, along x, at u-points
    dyC = np.full_like(PV, R * d)          # between centres, along y, at v-points
    dxG = R * np.cos(PV) * d               # cell south-face length
    dyG = np.full_like(PU, R * d)          # cell west-face length
    dxV = R * np.cos(PG) * d               # at corners
    dyU = np.full_like(PG, R * d)
    rAw = R ** 2 * d * (np.sin(PU + d / 2) - np.sin(PU - d / 2))
    rAs = R ** 2 * d * (np.sin(PV + d / 2) - np.sin(PV - d / 2))
    rAz = R ** 2 * d * (np.sin(PG + d / 2) - np.sin(PG - d / 2))

    def fit(a):  # u-/v-/corner-point arrays to (ny+1, nx+1)
        o = np.zeros((ny + 1, nx + 1)); o[:a.shape[0], :a.shape[1]] = a; return o

    fields = [pad(xC), pad(yC), pad(dxF), pad(dyF), pad(rA), xG, yG,
              fit(dxV), fit(dyU), fit(rAz), fit(dxC), fit(dyC), fit(rAw), fit(rAs), fit(dxG), fit(dyG)]
    os.makedirs(outdir, exist_ok=True)
    np.stack(fields).astype('>f8').tofile(os.path.join(outdir, 'Nares.mitgrid'))
    # grid angle (rotated x-axis vs local east) at centres, for vector rotation checks
    lonE, latE = to_geo(LC + 1e-6, PC, e1, e2, n)
    ang = np.arctan2(np.deg2rad(latE - yC), np.deg2rad(lonE - xC) * np.cos(np.deg2rad(yC)))
    meta = dict(nx=nx, ny=ny, dx_km=dx_km, A=A, B=B, half_width_km=HALF_W_KM,
                ext_km=[EXT_S_KM, EXT_N_KM], pole=n.tolist(), e1=e1.tolist(), e2=e2.tolist(),
                lon_range=[float(xG.min()), float(xG.max())], lat_range=[float(yG.min()), float(yG.max())],
                dx_range_m=[float(dxF.min()), float(dxF.max())],
                angle_deg_range=[float(np.rad2deg(ang.min())), float(np.rad2deg(ang.max()))])
    json.dump(meta, open(os.path.join(outdir, 'Nares_grid_meta.json'), 'w'), indent=1)
    print(json.dumps(meta, indent=1))


if __name__ == '__main__':
    main(sys.argv[1], *(float(v) for v in sys.argv[2:3]))
