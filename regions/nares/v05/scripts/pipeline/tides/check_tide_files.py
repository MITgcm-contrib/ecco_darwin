#!/usr/bin/env python3
"""Volume-budget and phase check of generated Nares tidal OBCS files (tides_plan.md step 2 post-check).

Reads OB{W,E,S,N}_{u,v}Tid{Am,Ph}.bin back exactly as obcs_init_fixed.F/obcs_add_tides.F use them
(real*8 big-endian, records = constituents, record length = boundary length; signal
u(t) = Am cos(w t - 2 pi Ph/T), t = myTime), forms the net inward tidal transport through the four
open edges with the model face thickness (sum_k hFac*drF at the u/v face obcs_add_tides.F writes:
west u at i=2, east u at i=1120, south v at j=2, north v at j=360, 1-based) times the face length,
and compares it with dV/dt of the FES2022 elevation tide over all wet cells, nodal-corrected to the
same model-time phase (Am = f A, phase = G - V0 - u). Also prints per-edge M2 maxima.

Pass criteria (plan): net/prism amplitude ratio within ~0.85-1.15 and phase difference < ~10 deg
for M2 and K1 (the pre-check from the FES atlases gave 1.10 / 3 deg and 1.10 / 9 deg).

Usage: check_tide_files.py <tide dir> <nares_ob_points.npz> <Nares_ncgrid.nc> <grid dir>
                           <fes2022 dir> [--start 2015-01-01 --nodal-date 2015-01-16]
Needs downscaling/python/tides on PYTHONPATH (imports gen_obcs_tides.nodal for f, u, V0).
"""
import argparse
import os
import re
import sys

import numpy as np
import xarray as xr
from scipy import ndimage
from scipy.interpolate import RegularGridInterpolator

NX, NY = 1120, 360


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('tide_dir'); ap.add_argument('points'); ap.add_argument('ncgrid')
    ap.add_argument('grid_dir'); ap.add_argument('fes2022')
    ap.add_argument('--start', default='2015-01-01'); ap.add_argument('--nodal-date', default='2015-01-16')
    a = ap.parse_args()
    sys.path.insert(0, os.path.expanduser('~/Documents/research/downscaling/python/tides'))
    from gen_obcs_tides import nodal

    snip = open(os.path.join(a.tide_dir, 'data.obcs_tides.snippet')).read()
    T = np.array([float(x) for x in re.search(r'OBCS_tidalPeriod *= *([^\n]*)', snip).group(1).split(',') if x.strip()])
    cons = re.search(r'# constituents ([A-Z0-9,]+);', snip).group(1).split(',')
    nc = len(cons); assert nc == T.size
    P = np.load(a.points)
    g = xr.open_dataset(a.ncgrid)
    drF = g.drF.values
    # face thickness and length at the faces obcs_add_tides.F writes (0-based indices)
    hW = lambda i: (g.HFacW[:, :, i].values * drF[:, None]).sum(0)          # (Ny,)
    hS = lambda j: (g.HFacS[:, j, :].values * drF[:, None]).sum(0)          # (Nx,)
    face = {'west': ('u', hW(1) * g.dyG.values[:, 1], +1),
            'east': ('u', hW(NX - 1) * g.dyG.values[:, NX - 1], -1),
            'south': ('v', hS(1) * g.dxG.values[1, :], +1),
            'north': ('v', hS(NY - 1) * g.dxG.values[NY - 1, :], -1)}
    X = {'west': 'W', 'east': 'E', 'south': 'S', 'north': 'N'}

    def rd(b, comp, kind, n):
        x = np.fromfile(os.path.join(a.tide_dir, f'OB{X[b]}_{comp}Tid{kind}.bin'), '>f8')
        assert x.size == nc * n, (b, comp, kind, x.size, nc * n)
        return x.reshape(nc, n)

    Q = np.zeros(nc, complex)
    for b, (comp, area, sgn) in face.items():
        n = P[f'{b}_lat'].size
        Am, Ph = rd(b, comp, 'Am', n), rd(b, comp, 'Ph', n)
        # zero where no OB (b_has False), as the model ignores those entries
        assert np.all(Am[:, ~P[f'{b}_has']] == 0)
        Z = Am * np.exp(-2j * np.pi * Ph / T[:, None])     # u(t) = Re(Z e^{i w t})
        Qb = sgn * (Z * area[None, :]).sum(1)
        Q += Qb
        k = cons.index('M2') if 'M2' in cons else 0
        print(f'{b:5s}: {int(P[f"{b}_has"].sum()):4d} OB | {comp}-normal {cons[k]} max {Am[k].max():.3f} m/s, '
              f'tangential max {rd(b, "v" if comp == "u" else "u", "Am", n)[k].max():.3f} m/s | '
              f'{cons[k]} transport amp {abs(Qb[k]) / 1e6:.2f} Sv')

    # FES2022 elevation prism in the same model-time phase convention
    bth = np.fromfile(f'{a.grid_dir}/Nares_bathymetry.bin', '>f8').reshape(NY, NX)
    shi = np.fromfile(f'{a.grid_dir}/Nares_shelficeTopo.bin', '>f8').reshape(NY, NX)
    wet = (bth < 0) & (shi > bth)
    yc, xc, rA = g.YC.values, g.XC.values % 360, g.rA.values
    lat_ref = float(np.mean(np.concatenate([P[f'{b}_lat'][P[f'{b}_has']] for b in face])))
    Tn, f, u, V0, _ = nodal(cons, a.start, a.nodal_date, lat_ref)
    assert np.allclose(Tn, T, rtol=1e-9), (Tn, T)
    print(f'\n{"con":4s} {"net in (Sv)":>12s} {"prism (Sv)":>11s} {"ratio":>6s} {"phase diff (deg)":>17s}')
    for k, c in enumerate(cons):
        d = xr.open_dataset(os.path.join(a.fes2022, f'{c.lower()}_fes2022.nc'))
        la, lo = d.lat.values, d.lon.values % 360
        sl, ol = (la >= 73) & (la <= 86), (lo >= 270) & (lo <= 325)
        A = d['amplitude'].values[np.ix_(sl, ol)].astype('f8') * 0.01
        G = d['phase'].values[np.ix_(sl, ol)].astype('f8')
        z = A * np.exp(1j * np.radians(G))
        bad = ~np.isfinite(z) | (A > 100)
        z = z[tuple(ndimage.distance_transform_edt(bad, return_distances=False, return_indices=True))]
        pts = np.c_[yc[wet], xc[wet]]
        zi = (RegularGridInterpolator((la[sl], lo[ol]), z.real)(pts)
              + 1j * RegularGridInterpolator((la[sl], lo[ol]), z.imag)(pts))
        # eta(t) = f A cos(w t + V0 + u - G) = Re(conj(Ze) e^{i w t}), Ze = f A exp(i(G - V0 - u))
        Ze = f[k] * zi * np.exp(-2j * np.pi * (V0[k] + u[k]))
        dVdt = 1j * (2 * np.pi / T[k]) * np.conj((Ze * rA[wet]).sum())
        r = abs(Q[k]) / abs(dVdt)
        dph = np.degrees(np.angle(Q[k] / dVdt))
        print(f'{c:4s} {abs(Q[k]) / 1e6:12.2f} {abs(dVdt) / 1e6:11.2f} {r:6.2f} {dph:17.1f}')


if __name__ == '__main__':
    main()
