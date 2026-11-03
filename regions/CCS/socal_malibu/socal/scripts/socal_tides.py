#!/usr/bin/env python3
"""Barotropic tidal OBCS for the socal 2-km nest from FES2014 currents (FES2022b has no currents).

1. Boundary points file for downscaling/python/tides/gen_obcs_tides.py from work/socal/socal_ncgrid.nc:
   west/east = columns 0 / nx-1 (Ny entries), south/north = rows 0 / ny-1 (Nx entries),
   lat/lon of the OB cell centers, AngleCS = 1, AngleSN = 0 (lon/lat grid), has = wet surface cell.
2. gen_obcs_tides.py -> OB{W,E,S,N}_{u,v}Tid{Am,Ph}.bin (pan_AO interface, records = constituents).
3. Convert to the OLDER pkg/obcs interface of the socal build (v05_V5r1 darwin3, Feb 2022):
   OB{N,S,E,W}amFile / OB{N,S,E,W}phFile, read by obcs_init_variables.F with READ_REC_{XZ,YZ}_RL
   as tidalComponents = 10 records of the boundary length; obcs_add_tides.F adds
   am*cos(2 pi (myTime - ph)/tidalPeriod) to the NORMAL velocity only (u at W/E, v at S/N), the same
   phase convention as gen_obcs_tides.py (ph = phase lag in seconds, phase zero at startDate).
   Unused records 9-10 are zero (tidalPeriod 0 there, so obcs_add_tides skips them).

Usage: python3 socal_tides.py <fes u dir> <fes v dir> <out dir> [start YYYY-MM-DD] [nodal date]
"""
import os
import subprocess
import sys

import numpy as np
import xarray as xr

HERE = os.path.dirname(os.path.abspath(__file__))
GRID = os.path.join(HERE, '..', 'work', 'socal', 'socal_ncgrid.nc')
GEN = os.path.expanduser('~/Documents/research/downscaling/python/tides/gen_obcs_tides.py')
NCOMP = 10                      # OBCS_PARAMS.h tidalComponents
CONS = ['M2', 'S2', 'N2', 'K2', 'K1', 'O1', 'P1', 'Q1']


def main():
    udir, vdir, out = sys.argv[1:4]
    start = sys.argv[4] if len(sys.argv) > 4 else '2025-01-01'
    nodal = sys.argv[5] if len(sys.argv) > 5 else '2025-04-01'     # midpoint of the 181-d run
    os.makedirs(out, exist_ok=True)
    g = xr.open_dataset(GRID)
    xc, yc, wet = g.XC.values, g.YC.values, g.HFacC.values[0] > 0
    ny, nx = xc.shape
    sl = {'west': (slice(None), 0), 'east': (slice(None), nx - 1), 'south': (0, slice(None)), 'north': (ny - 1, slice(None))}
    P = {}
    for b, s in sl.items():
        P[f'{b}_lat'], P[f'{b}_lon'], P[f'{b}_has'] = yc[s], xc[s] % 360.0, wet[s]
        P[f'{b}_cs'], P[f'{b}_sn'] = np.ones_like(yc[s]), np.zeros_like(yc[s])
        print(f'{b}: {wet[s].sum()} wet of {wet[s].size}')
    pts = os.path.join(out, 'ob_points_socal.npz')
    np.savez(pts, **P)
    subprocess.run([sys.executable, GEN, '--points', pts, '--u-dir', udir, '--v-dir', vdir, '--start', start,
                    '--nodal-date', nodal, '--out-dir', os.path.join(out, 'gen'), '--prec', '64'], check=True)
    per = [l for l in open(os.path.join(out, 'gen', 'data.obcs_tides.snippet')) if 'OBCS_tidalPeriod' in l][0]
    T = [float(x) for x in per.split('=')[1].strip().rstrip(',').split(',')]
    snip = [' useOBCStides = .TRUE.,', '# constituents %s; FES2014 currents; phase zero = startDate %s; nodal f,u at %s;'
            % (','.join(CONS), start, nodal), '# socal_tides.py (older pkg/obcs am/ph interface, normal velocity only)']
    snip += [f' tidalPeriod({k + 1}) = {t:.6f},' for k, t in enumerate(T)]
    for b, X, comp in (('west', 'W', 'u'), ('east', 'E', 'u'), ('south', 'S', 'v'), ('north', 'N', 'v')):
        n = P[f'{b}_lat'].size
        for kind, key in (('am', 'Am'), ('ph', 'Ph')):
            a = np.fromfile(os.path.join(out, 'gen', f'OB{X}_{comp}Tid{key}.bin'), '>f8').reshape(len(CONS), n)
            full = np.zeros((NCOMP, n))
            full[:len(CONS)] = a
            fn = f'OB{X}{kind}_tides.bin'
            full.astype('>f8').tofile(os.path.join(out, fn))
            snip.append(f" OB{X}{kind}File = '{fn}',")
        am = np.fromfile(os.path.join(out, f'OB{X}am_tides.bin'), '>f8').reshape(NCOMP, n)
        print(f'{b}: M2 normal-velocity amplitude median {np.median(am[0][P[f"{b}_has"]]):.4f}, '
              f'max {am[0].max():.4f} m/s; K1 median {np.median(am[4][P[f"{b}_has"]]):.4f} m/s')
    open(os.path.join(out, 'data.obcs_tides.snippet'), 'w').write('\n'.join(snip) + '\n')
    print('wrote', os.path.join(out, 'data.obcs_tides.snippet'))


if __name__ == '__main__':
    main()
