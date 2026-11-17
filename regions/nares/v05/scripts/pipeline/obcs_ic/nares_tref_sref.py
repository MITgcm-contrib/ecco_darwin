#!/usr/bin/env python3
"""tRef/sRef for the Nares run: horizontal mean of the starting T/S (pickup_THETA/SALT from
gen_pickups, 90 levels) over wet cells at each level (pan_AO method). Levels with no wet cell
repeat the deepest wet value. Prints namelist lines for input/data PARM01.
Usage: nares_tref_sref.py <forcings/pickups dir> <Nares_ncgrid.nc> [iter]"""
import sys, os, numpy as np, xarray as xr
d, g = sys.argv[1], sys.argv[2]; it = sys.argv[3] if len(sys.argv) > 3 else '0000604872'
NX, NY, NR = 1120, 360, 90
h = xr.open_dataset(g)['HFacC'].values               # (90, 360, 1120)
out = {}
for v in ('THETA', 'SALT'):
    a = np.fromfile(os.path.join(d, f'pickup_{v}.{it}.data'), '>f8').reshape(NR, NY, NX)
    prof = np.array([a[k][h[k] > 0].mean() if (h[k] > 0).any() else np.nan for k in range(NR)])
    last = np.where(np.isfinite(prof))[0].max(); prof[last + 1:] = prof[last]
    out[v] = prof
for nm, v in (('tRef', 'THETA'), ('sRef', 'SALT')):
    vals = ', '.join(f'{x:.4f}' for x in out[v])
    print(f' {nm} = {vals},')
np.savez(os.path.join(d, 'tref_sref.npz'), tRef=out['THETA'], sRef=out['SALT'])
