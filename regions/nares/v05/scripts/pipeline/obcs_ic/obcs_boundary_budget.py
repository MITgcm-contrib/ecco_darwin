#!/usr/bin/env python3
"""Volume budget of the Nares open boundaries, split by end of the strait.

Normal transport (Sv, positive = INTO the domain) through each open boundary from an OBCS directory,
time-mean over the records, with the north (j=Ny) and south (j=1) edges split into their Baffin-end
(i < ISPLIT) and Lincoln-end (i >= ISPLIT) segments. The implied along-strait throughflow is the net
inflow on the Lincoln end (east + Lincoln-end N/S segments), positive = southward (Lincoln -> Baffin),
before useOBCSbalance removes any global imbalance (the imbalance is printed too).
OBCS conventions: east u at face i=Nx-1 (0-based), west at i=1, north v at face j=Ny-1, south at j=1.
Usage: obcs_boundary_budget.py <nares dir> <OBCS subdir, e.g. OBCS_rot> [ISPLIT=560]
"""
import sys
import numpy as np
import xarray as xr

N, sub = sys.argv[1], sys.argv[2]
ISPLIT = int(sys.argv[3]) if len(sys.argv) > 3 else 560
d = f'{N}/forcings/{sub}'
g = xr.open_dataset(f'{N}/Nares_ncgrid.nc')
hW, hS, dyG, dxG, drF = g['HFacW'].values, g['HFacS'].values, g['dyG'].values, g['dxG'].values, g['drF'].values
NR = drF.size
NY, NX = g['XC'].shape


def tmean(f, n):
    a = np.memmap(f, '>f4', 'r').reshape(-1, NR, n)
    s = np.zeros((NR, n))
    for t in range(a.shape[0]):
        s += a[t]
    return s / a.shape[0], a.shape[0]


out = {}
uE, nt = tmean(f'{d}/UVEL_east.bin', NY); iE = min(NX - 1, hW.shape[-1] - 1)
out['east'] = -(uE * hW[:, :, iE] * dyG[:, iE][None] * drF[:, None]).sum() / 1e6        # +u leaves
uW, _ = tmean(f'{d}/UVEL_west.bin', NY)
out['west'] = (uW * hW[:, :, 1] * dyG[:, 1][None] * drF[:, None]).sum() / 1e6
vN, _ = tmean(f'{d}/VVEL_north.bin', NX); jN = min(NY - 1, hS.shape[1] - 1)
aN = -(vN * hS[:, jN, :] * dxG[jN, :][None] * drF[:, None])                             # +v leaves
vS, _ = tmean(f'{d}/VVEL_south.bin', NX)
aS = (vS * hS[:, 1, :] * dxG[1, :][None] * drF[:, None])
for nm, a in (('north', aN), ('south', aS)):
    out[f'{nm}_baffin'] = a[:, :ISPLIT].sum() / 1e6
    out[f'{nm}_lincoln'] = a[:, ISPLIT:].sum() / 1e6
print(f'{sub}: time-mean over {nt} records; + = into the domain (Sv)')
for k, v in out.items():
    print(f'  {k:15s} {v:+.3f}')
lin = out['east'] + out['north_lincoln'] + out['south_lincoln']
baf = out['west'] + out['north_baffin'] + out['south_baffin']
print(f'  net Lincoln end {lin:+.3f}   net Baffin end {baf:+.3f}   global imbalance {lin + baf:+.3f} (removed by useOBCSbalance)')
print(f'  implied southward throughflow ~ {(lin - baf) / 2:+.3f} Sv (balance-split estimate)')
