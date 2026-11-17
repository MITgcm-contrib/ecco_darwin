#!/usr/bin/env python3
"""OBCS index arrays for the Nares regional run (1120 x 360 rectangle, i runs S->N
along the strait). An edge cell is an open boundary only where it is wet
(bathy < 0 and shelficeTopo > bathy); elsewhere 0 = OB_indexUnset (obcs_readparms.F).
Model edges: W (i=1) = North Water / N Baffin, E (i=1120) = Lincoln Sea,
S (j=1) = Greenland side, N (j=360) = Ellesmere side.
Usage: python3 gen_ob_index.py <grid dir> <out namelist fragment>"""
import sys, numpy as np
g, out = sys.argv[1], sys.argv[2]
NX, NY = 1120, 360
b = np.fromfile(f'{g}/Nares_bathymetry.bin', '>f8').reshape(NY, NX)
s = np.fromfile(f'{g}/Nares_shelficeTopo.bin', '>f8').reshape(NY, NX)
wet = (b < 0) & (s > b)
arr = {'OB_Iwest':  np.where(wet[:, 0], 1, 0),      # one value per row j
       'OB_Ieast':  np.where(wet[:, -1], NX, 0),
       'OB_Jsouth': np.where(wet[0, :], 1, 0),      # one value per column i
       'OB_Jnorth': np.where(wet[-1, :], NY, 0)}
def fmt(a, per=15):
    return '\n'.join('   ' + ', '.join(f'{v:4d}' for v in a[k:k+per]) + ',' for k in range(0, a.size, per))
with open(out, 'w') as f:
    for k, a in arr.items():
        f.write(f' {k} =\n{fmt(a)}\n')
for k, a in arr.items():
    print(f'{k}: {a.size} entries, {np.count_nonzero(a)} set')
corners = [(0, 0), (0, NX-1), (NY-1, 0), (NY-1, NX-1)]
print('wet corners (on two boundaries):', [c for c in corners if wet[c]])
