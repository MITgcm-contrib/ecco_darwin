#!/usr/bin/env python3
"""Inputs for wad_estuary_3d: a 10 m deep channel along x (8 km, 40
columns of 200 m; |y - 3 km| < 600 m) shoaling to 3 m at its head, with
tidal flats on both sides rising to +1.3 m at the side walls. The first
1 km (next to the western open boundary) has no flats, so that the
Flather boundary (code/obcs_calc.F, 1.5 m M2 tide) is in deep water.
Salinity stratified 30 (z=0) to 34 (z=-10 m); f = 1e-4 1/s; KPP and
GM/Redi on. Model rest level DATUM = 3 m above sea level; initial state
at high water, at rest.
Also writes the input.ptr tracers (1: uniform; 2: dye over 3 < x < 5 km)
and the input.seiche initial surface (high water tilted by 1e-4 across the
estuary: a cross-channel seiche over the drying flats).

Usage:  python3 gendata.py [OUTDIR]
With no OUTDIR the files go to input/, input.ptr/ and input.seiche/ next
to this script. With OUTDIR they all go there. The grid can be refined for
runs outside testreport with WAD_NX (default 40), WAD_NY (30), WAD_NR (5):
the domain stays 8 x 6 km and 13 m deep (dx = 8000/NX, dz = 13/NR).
"""
import os
import sys

import numpy as np

NX = int(os.environ.get('WAD_NX', 40))
NY = int(os.environ.get('WAD_NY', 30))
NR = int(os.environ.get('WAD_NR', 5))
DX = 8000.0 / NX
DATUM = 3.0
x = (np.arange(NX) - 0.5) * DX            # i=0 is the open-boundary column
y = (np.arange(NY) + 0.5) * DX
X, Y = np.meshgrid(x, y)                  # shape (NY, NX)
zchan = -10.0 + 7.0 * np.clip(X, 0, None) / 8000.0
dy = np.abs(Y - 3000.0)
zflat = np.where(dy < 600.0, -99.0, -3.0 + 4.5 * (dy - 600.0) / 2400.0)
zb = np.maximum(zchan, zflat)
# open boundary in deep water: no flats in the first 1 km (estuary mouth)
zb = np.where(X < 1000.0, np.minimum(zchan, -10.0 + 7.0 * 1000.0 / 8000.0), zb)
ocean = np.ones((NY, NX), bool); ocean[:, -1] = False; ocean[-1, :] = False
r_low = np.where(ocean, zb - DATUM, 0.0)
eta = np.where(ocean, np.maximum(1.5 - DATUM, r_low + 0.05), 0.0)   # start at high water
# salinity: 30 at z=0 .. 34 at z=-10 (model z), per r* layer centre
DZ = 13.0 / NR
s = np.zeros((NR, NY, NX))
for k in range(NR):
    rc = -(k + 0.5) * DZ
    with np.errstate(invalid='ignore', divide='ignore'):
        z = eta + rc * (eta - r_low) / np.where(ocean, -r_low, 1.0)
    s[k] = np.where(ocean, 30.0 + 0.4 * np.clip(-(z + DATUM), 0, 10), 0.0)
oc3 = ocean[None].repeat(NR, 0)
ptr1 = np.where(oc3, 1.0, 0.0)
ptr2 = np.where(oc3 & ((X > 3000) & (X < 5000))[None], 1.0, 0.0)
# input.seiche: a cross-estuary seiche, the high-water surface tilted by
# 1e-4 across the estuary (+-0.29 m at the side walls; not at the open
# boundary column), released at once
tilt = 1e-4 * (Y - 3000.0)
tilt[:, 0] = 0.0
eta_seiche = np.where(ocean, np.maximum(eta + tilt, r_low + 0.05), 0.0)

if __name__ == '__main__':
    here = os.path.dirname(os.path.abspath(__file__))
    top = os.path.dirname(here)
    out = sys.argv[1] if len(sys.argv) > 1 else None

    def w(name, a, sub='input'):
        d = out or os.path.join(top, sub)
        a.astype('>f8').tofile(os.path.join(d, name))

    w('bathy.bin', r_low)
    w('eta.bin', eta)
    w('S.bin', s)
    w('ptr1.bin', ptr1, 'input.ptr')
    w('ptr2.bin', ptr2, 'input.ptr')
    w('eta_seiche.bin', eta_seiche, 'input.seiche')
    print(f'{NX}x{NY}x{NR}, dx={DX:g} m: bed range', zb[ocean].min(),
          zb[ocean].max(), ' dry at rest:',
          int(((eta - r_low) < 0.06)[ocean].sum()))
