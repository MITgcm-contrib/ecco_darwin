#!/usr/bin/env python3
"""Split gen_pickups' pickup_seaice.<iter>.data (6 x 360 x 1120 float64: siTICE, siAREA, siHEFF,
siHSNOW, siUICE, siVICE) into the single-field files data.seaice reads (AreaFile ... vIceFile).
Without these the model cold-starts with NO ice (pan_AO's 16 Sv day-1 freeze-up sink).
Usage: nares_seaice_ic.py <forcings/pickups dir> [iter]"""
import sys, os, numpy as np
d = sys.argv[1]; it = sys.argv[2] if len(sys.argv) > 2 else '0000604872'
NX, NY = 1120, 360
FLDS = ['siTICE', 'siAREA', 'siHEFF', 'siHSNOW', 'siUICE', 'siVICE']
RANGE = {'siAREA': (0, 1.0001), 'siHEFF': (0, 30), 'siHSNOW': (0, 10), 'siUICE': (-3, 3), 'siVICE': (-3, 3)}
a = np.fromfile(os.path.join(d, f'pickup_seaice.{it}.data'), '>f8')
assert a.size == 6 * NX * NY, a.size
a = a.reshape(6, NY, NX)
for k, f in enumerate(FLDS):
    if f == 'siTICE':
        continue                       # no cold-start namelist slot for TICE (as pan_AO)
    x = a[k]
    lo, hi = RANGE[f]; assert np.isfinite(x).all() and x.min() >= lo and x.max() <= hi, (f, x.min(), x.max())
    out = os.path.join(d, f'pickup_seaice_{f}.{it}.bin')
    x.astype('>f8').tofile(out)
    rb = np.fromfile(out, '>f8'); assert rb.size == NX * NY and np.array_equal(rb, x.ravel())
    print(f'{f}: min {x.min():.3f} max {x.max():.3f} mean {x.mean():.4f} -> {out}')
