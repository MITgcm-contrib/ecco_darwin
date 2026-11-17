#!/usr/bin/env python3
"""Expand Humboldt_Qsg_daily.npz into yearly iceplume runoffQsg files for the Nares run.
One file per year, runoffQsg_YYYY: ndays records of a 1120x360 field (big-endian
real*4, x fastest), m3/s at the plume cells, 0 elsewhere. Read by iceplume through
EXF with useExfYearlyFields (data.exf) and runoffQsgperiod = 86400,
runoffQsgstartdate = 2015-01-01 12:00 (data.iceplume): records are daily means
centred at noon. Because of that, the first half day of 2015 interpolates from
the last record of runoffQsg_2014, so a 2014 file is written too, holding the
2015-01-01 field (mid-winter, ~0) on every day.
Usage: python3 gen_runoffQsg_files.py <npz> <out dir> [case=baseline]"""
import sys, os, numpy as np
npz, out = sys.argv[1], sys.argv[2]
case = sys.argv[3] if len(sys.argv) > 3 else 'baseline'
NX, NY = 1120, 360
z = np.load(npz)
dates = z['dates']; cells = z[f'{case}_cells']; Q = z[f'{case}_Q']
assert Q.shape == (dates.size, cells.shape[0]) and np.isfinite(Q).all() and (Q >= 0).all()
os.makedirs(out, exist_ok=True)
pre = 'runoffQsg' if case == 'baseline' else f'runoffQsg_{case}'
years = sorted({d[:4] for d in dates})
def write(fn, rows):
    fld = np.zeros((NY, NX), '>f4')
    with open(fn, 'wb') as f:
        for q in rows:
            fld[:] = 0
            np.add.at(fld, (cells[:, 0], cells[:, 1]), q)   # cells are (j, i)
            fld.tofile(f)
    return len(rows)
first = int(years[0])
n14 = 366 if (first-1) % 4 == 0 else 365
write(os.path.join(out, f'{pre}_{first-1}'), [Q[0]] * n14)
for y in years:
    sel = np.array([d.startswith(y) for d in dates])
    fn = os.path.join(out, f'{pre}_{y}')
    n = write(fn, Q[sel])
    a = np.fromfile(fn, '>f4').reshape(n, NY, NX)
    assert a.shape[0] == (366 if int(y) % 4 == 0 else 365), (y, a.shape)
    err = abs(a.sum(axis=(1, 2)).astype(float) - Q[sel].sum(1)).max()
    assert err < 1e-3 * max(1.0, Q[sel].sum(1).max()), (y, err)
    print(f'{fn}: {n} days, mean {Q[sel].sum(1).mean():7.1f} max {Q[sel].sum(1).max():7.1f} m3/s, '
          f'readback err {err:.2e}, {os.path.getsize(fn)/1e6:.0f} MB')
print('done')
