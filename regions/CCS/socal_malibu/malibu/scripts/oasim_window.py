#!/usr/bin/env python3
"""Build 12-record OASIM files for a run that crosses a year boundary.

darwin3 24885b71 EXF accepts only fldPeriod = -12 (12 monthly means, cycled,
no year suffix) or > 0 for RT_Ed/RT_Es, so calendar-month OASIM files
(edsE?NN_YYYY, 12 records each) cannot be read with useExfYearlyFields.
This writes oasim_<tag>/edsE?NN with month m taken from year y1 for m in
months1 and from y0 otherwise, e.g. for 16 Nov 2024 - 1 Mar 2025:
Jan-Mar from 2025, all other months from 2024.

Usage: python3 oasim_window.py <OASIM dir> <out dir> <y0> <y1> <months from y1, e.g. 1,2,3>
"""
import os
import sys

import numpy as np

src, out, y0, y1, m1 = sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5]
m1 = [int(m) for m in m1.split(',')]
os.makedirs(out, exist_ok=True)
for e in ('Ed', 'Es'):
    for l in range(1, 14):
        a0 = np.fromfile(os.path.join(src, f'eds{e}{l:02d}_{y0}'), '>f4').reshape(12, -1)
        a1 = np.fromfile(os.path.join(src, f'eds{e}{l:02d}_{y1}'), '>f4').reshape(12, -1)
        o = a0.copy()
        for m in m1:
            o[m - 1] = a1[m - 1]
        o.astype('>f4').tofile(os.path.join(out, f'eds{e}{l:02d}'))
print('wrote', out, 'months from', y1, m1, 'others from', y0)
