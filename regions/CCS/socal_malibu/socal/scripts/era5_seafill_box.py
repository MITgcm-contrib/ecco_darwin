#!/usr/bin/env python3
"""Sea-filled ERA5 box for the socal 2-km nest's EXF forcing (2026-10-10).

Raw ERA5 (hzhang1 forcing/era5, 1280 x 640 Gaussian grid) is defined over land; EXF interpolates
bilinearly, so 604 of 1126 socal shelf cells (< 100 m) took part of their forcing from land cells,
which in Mar-Jun 2025 are +1.2 C warmer (daily max +4.5 C), 41% less windy and 0.6 g/kg drier than
the adjacent ocean cells (python/era5_box.py analysis). Same problem and fix as seafill_v5r1_exf.py.

Land mask: ERA5 cells whose Mar-Jun 2025 mean diurnal range of 2-m air temperature exceeds 3 C
(ocean cells 1-3 C, land up to 15 C; no ERA5 land-sea mask in the forcing directory).
Each land cell, every record and field, takes the value of the nearest ocean cell (distance in
grid-index space, scipy distance_transform_edt). Writes <out>/ERA5_<fld>_<year> (real*4, big-endian,
nlat x nlon box) for each year given, plus era5_box_grid.txt with the EXF_NML_04 lines.

Usage: era5_seafill_box.py <era5 dir> <data.exf> <out dir> lon0 lon1 lat0 lat1 year [year ...]
"""
import os
import re
import sys

import numpy as np
from scipy import ndimage

FIELDS = ('u10m', 'v10m', 'tmp2m_degC', 'spfh2m', 'dsw', 'dlw', 'rain', 'pres')
EXFNAME = {'u10m': 'uwind', 'v10m': 'vwind', 'tmp2m_degC': 'atemp', 'spfh2m': 'aqh', 'dsw': 'swdown',
           'dlw': 'lwdown', 'rain': 'precip', 'pres': 'apressure'}


def main():
    d, exf, out = sys.argv[1:4]
    lo0, lo1, la0, la1 = map(float, sys.argv[4:8])
    years = sys.argv[8:]
    s = open(exf).read()
    blk = re.search(r'uwind_lat_inc\s*=(.*?)uwind_nlon', s, re.S).group(1)
    inc = []
    for tok in re.findall(r'(?:(\d+)\*)?([0-9.]+)D0', blk):
        inc += [float(tok[1])] * (int(tok[0]) if tok[0] else 1)
    lat = -89.784885 + np.r_[0.0, np.cumsum(inc)]
    lon = 0.28125 * np.arange(1280) - 360.0
    j = np.where((lat >= la0) & (lat <= la1))[0]
    i = np.where((lon >= lo0) & (lon <= lo1))[0]
    sj, si = slice(j[0], j[-1] + 1), slice(i[0], i[-1] + 1)
    os.makedirs(out, exist_ok=True)
    t = np.memmap(f'{d}/ERA5_tmp2m_degC_2025', '>f4', 'r').reshape(-1, 640, 1280)[59 * 24:181 * 24, sj, si]
    t = np.array(t, 'f8').reshape(-1, 24, j.size, i.size)
    dr = (t.max(1) - t.min(1)).mean(0)
    land = dr > 3.0
    _, (jj, ii) = ndimage.distance_transform_edt(land, return_indices=True)
    print(f'box {j.size} x {i.size}: land {land.sum()}, ocean {(~land).sum()}; diurnal range {dr.min():.1f}-{dr.max():.1f} C')
    np.save(os.path.join(out, 'land_mask.npy'), land)
    for y in years:
        for f in FIELDS:
            a = np.array(np.memmap(f'{d}/ERA5_{f}_{y}', '>f4', 'r').reshape(-1, 640, 1280)[:, sj, si])
            a = a[:, jj, ii]                               # land cells <- nearest ocean cell
            a.astype('>f4').tofile(os.path.join(out, f'ERA5_{f}_{y}'))
            print(y, f, a.shape)
    dl = np.diff(lat[j])
    with open(os.path.join(out, 'era5_box_grid.txt'), 'w') as fh:
        for f in FIELDS:
            n = EXFNAME[f]
            fh.write(f' {n}_lon0 = {lon[i[0]]:.7f}D0,\n {n}_lon_inc = 0.2812500D0,\n'
                     f' {n}_lat0 = {lat[j[0]]:.7f}D0,\n {n}_lat_inc = ' + ','.join(f'{x:.7f}D0' for x in dl) + ',\n'
                     f' {n}_nlon = {i.size},\n {n}_nlat = {j.size},\n')
    print('ERA5_SEAFILL_DONE')


if __name__ == '__main__':
    main()
