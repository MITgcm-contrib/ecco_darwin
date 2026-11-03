#!/usr/bin/env python3
"""Initial erodible fine-sediment bed for pkg/sediment (sedBedFile, g/m^2).

Mud fraction: usSEABED DS-182 (Reid et al. 2006, https://pubs.usgs.gov/ds/2006/182/),
PAC_EXT lab-measured surface samples in the Santa Monica Bay box
(../work/sediment/usseabed_mud_smbay.csv, 451 samples; SCCWRP'94/'98, USGS
1997-98, USC Gorsline, NOAA). Each wet cell takes the median mud % of samples
within R_KM km; cells with < 3 samples take the median of their depth bin
(0-20, 20-50, 50-100, 100-500, > 500 m; 15, 52, 42, 75 % for the first four).
Erodible inventory: Palos Verdes shelf (Ferre et al. 2010, Cont. Shelf Res.,
doi:10.1016/j.csr.2010.01.011, p18 / Fig. 11d) mass remobilised at 0.32 Pa,
0.12-0.14 kg/m2 at 65 m (site B3, 65 % fines, Table 2) and 0.08-0.13 at 35 m.
Bed = M_REF * mud / MUD_REF, M_REF = 130 g/m2, MUD_REF = 65 %. The scaling with
the mud fraction is ASSUMED (no Santa Monica Bay erodibility measurements).

Usage: python3 make_sed_bed.py <out dir>   -> <out>/sed_bed_<nx>x<ny> (>f4, g/m^2)
"""
import json
import os
import sys

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
WORK = os.path.join(HERE, '..', 'work')
R_KM, M_REF, MUD_REF = 3.0, 130.0, 65.0
BINS = [0, 20, 50, 100, 500, 1e4]


def main():
    out = sys.argv[1]
    os.makedirs(out, exist_ok=True)
    g = json.load(open(os.path.join(WORK, 'grid_info.json')))
    nx, ny = g['nx'], g['ny']
    depth = -np.fromfile(os.path.join(WORK, g['bathy']), '>f4').reshape(ny, nx)
    s = pd.read_csv(os.path.join(WORK, 'sediment', 'usseabed_mud_smbay.csv')).dropna(subset=['mud_pct'])
    sd = s.dropna(subset=['depth_m'])
    binmed = [float(np.median(sd.mud_pct[(sd.depth_m >= a) & (sd.depth_m < b)])) for a, b in zip(BINS[:-1], BINS[1:])]
    X, Y = np.meshgrid(np.array(g['xc']), np.array(g['yc']))
    kx = 111.32 * np.cos(np.deg2rad(34.0))
    mud = np.full((ny, nx), np.nan)
    src = np.zeros((ny, nx), int)
    for j in range(ny):
        for i in range(nx):
            if depth[j, i] <= 0:
                continue
            d = np.hypot((s.lon.values - X[j, i]) * kx, (s.lat.values - Y[j, i]) * 111.32)
            near = s.mud_pct.values[d < R_KM]
            if near.size >= 3:
                mud[j, i], src[j, i] = np.median(near), 1
            else:
                mud[j, i] = binmed[np.searchsorted(BINS, depth[j, i], side='right') - 1]
    bed = np.where(depth > 0, M_REF * np.nan_to_num(mud) / MUD_REF, 0.0)
    bed.astype('>f4').tofile(os.path.join(out, f'sed_bed_{nx}x{ny}'))
    sh = (depth > 0) & (depth < 30)
    info = dict(depth_bin_median_mud_pct=dict(zip(['0-20', '20-50', '50-100', '100-500', '>500'], binmed)),
                frac_cells_from_samples=float(src[depth > 0].mean()),
                bed_g_m2_median={'<30 m': float(np.median(bed[sh])), '30-100 m': float(np.median(bed[(depth >= 30) & (depth < 100)])),
                                 '>100 m': float(np.median(bed[depth >= 100]))},
                M_REF_g_m2=M_REF, MUD_REF_pct=MUD_REF, R_km=R_KM)
    json.dump(info, open(os.path.join(out, 'sed_bed.json'), 'w'), indent=1)
    print(json.dumps(info))


if __name__ == '__main__':
    main()
