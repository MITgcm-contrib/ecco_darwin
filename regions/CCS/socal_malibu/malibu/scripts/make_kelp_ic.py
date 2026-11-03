#!/usr/bin/env python3
"""Kelp initial state from KelpWatch (Landsat canopy biomass) on the nest grid.

Replaces the placeholder kelp IC of make_ic_obcs.py (KELP_B0 = 100 gDW/m3 in
the bottom cell of every 5-25 m column) with the observed canopy for one quarter.

Data: SBC LTER, Bell, Cavanaugh & Siegel, "Time series of quarterly NetCDF
files of kelp biomass in the canopy from Landsat 5, 7 and 8, since 1984
(ongoing)", EDI knb-lter-sbc.74 rev. 33, doi:10.6073/pasta/1d4eed015c6d26d3fef87bbe782b2430;
box subset ../work/kelpwatch/kelpwatch_malibu_2024Q3_2025Q1.nc (biomass in
wet kg per 900 m2 pixel; Landsat calibration: SBC LTER Landsat protocol 2021,
Eq. 1).
Conversion (canopy wet -> whole-plant dry):
  DRY_WET = 0.094, mean dry:wet of blades (n = 1055, through 2016; stipes the
    same), SBC LTER Kelp NPP protocol p. 16;
  F_CANOPY = canopy / whole-plant mass, SBC LTER knb-lter-sbc.112.15
    (FSC_fraction_canopy, doi:10.6073/pasta/710844a10e2fe8873a4b49d97999d549):
    quarterly means Q1 0.277, Q2 0.336, Q3 0.263, Q4 0.191 (all site-dates: median
    0.258, 10-90% 0.016-0.52); calibration sites are 4.5-10.7 m deep, so the
    fraction at Malibu is UNCERTAIN.
  whole-plant dry per cell [g] = sum(pixel wet kg) x 1000 x DRY_WET / F_CANOPY
Each pixel goes to the nest cell containing it (pixels in model land cells go
to the nearest wet cell within 2 cells, else are dropped and counted). As in
make_ic_obcs.py, the stock sits in the bottom wet cell (kLowC) as Bmag
[gDW/m3] = mass / (cell area x drF x hFacC); Qmag = 20 mgN/gDW, Cmag = 0.3 B.

Usage: python3 make_kelp_ic.py <work dir with init/> <hFacC.data> [YYYYQn] [--write]
(without --write it only reports the comparison with the current IC)
"""
import json
import os
import sys

import numpy as np
import xarray as xr

HERE = os.path.dirname(os.path.abspath(__file__))
WORK = os.path.join(HERE, '..', 'work')
DRY_WET = 0.094
F_CANOPY = {1: 0.277, 2: 0.336, 3: 0.263, 4: 0.191}
KELP_Q0, KELP_CARBON = 20.0, 0.3        # as make_ic_obcs.py


def main():
    write = '--write' in sys.argv
    args = [a for a in sys.argv[1:] if a != '--write']
    wdir, hfile = args[:2]
    yq = args[2] if len(args) > 2 else '2024Q4'
    year, q = int(yq[:4]), int(yq[-1])
    g = json.load(open(os.path.join(WORK, 'grid_info.json')))
    nx, ny, nr = g['nx'], g['ny'], g['nr']
    xg, yg, yc = np.array(g['xg']), np.array(g['yg']), np.array(g['yc'])
    drf = np.array(g['delr'])
    hf = np.fromfile(hfile, '>f4').reshape(nr, ny, nx)
    wet = hf[0] > 0
    kbot = np.clip((hf > 0).sum(0) - 1, 0, nr - 1)
    area = (np.diff(xg)[None, :] * 111.32e3 * np.cos(np.deg2rad(yc))[:, None]) * \
        (np.diff(yg)[:, None] * 111.32e3)                               # m2
    kw = xr.open_dataset(os.path.join(WORK, 'kelpwatch', 'kelpwatch_malibu_2024Q3_2025Q1.nc'))
    k = int(np.where((kw.year.values == year) & (kw.quarter_of_year.values == q))[0][0])
    wetkg = np.nan_to_num(kw.biomass.values[:, k])
    lon, lat = kw.longitude.values, kw.latitude.values
    i = np.searchsorted(xg, lon) - 1
    j = np.searchsorted(yg, lat) - 1
    ok = (i >= 0) & (i < nx) & (j >= 0) & (j < ny) & (wetkg > 0)
    mass = np.zeros((ny, nx))                                           # canopy wet kg per cell
    dropped = 0.0
    wj, wi = np.where(wet)
    for jj, ii, m in zip(j[ok], i[ok], wetkg[ok]):
        if not wet[jj, ii]:
            dd = np.hypot(wj - jj, wi - ii)
            n = dd.argmin()
            if dd[n] > 2:
                dropped += m
                continue
            jj, ii = wj[n], wi[n]
        mass[jj, ii] += m
    dry = mass * 1e3 * DRY_WET / F_CANOPY[q]                            # whole-plant g dry per cell
    B = np.zeros((nr, ny, nx))
    jj, ii = np.where(dry > 0)
    kk = kbot[jj, ii]
    B[kk, jj, ii] = dry[jj, ii] / (area[jj, ii] * drf[kk] * hf[kk, jj, ii])
    # current (placeholder) IC for comparison
    tag = f'{nx}x{ny}x{nr}'
    B_old = np.fromfile(os.path.join(wdir, 'init', f'Bmag_{tag}'), '>f4').reshape(nr, ny, nx)
    vol = area[None] * drf[:, None, None] * hf
    tot_old, tot_new = (B_old * vol).sum() / 1e6, (B * vol).sum() / 1e6
    info = dict(quarter=yq, kelpwatch_pixels_with_kelp=int((wetkg > 0).sum()),
                canopy_wet_t=float(wetkg.sum() / 1e3), dropped_wet_t=float(dropped / 1e3),
                dry_wet=DRY_WET, f_canopy=F_CANOPY[q],
                cells_new=int((dry > 0).sum()), cells_old=int((B_old.sum(0) > 0).sum()),
                total_dry_t_new=float(tot_new), total_dry_t_old=float(tot_old),
                cell_dry_kg_m2_new_median=float(np.median(dry[dry > 0] / 1e3 / area[dry > 0])),
                cell_dry_kg_m2_new_max=float((dry / 1e3 / area).max()))
    print(json.dumps(info, indent=1))
    if write:
        Q = np.where(B > 0, KELP_Q0, 0.0)
        C = KELP_CARBON * B
        for name, arr in (('Bmag', B), ('Cmag', C), ('Qmag', Q)):
            arr.astype('>f4').tofile(os.path.join(wdir, 'init', f'{name}_{tag}'))
        json.dump(info, open(os.path.join(wdir, 'kelp_ic.json'), 'w'), indent=1)
        print('wrote Bmag/Cmag/Qmag in', os.path.join(wdir, 'init'))


if __name__ == '__main__':
    main()
