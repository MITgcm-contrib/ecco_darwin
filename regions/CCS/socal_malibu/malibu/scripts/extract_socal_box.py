#!/usr/bin/env python3
"""Cut the Malibu nest box (+ margin) out of the socal 2-km nest's daily 3-D output.

Source: socal regional_run/run_prod/diags/daily/state_3d.<iter>.data (THETA, SALT, UVEL, VVEL,
WVEL; 168 x 101 x 60, daily means, deltaT 120 s, startDate 2025-01-01) and the run's grid
(socal_ncgrid.nc fields written to grid_socal.npz by the caller's copy, see below).
Out (same layout as extract_v5r1_box.py, so make_ic_obcs.py can read it):
  <out>/{THETA,SALT,uVel_C,vVel_C}.bin  big-endian float32 (nrec, nz, ny, nx), daily means;
  U and V averaged from the C-grid faces to cell centers (east/north: the socal grid is lon/lat)
  <out>/hFacC.bin, <out>/box.json (lon, lat, RC, DRF, daily_iters, daily_mid)

Monthly BGC (2026-10-10): diags/monthly/ptracers_3d.<iter>.data (31 darwin v05 tracers, calendar-month
means written at month end) -> <out>/<name>.bin with make_ic_obcs.py's PARENT_TRACERS names
(DIC ... O2, c1-c7, Chl1-Chl5), box.json monthly_iters / monthly_mid (mid-month of each mean).

Usage: python3 extract_socal_box.py <socal run dir> <grid npz> <out dir> [lon0 lon1 lat0 lat1]
"""
import datetime as dt
import glob
import json
import os
import sys

import numpy as np

NX, NY, NR, DT = 168, 101, 60, 120.0
T0 = dt.datetime(2025, 1, 1)


def main():
    run, gfile, out = sys.argv[1:4]
    lon0, lon1, lat0, lat1 = map(float, sys.argv[4:8]) if len(sys.argv) > 7 else (-119.25, -118.30, 33.65, 34.25)
    g = np.load(gfile)
    xc, yc, hf, drf = g['XC'][0], g['YC'][:, 0], g['HFacC'], g['drF']
    i = np.where((xc >= lon0) & (xc <= lon1))[0]
    j = np.where((yc >= lat0) & (yc <= lat1))[0]
    i0, i1, j0, j1 = i[0], i[-1] + 1, j[0], j[-1] + 1
    os.makedirs(out, exist_ok=True)
    fs = sorted(glob.glob(os.path.join(run, 'diags', 'daily', 'state_3d.*.data')))
    n = len(fs)
    outs = {k: np.memmap(os.path.join(out, f'{k}.bin'), '>f4', 'w+', shape=(n, NR, j1 - j0, i1 - i0))
            for k in ('THETA', 'SALT', 'uVel_C', 'vVel_C')}
    iters, mids = [], []
    for t, f in enumerate(fs):
        a = np.fromfile(f, '>f4').reshape(5, NR, NY, NX)
        outs['THETA'][t] = a[0, :, j0:j1, i0:i1]
        outs['SALT'][t] = a[1, :, j0:j1, i0:i1]
        u, v = a[2], a[3]
        uc = 0.5 * (u[:, :, :-1] + u[:, :, 1:])                  # west faces -> centers (i = 0..NX-2)
        vc = 0.5 * (v[:, :-1, :] + v[:, 1:, :])
        uc = np.concatenate([uc, u[:, :, -1:]], axis=2)            # last column: its west face
        vc = np.concatenate([vc, v[:, -1:, :]], axis=1)
        outs['uVel_C'][t] = uc[:, j0:j1, i0:i1]
        outs['vVel_C'][t] = vc[:, j0:j1, i0:i1]
        it = int(f.split('.')[-2])
        iters.append(f'{it:010d}')
        mids.append((T0 + dt.timedelta(seconds=it * DT - 43200)).strftime('%Y-%m-%dT%H:%M:%S'))
    for m in outs.values():
        m.flush()
    names = ['DIC', 'NO3', 'NO2', 'NH4', 'PO4', 'FeT', 'SiO2', 'DOC', 'DON', 'DOP', 'DOFe', 'POC', 'PON', 'POP',
             'POFe', 'POSi', 'PIC', 'ALK', 'O2', 'c1', 'c2', 'c3', 'c4', 'c5', 'c6', 'c7',
             'Chl1', 'Chl2', 'Chl3', 'Chl4', 'Chl5']
    fm = sorted(glob.glob(os.path.join(run, 'diags', 'monthly', 'ptracers_3d.*.data')))
    miters, mmid = [], []
    if fm:
        mo = {k: np.memmap(os.path.join(out, f'{k}.bin'), '>f4', 'w+', shape=(len(fm), NR, j1 - j0, i1 - i0))
              for k in names}
        prev = T0
        for t, f in enumerate(fm):
            a = np.fromfile(f, '>f4').reshape(31, NR, NY, NX)
            for q, k in enumerate(names):
                mo[k][t] = a[q, :, j0:j1, i0:i1]
            it = int(f.split('.')[-2])
            end = T0 + dt.timedelta(seconds=it * DT)
            miters.append(f'{it:010d}')
            mmid.append((prev + (end - prev) / 2).strftime('%Y-%m-%dT%H:%M:%S'))
            prev = end
        for m in mo.values():
            m.flush()
        print(f'monthly BGC: {len(fm)} means, mid {mmid}')
    hf[:, j0:j1, i0:i1].astype('>f4').tofile(os.path.join(out, 'hFacC.bin'))
    zf = -np.r_[0.0, np.cumsum(drf)]
    json.dump(dict(lon=xc[i0:i1].tolist(), lat=yc[j0:j1].tolist(), RC=(0.5 * (zf[1:] + zf[:-1])).tolist(),
                   DRF=drf.tolist(), daily_iters=iters, daily_mid=mids, monthly_iters=miters, monthly_mid=mmid, source=run,
                   index=dict(i0=int(i0), i1=int(i1), j0=int(j0), j1=int(j1))),
              open(os.path.join(out, 'box.json'), 'w'), indent=1)
    print(f'{n} days {mids[0]} .. {mids[-1]}; box i {i0}-{i1 - 1}, j {j0}-{j1 - 1} '
          f'({i1 - i0} x {j1 - j0}); THETA surface day 1 '
          f'{outs["THETA"][0, 0][hf[0, j0:j1, i0:i1] > 0].min():.2f}-{outs["THETA"][0, 0][hf[0, j0:j1, i0:i1] > 0].max():.2f} C')
    print('EXTRACT_DONE')


if __name__ == '__main__':
    main()
