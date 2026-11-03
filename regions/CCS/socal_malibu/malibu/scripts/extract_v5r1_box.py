#!/usr/bin/env python3
"""Cut a Southern California Bight box out of the global ECCO-Darwin v05 LLC270
(ECCO V5r1) run on Pleiades, for the malibu_v01 nest in PACE years (2024-25).

Run on a Pleiades compute node (not pfe), e.g. from extract_v5r1_box.pbs.

  python3 extract_v5r1_box.py <out dir> [first month YYYY-MM] [last month YYYY-MM] [years]

Reads (big-endian float32, LLC270 compact format, 270 x 3510 per level):
  run2/diags/monthly/<FLD>.<iter>.data   monthly means, 50 levels
     THETA, SALTanom (= SALT - 35), uVel_C / vVel_C (east / north velocity at
     cell centres, UE_VEL_C / VN_VEL_C), and the 31 darwin tracers
  run2/iter70/EXF<f>_6hourly_<YYYY>      V5r1 atmosphere, 6-hourly
  run2/llc270_Mahowald_2009_soluble_iron_dust.bin   12 monthly records
  run_griddump/{XC,YC,RC,DRF,hFacC}.data
Writes <out>/<name>.bin as (nrec, [nz,] ny, nx) big-endian float32 on a
lon/lat box (ascending lon and lat; the LLC270 face is rotated, so the cut is
re-oriented), plus box.json with lon/lat/levels/record times.
Which face: the box must lie inside one LLC face (checked).
"""
import glob
import json
import os
import sys

import numpy as np

RUN = '/nobackupp19/dcarrol2/v05_V5r1/darwin3/run2'
GRID = '/nobackupp19/dcarrol2/v05_V5r1/darwin3/run_griddump'
N, NZ = 270, 50
BOX = dict(lon=(-122.0, -116.0), lat=(32.0, 36.0))
T0 = np.datetime64('1992-01-01T00:00:00')
DT = 1200.0
TRACERS = ['DIC', 'NO3', 'NO2', 'NH4', 'PO4', 'FeT', 'SiO2', 'DOC', 'DON', 'DOP', 'DOFe',
           'POC', 'PON', 'POP', 'POFe', 'POSi', 'PIC', 'ALK', 'O2',
           'c1', 'c2', 'c3', 'c4', 'c5', 'c6', 'c7', 'Chl1', 'Chl2', 'Chl3', 'Chl4', 'Chl5']
EXF = ['atemp', 'aqh', 'preci', 'uwind', 'vwind', 'swdn', 'lwdn']

# compact LLC layout: faces 1-2 (3N x N, stored j-major as 3N rows of N),
# face 3 (N x N), faces 4-5 (N rows of 3N); offsets in cells of one level
FACES = {1: (0, 3 * N, N), 2: (3 * N * N, 3 * N, N), 3: (6 * N * N, N, N),
         4: (7 * N * N, N, 3 * N), 5: (10 * N * N, N, 3 * N)}
LEV = 13 * N * N


def face(flat, f):
    off, ny, nx = FACES[f]
    return flat[off:off + ny * nx].reshape(ny, nx)


def locate():
    xc = np.fromfile(os.path.join(GRID, 'XC.data'), '>f4')
    yc = np.fromfile(os.path.join(GRID, 'YC.data'), '>f4')
    for f in FACES:
        X, Y = face(xc, f), face(yc, f)
        m = (X >= BOX['lon'][0]) & (X <= BOX['lon'][1]) & (Y >= BOX['lat'][0]) & (Y <= BOX['lat'][1])
        if m.any():
            jj, ii = np.where(m)
            sl = (slice(jj.min(), jj.max() + 1), slice(ii.min(), ii.max() + 1))
            Xs, Ys = X[sl], Y[sl]
            # orient to (lat, lon) ascending
            if np.ptp(Xs[:, 0]) > np.ptp(Xs[0, :]):       # lon varies along axis 0: transpose
                tr = True
                Xs, Ys = Xs.T, Ys.T
            else:
                tr = False
            fy = Ys[1, 0] < Ys[0, 0]
            fx = Xs[0, 1] < Xs[0, 0]
            Xs = Xs[::-1] if fy else Xs
            Xs = Xs[:, ::-1] if fx else Xs
            Ys = Ys[::-1] if fy else Ys
            Ys = Ys[:, ::-1] if fx else Ys
            assert np.allclose(Xs, Xs[0:1, :], atol=1e-3) and np.allclose(Ys, Ys[:, 0:1], atol=1e-3), \
                'box is not a lon/lat rectangle on this face'
            return dict(face=f, sl=sl, tr=tr, fy=fy, fx=fx, lon=Xs[0].astype(float), lat=Ys[:, 0].astype(float))
    raise SystemExit('box not found')


def cut(level_flat, L):
    a = (level_flat if level_flat.ndim == 2 else face(level_flat, L['face']))[L['sl']]
    if L['tr']:
        a = a.T
    if L['fy']:
        a = a[::-1]
    if L['fx']:
        a = a[:, ::-1]
    return a


FACE_OF_BOX = [None]


def read_level(fn, rec, k, nz):
    """One level; only the box's face is read once the box is located."""
    mm = np.memmap(fn, '>f4', 'r')
    off = (rec * nz + k) * LEV
    f = FACE_OF_BOX[0]
    if f is None:
        return np.asarray(mm[off:off + LEV])
    o, ny, nx = FACES[f]
    return np.asarray(mm[off + o:off + o + ny * nx]).reshape(ny, nx)


def meta_times(fn):
    txt = open(fn[:-5] + '.meta').read()
    a, b = txt.split('timeInterval = [')[1].split(']')[0].split()
    return float(a), float(b)


def main():
    out = sys.argv[1]
    m0 = np.datetime64(sys.argv[2] if len(sys.argv) > 2 else '2024-10')
    m1 = np.datetime64(sys.argv[3] if len(sys.argv) > 3 else '2025-03')
    years = [int(y) for y in (sys.argv[4].split(',') if len(sys.argv) > 4 else ['2024', '2025'])]
    os.makedirs(out, exist_ok=True)
    L = locate()
    FACE_OF_BOX[0] = L['face']
    print('face', L['face'], 'lon', L['lon'][0], L['lon'][-1], L['lon'].size,
          'lat', L['lat'][0], L['lat'][-1], L['lat'].size, flush=True)
    info = dict(lon=L['lon'].tolist(), lat=L['lat'].tolist(), face=L['face'])
    for g, nz in [('hFacC', NZ), ('XC', 1), ('YC', 1)]:
        a = np.stack([cut(read_level(os.path.join(GRID, g + '.data'), 0, k, nz), L) for k in range(nz)])
        a.astype('>f4').tofile(os.path.join(out, g + '.bin'))
    info['RC'] = np.fromfile(os.path.join(GRID, 'RC.data'), '>f4').astype(float).tolist()
    info['DRF'] = np.fromfile(os.path.join(GRID, 'DRF.data'), '>f4').astype(float).tolist()

    # ---- monthly 3-D means: pick files whose averaging interval centre is in [m0, m1]
    ref = sorted(glob.glob(os.path.join(RUN, 'diags', 'monthly', 'THETA.*.data')))
    sel, mids = [], []
    for fn in ref:
        a, b = meta_times(fn)
        mid = T0 + np.timedelta64(int((a + b) / 2), 's')
        if m0 <= mid.astype('datetime64[M]') <= m1:
            sel.append(fn.split('.')[-2])
            mids.append(str(mid))
    info['monthly_iters'], info['monthly_mid'] = sel, mids
    print('months', mids, flush=True)
    for fld in ['THETA', 'SALTanom', 'uVel_C', 'vVel_C'] + TRACERS:
        recs = []
        for it in sel:
            fn = os.path.join(RUN, 'diags', 'monthly', f'{fld}.{it}.data')
            recs.append(np.stack([cut(read_level(fn, 0, k, NZ), L) for k in range(NZ)]))
        a = np.stack(recs)
        if fld == 'SALTanom':
            a = np.where(a != 0, a + 35.0, 0.0)
            fld = 'SALT'
        a.astype('>f4').tofile(os.path.join(out, f'{fld}.bin'))
        print(fld, a.shape, float(np.nanmin(a)), float(np.nanmax(a)), flush=True)

    # ---- 6-hourly atmosphere (whole years) and monthly iron dust
    for y in years:
        for f in EXF:
            fn = os.path.join(RUN, 'iter70', f'EXF{f}_6hourly_{y}')
            if not os.path.exists(fn):
                print('missing', fn, flush=True)
                continue
            nrec = os.path.getsize(fn) // (LEV * 4)
            a = np.stack([cut(read_level(fn, r, 0, 1), L) for r in range(nrec)])
            a.astype('>f4').tofile(os.path.join(out, f'EXF{f}_6hourly_{y}'))
            print(f, y, a.shape, flush=True)
    fn = os.path.join(RUN, 'llc270_Mahowald_2009_soluble_iron_dust.bin')
    nrec = os.path.getsize(fn) // (LEV * 4)
    np.stack([cut(read_level(fn, r, 0, 1), L) for r in range(nrec)]).astype('>f4').tofile(
        os.path.join(out, 'iron_dust.bin'))
    json.dump(info, open(os.path.join(out, 'box.json'), 'w'), indent=1)
    print('done', out, flush=True)


if __name__ == '__main__':
    main()
