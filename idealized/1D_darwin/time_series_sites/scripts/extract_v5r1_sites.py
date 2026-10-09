#!/usr/bin/env python3
"""
Extract 1-D column inputs and validation fields at ocean time-series sites
from the ECCO-Darwin v05 llc270 v5r1 solution (Pleiades).

Read-only on the v5r1 inputs/run; writes to OUT/<site>/.
Python 3.6 + numpy 1.14 compatible (pfe system python).

Per site it writes (big-endian float32, 1x1 horizontal):
  bathy_1x1                 depth (negative, m)
  diffkr_1x1x50             xx_diffkr.effective.V5alpha.iter50 column
  T_ini_1x1x50, S_ini_1x1x50  from pickup.0000000002 (fill below bottom)
  ptracer01_ini..ptracer31_ini from pickup_ptracers.0000000002
  EXF/EXF<fld>_6hourly_YYYY  1991-2025, 6-hourly (as iter70)
  iron_dust_1x1_mon         Mahowald 2009 soluble Fe dust, 12 months
  apCO2/apCO2_YYYY          daily NOAA MBL apCO2 interpolated to site lat
  relax_T_1x1x50_mon, relax_S_1x1x50_mon  monthly THETA/SALT, padded
                            (Jan1992, Jan1992..Dec2025, Dec2025) = 410 recs
  v5r1_monthly_<var>.bin, v5r1_daily_<var>.bin  validation time series
  v5r1_monthly_steps.txt, v5r1_daily_steps.txt  timeStepNumber of records
  site_info.txt
"""
import os, sys, glob, math
import numpy as np
from multiprocessing import Pool

RUN = '/nobackup/dcarrol2/v05_V5r1/darwin3/run2'
INP = '/nobackup/hzhang1/pub/llc270_FWD/v5r1'
OUT = sys.argv[1] if len(sys.argv) > 1 else '/nobackup/dcarrol2/1-D/sites'
NPROC = int(os.environ.get('NPROC', '24'))

SITES = [  # name, lat, lon (deg E, -180..180), nominal depth (m) for sanity check
    ('HOT',    22.75, -158.00, 4800.),
    ('BATS',   31.67,  -64.17, 4600.),
    ('HydroS', 32.17,  -64.50, 3200.),
    ('PAPA',   50.10, -144.90, 4200.),
    ('PAP',    49.00,  -16.50, 4800.),
]

NX = 270
NR = 50
N2 = NX * NX * 13                       # compact 270 x 3510
FACES = [(270, 810), (270, 810), (270, 270), (810, 270), (810, 270)]  # (nx, ny)
DRF = np.array([
    10.00, 10.00, 10.00, 10.00, 10.00, 10.00, 10.00, 10.01,
    10.03, 10.11, 10.32, 10.80, 11.76, 13.42, 16.04, 19.82, 24.85,
    31.10, 38.42, 46.50, 55.00, 63.50, 71.58, 78.90, 85.15, 90.18,
    93.96, 96.58, 98.25, 99.25, 100.01, 101.33, 104.56, 111.33, 122.83,
    139.09, 158.94, 180.83, 203.55, 226.50, 249.50, 272.50, 295.50, 318.50,
    341.50, 364.50, 387.50, 410.50, 433.50, 456.50])
RF = np.concatenate([[0.], -np.cumsum(DRF)])
EXF_FLDS = ['atemp', 'aqh', 'preci', 'uwind', 'vwind', 'swdn', 'lwdn']
YEARS = list(range(1991, 2026))


def wbin(path, a):
    d = os.path.dirname(path)
    if d and not os.path.isdir(d):
        os.makedirs(d)
    np.asarray(a, dtype='>f4').tofile(path)


def read_grid():
    """xC, yC on the compact (flat) llc270 index from tile00N.mitgrid."""
    xc = np.zeros(N2); yc = np.zeros(N2)
    off = 0
    for f, (nx, ny) in enumerate(FACES):
        a = np.fromfile(os.path.join(RUN, 'tile%03d.mitgrid' % (f + 1)), '>f8')
        a = a.reshape(-1, ny + 1, nx + 1)
        xc[off:off + nx * ny] = a[0, :ny, :nx].ravel()
        yc[off:off + nx * ny] = a[1, :ny, :nx].ravel()
        off += nx * ny
    return xc, yc


def flat2face(p):
    off = 0
    for f, (nx, ny) in enumerate(FACES):
        if p < off + nx * ny:
            q = p - off
            return f + 1, q % nx + 1, q // nx + 1   # face, i, j (1-based)
        off += nx * ny


def find_sites(xc, yc, depth):
    out = []
    lat0 = np.radians(yc); lon0 = np.radians(xc)
    for name, lat, lon, zn in SITES:
        la, lo = math.radians(lat), math.radians(lon)
        d = np.arccos(np.clip(np.sin(la) * np.sin(lat0) +
                              np.cos(la) * np.cos(lat0) * np.cos(lon0 - lo), -1, 1))
        d[depth >= 0] = 1e9                 # wet points only
        p = int(np.argmin(d))
        out.append(dict(name=name, lat=lat, lon=lon, znom=zn, p=p,
                        xc=xc[p], yc=yc[p], dist_km=d[p] * 6371.,
                        depth=-depth[p], fij=flat2face(p)))
    return out


def col(path, prec, ps, nrec_off, nlev=NR):
    """Read columns at flat points ps; records nrec_off..nrec_off+nlev-1."""
    m = np.memmap(path, dtype=prec, mode='r')
    idx = (nrec_off + np.arange(nlev))[None, :] * N2 + np.asarray(ps)[:, None]
    return np.array(m[idx.ravel()]).reshape(len(ps), nlev).astype('f8')


def fill_below(prof, nwet):
    p = prof.copy()
    if nwet < NR:
        p[nwet:] = p[nwet - 1]
    return p


def meta_info(path):
    txt = open(path).read()
    nd = int(txt.split('nDims = [')[1].split(']')[0])
    step = int(txt.split('timeStepNumber = [')[1].split(']')[0])
    return nd, step


def task_diag(args):
    path, ps, nlev = args
    return col(path, '>f4', ps, 0, nlev)


def task_exf(args):
    path, ps = args
    n = os.path.getsize(path) // (4 * N2)
    m = np.memmap(path, dtype='>f4', mode='r')
    idx = (np.arange(n)[None, :] * N2 + np.asarray(ps)[:, None])
    return np.array(m[idx.ravel()]).reshape(len(ps), n)


def main():
    pool = Pool(NPROC)
    xc, yc = read_grid()
    depth = np.fromfile(os.path.join(INP, 'input_bin', 'bathy270_filled_noCaspian_r4'), '>f4').astype('f8')
    sites = find_sites(xc, yc, depth)
    ps = [s['p'] for s in sites]
    for s in sites:
        s['nwet'] = int(np.sum(RF[:-1] > -s['depth']))
        s['dir'] = os.path.join(OUT, s['name'])
        if not os.path.isdir(s['dir']):
            os.makedirs(s['dir'])
        print('%-7s lat %7.3f lon %8.3f -> face %d i %d j %d  xc %8.3f yc %7.3f  %.1f km  '
              'depth %.0f m (nominal %.0f), nwet %d' % (
                  s['name'], s['lat'], s['lon'], s['fij'][0], s['fij'][1], s['fij'][2],
                  s['xc'], s['yc'], s['dist_km'], s['depth'], s['znom'], s['nwet']), flush=True)

    # ---- static fields and initial conditions
    dkr = col(os.path.join(INP, 'input_bin', 'xx_diffkr.effective.V5alpha.iter50'), '>f4', ps, 0)
    pk = os.path.join(RUN, 'pickup.0000000002.data')     # Uvel Vvel Theta Salt ... (3-D x50 each)
    T0 = col(pk, '>f8', ps, 2 * NR); S0 = col(pk, '>f8', ps, 3 * NR)
    ptr = [col(os.path.join(INP, 'input_darwin_bin', 'pickup_ptracers.0000000002.data'),
               '>f8', ps, it * NR) for it in range(31)]
    fe = col(os.path.join(INP, 'input_darwin_bin', 'llc270_Mahowald_2009_soluble_iron_dust.bin'),
             '>f4', ps, 0, 12)
    for n, s in enumerate(sites):
        d, nw = s['dir'], s['nwet']
        wbin(os.path.join(d, 'bathy_1x1'), [-s['depth']])
        wbin(os.path.join(d, 'diffkr_1x1x50'), fill_below(dkr[n], nw))
        wbin(os.path.join(d, 'T_ini_1x1x50'), fill_below(T0[n], nw))
        wbin(os.path.join(d, 'S_ini_1x1x50'), fill_below(S0[n], nw))
        for it in range(31):
            wbin(os.path.join(d, 'ptracer%02d_ini' % (it + 1)), fill_below(ptr[it][n], nw))
        wbin(os.path.join(d, 'iron_dust_1x1_mon'), fe[n])
        with open(os.path.join(d, 'site_info.txt'), 'w') as f:
            for k in ['name', 'lat', 'lon', 'znom', 'xc', 'yc', 'dist_km', 'depth', 'nwet', 'fij', 'p']:
                f.write('%s = %s\n' % (k, s[k]))
    print('static + IC done', flush=True)

    # ---- apCO2 (2 x 256 lat grid, zonally uniform), daily yearly files
    inc = np.array([0.6958694, 0.6999817, 0.7009048, 0.7012634, 0.7014313] + [0.7017418] * 245 +
                   [0.7014313, 0.7012634, 0.7009048, 0.6999817, 0.6958694])
    lats = -89.4628220 + np.concatenate([[0.], np.cumsum(inc)])
    for y in YEARS + [2026]:          # 2026: daily interpolation on 31 Dec 2025 needs 1 Jan 2026
        a = np.fromfile(os.path.join(RUN, 'apCO2_%d' % y), '>f4').reshape(-1, 256, 2)[:, :, 0]
        for s in sites:
            v = np.array([np.interp(s['yc'], lats, r) for r in a])
            wbin(os.path.join(s['dir'], 'apCO2', 'apCO2_%d' % y), v)
    print('apCO2 done', flush=True)

    # ---- EXF 6-hourly yearly
    jobs = [(os.path.join(INP, 'input_bin', 'iter70', 'EXF%s_6hourly_%d' % (f, y)), ps)
            for f in EXF_FLDS for y in YEARS]
    res = pool.map(task_exf, jobs, chunksize=1)
    for (path, _), r in zip(jobs, res):
        for n, s in enumerate(sites):
            wbin(os.path.join(s['dir'], 'EXF', os.path.basename(path)), r[n])
    print('EXF done', flush=True)

    # ---- monthly and daily diagnostics
    for freq in ['monthly', 'daily']:
        ddir = os.path.join(RUN, 'diags', freq)
        names = sorted(set(os.path.basename(f).split('.')[0] for f in glob.glob(ddir + '/*.meta')))
        for v in names:
            metas = sorted(glob.glob(os.path.join(ddir, v + '.*.meta')))
            nd, _ = meta_info(metas[0])
            nlev = NR if nd == 3 else 1
            steps = [int(os.path.basename(m).split('.')[1]) for m in metas]
            res = pool.map(task_diag, [(m[:-5] + '.data', ps, nlev) for m in metas], chunksize=8)
            arr = np.stack(res, axis=1)            # (nsite, ntime, nlev)
            for n, s in enumerate(sites):
                wbin(os.path.join(s['dir'], 'v5r1_%s' % freq, v + '.bin'), arr[n])
                np.savetxt(os.path.join(s['dir'], 'v5r1_%s' % freq, 'steps.txt'), steps, fmt='%d')
            print('%s %s: %d records, nlev %d' % (freq, v, len(steps), nlev), flush=True)

    # ---- relaxation files from monthly THETA / SALTanom (+35)
    for n, s in enumerate(sites):
        d = os.path.join(s['dir'], 'v5r1_monthly')
        T = np.fromfile(os.path.join(d, 'THETA.bin'), '>f4').reshape(-1, NR)
        S = np.fromfile(os.path.join(d, 'SALTanom.bin'), '>f4').reshape(-1, NR) + 35.
        for nm, a in [('T', T), ('S', S)]:
            a = np.array([fill_below(r, s['nwet']) for r in a])
            a = np.concatenate([a[:1], a, a[-1:]], axis=0)
            wbin(os.path.join(s['dir'], 'relax_%s_1x1x50_mon' % nm), a)
    print('relax done', flush=True)


if __name__ == '__main__':
    main()
