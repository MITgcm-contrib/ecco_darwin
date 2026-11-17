#!/usr/bin/env python3
"""gen_obcs_tides.py -- barotropic tidal OBCS files (pkg/obcs, ALLOW_OBCS_TIDES) from
FES2014 depth-mean currents, for any regional child grid (2026-10-08).

Region-agnostic: the grid enters only through a points file (--points, .npz) holding,
for each open boundary b in west/east/south/north present,
    b_lat, b_lon   (n,)  lat/lon of the OB cell centre of every entry along the boundary
    b_cs, b_sn     (n,)  CHILD-grid AngleCS/AngleSN at that cell
    b_has          (n,)  True where the entry has an open boundary
n is the full boundary length (Ny for west/east, Nx for south/north) as MITgcm reads it
with MDS_READ_SEC_XZ (globalFiles). Make it with a small per-region extractor
(e.g. pan_AO/prep/extract_ob_points.py).

What MITgcm does with the files (darwin3 ckpt69k pkg/obcs, obcs_init_fixed.F:535-575 and
obcs_add_tides.F):
    Ph file = phase lag in SECONDS; locPh = Ph*2*pi/T; Cs = Am*cos(locPh), Sn = Am*sin(locPh)
    u_tide(t) = Cs*cos(w t) + Sn*sin(w t) = Am*cos(w t - 2 pi Ph/T),  t = myTime, w = 2 pi/T
and the same u_tide is added at every level of OB?u and v_tide to OB?v (both normal and
tangential, grid-x / grid-y components). Phase zero is therefore myTime = 0, i.e. the
run's startDate (myTime stays continuous across pickups, so one file serves all segments).

Atlas convention (FES2014, Greenwich phase lag G, UTide/t_tide):
    u(t) = f A cos(V(t) + u - G),  V(t) = V0 + w (t - t0)
so with t0 = startDate:   Am = f A,   2 pi Ph/T = G - V0 - u   (mod 2 pi).
V0 is evaluated at startDate exactly; f and u (18.6-yr nodal modulation) at --nodal-date
(default: startDate + 182 d). f,u drift with the 18.6-yr nodal cycle: tested against UTide's
hourly nodal correction, holding f,u fixed for a year gives errors (fraction of amplitude) of
M2 1.4 %, K1 3.9 %, O1 7 %, K2 10 % after 360 d. obcs_init_fixed reads the files at every
(re)start, so REGENERATE PER CHAIN SEGMENT with --nodal-date = segment midpoint (keep --start =
the run's startDate: V0 stays referenced to myTime = 0).

Rotation into the child grid acts on the complex harmonic constants (linear, real angles):
    Z_x =  CS Z_E + SN Z_N,   Z_y = -SN Z_E + CS Z_N,   Z = A exp(i G)
(the same rotation as rotate_obcs_velocity.py; the AngleCS/SN must be the CHILD grid's).

Periods come from UTide's constituent frequencies (cycles/hour), T = 3600/f exactly, and
are printed as the OBCS_tidalPeriod namelist line; OBCS_tideCompSize = 10 caps the count.

Atlas files: --u-dir / --v-dir with one netCDF per constituent named <c>*.nc or <c>_*.nc
(lower case, e.g. m2.nc, m2_45N.nc), variables Ua/Ug (Va/Vg), cm/s and degrees, on a
regular lat/lon grid (may be lat-cropped; longitude must cover the boundary points).

    python gen_obcs_tides.py --points ob_points.npz --u-dir FES/eastward_velocity \
        --v-dir FES/northward_velocity --start 2020-01-01 --out-dir OUT --prec 64 [--selftest]

Writes OB{W,E,S,N}_{u,v}Tid{Am,Ph}.bin (big-endian, --prec 32/64 = readBinaryPrec),
records = constituents in --cons order, and data.obcs_tides.snippet.
"""
import argparse
import os
import re
import sys

import numpy as np
import xarray as xr
from scipy import ndimage
from scipy.interpolate import RegularGridInterpolator

CONS = ['M2', 'S2', 'N2', 'K2', 'K1', 'O1', 'P1', 'Q1']
BND = {'west': 'W', 'east': 'E', 'south': 'S', 'north': 'N'}


# ---------------------------------------------------------------- astronomy (UTide)
def _utide_time(date):
    from utide._time_conversion import _normalize_time
    return float(_normalize_time(np.array([0.0]), date)[0])


def nodal(cons, start, nodal_date, lat):
    """T (s), f, u (cycles, at nodal_date), V0 (cycles, at start) per constituent"""
    from utide import harmonics, ut_constants
    names = [n.strip() for n in ut_constants.const.name]
    lind = np.array([names.index(c) for c in cons])
    freq = np.asarray(ut_constants.const.freq)[lind]          # cycles/hour
    t0, tn = _utide_time(start), _utide_time(nodal_date)
    F, U, _ = harmonics.FUV(np.array([tn]), tn, lind, lat, [0, 0, 0, 0])
    _, _, V = harmonics.FUV(np.array([t0]), t0, lind, lat, [0, 0, 0, 0])
    return 3600.0 / freq, F[0], U[0], V[0], freq


# ---------------------------------------------------------------- atlas
def _atlas_file(d, c):
    pat = re.compile(rf'^{c.lower()}(\.|_).*nc$')
    hits = sorted(f for f in os.listdir(d) if pat.match(f))
    if not hits:
        raise SystemExit(f'no {c} file in {d}')
    return os.path.join(d, hits[0])


def atlas_at(d, comp, c, lat, lon):
    """complex constants A exp(iG) (m/s) at points; land filled from nearest sea cell"""
    ds = xr.open_dataset(_atlas_file(d, c))
    a = ds[f'{comp}a'].values.astype('f8')
    g = ds[f'{comp}g'].values.astype('f8')
    la, lo = ds['lat'].values.astype('f8'), ds['lon'].values.astype('f8') % 360
    ds.close()
    z = a * 0.01 * np.exp(1j * np.radians(g))                   # cm/s -> m/s
    bad = ~np.isfinite(z) | (np.abs(a) > 1e5)
    if bad.all():
        raise SystemExit(f'{c} {comp}: atlas all missing')
    if bad.any():
        idx = ndimage.distance_transform_edt(bad, return_distances=False, return_indices=True)
        z = z[tuple(idx)]
    o = np.argsort(lo)
    lo, z = lo[o], z[:, o]
    # periodic in longitude if the atlas is global in lon
    if lo[-1] - lo[0] > 359.0 - 2 * np.diff(lo).max():
        lo = np.r_[lo[-1] - 360, lo, lo[0] + 360]
        z = np.c_[z[:, -1:], z, z[:, :1]]
    plon = lon % 360
    if lat.min() < la.min() or lat.max() > la.max():
        raise SystemExit(f'{c}: points lat {lat.min():.2f}..{lat.max():.2f} outside atlas '
                         f'{la.min():.2f}..{la.max():.2f}')
    pts = np.c_[lat, plon]
    f = lambda w: RegularGridInterpolator((la, lo), w, bounds_error=True)(pts)
    return f(z.real) + 1j * f(z.imag)


# ---------------------------------------------------------------- self-test
def selftest(cons, start, nodal_date, lat):
    """MITgcm form Am cos(w t - 2 pi Ph/T) must equal UTide's nodal-corrected
    reconstruction of (A, G) over 2 days from start (checks the phase convention; the
    nodal f,u drift over 2 d is < 5e-4 of the amplitude)."""
    import utide
    T, f, u, V0, freq = nodal(cons, start, start, lat)          # f,u at start for the test
    t = np.arange(0, 2 * 24 * 6) / 144.0                          # 2 days, 10-min steps
    tt = t + 0.0
    rng = np.random.default_rng(1)
    A = rng.uniform(0.05, 0.5, len(cons))
    G = rng.uniform(0, 360, len(cons))
    s = np.cos(2 * np.pi * np.arange(0, 30 * 24) / 24.0 / 0.5)
    tfit = np.arange(0, 30 * 24) / 24.0
    coef = utide.solve(tfit, s, lat=lat, constit=cons, method='ols', conf_int='none',
                       nodal=True, trend=False, verbose=False, epoch=start)
    order = [cons.index(n) for n in coef['name']]
    coef['A'] = A[order]
    coef['g'] = G[order]
    coef['mean'] = 0.0
    coef['slope'] = 0.0
    ref = utide.reconstruct(tt, coef, verbose=False, min_SNR=0, min_PE=0, epoch=start).h
    Am = f * A
    Ph = (np.mod(G - 360.0 * (V0 + u), 360.0) / 360.0) * T
    ts = t * 86400.0
    mit = sum(Am[k] * np.cos(2 * np.pi * ts / T[k] - 2 * np.pi * Ph[k] / T[k])
              for k in range(len(cons)))
    err = np.abs(mit - ref).max() / np.abs(ref).max()
    print(f'SELFTEST max |MITgcm - UTide| / max|UTide| over 2 d = {err:.2e}')
    if err > 1e-3:
        raise SystemExit('SELFTEST FAIL')
    print('SELFTEST PASS')


# ---------------------------------------------------------------- main
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--points', required=True)
    ap.add_argument('--u-dir', required=True)
    ap.add_argument('--v-dir', required=True)
    ap.add_argument('--start', required=True, help='run startDate, YYYY-MM-DD[THH:MM]')
    ap.add_argument('--nodal-date', default=None)
    ap.add_argument('--cons', default=','.join(CONS))
    ap.add_argument('--out-dir', required=True)
    ap.add_argument('--prec', type=int, default=64, choices=(32, 64))
    ap.add_argument('--selftest', action='store_true', help='only run the phase self-test')
    a = ap.parse_args()
    cons = [c.strip().upper() for c in a.cons.split(',')]
    if len(cons) > 10:
        raise SystemExit('OBCS_tideCompSize = 10')
    P = np.load(a.points)
    bnds = [b for b in BND if f'{b}_lat' in P.files]
    lat_ref = float(np.mean(np.concatenate([P[f'{b}_lat'][P[f'{b}_has']] for b in bnds])))
    nod = a.nodal_date or str(np.datetime64(a.start[:10]) + np.timedelta64(182, 'D'))
    selftest(cons, a.start, a.start, lat_ref)
    if a.selftest:
        return
    T, f, u, V0, _ = nodal(cons, a.start, nod, lat_ref)
    print(f'start {a.start}  nodal f,u at {nod}  lat_ref {lat_ref:.2f}')
    for k, c in enumerate(cons):
        print(f'  {c:3s} T = {T[k]:.6f} s  f = {f[k]:.4f}  u = {360 * u[k]:7.2f} deg  '
              f'V0 = {360 * (V0[k] % 1):7.2f} deg')
    os.makedirs(a.out_dir, exist_ok=True)
    dt = '>f8' if a.prec == 64 else '>f4'
    snippet = []
    for b in bnds:
        lat, lon = P[f'{b}_lat'].astype('f8'), P[f'{b}_lon'].astype('f8')
        cs, sn, has = P[f'{b}_cs'].astype('f8'), P[f'{b}_sn'].astype('f8'), P[f'{b}_has']
        r = np.hypot(cs, sn)
        r[r == 0] = 1.0
        cs, sn = cs / r, sn / r
        n = lat.size
        out = {k: np.zeros((len(cons), n)) for k in ('uAm', 'uPh', 'vAm', 'vPh')}
        for k, c in enumerate(cons):
            ZE = np.zeros(n, complex)
            ZN = np.zeros(n, complex)
            ZE[has] = atlas_at(a.u_dir, 'U', c, lat[has], lon[has])
            ZN[has] = atlas_at(a.v_dir, 'V', c, lat[has], lon[has])
            Zx = cs * ZE + sn * ZN
            Zy = -sn * ZE + cs * ZN
            for comp, Z in (('u', Zx), ('v', Zy)):
                A = np.abs(Z)
                G = np.degrees(np.angle(Z))
                ph = np.mod(G - 360.0 * (V0[k] + u[k]), 360.0) / 360.0 * T[k]
                out[f'{comp}Am'][k] = np.where(has, f[k] * A, 0.0)
                out[f'{comp}Ph'][k] = np.where(has, ph, 0.0)
        X = BND[b]
        for key, arr in out.items():
            comp, kind = key[0], key[1:]
            fn = f'OB{X}_{comp}Tid{kind}.bin'
            arr.astype(dt).tofile(os.path.join(a.out_dir, fn))
            snippet.append(f" OB{X}_{comp}Tid{kind}File = '{fn}',")
        m2 = cons.index('M2') if 'M2' in cons else 0
        print(f'{b:5s}: {has.sum():5d} OB entries | {cons[m2]} |u| median {np.median(out["uAm"][m2][has]):.3f} '
              f'max {out["uAm"][m2][has].max():.3f} m/s, |v| median {np.median(out["vAm"][m2][has]):.3f} '
              f'max {out["vAm"][m2][has].max():.3f} m/s')
    per = ', '.join(f'{x:.6f}' for x in T)
    head = [' useOBCStides = .TRUE.,', f' OBCS_tidalPeriod = {per},',
            f'# constituents {",".join(cons)}; phase zero = startDate {a.start}; '
            f'nodal f,u at {nod}; gen_obcs_tides.py']
    with open(os.path.join(a.out_dir, 'data.obcs_tides.snippet'), 'w') as fh:
        fh.write('\n'.join(head + snippet) + '\n')
    print('wrote', os.path.join(a.out_dir, 'data.obcs_tides.snippet'))
    print('GEN_TIDES_DONE')


if __name__ == '__main__':
    sys.exit(main())
