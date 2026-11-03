#!/usr/bin/env python3
"""Initial and open-boundary conditions for the Malibu nest from the CCS parent.

Parent: ecco_darwin regions/CCS v05/v05_kelp, 80x100x48, xgOrigin=-131.66667,
dx=1/3 deg, ygOrigin=19.8998, delYFile (100 values), 48 levels.

Two parent sources:
  --mode init    one snapshot per field (CCS run_template/init/*_80x100x48.16-Jan-1992).
                 IC from it; OBCS = the same snapshot repeated (2 records), i.e.
                 steady boundaries for a short test.
  --mode series  per-field files of nrec parent snapshots (e.g. monthly or daily
                 output of the CCS run, converted to UVEL/VVEL, not *MASS), named
                 <dir>/<FIELD>.bin; IC from record 1, OBCS from all records.

Method: fill parent land by nearest wet value on each level, bilinear in
lon/lat, linear in depth (constant above the first / below the last parent
level). Normal velocity at each open boundary has its net transport removed
(the model also does this with useOBCSbalance). Kelp tracers (Bmag, Cmag,
Qmag) start in a depth-defined kelp mask; TSS starts at zero. CDOM (darwin
DARWIN_ALLOW_CDOM, P units) is not in the parent: IC and OBCS are a fixed
fraction of the parent DOP (see CDOM_* below), taken out of DOP/DON/DOC/DOFe
with darwin's CDOM stoichiometry so P, N, C and Fe are unchanged.

Writes big-endian float32 to ../work/{init,obcs}/ and prints the PTRACERS /
OBCS file lists. Usage:
  python3 make_ic_obcs.py --mode init --parent ~/Library/CloudStorage/Box-Box/KelpProject/CCS/run_template

--mode socal (2026-10-10): physics only (THETA, SALT, U, V) from the socal 2-km nest's daily means,
cut by extract_socal_box.py (box.json, hFacC.bin, {THETA,SALT,uVel_C,vVel_C}.bin; U/V east/north at
cell centers). IC = the daily mean of --start; OBCS = one record per day, centered at 12:00, with the
first record repeated at --start - 12 h, so the OBCS startdate is (start - 1 d) 12:00 and period
86400 s (no record 0 needed at myTime = 0). If the box has monthly BGC (extract_socal_box.py, monthly_mid),
the 31 darwin tracers are interpolated linearly in time between the monthly means (constant before the
first / after the last mid-month) to the same daily records, then CDOM, kelp and TSS as in the other
modes; otherwise only the physics is written.
  python3 make_ic_obcs.py --mode socal --parent ../work/socal_box_2025 --start 2025-01-01 --days 180 --out ../work/socal2025
"""
import argparse
import json
import os

import numpy as np
from scipy import ndimage
from scipy.interpolate import RegularGridInterpolator

HERE = os.path.dirname(os.path.abspath(__file__))
WORK = os.path.join(HERE, '..', 'work')
OUT = [WORK]      # output dir (init/, obcs/, obcs_files.json); --out
PAD_END = [False]  # socal mode: repeat the last daily record once at the end (--pad-end)

# ---- parent (CCS v05) ------------------------------------------------------
PNX, PNY, PNR = 80, 100, 48
PXG0, PDX, PYG0 = -131.66667, 1.0 / 3.0, 19.8998
PDELR = np.array([
    10.00, 10.00, 10.00, 10.00, 10.00, 10.00, 10.00, 10.01,
    10.03, 10.11, 10.32, 10.80, 11.76, 13.42, 16.04, 19.82, 24.85,
    31.10, 38.42, 46.50, 55.00, 63.50, 71.58, 78.90, 85.15, 90.18,
    93.96, 96.58, 98.25, 99.25, 100.01, 101.33, 104.56, 111.33, 122.83,
    139.09, 158.94, 180.83, 203.55, 226.50, 249.50, 272.50, 295.50, 318.50,
    341.50, 364.50, 387.50, 410.50])

# parent tracers in parent order; the nest order (data.ptracers, darwin
# DARWIN_INDICES.h with DARWIN_ALLOW_CDOM) is NEST_TRACERS below:
# 1-19 parent, 20 CDOM, 21-32 parent plankton, 33-35 kelp, 36 TSS
PARENT_TRACERS = ['DIC', 'NO3', 'NO2', 'NH4', 'PO4', 'FeT', 'SiO2', 'DOC', 'DON',
                  'DOP', 'DOFe', 'POC', 'PON', 'POP', 'POFe', 'POSi', 'PIC', 'ALK',
                  'O2', 'c1', 'c2', 'c3', 'c4', 'c5', 'c6', 'c7',
                  'Chl1', 'Chl2', 'Chl3', 'Chl4', 'Chl5']
LOCAL_TRACERS = ['Bmag', 'Cmag', 'Qmag', 'TSS']
NEST_TRACERS = PARENT_TRACERS[:19] + ['CDOM'] + PARENT_TRACERS[19:] + LOCAL_TRACERS

# ---- CDOM from DOP (option 1, 2026-10-06) ---------------------------------
# CDOM = f DOP with f fixed so that offshore surface CDOM absorption is
# a_CDOM(443) = CDOM_A443 for the median surface DOP of the deep (> 500 m) part
# of the nest. darwin: a_CDOM(l) = CDOMcoeff CDOM exp(-Sdom (l - 450)),
# CDOMcoeff = 100 m-1/(mmol P m-3), Sdom = 0.014 nm-1 (darwin3 24885b71 code
# defaults, darwin_read_params.F; not the Dutkiewicz et al. 2015 values).
# CDOM_A443 = 0.045 m-1: midpoint of the measured Santa Monica Bay surface
# a_g(443) on the dry-weather dates of 5 and 11 Nov 2015, 0.028-0.064 m-1
# (Sept-Nov 2015 overall 0.015-0.32; Harringmeyer et al. 2021, Front. Environ.
# Sci. 9:647966, Table 1). Was 0.02 (low end) until 2026-10-06. The parent
# (CCS v05 / ECCO-Darwin v05 V5r1) has no CDOM tracer, so this sets the IC and
# the open-boundary CDOM. PACE is kept as an independent check, not used to set f.
CDOM_A443 = 0.045
CDOM_COEFF, CDOM_SDOM = 100.0, 0.014
R_NP_CDOM, R_CP_CDOM, R_FEP_CDOM = 16.0, 120.0, 1.0e-3   # darwin3 defaults

# ---- kelp initial state (giant kelp, G-MACMODS parameters in data.darwin) --
KELP_ZMIN, KELP_ZMAX = 5.0, 25.0   # bed depth range [m]
KELP_B0 = 100.0    # initial biomass in the bottom cell [gDW/m3]; see note below
KELP_Q0 = 20.0     # N quota [mgN/gDW] (Qmin 10.18, Qmax 54)
KELP_CARBON = 0.3  # mp_spp_carbon [gC/gDW]
TSS0 = 0.0         # initial TSS [g/m3]: none; sediment only from Malibu Creek and
                   # resuspension (a uniform 1 g/m3 made Rrs 3-4x too high vs PACE)
OB_EDGE_GAP = 3    # no kelp within this many cells of an open boundary
HFACMIN, HFACMINDR = 0.3, 1.0   # must match hFacMin / hFacMinDr in input/data


def parent_grid(parent_dir):
    dely = np.fromfile(os.path.join(parent_dir, 'delYFile'), '>f4').astype(float)
    assert dely.size == PNY, dely.size
    yg = PYG0 + np.r_[0.0, np.cumsum(dely)]
    xg = PXG0 + PDX * np.arange(PNX + 1)
    zf = -np.r_[0.0, np.cumsum(PDELR)]
    return dict(xg=xg, yg=yg, xc=0.5 * (xg[1:] + xg[:-1]), yc=0.5 * (yg[1:] + yg[:-1]),
                zc=0.5 * (zf[1:] + zf[:-1]))


def nest_grid():
    with open(os.path.join(WORK, 'grid_info.json')) as f:
        g = json.load(f)
    for k in ('xc', 'yc', 'xg', 'yg', 'delr'):
        g[k] = np.array(g[k])
    zf = -np.r_[0.0, np.cumsum(g['delr'])]
    g['zc'] = 0.5 * (zf[1:] + zf[:-1])
    g['depth'] = -np.fromfile(os.path.join(WORK, g['bathy']), '>f4').reshape(g['ny'], g['nx'])
    g['dz'] = g['delr']
    return g


def fill_land(f3, wet3):
    """Nearest-wet fill on each level; levels with no wet copy the level above."""
    out = np.array(f3, dtype=float)
    for k in range(out.shape[0]):
        w = wet3[k]
        if w.any():
            _, (jj, ii) = ndimage.distance_transform_edt(~w, return_indices=True)
            out[k] = out[k][jj, ii]
        elif k > 0:
            out[k] = out[k - 1]
    return out


def to_points(f3, px, py, pz, x, y, z):
    """f3 (PNR, PNY, PNX) on parent (pz, py, px) -> values at points x, y (1-D,
    same length) and depths z (1-D): returns (len(z), len(x))."""
    hz = np.empty((f3.shape[0], x.size))
    for k in range(f3.shape[0]):
        hz[k] = RegularGridInterpolator((py, px), f3[k], bounds_error=False,
                                        fill_value=None)(np.c_[y, x])
    # pz decreasing (negative down): interpolate on -pz
    out = np.empty((z.size, x.size))
    for n in range(x.size):
        out[:, n] = np.interp(-z, -pz, hz[:, n])
    return out


def to_volume(f3, px, py, pz, xs, ys, z):
    X, Y = np.meshgrid(xs, ys)
    v = to_points(f3, px, py, pz, X.ravel(), Y.ravel(), z)
    return v.reshape(z.size, ys.size, xs.size)


def read_parent(path, nrec=None):
    a = np.fromfile(path, '>f4')
    n = a.size // (PNX * PNY * PNR)
    a = a[:n * PNX * PNY * PNR].reshape(n, PNR, PNY, PNX)
    return a if nrec is None else a[:nrec]


def boundary_points(g, side, kind):
    """(x, y) of OB points. kind: 't' tracer/tangential row, 'u' / 'v' faces.
    MITgcm OBCS: W normal u at i=OB_Iwest+1 face, E at OB_Ieast face,
    S normal v at j=OB_Jsouth+1 face, N at OB_Jnorth face; others at cell centres."""
    xc, yc, xg, yg = g['xc'], g['yc'], g['xg'], g['yg']
    if side in ('west', 'east'):
        i = 0 if side == 'west' else -1
        if kind == 'u':
            x0 = xg[1] if side == 'west' else xg[-2]
            return np.full(yc.size, x0), yc
        if kind == 'v':
            return np.full(yc.size, xc[i]), yg[:-1]
        return np.full(yc.size, xc[i]), yc
    j = 0 if side == 'south' else -1
    if kind == 'v':
        y0 = yg[1] if side == 'south' else yg[-2]
        return xc, np.full(xc.size, y0)
    if kind == 'u':
        return xg[:-1], np.full(xc.size, yc[j])
    return xc, np.full(xc.size, yc[j])


def hfac_c(g):
    """hFacC as MITgcm ini_masks_etc computes it (no shelfice): a partial
    bottom cell thinner than half of hFacMnSz is dropped (land), otherwise
    it is at least hFacMnSz."""
    dz = np.asarray(g['dz'], float)
    zf = np.r_[0.0, np.cumsum(dz)]
    d = g['depth']
    h = np.clip((d[None] - zf[:-1, None, None]) / dz[:, None, None], 0.0, 1.0)
    mn = np.maximum(HFACMIN, np.minimum(HFACMINDR / dz, 1.0))[:, None, None]
    return np.where(h < 0.5 * mn, 0.0, np.maximum(h, mn) * (h > 0))


def boundary_wet(g, side):
    d = g['depth']
    line = {'west': d[:, 0], 'east': d[:, -1], 'south': d[0], 'north': d[-1]}[side]
    zf_top = np.r_[0.0, np.cumsum(g['dz'])][:-1]
    return (line[None, :] > zf_top[:, None])          # (nr, n) crude wet mask


def balance(ob, g, side):
    """Remove net transport through one boundary (each record)."""
    wet = boundary_wet(g, side)
    ds = (np.array(g['dely']) * 111.2e3 if side in ('west', 'east')
          else np.full(g['nx'], 111.2e3 * np.cos(np.deg2rad(34.0)) * g['dlon']))
    area = wet * g['dz'][:, None] * ds[None, :]
    for t in range(ob.shape[0]):
        ob[t] -= (ob[t] * area).sum() / area.sum() * wet
    return ob


def v5r1_grid(d):
    """Parent grid of the V5r1 LLC270 box (extract_v5r1_box.py): 1/3 deg lon,
    LLC270 latitudes, 50 levels; wet = hFacC > 0."""
    box = json.load(open(os.path.join(d, 'box.json')))
    lon, lat = np.array(box['lon']), np.array(box['lat'])
    rc = np.array(box['RC'])
    hf = np.fromfile(os.path.join(d, 'hFacC.bin'), '>f4').reshape(rc.size, lat.size, lon.size)
    return dict(xc=lon, yc=lat, zc=rc, wet=hf > 0), box


def v5r1_series(d, box, name, start, days):
    """Monthly means -> records at start + 30 d n (n = 0..), linear in time
    between month centres (constant beyond the first/last month)."""
    nz, ny, nx = len(box['RC']), len(box['lat']), len(box['lon'])
    F = np.fromfile(os.path.join(d, f'{name}.bin'), '>f4').reshape(-1, nz, ny, nx).astype(float)
    # run2 'monthly' files are calendar-month means written at month end; the
    # meta timeInterval start is not reliable, so use the calendar mid-month of
    # the month that the record's interval centre falls in
    mon = np.array([np.datetime64(t[:7], 'M') for t in box['monthly_mid']])
    tm = (mon.astype('datetime64[s]') + ((mon + np.timedelta64(1, "M")).astype('datetime64[s]')
                                         - mon.astype('datetime64[s]')) / 2)
    t0 = np.datetime64(start, 's')
    nrec = int(np.ceil(days / 30.0)) + 2
    tt = t0 + np.arange(nrec) * np.timedelta64(30 * 86400, 's')
    x = (tm - t0).astype(float)
    out = np.empty((nrec, nz, ny, nx))
    for n, t in enumerate((tt - t0).astype(float)):
        j = np.clip(np.searchsorted(x, t) - 1, 0, len(x) - 2)
        w = np.clip((t - x[j]) / (x[j + 1] - x[j]), 0.0, 1.0)
        out[n] = (1 - w) * F[j] + w * F[j + 1]
    return out


def derive_cdom(g, tag, sides, obfiles):
    """CDOM = f DOP for the IC and OBCS; remove the same P, N, C, Fe from
    DOP, DON, DOC, DOFe (clipped at zero)."""
    rd = lambda fn: np.fromfile(fn, '>f4').astype(float)
    ic = lambda n: os.path.join(OUT[0], 'init', f'{n}_{tag}')
    dop = rd(ic('DOP')).reshape(g['nr'], g['ny'], g['nx'])
    deep = (g['depth'] > 500) & (dop[0] > 0)
    dop_ref = float(np.median(dop[0][deep]))
    f = CDOM_A443 / (CDOM_COEFF * dop_ref * np.exp(-CDOM_SDOM * (443 - 450)))
    ratios = {'DOP': 1.0, 'DON': R_NP_CDOM, 'DOC': R_CP_CDOM, 'DOFe': R_FEP_CDOM}

    def split(paths_in):
        cd = f * rd(paths_in['DOP'])
        for n, rr in ratios.items():
            v = rd(paths_in[n])
            np.maximum(v - rr * cd, 0.0).astype('>f4').tofile(paths_in[n])
        return cd

    split_ic = {n: ic(n) for n in ratios}
    cd = split(split_ic)
    cd.astype('>f4').tofile(ic('CDOM'))
    files = {}
    for side in sides:
        ob = lambda n: os.path.join(OUT[0], 'obcs', obfiles[n][side])
        cd_ob = split({n: ob(n) for n in ratios})
        fn = f'malibu_{tag}_CDOM_{side.capitalize()}'
        cd_ob.astype('>f4').tofile(os.path.join(OUT[0], 'obcs', fn))
        files[side] = fn
    obfiles['CDOM'] = files
    print(f'  CDOM  = {f:.3e} x DOP (surface DOP {dop_ref:.3f} mmol P/m3 offshore -> '
          f'a_CDOM(443) {CDOM_A443} m-1); ic [{cd.min():.3g}, {cd.max():.3g}]')
    json.dump(dict(f=f, dop_ref=dop_ref, a443=CDOM_A443),
              open(os.path.join(OUT[0], 'cdom_ic.json'), 'w'), indent=1)


def socal_grid(d):
    """socal box (extract_socal_box.py): 1/48 deg lon/lat, 60 levels; wet = hFacC > 0."""
    box = json.load(open(os.path.join(d, 'box.json')))
    lon, lat, rc = np.array(box['lon']), np.array(box['lat']), np.array(box['RC'])
    hf = np.fromfile(os.path.join(d, 'hFacC.bin'), '>f4').reshape(rc.size, lat.size, lon.size)
    return dict(xc=lon, yc=lat, zc=rc, wet=hf > 0), box


def socal_series(d, box, name, start, days):
    """Daily means from --start for --days, with the first one repeated in front (see docstring)."""
    nz, ny, nx = len(box['RC']), len(box['lat']), len(box['lon'])
    F = np.fromfile(os.path.join(d, f'{name}.bin'), '>f4').reshape(-1, nz, ny, nx)
    k0 = [m[:10] for m in box['daily_mid']].index(start)
    n = int(days)
    assert k0 + n <= F.shape[0], f'socal box has {F.shape[0] - k0} days from {start}, asked {n}'
    F = np.concatenate([F[k0:k0 + 1], F[k0:k0 + n]] + ([F[k0 + n - 1:k0 + n]] if PAD_END[0] else []))
    return F.astype(float)


def socal_monthly_series(d, box, name, start, days):
    """Monthly means (mid-month times) -> the daily OBCS record times of socal_series (+ the IC at start)."""
    nz, ny, nx = len(box['RC']), len(box['lat']), len(box['lon'])
    F = np.fromfile(os.path.join(d, f'{name}.bin'), '>f4').reshape(-1, nz, ny, nx).astype(float)
    tm = np.array([np.datetime64(t, 's') for t in box['monthly_mid']])
    t0 = np.datetime64(start, 's')
    tt = np.r_[t0 - np.timedelta64(12 * 3600, 's'),
               t0 + np.timedelta64(12 * 3600, 's') + np.arange(int(days) + int(PAD_END[0])) * np.timedelta64(86400, 's')]
    x = (tm - t0).astype(float)
    out = np.empty((tt.size, nz, ny, nx))
    for n, t in enumerate((tt - t0).astype(float)):
        if x.size == 1 or t <= x[0]:
            out[n] = F[0]
        elif t >= x[-1]:
            out[n] = F[-1]
        else:
            j = np.searchsorted(x, t) - 1
            w = (t - x[j]) / (x[j + 1] - x[j])
            out[n] = (1 - w) * F[j] + w * F[j + 1]
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--mode', choices=('init', 'series', 'v5r1', 'socal'), default='init')
    ap.add_argument('--parent', required=True,
                    help='run_template dir (init mode) or dir of <FIELD>.bin series')
    ap.add_argument('--grid-parent', default=None,
                    help='dir with parent delYFile (default: --parent)')
    ap.add_argument('--nrec', type=int, default=None)
    ap.add_argument('--start', default=None,
                    help='v5r1: nest start date YYYY-MM-DD (IC time, first OBCS record)')
    ap.add_argument('--days', type=float, default=105.,
                    help='v5r1: run length [d]; OBCS records every 30 d to cover it')
    ap.add_argument('--out', default=None, help='output dir (default ../work)')
    ap.add_argument('--pad-end', action='store_true',
                    help='socal: one more daily record at the end (copy of the last day), so a run of --days '
                         'never asks for a record past the file: run_2025h1 (181 records, 180 d) stopped at '
                         'day 179.5 reading record 182: exf_getffieldrec.F interpolates between records count0 = '
                         'INT(t/period)+1 and count0+1 (t from the 12:00 startdate), so after 12:00 of the last '
                         'day it needs one record past the last daily mean')
    a = ap.parse_args()
    PAD_END[0] = a.pad_end
    if a.out:
        OUT[0] = a.out

    g = nest_grid()
    if a.mode == 'v5r1':
        p, box = v5r1_grid(a.parent)
    elif a.mode == 'socal':
        p, box = socal_grid(a.parent)
    else:
        p = parent_grid(a.grid_parent or a.parent)
    os.makedirs(os.path.join(OUT[0], 'init'), exist_ok=True)
    os.makedirs(os.path.join(OUT[0], 'obcs'), exist_ok=True)
    tag = f"{g['nx']}x{g['ny']}x{g['nr']}"
    sides = [s for s in ('west', 'south', 'east') if g['open'][s]]

    def src(name):
        if a.mode == 'v5r1':
            return v5r1_series(a.parent, box, name, a.start, a.days)
        if a.mode == 'socal':
            if name in ('THETA', 'SALT', 'uVel_C', 'vVel_C'):
                return socal_series(a.parent, box, name, a.start, a.days)
            return socal_monthly_series(a.parent, box, name, a.start, a.days)
        if a.mode == 'init':
            fn = os.path.join(a.parent, 'init', f'{name}_80x100x48.16-Jan-1992')
            return read_parent(fn)                      # (1, nr, ny, nx)
        return read_parent(os.path.join(a.parent, f'{name}.bin'), a.nrec)

    S = src('SALT')
    wetT = S[0] > 0
    wetU = wetT & np.roll(wetT, 1, axis=2)              # parent U at west faces
    wetV = wetT & np.roll(wetT, 1, axis=1)
    pz = p['zc']
    if a.mode in ('v5r1', 'socal'):   # east/north velocities at cell centres
        wetT = p['wet']
        wetU = wetV = wetT

    def do_field(name, wet, px, py, kind_ic, out_name, clip=False):
        F = src(name)
        if clip:   # parent BGC ICs hold small negatives (e.g. DOC, c1-c7, Chl)
            F = np.maximum(F, 0.0)
        nrec = F.shape[0]
        recs = [fill_land(F[t], wet) for t in range(nrec)]
        # IC
        xs = {'t': g['xc'], 'u': g['xg'][:-1], 'v': g['xc']}[kind_ic]
        ys = {'t': g['yc'], 'u': g['yc'], 'v': g['yg'][:-1]}[kind_ic]
        ic = to_volume(recs[1 if a.mode == 'socal' else 0], px, py, pz, xs, ys, g['zc'])
        ic.astype('>f4').tofile(os.path.join(OUT[0], 'init', f'{out_name}_{tag}'))
        # OBCS (steady test: repeat the snapshot)
        kind_ob = {'t': 't', 'u': 'u', 'v': 'v'}[kind_ic]
        files = {}
        for side in sides:
            x, y = boundary_points(g, side, kind_ob)
            ob = np.stack([to_points(r, px, py, pz, x, y, g['zc']) for r in recs])
            if (kind_ob == 'u' and side in ('west', 'east')) or \
               (kind_ob == 'v' and side in ('south', 'north')):
                ob = balance(ob, g, side)
            if nrec == 1:
                ob = np.concatenate([ob, ob])
            fn = f'malibu_{tag}_{out_name}_{side.capitalize()}'
            ob.astype('>f4').tofile(os.path.join(OUT[0], 'obcs', fn))
            files[side] = fn
        print(f'  {out_name:5s} ic [{ic.min():.4g}, {ic.max():.4g}]  obcs recs {max(nrec, 2)}')
        return files

    obfiles = {}
    obfiles['THETA'] = do_field('THETA', wetT, p['xc'], p['yc'], 't', 'THETA')
    obfiles['SALT'] = do_field('SALT', wetT, p['xc'], p['yc'], 't', 'SALT')
    if a.mode in ('v5r1', 'socal'):
        obfiles['U'] = do_field('uVel_C', wetU, p['xc'], p['yc'], 'u', 'U')
        obfiles['V'] = do_field('vVel_C', wetV, p['xc'], p['yc'], 'v', 'V')
    else:
        obfiles['U'] = do_field('U', wetU, p['xg'][:-1], p['yc'], 'u', 'U')
        obfiles['V'] = do_field('V', wetV, p['xc'], p['yg'][:-1], 'v', 'V')
    if a.mode == 'socal' and not box.get('monthly_mid'):
        np.zeros((g['ny'], g['nx'])).astype('>f4').tofile(os.path.join(OUT[0], 'init', f"ETAN_{g['nx']}x{g['ny']}"))
        st = np.datetime64(a.start, 's') - np.timedelta64(12 * 3600, 's')
        with open(os.path.join(OUT[0], 'obcs_files.json'), 'w') as f:
            json.dump(dict(tag=tag, sides=sides, files=obfiles, source=a.parent, physics_only=True,
                           obcs_startdate=str(st), obcs_period=86400.0, nrec=int(a.days) + 1 + int(a.pad_end)), f, indent=1)
        print('physics only; OBCS startdate', st, 'period 86400 s; wrote', os.path.join(OUT[0], 'obcs_files.json'))
        return
    for name in PARENT_TRACERS:
        obfiles[name] = do_field(name, wetT, p['xc'], p['yc'], 't', name, clip=True)
    derive_cdom(g, tag, sides, obfiles)
    if a.mode == 'socal':
        print(f'socal mode: BGC from monthly means {box["monthly_mid"]}; OBCS startdate',
              np.datetime64(a.start, 's') - np.timedelta64(12 * 3600, 's'), 'period 86400 s')
    eta = np.zeros((g['ny'], g['nx']))
    eta.astype('>f4').tofile(os.path.join(OUT[0], 'init', f"ETAN_{g['nx']}x{g['ny']}"))

    # ---- kelp + TSS initial fields -----------------------------------------
    # The kelp model holds Bmag per cell volume and grows from the bottom
    # (biomass_down = bottom cell), so the initial stock sits in the bottom
    # wet cell of each kelp column. KELP_B0 is a placeholder: replace it with
    # KelpWatch (Landsat) canopy biomass mapped to the grid before science runs.
    d = g['depth']
    mask = (d >= KELP_ZMIN) & (d <= KELP_ZMAX)
    gap = OB_EDGE_GAP
    edge = np.zeros_like(mask)
    if g['open']['west']:
        edge[:, :gap] = True
    if g['open']['east']:
        edge[:, -gap:] = True
    if g['open']['south']:
        edge[:gap, :] = True
    mask &= ~edge
    zf = np.r_[0.0, np.cumsum(g['dz'])]
    # bottom wet cell = MITgcm kLowC (a raw-depth bottom cell can be land)
    hf = hfac_c(g)
    kbot = np.clip((hf > 0).sum(axis=0) - 1, 0, g['nr'] - 1)
    mask &= hf[0] > 0
    B = np.zeros((g['nr'], g['ny'], g['nx']))
    jj, ii = np.where(mask)
    B[kbot[jj, ii], jj, ii] = KELP_B0
    Q = np.where(B > 0, KELP_Q0, 0.0)
    C = KELP_CARBON * B
    T = np.where(hf > 0, TSS0, 0.0)
    for name, arr in zip(LOCAL_TRACERS, (B, C, Q, T)):
        arr.astype('>f4').tofile(os.path.join(OUT[0], 'init', f'{name}_{tag}'))
    print(f'  kelp columns {mask.sum()} (depth {KELP_ZMIN}-{KELP_ZMAX} m, '
          f'{gap}-cell gap from open boundaries)')

    with open(os.path.join(OUT[0], 'obcs_files.json'), 'w') as f:
        json.dump(dict(tag=tag, sides=sides, files=obfiles,
                       tracers=NEST_TRACERS), f, indent=1)
    print('wrote', os.path.join(OUT[0], 'obcs_files.json'))


if __name__ == '__main__':
    main()
