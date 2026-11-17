#!/usr/bin/env python3
"""gen_humboldt_plumes.py -- iceplume inputs for Humboldt Glacier (Sermersuaq) in the
Nares Strait grid, from Mankoff et al. Greenland liquid water discharge v8.0
(doi:10.22008/FK2/XKQVL7, ice-margin outlets, daily, MAR).

Decisions (model_setup/nares_strait/DECISIONS.md, 2026-10-06):
  * baseline: the 59 marine-terminating Humboldt outlets are lumped into the 5 that
    carry the most discharge (95% of peak discharge); every other outlet's daily
    discharge is added to its nearest of the five, so the daily total is conserved.
    The five become cone plumes (plumeMask +-3); all other Humboldt front cells are
    melt-only (+-1).
  * sensitivity 'sheet': the daily total spread along all Humboldt front cells in
    proportion to front length, as sheet plumes (+-2).
  * sensitivity 'nodischarge': front cells melt-only, no discharge.
  * MAR v8 ends 2025-06-11: 2025-06-12..2025-12-31 is filled per outlet with the
    2015-2024 mean for that day of year.
Mask sign (iceplume_calc.F): > 0 -> front normal along grid x (dxG), < 0 -> along y.

Outputs (in <out>/): Humboldt_plumes_{baseline,sheet,nodischarge}_mask.bin and
_length.bin (real*8, ny x nx), Humboldt_plume_cells.csv, Humboldt_Qsg_daily.npz
(dates, per-cell daily discharge for each case) -- expanded to yearly 2-D
runoffQsg files on pfe when the run is staged (~590 MB/yr as real*4).

Usage: python3 gen_humboldt_plumes.py <grid dir> <mankoff v8 dir> <out dir>
"""
import os
import sys
import json
import datetime as dt
import numpy as np
import pandas as pd
import netCDF4 as nc
from scipy import ndimage

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import gen_nares_mitgrid as G

YEARS = range(2015, 2026)
NTOP = 5
SEARCH = 3      # cells: move each lumped plume to the deepest front cell within this radius
END = dt.date(2025, 12, 31)


def main(gdir, mdir, out):
    m = json.load(open(os.path.join(gdir, 'Nares_grid_meta.json')))
    nx, ny = m['nx'], m['ny']
    g = np.fromfile(os.path.join(gdir, 'Nares.mitgrid'), '>f8').reshape(16, ny + 1, nx + 1)
    xC, yC = g[0][:-1, :-1], g[1][:-1, :-1]
    dxG, dyG = g[14][:ny, :nx], g[15][:ny, :nx]
    wet = np.fromfile(os.path.join(gdir, 'Nares_bathymetry.bin'), '>f8').reshape(ny, nx) < 0
    depth = -np.fromfile(os.path.join(gdir, 'Nares_bathymetry.bin'), '>f8').reshape(ny, nx)
    e1, e2, n = (np.array(m[k]) for k in ('e1', 'e2', 'pole'))
    lamG, _ = G.to_rot(g[5][0, :], g[6][0, :], e1, e2, n)
    _, phiG = G.to_rot(g[5][:, 0], g[6][:, 0], e1, e2, n)
    d = np.mean(np.diff(lamG))

    o = pd.read_csv(os.path.join(mdir, 'ice', 'outlets.csv'))
    lam, phi = G.to_rot(o.lon.values, o.lat.values, e1, e2, n)
    inside = (lam > lamG[0]) & (lam < lamG[-1]) & (phi > phiG[0]) & (phi < phiG[-1])
    h = o[(o.B2015_name == 'Sermersuaq') & (o.elev < 0) & inside].copy()
    lam, phi = G.to_rot(h.lon.values, h.lat.values, e1, e2, n)
    h['ci'] = np.floor((lam - lamG[0]) / d).astype(int)
    h['ri'] = np.floor((phi - phiG[0]) / d).astype(int)

    # front cells: wet cells sharing an edge with land, within 8 km of a Humboldt outlet
    land = ~wet
    nbrE = np.zeros_like(land); nbrE[:, :-1] = land[:, 1:]
    nbrW = np.zeros_like(land); nbrW[:, 1:] = land[:, :-1]
    nbrN = np.zeros_like(land); nbrN[:-1, :] = land[1:, :]
    nbrS = np.zeros_like(land); nbrS[1:, :] = land[:-1, :]
    coast = wet & (nbrE | nbrW | nbrN | nbrS)
    near = np.zeros_like(wet)
    near[h.ri.clip(0, ny - 1), h.ci.clip(0, nx - 1)] = True
    near = ndimage.binary_dilation(near, iterations=8)
    front = coast & near
    fr, fc = np.where(front)
    # each outlet -> nearest front cell (Mankoff outlets lie on the ice margin, i.e. on land here)
    h['fr'] = [fr[np.argmin((fr - r) ** 2 + (fc - c) ** 2)] for r, c in zip(h.ri, h.ci)]
    h['fc'] = [fc[np.argmin((fr - r) ** 2 + (fc - c) ** 2)] for r, c in zip(h.ri, h.ci)]

    # daily discharge per outlet, 2015-01-01 .. END
    ids = h['cat'].values
    dates, Q = [], []
    for y in YEARS:
        f = os.path.join(mdir, 'ice', 'discharge', f'MAR_{y}.nc')
        ds = nc.Dataset(f)
        st = ds['station'][:]
        q = np.ma.getdata(ds['discharge'][:])
        if q.shape[0] == st.size:
            q = q.T
        col = {s: k for k, s in enumerate(st)}
        Q.append(np.stack([q[:, col[i]] for i in ids], axis=1))
        t = ds['time']
        dates += [x.date() if hasattr(x, 'date') else x for x in
                  nc.num2date(t[:], t.units, getattr(t, 'calendar', 'standard'),
                              only_use_cftime_datetimes=False)]
    Q = np.concatenate(Q).astype(float)
    dates = np.array([dt.date(x.year, x.month, x.day) for x in dates])
    last = dates[-1]
    # climatological fill to END
    doy = np.array([x.timetuple().tm_yday for x in dates])
    clim_rows = dates < dt.date(2025, 1, 1)
    fill_dates = [last + dt.timedelta(days=k) for k in range(1, (END - last).days + 1)]
    fill = np.stack([Q[clim_rows & (doy == x.timetuple().tm_yday)].mean(0) for x in fill_dates])
    Q = np.concatenate([Q, fill]); dates = np.concatenate([dates, np.array(fill_dates)])
    Q[~np.isfinite(Q)] = 0.0

    # baseline lumping: top NTOP by total discharge, others -> nearest top outlet (in cells)
    tot = Q.sum(0)
    top = np.argsort(-tot)[:NTOP]
    assign = np.empty(len(ids), int)
    for k in range(len(ids)):
        dd = (h.fr.values[top] - h.fr.values[k]) ** 2 + (h.fc.values[top] - h.fc.values[k]) ** 2
        assign[k] = top[np.argmin(dd)]
    # 1-km cell-mean depths at the front are much shallower than the outlets' BedMachine
    # depths; put each plume in the deepest front cell within SEARCH cells of its outlet
    for t in top:
        r0, c0 = int(h.fr.values[t]), int(h.fc.values[t])
        cand = [(depth[r, c], r, c) for r, c in zip(fr, fc) if abs(r - r0) <= SEARCH and abs(c - c0) <= SEARCH]
        _, r1, c1 = max(cand)
        h.loc[h.index[t], 'fr'] = r1; h.loc[h.index[t], 'fc'] = c1
    plume_cells = sorted({(int(h.fr.values[t]), int(h.fc.values[t])) for t in top})

    def orient(r, c):
        # land to the E/W of the cell -> front faces x (normal along x, mask > 0); else y (< 0)
        return 1 if (nbrE[r, c] or nbrW[r, c]) else -1

    cases = {}
    # baseline
    mask = np.zeros((ny, nx)); length = np.zeros((ny, nx))
    for r, c in zip(fr, fc):
        mask[r, c] = 1 * orient(r, c)
        length[r, c] = dyG[r, c] if orient(r, c) > 0 else dxG[r, c]
    qcell_base = {}
    for t in top:
        r, c = int(h.fr.values[t]), int(h.fc.values[t])
        mask[r, c] = 3 * orient(r, c)
        qcell_base[(r, c)] = qcell_base.get((r, c), 0) + Q[:, assign == t].sum(1)
    cases['baseline'] = (mask.copy(), length.copy(), qcell_base)
    # sheet: daily total spread over all front cells by front length
    ms = np.where(mask != 0, 2 * np.sign(mask), 0.0)
    w = length[fr, fc] / length[fr, fc].sum()
    qsheet = {(int(r), int(c)): Q.sum(1) * wk for r, c, wk in zip(fr, fc, w)}
    cases['sheet'] = (ms, length.copy(), qsheet)
    # no discharge
    cases['nodischarge'] = (np.where(mask != 0, np.sign(mask), 0.0), length.copy(), {})

    os.makedirs(out, exist_ok=True)
    npz = {'dates': np.array([x.isoformat() for x in dates])}
    for name, (mk, ln, qd) in cases.items():
        mk.astype('>f8').tofile(os.path.join(out, f'Humboldt_plumes_{name}_mask.bin'))
        ln.astype('>f8').tofile(os.path.join(out, f'Humboldt_plumes_{name}_length.bin'))
        cells = np.array(sorted(qd)) if qd else np.zeros((0, 2), int)
        npz[f'{name}_cells'] = cells
        npz[f'{name}_Q'] = np.stack([qd[tuple(cc)] for cc in cells], 1) if qd else np.zeros((len(dates), 0))
    np.savez_compressed(os.path.join(out, 'Humboldt_Qsg_daily.npz'), **npz)
    rows = []
    for k, t in enumerate(top):
        r, c = int(h.fr.values[t]), int(h.fc.values[t])
        rows.append(dict(plume=k + 1, j=r, i=c, lat=float(yC[r, c]), lon=float(xC[r, c]),
                         model_depth_m=float(depth[r, c]), outlet_elev_m=float(h.elev.values[t]),
                         n_outlets_lumped=int((assign == t).sum()),
                         mean_Q=float(Q[:, assign == t].sum(1).mean()),
                         max_Q=float(Q[:, assign == t].sum(1).max())))
    pd.DataFrame(rows).to_csv(os.path.join(out, 'Humboldt_plume_cells.csv'), index=False)
    summary = dict(outlets=int(len(ids)), front_cells=int(front.sum()), plume_cells=len(plume_cells),
                   days=int(len(dates)), mar_last_day=last.isoformat(), filled_days=len(fill_dates),
                   total_mean_Q=float(Q.sum(1).mean()), total_max_Q=float(Q.sum(1).max()),
                   conserved=bool(np.allclose(sum(v for v in qcell_base.values()), Q.sum(1))))
    json.dump(summary, open(os.path.join(out, 'Humboldt_plumes_summary.json'), 'w'), indent=1)
    print(json.dumps(summary, indent=1))
    print(pd.DataFrame(rows).round(1).to_string(index=False))


if __name__ == '__main__':
    main(*sys.argv[1:4])
