#!/usr/bin/env python3
"""Temperature bias of the SoCal 2-km nest's parent-derived IC and OBCS vs CalCOFI CTD casts.

malibu_v01 / socal warm-bias fix (2026-10-19). Parent = LLC270 V5r1 (ECCO-Darwin v05).
1. IC: parent pickup THETA at iter 867888 (2025-01-01, before socal_ic_murfix.py) at every
   CalCOFI 2501RL cast inside the SoCal box (Jan 4-15 2025): model minus CTD TempAve on the
   SoCal RC levels; split by bottom depth (shelf < 500 m, offshore >= 500 m).
2. OBCS: for casts of cruises 2501RL and 2504SH within --near km of an open boundary, the
   hourly OBCS THETA record at the cast's UTC hour at the nearest boundary point, minus CTD.
   Records are read straight from the files on Pleiades with dd (no copy of whole files).
Output: json with per-cast profiles and the medians by depth level (all, shelf, offshore).

usage: socal_ctd_bias.py <dir with XC/YC/RC/hFacC.data + pickup_THETA.0000867888.data> <out.json>
"""
import json, os, subprocess, sys
import numpy as np
import pandas as pd

NX, NY, NR = 168, 101, 60
OB = '/nobackup/dcarrol2/downscaling/socal/forcings/OBCS'
HERE = os.path.dirname(os.path.abspath(__file__))
OBS = os.path.join(HERE, '..', '..', '..', 'raw_data', 'obs_2024-25')
NEAR = 60.0     # km


def ob_record(side, k):
    n = NY if side in ('west', 'east') else NX
    size = NR * n * 4
    cmd = 'dd if=%s/THETA_%s.bin bs=%d skip=%d count=1 status=none' % (OB, side, size, k)
    raw = subprocess.run(['ssh', 'pfe', cmd], capture_output=True, check=True).stdout
    return np.frombuffer(raw, '>f4').reshape(NR, n).astype(float)


def main():
    S, out = sys.argv[1:3]
    xc = np.fromfile(os.path.join(S, 'XC.data'), '>f4').reshape(NY, NX).astype(float)
    yc = np.fromfile(os.path.join(S, 'YC.data'), '>f4').reshape(NY, NX).astype(float)
    rc = -np.fromfile(os.path.join(S, 'RC.data'), '>f4').astype(float)
    h = np.fromfile(os.path.join(S, 'hFacC.data'), '>f4').reshape(NR, NY, NX)
    ic = np.fromfile(os.path.join(S, 'pickup_THETA.0000867888.data'), '>f8').reshape(NR, NY, NX)
    coslat = np.cos(np.radians(33.5))
    sides = {'west': (yc[:, 0], xc[:, 0], h[:, :, 0]), 'east': (yc[:, -1], xc[:, -1], h[:, :, -1]),
             'south': (yc[0, :], xc[0, :], h[:, 0, :]), 'north': (yc[-1, :], xc[-1, :], h[:, -1, :])}
    c = pd.concat([pd.read_csv(os.path.join(OBS, f), low_memory=False) for f in
                   ('calcofi_prelim_CTD_downcast_SCB_2411SR_2501RL_2502RL.csv',
                    'calcofi_prelim_CTD_downcast_SCB_2504SH.csv')])
    c = c[c.Lat_Dec.between(yc.min(), yc.max()) & c.Lon_Dec.between(xc.min(), xc.max())]
    res = dict(levels_m=rc.tolist(), ic={}, obcs={})
    t0 = pd.Timestamp('2025-01-01 00:00')
    for cid, g in c.groupby('Cast_ID'):
        cr = cid[:4]
        if cr not in ('2501', '2504'):
            continue
        g = g.sort_values('Depth')
        oz, ot = g.Depth.values.astype(float), g.TempAve.values.astype(float)
        ok = np.isfinite(ot)
        oz, ot = oz[ok], ot[ok]
        lat, lon = g.Lat_Dec.iloc[0], g.Lon_Dec.iloc[0]
        obs = np.where(rc <= oz.max(), np.interp(rc, oz, ot), np.nan)
        obs[rc < oz.min()] = ot[0]
        if cr == '2501':
            d2 = (yc - lat) ** 2 + ((xc - lon) * coslat) ** 2
            d2 = np.where(h[0] > 0, d2, np.inf)
            j, i = np.unravel_index(d2.argmin(), d2.shape)
            m = ic[:, j, i]
            wet = h[:, j, i] > 0
            bottom = float(rc[wet].max())
            res['ic'][cid] = dict(lat=lat, lon=lon, bottom=bottom,
                                  bias=[float(a - b) if w and np.isfinite(b) else None
                                        for a, b, w in zip(m, obs, wet)])
        # nearest boundary point
        best = None
        for s, (la, lo, hh) in sides.items():
            dk = np.hypot((la - lat) * 111.0, (lo - lon) * 111.0 * coslat)
            dk = np.where(hh[0] > 0, dk, np.inf)
            p = int(dk.argmin())
            if best is None or dk[p] < best[1]:
                best = (s, float(dk[p]), p)
        if best[1] > NEAR:
            continue
        s, dist, p = best
        tu = pd.to_datetime(g.Date_Time_UTC.iloc[0])
        k = int((tu - t0).total_seconds() // 3600)
        rec = ob_record(s, k)
        m = rec[:, p]
        wet = sides[s][2][:, p] > 0
        res['obcs'][cid] = dict(side=s, dist_km=dist, record=k, bottom=float(rc[wet].max()),
                                bias=[float(a - b) if w and np.isfinite(b) else None
                                      for a, b, w in zip(m, obs, wet)])
        print(cid, s, '%.0f km' % dist, 'rec', k)
    for part in ('ic', 'obcs'):
        for name, sel in (('all', lambda d: True), ('shelf', lambda d: d['bottom'] < 500),
                          ('offshore', lambda d: d['bottom'] >= 500)):
            for cr in (('2501', '2504') if part == 'obcs' else ('2501',)):
                rows = [np.array([np.nan if x is None else x for x in d['bias']])
                        for cid, d in res[part].items() if cid.startswith(cr) and sel(d)]
                if not rows:
                    continue
                a = np.array(rows)
                med = np.nanmedian(a, 0)
                n = np.isfinite(a).sum(0)
                res['%s_median_%s_%s' % (part, cr, name)] = [None if k < 3 else float(v) for v, k in zip(med, n)]
                res['%s_n_%s_%s' % (part, cr, name)] = n.tolist()
                print(part, cr, name, 'n=%d' % len(rows), ' '.join(
                    '%.0fm:%+.2f' % (z, v) for z, v, k in zip(rc, med, n) if k >= 3 and z <= 150 and int(z) % 10 in (1, 5, 0) and z < 160))
    json.dump(res, open(out, 'w'), indent=1)


if __name__ == '__main__':
    main()
