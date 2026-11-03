#!/usr/bin/env python3
"""PACE-scaled boundary phytoplankton for the Malibu nest, Feb-Apr 2025 (2026-10-23, March bloom test).

The v4 boundaries (work/malibu2025_v4/obcs, CalCOFI-based, calcofi_bgc.py) hold surface total Chl of
0.9-1.6 mg m-3 through March, while PACE OCI daily L3 Chl shows a regional bloom on Mar ~8-23
(shelf 4-24, basin 2-12 mg m-3). The 0-30 m budget shows that the nearshore loss is horizontal
export (residence ~4-8 d), so the nest can only show the bloom if the boundaries bring it in.

Factor per boundary and day: f(t) = PACE(t) / PACE_ref, with
  PACE(t)  = median chlor_a of the PACE pixels within STRIP deg inside the boundary line (and inside the
             nest box), days with < MINPIX valid pixels treated as missing, gaps filled by linear
             interpolation in time, then a 5-day running median;
  PACE_ref = median of PACE(t) over Feb 1 - Mar 5 (before the bloom).
f is clipped to [FMIN, FMAX], set to 1 before Feb 1 and after Apr 30, and ramped linearly to 1 over the
first and last 5 days of Feb 1 - Apr 30. The relative change keeps the CalCOFI calibration of the
boundaries and does not import a constant PACE bias in coastal water. Phytoplankton C (c1-c5) and Chl
(Chl1-Chl5) are multiplied by f at all depths; zooplankton (c6, c7) and nutrients are unchanged.

usage: make_pace_obcs_phyto.py <v4 obcs dir (local copy of c1-c5, Chl1-Chl5 files)> <out dir> <json out>
"""
import datetime as dt
import glob
import json
import os
import sys
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
WORK = os.path.join(HERE, '..', 'work')
NR, NY, NX, NREC = 50, 72, 128, 182
STRIP, MINPIX, FMIN, FMAX = 0.09, 3, 0.5, 15.0
T0 = dt.date(2025, 1, 1)
D1, D2, RAMP = dt.date(2025, 2, 1), dt.date(2025, 4, 30), 5


def runmed(x, w=5):
    h = w // 2
    return np.array([np.nanmedian(x[max(0, i - h):i + h + 1]) for i in range(len(x))])


def main():
    src, out, jout = sys.argv[1:4]
    g = json.load(open(os.path.join(WORK, 'grid_info.json')))
    xc, yc = np.array(g['xc']), np.array(g['yc'])
    x0, x1, y0 = xc.min(), xc.max(), yc.min()
    days = [T0 + dt.timedelta(k) for k in range(NREC)]
    sides = {'West': lambda la, lo: (lo >= x0) & (lo <= x0 + STRIP),
             'East': lambda la, lo: (lo <= x1) & (lo >= x1 - STRIP),
             'South': lambda la, lo: (la >= y0) & (la <= y0 + STRIP) & (lo >= x0) & (lo <= x1)}
    raw = {s: np.full(NREC, np.nan) for s in sides}
    npix = {s: np.zeros(NREC, int) for s in sides}
    for k, d in enumerate(days):
        f = os.path.join(WORK, 'pace', 'pace_bgc_%s_smbay.npz' % d.strftime('%Y%m%d'))
        if not os.path.exists(f):
            continue
        p = np.load(f)
        LA, LO = np.meshgrid(p['lat'], p['lon'], indexing='ij')
        inbox = (LA >= y0) & (LA <= yc.max()) & np.isfinite(p['chlor_a'])
        for s, sel in sides.items():
            m = inbox & sel(LA, LO)
            npix[s][k] = m.sum()
            if m.sum() >= MINPIX:
                raw[s][k] = np.median(p['chlor_a'][m])
    res = dict(strip_deg=STRIP, minpix=MINPIX, clip=[FMIN, FMAX], ref='median Feb 1 - Mar 5', sides={})
    fac = {}
    kd1, kd2 = (D1 - T0).days, (D2 - T0).days
    for s in sides:
        r = raw[s].copy()
        ok = np.isfinite(r)
        win = np.arange(NREC)
        r[kd1:kd2 + 1] = np.interp(win[kd1:kd2 + 1], win[ok], r[ok]) if ok.sum() >= 2 else np.nan
        sm = runmed(r)
        ref = np.nanmedian(sm[kd1:(dt.date(2025, 3, 5) - T0).days + 1])
        f = np.ones(NREC)
        f[kd1:kd2 + 1] = np.clip(sm[kd1:kd2 + 1] / ref, FMIN, FMAX)
        for i in range(RAMP):                                           # ramp in and out
            w = (i + 1) / (RAMP + 1)
            f[kd1 + i] = 1 + w * (f[kd1 + i] - 1)
            f[kd2 - i] = 1 + w * (f[kd2 - i] - 1)
        fac[s] = f
        res['sides'][s] = dict(ref_mgm3=float(ref), days_with_pace=int(ok[kd1:kd2 + 1].sum()),
                               factor_max=float(f.max()), factor_by_10d={str(days[d]): round(float(f[d]), 2)
                                                                        for d in range(kd1, kd2 + 1, 10)},
                               raw=[None if not np.isfinite(v) else round(float(v), 3) for v in raw[s]],
                               npix=npix[s].tolist(), factor=[round(float(v), 4) for v in f])
        print(s, 'ref %.2f mg m-3, %d PACE days, factor max %.1f' % (ref, ok[kd1:kd2 + 1].sum(), f.max()),
              ' '.join('%s:%.1f' % (days[d].strftime('%m-%d'), f[d]) for d in range(kd1, kd2 + 1, 7)))
    os.makedirs(out, exist_ok=True)
    for s, n in (('West', NY), ('East', NY), ('South', NX)):
        for v in ['c%d' % k for k in range(1, 6)] + ['Chl%d' % k for k in range(1, 6)]:
            fn = 'malibu_128x72x50_%s_%s' % (v, s)
            a = np.fromfile(os.path.join(src, fn), '>f4').reshape(NREC, NR, n)
            (a * fac[s][:, None, None]).astype('>f4').tofile(os.path.join(out, fn))
    json.dump(res, open(jout, 'w'), indent=1)


if __name__ == '__main__':
    main()
