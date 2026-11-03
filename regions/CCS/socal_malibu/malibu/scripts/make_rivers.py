#!/usr/bin/env python3
"""River point sources for the Malibu nest: freshwater (EXF runoff), suspended
sediment (darwin TSSrunoff -> TSS tracer 36, handled by pkg/sediment) and
nutrients / carbon (darwin NO3, NH4, IP(=PO4), DSi, DOC, DON, DOP, POC, PON,
POP, DIC runoff fluxes), with an optional wildfire effect (Woolsey 2018).

For each river in RIVERS the mouth is put in the nearest wet surface cell
(model hFacC rules) next to land, and daily discharge Q [m3/s] is turned into
  runoff  = Q / A_cell                       [m/s]     (EXF runofffile)
  TSSflux = Q * SSC(Q) / A_cell              [g/m2/s]  (darwin TSSrunofffile)
with SSC [mg/L = g/m3] from a rating curve SSC = a * Q**b.
These replace the parent's interpolated JRA55-do runoff and TSS_CCS_cutout
(which smear river input over 1/3-deg cells) inside the nest.

Discharge sources (per river, see RIVERS):
  f130r             DEFAULT. Observed LA County F130-R daily means (Oct 1991 -
                    Sep 2026, f130r_daily.py), carried to the mouth with the
                    fit to USGS 11105510 in work/rivers/F130R_to_mouth.json.
  synthetic         Rainfall-runoff model driven by the parent model's own rain
                    (malibu_precip.py, synth_runoff.py; fitted to USGS 11105510
                    2008-12, NSE 0.77); fallback for years outside F130-R.
  usgs:<site>       USGS daily values (waterservices.usgs.gov)
  monthly:<table>   monthly means (cfs) from the LA County Hydrologic Reports,
                    constant within each month
  csv:<file>        your own daily file: date (YYYY-MM-DD), Q (m3/s)

Writes yearly files (useExfYearlyFields = T) to ../work/rivers/:
  rivers_runoff_<nx>x<ny>_<YYYY>, rivers_TSSflux_<nx>x<ny>_<YYYY>
  (one record per day, big-endian float32), and rivers_<YYYY>.csv (Q, SSC).

Usage: python3 make_rivers.py 1992 [1993 ...]
"""
import csv
import datetime as dt
import json
import os
import sys
import subprocess

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
WORK = os.path.join(HERE, '..', 'work')
sys.path.insert(0, HERE)
from make_ic_obcs import hfac_c, nest_grid  # noqa: E402

CFS = 0.0283168466  # m3/s per cfs

# Malibu Creek below Cold Creek, LA County station F130-R (104.96 mi2), monthly
# mean discharge [cfs] from the LA County DPW Hydrologic Reports, p. C19:
# WY1991-92 and WY1992-93 (Dec 1992: "no records", filled with Nov 1992).
F130R_MONTHLY_CFS = {
    (1991, 10): 2.5, (1991, 11): 3.3, (1991, 12): 36.8,
    (1992, 1): 46.9, (1992, 2): 669.0, (1992, 3): 252.0, (1992, 4): 58.7,
    (1992, 5): 25.9, (1992, 6): 29.2, (1992, 7): 8.5, (1992, 8): 3.1,
    (1992, 9): 5.5, (1992, 10): 8.6, (1992, 11): 9.3, (1992, 12): 9.3,
    (1993, 1): 440.0, (1993, 2): 187.0, (1993, 3): 68.0, (1993, 4): 45.0,
    (1993, 5): 24.0, (1993, 6): 14.1, (1993, 7): 4.7, (1993, 8): 6.3,
}

RIVERS = [
    dict(name='Malibu Creek',
         # Malibu Lagoon outlet at Surfrider Beach (GNIS mouth 34.03167 N,
         # 118.67972 W, ~30 m away)
         lon=-118.6800, lat=34.0320,
         # OBSERVED daily discharge (default): LA County Public Works F130-R
         # (Malibu Creek below Cold Creek), Oct 1991 - Sep 2026, emailed
         # 2026-10-06 (raw_data/F130-R; python/f130r_daily.py -> work/rivers/
         # F130R_daily.csv, scaled x109.2/105.0 by drainage area, 262 missing days
         # interpolated incl. 182 days of 1998; data from Oct 2025 provisional).
         # Delivered to the mouth as Q = max(s Q_F130 - L, 0), s = 1.23,
         # L = 0.0925 m3/s, fitted to USGS 11105510 (Malibu Creek at Malibu,
         # Dec 2007 - Oct 2012): NSE 0.96, zero-flow days 25.8% vs 26.3% observed
         # (work/rivers/F130R_to_mouth.json; s is empirical, larger than the
         # 4% area difference). 'synthetic' (rainfall-runoff model) only for
         # years outside the record.
         discharge={'default': 'f130r'},
         # PLACEHOLDER rating curve: exponent b = 1.5 (ASSUMED, no Malibu daily
         # curve found; Warrick et al. 2015 J. Hydrol. Reg. Stud. 4:516 show SSC
         # in small Santa Barbara-coast basins is only weakly tied to Q), a
         # anchored to the geometric mean of the two USGS samples at 11105500
         # (1986-02-01: 2.9 m3/s, 114 mg/L; 1987-03-06: 2.1 m3/s, 21 mg/L).
         ssc_a=12.6, ssc_b=1.5,
         # base chemistry [mmol/m3], flow-independent, wet-season values where
         # available (the load is carried by winter flows):
         #  NO3 2.60 mg N/L: median NO3-N, LACDPW mass-emission station on Malibu
         #    Creek, discharge period Nov 16 - Apr 14, n = 43 (2013 USEPA Malibu
         #    Creek & Lagoon TMDL, Table 7-12)
         #  NH4 0.26 mg N/L: median of the 4 pre-fire storm events with a
         #    detect (Jan-Mar 2018), LA County Public Works station S02, project
         #    LACFCD_ME02 (raw_data/LA_public_works/Chem_102026_0925AM.xlsx,
         #    received 2026-10-20). It replaces 0.06 mg N/L (Heal the Bay MC-1,
         #    1998-2010, mostly dry weather), which was 4x lower.
         #    S02 pre-fire storm NO3+NO2 median 2.8 mg N/L supports NO3 = 2.60.
         #  PO4 0.46 mg P/L: median at Heal the Bay Stream Team MC-1 (lower
         #    Malibu Creek), 1998-2010, n = 117 (TMDL Table 7-7). S02 has total P only.
         #  DON 0.49 mg N/L, DOP 0.026 mg P/L: TDN - DIN and TDP - SRP below
         #    Tapia, Aug 2002 (SCCWRP Tech. Rep. 412, Busse et al. 2003, Table 2;
         #    dry season)
         #  DOC 3.4 mg C/L: median Malibu Creek main stem, May-Jul 2009-2020
         #    (CEDEN/SWAMP; baseflow, no storm DOC found)
         #  DSi 380 uM: median SiO2 23 mg/L, USGS 11105500, 1982-88 (n = 15, Water
         #    Quality Portal)
         #  DIC 3500 uM: median alkalinity 166 mg/L CaCO3 (3.3 meq/L), pH 7.6-7.9
         #    at USGS 11105500, 1982-88; DIC = ALK (1 + [H+]/K1) = 3.45 mM at
         #    pH 7.75, pK1 6.42 (freshwater, ~15 C)
         conc=dict(NO3=2.60e3 / 14.007, NH4=0.26e3 / 14.007, IP=0.46e3 / 30.974,
                   DON=0.49e3 / 14.007, DOP=0.026e3 / 30.974, DOC=3.4e3 / 12.011,
                   DSi=380.0, DIC=3500.0),
         # particulate organic matter on the sediment (PLACEHOLDER, no Malibu
         # data): 1.5% OC by mass, C:N 10, C:P 200 (molar)
         poc_frac=0.015, pom_cn=10.0, pom_cp=200.0,
         # Woolsey Fire (8 Nov 2018) burned 52.0% of the basin above USGS
         # 11105510 (our GIS: NIFC perimeter x USGS NLDI basin, 147.5 of 283.5
         # km2; no published value found). Factors F for a fully burned basin,
         # scaled as 1 + fb (F - 1) exp(-t/tau):
         #  runoff 2.3, tau 3.5 yr: ONLY used for synthetic discharge (observed
         #    F130-R already contains the fire). Fit to yearly runoff-ratio
         #    increases after the 2003 Old Fire (Devil Canyon 2.5x -> 1.4x and
         #    City Creek 1.6x -> 1.2x, WY2004-08; Kinoshita 2012 UCLA PhD Table
         #    4-1 = Kinoshita & Hogue 2011, Catena 87:240)
         #  SSC 10, NH4 32, NO3 5.5, DON 2.8, PO4 2.2: burned/unburned ANNUAL
         #    EXPORT ratios, first post-fire year, 2004 Gaviota Fire, chaparral,
         #    where annual runoff was similar (Coombs & Melack 2013, Hydrol.
         #    Process. 27:3842, abstract), applied here as concentration factors.
         #    CONFLICT for SSC: Malibu WY2019 sediment load only 1.3x non-fire
         #    years from a rating curve (Lopez et al. 2025, Earth Space Sci.
         #    12:e2024EA003843), which cannot show a post-fire SSC shift.
         #  tau: NH4 0.25 yr (export mostly in the first three storms, Coombs &
         #    Melack 2013); DON/PO4 1.5 yr are ASSUMED (no published recovery
         #    time found).
         #  Check vs LA County S02 storm samples (2026-10-20; grab samples vs
         #    F130-R daily-mean Q, so each event ratio is uncertain by ~2x):
         #    SSC / rating curve = 0.9-1.2 before the fire, 2.4-7.7 in Jan 2019
         #    (model 5.2-5.3), 0.16-0.89 in Nov 2019-Mar 2020 (old model 3.4-3.8).
         #    So SSC tau 2 yr -> 0.5 yr. Post-fire NO3+NO2 0.9-2.0 mg N/L, with no
         #    rise above pre-fire (median 2.8), so the NO3 factor (5.5, 1.5 yr) is
         #    removed. Post-fire NH4 0.67-1.05 mg N/L (~3x pre-fire); the NH4
         #    factor is not changed in this update.
         fires=[dict(name='Woolsey', date=dt.date(2018, 11, 8), fb=0.520,
                     F=dict(runoff=(2.3, 3.5), SSC=(10.0, 0.5), NH4=(32.0, 0.25),
                            DON=(2.8, 1.5), IP=(2.2, 1.5)))],
         ncell=1),
]
TABLES = {'F130R': F130R_MONTHLY_CFS}
# darwin runoff fields written per river (IP is used as PO4: R_DIP_IP_runoff = 1)
BGC = ['NO3', 'NH4', 'IP', 'DSi', 'DOC', 'DON', 'DOP', 'POC', 'PON', 'POP', 'DIC']


def cell_areas(g):
    R = 6371.0e3
    dx = np.deg2rad(g['dlon']) * R * np.cos(np.deg2rad(np.array(g['yc'])))
    dy = np.deg2rad(np.array(g['dely'])) * R
    return (dy * dx)[:, None] * np.ones(g['nx'])[None, :]


def mouth_cells(g, hf0, lon, lat, n):
    """n nearest wet surface cells that touch land."""
    wet = hf0 > 0
    land = ~wet
    coast = wet & (np.roll(land, 1, 0) | np.roll(land, -1, 0) |
                   np.roll(land, 1, 1) | np.roll(land, -1, 1))
    X, Y = np.meshgrid(g['xc'], g['yc'])
    dist = np.hypot((X - lon) * np.cos(np.deg2rad(lat)), Y - lat)
    dist[~coast] = np.inf
    idx = np.argsort(dist.ravel())[:n]
    return [np.unravel_index(i, dist.shape) for i in idx], dist.ravel()[idx] * 111.2


def daily_q(spec, year, r=None):
    days = [dt.date(year, 1, 1) + dt.timedelta(d) for d in
            range((dt.date(year + 1, 1, 1) - dt.date(year, 1, 1)).days)]
    if spec == 'synthetic':
        return days, synthetic_q(year, r)
    if spec == 'f130r':
        fn = os.path.join(WORK, 'rivers', 'F130R_daily.csv')
        fit = json.load(open(os.path.join(WORK, 'rivers', 'F130R_to_mouth.json')))
        q = {row[0]: float(row[1]) for row in csv.reader(open(fn)) if row and row[0][:4].isdigit()}
        if not all(d.isoformat() in q for d in days):
            print(f'  F130-R does not cover {year}: synthetic discharge')
            return days, synthetic_q(year, r)
        qf = np.array([q[d.isoformat()] for d in days])
        return days, np.maximum(fit['scale'] * qf - fit['loss_m3s'], 0.0)
    kind, arg = spec.split(':', 1)
    if kind == 'monthly':
        tab = TABLES[arg]
        return days, np.array([tab[(d.year, d.month)] * CFS for d in days])
    if kind == 'usgs':
        url = ('https://waterservices.usgs.gov/nwis/dv/?format=rdb&sites=' + arg +
               f'&parameterCd=00060&startDT={year}-01-01&endDT={year}-12-31')
        txt = subprocess.run(['curl', '-s', '-m', '120', url], capture_output=True,
                             text=True, check=True).stdout
        q = {}
        for line in txt.splitlines():
            p = line.split('\t')
            if len(p) > 3 and p[0] == 'USGS':
                try:
                    q[p[2]] = float(p[3]) * CFS
                except ValueError:
                    pass
        vals = np.array([q.get(d.isoformat(), np.nan) for d in days])
        if np.isnan(vals).any():
            print(f'  {arg} {year}: {np.isnan(vals).sum()} missing days filled by interpolation')
            ok = ~np.isnan(vals)
            vals = np.interp(np.arange(len(vals)), np.where(ok)[0], vals[ok])
        return days, vals
    if kind == 'csv':
        q = {r[0]: float(r[1]) for r in csv.reader(open(arg)) if r and r[0][:4].isdigit()}
        return days, np.array([q[d.isoformat()] for d in days])
    raise ValueError(spec)


def fire_mult(r, days, key):
    """Multiplier for one term from the river's fires (1 outside fires)."""
    m = np.ones(len(days))
    for f in r.get('fires', []):
        if key not in f['F']:
            continue
        F, tau = f['F'][key]
        t = np.array([(d - f['date']).days / 365.25 for d in days])
        m *= np.where(t >= 0, 1.0 + f['fb'] * (F - 1.0) * np.exp(-np.clip(t, 0, None) / tau), 1.0)
    return m


def synthetic_q(year, r=None):
    """Fitted rainfall-runoff model; spun up on the previous year's rain
    (or on the same year when that is not cached)."""
    import synth_runoff as sr
    par = json.load(open(os.path.join(WORK, 'rivers', 'synth_runoff_params.json')))['params']
    par = [par[k] for k in ('Smax_mm', 'c', 'p', 'kf_days', 'ks_days', 'phi_fast',
                            'ET0_mm_day', 'Qb_m3s', 'Qloss_m3s')]
    prev = year - 1 if os.path.exists(os.path.join(WORK, 'rivers', f'precip_malibu_{sr.RAIN}{year - 1}.npy')) else year
    days, p = sr.precip([prev, year])
    q = sr.simulate(par, days, p, cmul=fire_mult(r, days, 'runoff') if r else None)
    n = (dt.date(year + 1, 1, 1) - dt.date(year, 1, 1)).days
    return q[-n:]


def main():
    years = [int(y) for y in sys.argv[1:]] or [1992]
    g = nest_grid()
    hf = hfac_c(g)
    area = cell_areas(g)
    out = os.path.join(WORK, 'rivers')
    os.makedirs(out, exist_ok=True)
    nx, ny = g['nx'], g['ny']
    for year in years:
        ndays = (dt.date(year + 1, 1, 1) - dt.date(year, 1, 1)).days
        runoff = np.zeros((ndays, ny, nx))
        tss = np.zeros((ndays, ny, nx))
        bgc = {k: np.zeros((ndays, ny, nx)) for k in BGC}
        rows = []
        for r in RIVERS:
            spec = r['discharge'].get(year, r['discharge']['default'])
            days, q = daily_q(spec, year, r)
            ssc = (r['ssc_a'] * fire_mult(r, days, 'SSC')
                   * np.maximum(q, 0.0) ** r['ssc_b'])                     # g/m3
            cc = {k: v * fire_mult(r, days, k) for k, v in r['conc'].items()}
            cc['POC'] = r['poc_frac'] * ssc / 12.011 * 1e3                 # mmol/m3
            cc['PON'] = cc['POC'] / r['pom_cn']
            cc['POP'] = cc['POC'] / r['pom_cp']
            cells, dkm = mouth_cells(g, hf[0], r['lon'], r['lat'], r['ncell'])
            atot = sum(area[j, i] for j, i in cells)
            for j, i in cells:
                runoff[:, j, i] += q / atot
                tss[:, j, i] += q * ssc / atot
                for k in BGC:
                    bgc[k][:, j, i] += q * cc[k] / atot                    # mmol/m2/s
            print(f"{r['name']} {year}: {spec}, mouth cell(s) {[(int(j), int(i)) for j, i in cells]} "
                  f"({dkm.max():.2f} km from mouth), Q mean {q.mean():.2f} max {q.max():.1f} m3/s, "
                  f"SSC max {ssc.max():.0f} mg/L, sediment {np.sum(q * ssc) * 86400 / 1e6:.0f} t/yr, "
                  f"NO3 {np.sum(q * cc['NO3']) * 86400 * 14.007e-9:.1f} t N/yr")
            rows += [(d.isoformat(), r['name'], f'{qq:.4f}', f'{s:.2f}') for d, qq, s in zip(days, q, ssc)]
        runoff.astype('>f4').tofile(os.path.join(out, f'rivers_runoff_{nx}x{ny}_{year}'))
        tss.astype('>f4').tofile(os.path.join(out, f'rivers_TSSflux_{nx}x{ny}_{year}'))
        for k in BGC:
            bgc[k].astype('>f4').tofile(os.path.join(out, f'rivers_{k}flux_{nx}x{ny}_{year}'))
        with open(os.path.join(out, f'rivers_{year}.csv'), 'w', newline='') as f:
            w = csv.writer(f)
            w.writerow(['date', 'river', 'Q_m3s', 'SSC_mgL'])
            w.writerows(rows)
    with open(os.path.join(out, 'rivers.json'), 'w') as f:
        json.dump([{k: v for k, v in r.items() if k not in ('discharge', 'fires')} for r in RIVERS], f, indent=1)


if __name__ == '__main__':
    main()
