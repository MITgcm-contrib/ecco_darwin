#!/usr/bin/env python3
"""
G-POM supplement (Fumenia et al. 2026, Sci Data, doi:10.1038/s41597-026-08413-7; data SEANOE doi:10.17882/113424,
CC-BY) -> obs/<SITE>/<SITE>_obs_gpom.csv, read by scripts/obs_utils.load_obs next to <SITE>_obs.csv (2026-10-18).

  python3 obs/scripts/std_gpom.py [obs/gpom/G-POM_revised.csv]

Selection: stations within 50 km of the station position; glass-fiber GF/F filters (Flag_filter_type_POC_PON = 1);
POC/PON molar-ratio flag 1 (valid) or 4 (not computable); replicate CV flag not "> 5%". Acidification is NOT required:
HOT reports total particulate C without acidification for its whole record (as in obs/HOT, "HOT PC"). PON = CHN, else
chemical wet oxidation; POP = ash hydrolysis, else wet oxidation. At sites that already have POC/PON/POP from the same
program (HOT, BATS: BCO-DMO 3773/3782), days already present for that variable are dropped, so only new cruises are
added. Units are kept as reported (mg m-3); obs_utils converts to mmol m-3.
"""
import os, sys
import numpy as np
import pandas as pd

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..')
SRC = sys.argv[1] if len(sys.argv) > 1 else os.path.join(ROOT, 'gpom', 'G-POM_revised.csv')
SITES = {'HOT': (22.75, -158.00), 'BATS': (31.67, -64.17), 'HydroS': (32.17, -64.50), 'Papa': (50.10, -144.90),
         'PAP': (49.00, -16.50), 'SOTS': (-46.94, 142.28), 'K2': (47.00, 160.00), 'MUNIDA': (-45.83, 171.53),
         'OWSM': (66.00, 2.00), 'CARIACO': (10.50, -64.67), 'DYFAMED': (43.4167, 7.8667)}
RKM = 50.


def dist(lat, lon, a, b):
    p = np.radians
    return 2 * 6371. * np.arcsin(np.sqrt(np.sin(p(lat - a) / 2) ** 2 +
                                        np.cos(p(a)) * np.cos(p(lat)) * np.sin(p(lon - b) / 2) ** 2))


g = pd.read_csv(SRC, low_memory=False)
g['time'] = pd.to_datetime(g['Date'].astype(str), format='%Y%m%d', errors='coerce') + \
    pd.to_timedelta(g['HH'].fillna(0), unit='h') + pd.to_timedelta(g['mm'].fillna(0), unit='m')
ok = (g['Flag_filter_type_POC_PON'] == 1) & g['QF_POC_PON_molar_ratio'].isin([1, 4]) & \
     (g['QF_CV_POC_CHN'] != 2) & (g['QF_CV_PON_CHN'] != 2)
g = g[ok & g['time'].notna()]
val = {'POC': g['POC_CHN (mg.m3)'],
       'PON': g['PON_CHN (mg.m3)'].fillna(g['PON_CWO (mg.m3)']),
       'POP': g['POP_ash_hydrolysis (mg.m3)'].fillna(g['POP_CWO (mg.m3)'])}
for s, (a, b) in SITES.items():
    near = dist(g['Latitude'].values, g['Longitude'].values, a, b) <= RKM
    f = os.path.join(ROOT, s, '%s_obs.csv' % s)
    if not near.any() or not os.path.exists(f):
        continue
    have = pd.read_csv(f, usecols=['time', 'variable'], low_memory=False)
    have = have[have['variable'].isin(val)]
    have['day'] = pd.to_datetime(have['time'], errors='coerce', utc=True).dt.tz_localize(None).dt.normalize()
    rows = []
    for v, x in val.items():
        k = near & x.notna().values
        d = pd.DataFrame({'time': g['time'][k], 'depth_m': g['Depth'][k], 'variable': v, 'value': x[k],
                          'units': 'mg/m3 (G-POM)', 'source': 'gpom_' + g['Cruise'][k].astype(str), 'qc': 2})
        old = set(have.loc[have['variable'] == v, 'day'])
        d = d[~d['time'].dt.normalize().isin(old)]
        d['time'] = d['time'].dt.strftime('%Y-%m-%dT%H:%M:%SZ')          # same ISO form as <SITE>_obs.csv
        rows.append(d)
    out = pd.concat(rows)
    if len(out):
        out.to_csv(os.path.join(ROOT, s, '%s_obs_gpom.csv' % s), index=False)
    print('%-8s %s' % (s, out.groupby('variable').size().to_dict()))
