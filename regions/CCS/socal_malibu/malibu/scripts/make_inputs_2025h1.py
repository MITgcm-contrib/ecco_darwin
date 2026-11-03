#!/usr/bin/env python3
"""Namelists for the Jan-Jun 2025 Malibu run (socal-nested), from the S2/S3 namelists (2026-10-10).

Starts from input_pace2024_s2 (S3 = S2 namelists + the A23 build) and changes only what the new
window and forcing need; tracers, traits, kelp, optics and sediment settings stay as in S3:
- data.cal: startDate 2025-01-01; data: nTimeSteps for --days at deltaT 90 s.
- data.exf: atmosphere from the sea-filled ERA5 box of the socal nest (era5_seafill/, hourly,
  era5_seafill_box.py; land cells <- nearest ocean cell) instead of the V5r1 6-hourly fields;
  EXF_NML_03 as in socal (degC offset, -1 for the ERA5 radiation sign, precip x1); EXF_NML_04
  grid lines from era5_box_grid.txt; waves grid from the new waves.json; OBCS startdate
  2024-12-31 12:00 and period 86400 s (daily socal-mode records, make_ic_obcs.py --mode socal).
- data.sediment: sedWaveNrec from the new waves.json.

Usage: python3 make_inputs_2025h1.py <out dir> <era5_box_grid.txt> <waves dir> [--days 180]
"""
import json
import os
import re
import shutil
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
SRC = os.path.join(HERE, '..', 'input_pace2024_s2')
DELTAT = 90.0
ERA5 = {'atemp': 'tmp2m_degC', 'aqh': 'spfh2m', 'precip': 'rain', 'uwind': 'u10m', 'vwind': 'v10m',
        'swdown': 'dsw', 'lwdown': 'dlw'}
V5 = {'atemp': 'EXFatemp', 'aqh': 'EXFaqh', 'precip': 'EXFpreci', 'uwind': 'EXFuwind', 'vwind': 'EXFvwind',
      'swdown': 'EXFswdn', 'lwdown': 'EXFlwdn'}


def sub1(s, pat, rep, flags=0):
    s2, n = re.subn(pat, rep, s, count=1, flags=flags)
    assert n == 1, pat
    return s2


def main():
    out, gridtxt, wdir = sys.argv[1:4]
    days = float(sys.argv[sys.argv.index('--days') + 1]) if '--days' in sys.argv else 180.0
    if os.path.exists(out):
        sys.exit(f'{out} exists')
    shutil.copytree(SRC, out)
    p = lambda f: os.path.join(out, f)

    s = open(p('data.cal')).read()
    s = sub1(s, r'startDate_1=\d+', 'startDate_1=20250101')
    s = s.replace('PACE-year case: ICs from ECCO-Darwin v05 V5r1', 'Jan-Jun 2025 case: ICs/OBCS from the socal 2-km nest')
    open(p('data.cal'), 'w').write(s)

    s = open(p('data')).read()
    s = sub1(s, r' nTimeSteps=\d+,', f' nTimeSteps={int(round(days * 86400 / DELTAT))},')
    open(p('data'), 'w').write(s)

    s = open(p('data.exf')).read()
    for f, e in ERA5.items():
        s = sub1(s, rf"  {f}file\s*= 'v5r1/{V5[f]}_6hourly',", f"  {f}file = 'era5_seafill/ERA5_{e}',")
        s = sub1(s, rf"  {f}period\s*= 21600.0,", f"  {f}period = 3600.0,")
        if f in ('precip', 'swdown', 'lwdown'):   # ERA5 hourly accumulations: centered at hh:30
            s = sub1(s, rf"  {f}startdate2\s*= 030000,", f"  {f}startdate2 = 003000,")
    s = s.replace(' &EXF_NML_02\n', ' &EXF_NML_02\n# atmosphere: sea-filled ERA5 hourly box of the socal 2-km nest (era5_seafill_box.py)\n', 1)
    s = sub1(s, r' &EXF_NML_03\n(?:[^/]*?\n)? /',
             ' &EXF_NML_03\n# ERA5 (hzhang1 forcing/era5) conventions, as in socal regional_run data.exf\n'
             '  exf_offset_atemp  = 273.15,\n  exf_inscal_swdown = -1.0,\n  exf_inscal_lwdown = -1.0,\n'
             '  exf_inscal_precip = 1.0,\n /', flags=re.S)
    grid = open(gridtxt).read()
    wv = json.load(open(os.path.join(wdir, 'waves.json')))
    lat_w = np.array(wv['lat'])
    wave = ''.join(f"  {f}_lon0 = {wv['lon0']:.4f},\n  {f}_lon_inc = {wv['dlon']:.4f},\n"
                   f"  {f}_lat0 = {lat_w[0]:.4f},\n  {f}_lat_inc = {lat_w.size - 1}*{np.diff(lat_w).mean():.4f},\n"
                   f"  {f}_nlon = {wv['nx']},\n  {f}_nlat = {wv['ny']},\n  {f}_interpMethod = 1,\n"
                   for f in ('WvHeight', 'WvPeriod'))
    grid = ''.join(l + '\n' for l in grid.splitlines() if not l.strip().startswith('apressure'))
    s = sub1(s, r' &EXF_NML_04\n.*?\n /', lambda m: ' &EXF_NML_04\n# ERA5 box grid (era5_box_grid.txt)\n' + grid
             + wave + '  runoff_interpMethod = 0,\n /', flags=re.S)
    obcs = ' &EXF_NML_OBCS\n# daily socal-mode records: record 1 = 2024-12-31 12:00 (copy of the Jan 1 mean)\n' + ''.join(
        f'  obcs{c}startdate1 = 20241231,\n  obcs{c}startdate2 = 120000,\n  obcs{c}period = 86400.0,\n'
        for c in ('W', 'S', 'E')) + ' /'
    s = sub1(s, r' &EXF_NML_OBCS.*?\n /', lambda m: obcs, flags=re.S)
    assert 'v5r1/EXF' not in s
    open(p('data.exf'), 'w').write(s)

    s = open(p('data.sediment')).read()
    s = sub1(s, r' sedWaveNrec = \d+,', f" sedWaveNrec = {wv['nrec_sed']},")
    s = re.sub(r'hourly from[^\n]*', f"hourly from {wv['start'][:10]} (waves_2025h1, make_waves.py)", s, count=1)
    open(p('data.sediment'), 'w').write(s)
    print(f'wrote {out}: start 2025-01-01, {days:g} d, ERA5 sea-filled atmosphere, waves nrec {wv["nrec_sed"]}')


if __name__ == '__main__':
    main()
