#!/usr/bin/env python3
"""LR17 inputs for runs outside gendata.py's built-in download window, without editing gendata.py.

gendata.py fixes its input files (TAG = '2022-02-24_2022-03-18', RIVER_FILE for Feb-Mar 2022). This wrapper
imports it, points TAG and RIVER_FILE at another download (made with LR17_ibi.py, LR17_era5.py and
LR17_rivers.py for the same START END), and calls gendata.main() with the v2 defaults (IBI tide, no sed).
Note: gendata takes IBI mean sea level as the time mean of zos over the whole download, so a longer
download shifts MSL slightly relative to v2.

Usage: [LR17_S_OFFSET=1.3] python3 LR17_gendata_somlit.py RUN_START DAYS OUTDIR [DL_START DL_END]
                                 [--fes-boundary[=eta] [--band-download=START_END]] [--zb=FILE] [--sed]
  --fes-boundary, --band-download, --zb and --sed are passed to gendata.main (v3 baseline: --fes-boundary;
  v4 baseline: --fes-boundary=eta --zb=data/LR17/work/LR17_zb_v4.bin).
  e.g. python3 LR17_gendata_somlit.py 2022-02-23 43 ../LR17_model/input/run_2022-02-23_somlit 2022-02-20 2022-04-10
LR17_S_OFFSET_FILE (CSV date,offset) gives a daily, time-varying offset instead (takes precedence).
LR17_S_OFFSET (default 0) is added to IBI salinity everywhere before gendata uses it, so it reaches both the
initial S and the open-boundary S. IBI is 1.2-2.3 too fresh vs the SOMLIT Antioche CTD in every month of
2012-2026 (-1.35 in Feb-Apr 2022) and 1.68 too fresh vs the CocoriCO2 ROCHL mooring in Mar 2022.
LR17_S_OFFSET_IC (e.g. 1.3) is added to the initial S only (wet cells of S.bin), leaving the open-boundary S as IBI:
the 43-day runs showed the pertuis flushes IBI's fresh initial state within ~3 weeks, while an offset on the open
boundaries makes late-March salinity 0.7-1.8 too high (run_2022-02-23_somlit_Sic, Oct 2026).
"""
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import gendata  # noqa: E402
import xarray as xr  # noqa: E402

FES = next(('eta' if a == '--fes-boundary=eta' else 'all' for a in sys.argv
            if a in ('--fes-boundary', '--fes-boundary=eta')), False)
SED = '--sed' in sys.argv
BAND = next((a.split('=', 1)[1] for a in sys.argv if a.startswith('--band-download=')), None)
ZB = next((os.path.abspath(a.split('=', 1)[1]) for a in sys.argv if a.startswith('--zb=')), None)
argv = [a for a in sys.argv if not a.startswith('--')]
start, days, out = argv[1], float(argv[2]), os.path.abspath(argv[3])
s, e = (argv[4], argv[5]) if len(argv) > 5 else ('2022-02-20', '2022-04-10')
gendata.TAG = f'{s}_{e}'
gendata.RIVER_FILE = os.path.join(gendata.D, 'hubeau', f'LR17_rivers_{s}_{e}.json')
for f in [os.path.join(gendata.D, 'ibi', f'ibi_{k}_{gendata.TAG}.nc') for k in ('ssh', 'cur', 'cur3d', 'temp3d', 'sal')] + \
         [os.path.join(gendata.D, 'era5', f'era5_{gendata.TAG}_{k}.nc') for k in ('instant', 'accum')] + [gendata.RIVER_FILE]:
    if not os.path.exists(f):
        sys.exit('missing ' + f)
S_OFF = float(os.environ.get('LR17_S_OFFSET', '0'))
S_OFF_FILE = os.environ.get('LR17_S_OFFSET_FILE', '')      # CSV date,offset: daily offset (e.g. 7-day obs SSS - IBI)
if S_OFF_FILE:
    import pandas as pd
    _off = pd.read_csv(S_OFF_FILE, parse_dates=['date']).set_index('date')['offset']
    _ibi0 = gendata.ibi

    def ibi_offset_t(key):
        ds = _ibi0(key)
        if key == 'sal':
            days = pd.DatetimeIndex(ds.time.values).normalize()
            off = _off.reindex(days).interpolate(limit_direction='both').values
            ds['so'] = ds.so + xr.DataArray(off, dims='time', coords={'time': ds.time})
            print(f'IBI salinity offset from {os.path.basename(S_OFF_FILE)}: {off.min():+.2f}..{off.max():+.2f}')
        return ds
    gendata.ibi = ibi_offset_t
elif S_OFF:
    _ibi = gendata.ibi

    def ibi_offset(key):
        ds = _ibi(key)
        if key == 'sal':
            ds['so'] = ds.so + S_OFF
        return ds
    gendata.ibi = ibi_offset
    print(f'IBI salinity offset {S_OFF:+.2f}')
print(f'gendata with TAG={gendata.TAG}, rivers {os.path.basename(gendata.RIVER_FILE)}')
gendata.main(start, days, out, SED, False, False, None, FES, BAND, ZB)
S_OFF_IC = float(os.environ.get('LR17_S_OFFSET_IC', '0'))
if S_OFF_IC:
    import numpy as np
    fn = os.path.join(out, 'S.bin')
    S0 = np.fromfile(fn, '>f4')
    np.where(S0 != 0, S0 + np.float32(S_OFF_IC), S0).astype('>f4').tofile(fn)
    print(f'initial salinity offset {S_OFF_IC:+.2f} (open boundaries unchanged)')
