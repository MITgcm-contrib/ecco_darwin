#!/usr/bin/env python3
"""Download ERA5 hourly single-level fields over the LR17 domain (CDS,
reanalysis-era5-single-levels) for the prototype's atmospheric forcing.

Usage:  python3 LR17_era5.py START END [OUTDIR]   (dates YYYY-MM-DD)
Box: lon -2.25 .. -0.25, lat 45.25 .. 47.0 (ERA5 0.25 deg; two cells of
margin for EXF interpolation). Variables: 10 m wind, 2 m temperature and
dewpoint, mean sea-level pressure, downward short- and long-wave radiation
(accumulated, J/m2 per hour) and total precipitation (m per hour).
Writes OUTDIR/era5_START_END_instant.nc (wind, temperature, dewpoint,
pressure) and ..._accum.nc (radiation, precipitation); the CDS delivers the
two as one zip. Needs ~/.cdsapirc.
"""
import datetime as dt
import os
import sys
import zipfile

import cdsapi

TOP = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
VARS = ['10m_u_component_of_wind', '10m_v_component_of_wind', '2m_temperature',
        '2m_dewpoint_temperature', 'mean_sea_level_pressure',
        'surface_solar_radiation_downwards', 'surface_thermal_radiation_downwards',
        'total_precipitation']


def main(start, end, out):
    os.makedirs(out, exist_ok=True)
    s, e = dt.date.fromisoformat(start), dt.date.fromisoformat(end)
    days = [s + dt.timedelta(days=k) for k in range((e - s).days + 1)]
    fn = os.path.join(out, f'era5_{start}_{end}.zip')
    req = {
        'product_type': ['reanalysis'], 'variable': VARS,
        'year': sorted({f'{d.year}' for d in days}),
        'month': sorted({f'{d.month:02d}' for d in days}),
        'day': sorted({f'{d.day:02d}' for d in days}),
        'time': [f'{h:02d}:00' for h in range(24)],
        'area': [47.0, -2.25, 45.25, -0.25],       # N, W, S, E
        'data_format': 'netcdf', 'download_format': 'unarchived',
    }
    cdsapi.Client().retrieve('reanalysis-era5-single-levels', req, fn)
    with zipfile.ZipFile(fn) as z:
        for n in z.namelist():
            kind = 'accum' if 'accum' in n else 'instant'
            dst = os.path.join(out, f'era5_{start}_{end}_{kind}.nc')
            open(dst, 'wb').write(z.read(n))
            print('wrote', dst, f'{os.path.getsize(dst) / 1e6:.1f} MB')
    os.remove(fn)


if __name__ == '__main__':
    main(sys.argv[1], sys.argv[2],
         sys.argv[3] if len(sys.argv) > 3 else os.path.join(TOP, 'data', 'LR17', 'era5'))
