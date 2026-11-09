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
two as one zip. Needs ~/.cdsapirc. A span over several months is requested month by month
(one CDS request crosses years, months and days, so it would fetch far more than the span) and
joined, which keeps the request under the CDS size limit.
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


def by_month(start, end, out):
    import xarray as xr
    s, e = dt.date.fromisoformat(start), dt.date.fromisoformat(end)
    parts, a = [], s
    while a <= e:
        b = min((a.replace(day=1) + dt.timedelta(days=32)).replace(day=1) - dt.timedelta(days=1), e)
        parts.append((a.isoformat(), b.isoformat()))
        if not os.path.exists(os.path.join(out, f'era5_{parts[-1][0]}_{parts[-1][1]}_accum.nc')):
            main(*parts[-1], out)                    # finished months are kept, so a rerun resumes
        a = b + dt.timedelta(days=1)
    for kind in ('instant', 'accum'):
        fs = [os.path.join(out, f'era5_{p}_{q}_{kind}.nc') for p, q in parts]
        ds = xr.concat([xr.open_dataset(f) for f in fs], 'valid_time')
        ds.to_netcdf(os.path.join(out, f'era5_{start}_{end}_{kind}.nc'))
        print('joined', kind, dict(ds.sizes))
    for p, q in parts:
        for kind in ('instant', 'accum'):
            os.remove(os.path.join(out, f'era5_{p}_{q}_{kind}.nc'))


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
    d = sys.argv[3] if len(sys.argv) > 3 else os.path.join(TOP, 'data', 'LR17', 'era5')
    a, b = dt.date.fromisoformat(sys.argv[1]), dt.date.fromisoformat(sys.argv[2])
    ((by_month if (b - a).days > 62 else main))(sys.argv[1], sys.argv[2], d)
