"""Station / river-mouth datasets -> SourceSet.

Many discharge products are already a list of rivers with a mouth position and a time
series (Dai & Trenberth, ArcticGRO, national gauge networks, hydrological model output at
outlets). Two generic readers:

from_csv(outlets_csv, series_csv)
    outlets_csv : one row per river, columns  name, lon, lat  [, kind]
    series_csv  : first column a date (YYYY-MM-DD or YYYY-MM), then one column per river
                  NAME (m3 s-1). Rows are records (daily, monthly, ...).
    temp_csv    : optional, same layout as series_csv, water temperature in deg C
                  (blank = unknown; fill with readers.temperature.fill_missing).

from_netcdf(path, q_var, lon_var, lat_var, ...)
    a netCDF with discharge (time, station) or (station, time) and station lon/lat.

Gauges are usually upstream of the mouth: if a station is far inland, pass the MOUTH
position (and scale the discharge to the mouth if you know the ungauged fraction).
"""
import csv

import numpy as np

from ..sources import SourceSet


def _parse_date(s):
    s = s.strip()
    if len(s) == 7:            # YYYY-MM
        s = s + "-01"
    return np.datetime64(s[:10], "D")


def _read_wide(path, names, date_col=None):
    with open(path, newline="") as fh:
        rd = csv.reader(fh)
        head = [h.strip() for h in next(rd)]
        data = [r for r in rd if r and r[0].strip()]
    dc = 0 if date_col is None else head.index(date_col)
    missing = [n for n in names if n not in head]
    if missing:
        raise ValueError(f"{path} has no column for: {missing}")
    col = [head.index(n) for n in names]
    time = np.array([_parse_date(r[dc]) for r in data])
    val = np.array([[float(r[c]) if r[c].strip() not in ("", "nan", "NaN", "NA") else np.nan
                     for c in col] for r in data])
    return time, val


def from_csv(outlets_csv, series_csv, scale=1.0, date_col=None, temp_csv=None):
    """SourceSet from an outlet table and a wide time-series table (see module doc)."""
    with open(outlets_csv, newline="") as fh:
        rows = list(csv.DictReader(fh))
    names = [r["name"].strip() for r in rows]
    lon = np.array([float(r["lon"]) for r in rows])
    lat = np.array([float(r["lat"]) for r in rows])
    kind = [r.get("kind", "station") or "station" for r in rows]
    time, q = _read_wide(series_csv, names, date_col)
    q = q * scale
    temp = None
    if temp_csv:
        tt, temp = _read_wide(temp_csv, names, date_col)
        if not np.array_equal(tt, time):
            raise ValueError(f"{temp_csv}: dates differ from {series_csv}")
    return SourceSet(lon, lat, q, time, name=names, kind=kind, temp=temp,
                     attrs=dict(dataset="table", outlets=outlets_csv, series=series_csv))


def from_netcdf(path, q_var, lon_var, lat_var, time_var="time", name_var=None,
                station_dim=None, scale=1.0):
    """SourceSet from a station netCDF. q_var may be (time, station) or (station, time)."""
    from netCDF4 import Dataset, num2date
    with Dataset(path) as ds:
        v = ds[q_var]
        q = np.ma.filled(v[:], np.nan).astype(float)
        sdim = station_dim or ds[lon_var].dimensions[0]
        if v.dimensions[0] == sdim:
            q = q.T
        lon = np.asarray(ds[lon_var][:], float)
        lat = np.asarray(ds[lat_var][:], float)
        tv = ds[time_var]
        t = num2date(tv[:], tv.units, getattr(tv, "calendar", "standard"),
                     only_use_cftime_datetimes=False, only_use_python_datetimes=True)
        time = np.array([np.datetime64(x, "D") for x in t])
        names = None
        if name_var:
            raw = ds[name_var][:]
            if getattr(raw, "dtype", None) is not None and raw.dtype.kind == "S":
                from netCDF4 import chartostring
                names = chartostring(raw)
            else:
                names = np.asarray(raw).astype(str)
    return SourceSet(lon, lat, q * scale, time, name=names,
                     kind=np.full(lon.size, "station", dtype=object),
                     attrs=dict(dataset="station netCDF", file=path, variable=q_var))
