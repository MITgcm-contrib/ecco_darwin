"""Gridded runoff / discharge fields -> SourceSet.

For products that already put river water on a grid at the coast -- e.g. JRA55-do
`friver` (kg m-2 s-1 on a 0.25 deg grid, non-zero at river-mouth cells), a GloFAS-style
"runoff map", or another model's coastal runoff. Every cell that carries water on ANY
record becomes a source at its cell centre; its discharge is value x area (flux units)
or the value itself (m3 s-1).

Do NOT feed a full routing-model discharge field (every land cell, e.g. raw GloFAS
avg_dis) to this reader: each river would be counted once per cell along its course.
Use readers.glofas (outlets) for that.
"""
import numpy as np

from ..geometry import regular_cell_area, wrap_lon
from ..sources import SourceSet

_FLUX = {"kg m-2 s-1": 1.0e-3, "kg/m2/s": 1.0e-3, "kg m**-2 s**-1": 1.0e-3,
         "m s-1": 1.0, "m/s": 1.0, "mm/day": 1.0e-3 / 86400.0, "mm d-1": 1.0e-3 / 86400.0}
_VOLUME = {"m3 s-1": 1.0, "m3/s": 1.0, "m**3 s**-1": 1.0}


def from_field(lon, lat, field, time, units, area=None, min_value=0.0, name_prefix="cell"):
    """lon, lat : 1-D (regular grid) or 2-D cell centres; field (nt, ny, nx).
    units: one of the flux units (kg m-2 s-1, m s-1, mm/day: converted with the cell area,
    water density 1000 kg m-3) or volume units (m3 s-1)."""
    field = np.asarray(field, dtype=float)
    if field.ndim == 2:
        field = field[None]
    lon = np.asarray(lon, float)
    lat = np.asarray(lat, float)
    if lon.ndim == 1 and lat.ndim == 1:
        if area is None:
            area = regular_cell_area(lon, lat)
        lon, lat = np.meshgrid(lon, lat)
    u = units.strip()
    if u in _FLUX:
        if area is None:
            raise ValueError("a flux field needs cell areas (pass area=)")
        fac = np.asarray(area, float) * _FLUX[u]
    elif u in _VOLUME:
        fac = np.ones(lon.shape) * _VOLUME[u]
    else:
        raise ValueError(f"unknown units {units!r}; known: {sorted(_FLUX) + sorted(_VOLUME)}")
    clean = np.where(np.isfinite(field), field, 0.0)
    carry = (clean > min_value).any(axis=0)
    jj, ii = np.nonzero(carry)
    q = clean[:, jj, ii] * fac[jj, ii]
    names = [f"{name_prefix}_{j}_{i}" for j, i in zip(jj, ii)]
    return SourceSet(wrap_lon(lon[jj, ii]), lat[jj, ii], q, np.asarray(time), name=names,
                     kind=np.full(jj.size, "gridded", dtype=object),
                     attrs=dict(dataset="gridded field", units=units))


def from_netcdf(path, var, lon="lon", lat="lat", time="time", units=None, area_var=None,
                records=None, min_value=0.0):
    """Read a gridded field from netCDF (units from the file's attribute unless given)."""
    from netCDF4 import Dataset, num2date
    with Dataset(path) as ds:
        v = ds[var]
        sl = slice(None) if records is None else records
        f = np.ma.filled(v[sl], np.nan).astype(float)
        LON = np.asarray(ds[lon][:], float)
        LAT = np.asarray(ds[lat][:], float)
        tv = ds[time]
        t = num2date(tv[sl], tv.units, getattr(tv, "calendar", "standard"),
                     only_use_cftime_datetimes=False, only_use_python_datetimes=True)
        tt = np.array([np.datetime64(x, "s") for x in np.atleast_1d(t)])
        A = np.asarray(ds[area_var][:], float) if area_var else None
        u = units or getattr(v, "units", None)
    if u is None:
        raise ValueError(f"{var} has no units attribute; pass units=")
    s = from_field(LON, LAT, f, tt, u, area=A, min_value=min_value)
    s.attrs.update(file=path, variable=var)
    return s
