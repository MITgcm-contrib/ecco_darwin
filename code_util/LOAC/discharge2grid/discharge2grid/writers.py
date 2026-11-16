"""Output writers.

write_mitgcm_runoff  dense big-endian float32 records (nt, *grid.shape), one file per year:
                     what MITgcm/EXF reads (runoffFile, useExfYearlyFields). Units m s-1.

write_sparse_netcdf  only the cells that ever receive water ("points"), so the file holds
                     no zeros outside the river footprints:
                         runoff(time, point)  m s-1      the forcing itself
                         point_index(point)          flat (C-order) index into the grid
                         <dim>(point)                one 0-based index per grid dim,
                                                     e.g. tile, j, i for an LLC file
                         lon, lat, area (point)
                     plus, for traceability, the sources and the sparse operator:
                         source_lon/lat/name/kind/status/snap_km/dest_point (source)
                         discharge(time, source)     m3 s-1 (optional, --with-discharge)
                         link_point, link_source, link_weight (link)  -> runoff = W q
                     and, when the sources carry temperature:
                         temperature(time, point)    deg C, flux-weighted; NaN where no water
                         source_temperature(time, source)
                     Reconstruct a dense field with read_sparse_netcdf(...).to_dense(k).

write_mitgcm_temperature  dense big-endian float32 river temperature (deg C), 0 where no
                     water, the companion of the runoff file (EXF runoftempfile).
"""
import datetime

import numpy as np

from .mapping import STATUS_NAMES


def write_mitgcm_runoff(path, delivery, dtype=">f4", records=None):
    """Dense runoff (m s-1), all records (or the given record indices), grid shape per record.
    Returns the number of records written."""
    recs = range(delivery.sources.nt) if records is None else records
    n = 0
    with open(path, "wb") as fh:
        for k in recs:
            fh.write(delivery.runoff_dense(k).astype(dtype).tobytes())
            n += 1
    return n


def write_mitgcm_temperature(path, delivery, dtype=">f4", records=None, fill=0.0):
    """Dense river temperature (deg C), flux-weighted per cell; `fill` (0, as in the
    production files) where a cell gets no water or only water of unknown temperature."""
    recs = range(delivery.sources.nt) if records is None else records
    n = 0
    with open(path, "wb") as fh:
        for k in recs:
            fh.write(delivery.temperature_dense(k, fill=fill).astype(dtype).tobytes())
            n += 1
    return n


def _time_values(time):
    """datetime64 array -> (float days since 1970-01-01, units string)."""
    t = np.asarray(time).astype("datetime64[s]").astype("int64") / 86400.0
    return t, "days since 1970-01-01 00:00:00"


def write_sparse_netcdf(path, delivery, title="river runoff", with_discharge=True,
                        with_operator=True, attrs=None, complevel=4, float_type="f4"):
    """Point-format runoff file (see module docstring). Returns path."""
    from netCDF4 import Dataset
    g = delivery.mapper.grid
    s, sn = delivery.sources, delivery.snap
    npt, ns = delivery.n_point, s.n
    with Dataset(path, "w", format="NETCDF4") as ds:
        ds.title = title
        ds.Conventions = "CF-1.8 (point list on a model grid; see description)"
        ds.description = ("Runoff on the model cells that receive river water. runoff(time, point) "
                          "is the volume flux per unit cell area (m s-1); sum(runoff*area) over "
                          "points equals the discharge of the kept sources. point_index is the "
                          "flat C-order index into the grid array of shape grid_shape; the "
                          "per-dimension indices are 0-based.")
        ds.grid_name = g.name
        ds.grid_shape = np.asarray(g.shape, dtype=np.int64)
        ds.grid_dims = " ".join(g.dims)
        for k, v in delivery.mapper.describe().items():
            setattr(ds, k, v)
        ds.history = f"{datetime.datetime.now(datetime.timezone.utc):%Y-%m-%d %H:%M} UTC discharge2grid"
        for k, v in (attrs or {}).items():
            setattr(ds, k, v)
        for k, v in s.attrs.items():
            setattr(ds, f"source_{k}", str(v))

        ds.createDimension("time", None)
        ds.createDimension("point", npt)
        ds.createDimension("source", ns)
        tv, tu = _time_values(s.time)
        t = ds.createVariable("time", "f8", ("time",))
        t.units, t.calendar, t.standard_name = tu, "proleptic_gregorian", "time"
        t[:] = tv

        pi = ds.createVariable("point_index", "i8", ("point",))
        pi.long_name = "flat C-order index into the model grid array (grid_shape)"
        pi[:] = delivery.points
        for name, ix in zip(g.dims, g.unravel(delivery.points)):
            v = ds.createVariable(name, "i4", ("point",))
            v.long_name = f"0-based index along grid dimension {name!r}"
            v[:] = ix
        for name, arr, units in (("lon", g.lon.ravel()[delivery.points], "degrees_east"),
                                 ("lat", g.lat.ravel()[delivery.points], "degrees_north"),
                                 ("area", delivery.point_area, "m2")):
            v = ds.createVariable(name, "f8", ("point",))
            v.units = units
            v[:] = arr

        r = ds.createVariable("runoff", float_type, ("time", "point"), zlib=complevel > 0,
                              complevel=complevel, chunksizes=(1, max(npt, 1)))
        r.units = "m s-1"
        r.long_name = "river runoff as volume flux per unit cell area"
        r.coordinates = "lon lat"
        for k in range(s.nt):
            r[k, :] = delivery.runoff_points(k)
        if s.has_temp:
            tv = ds.createVariable("temperature", float_type, ("time", "point"), zlib=complevel > 0,
                                   complevel=complevel, chunksizes=(1, max(npt, 1)),
                                   fill_value=np.float32(np.nan))
            tv.units = "degree_Celsius"
            tv.long_name = ("river water temperature, flux-weighted over the water reaching the "
                            "cell; NaN where the cell gets no water (or only water of unknown "
                            "temperature) on that record")
            tv.coordinates = "lon lat"
            for k in range(s.nt):
                tv[k, :] = delivery.temperature_points(k, fill=np.nan)

        # sources
        for name, arr, typ, units in (("source_lon", s.lon, "f8", "degrees_east"),
                                      ("source_lat", s.lat, "f8", "degrees_north"),
                                      ("source_snap_km", sn.snap_km, "f4", "km")):
            v = ds.createVariable(name, typ, ("source",))
            v.units = units
            v[:] = arr
        st = ds.createVariable("source_status", "i1", ("source",))
        st.flag_values = np.array(sorted(STATUS_NAMES), dtype=np.int8)
        st.flag_meanings = " ".join(STATUS_NAMES[k] for k in sorted(STATUS_NAMES))
        st[:] = sn.status
        # destination cell of each kept source, as an index into "point" (-1 if dropped)
        pos = {int(w): j for j, w in enumerate(delivery.point_wet)}
        dp = ds.createVariable("source_dest_point", "i4", ("source",))
        dp.long_name = "point index of the cell each kept source was snapped to (-1 = dropped)"
        dp[:] = np.array([pos.get(int(d), -1) if d >= 0 else -1 for d in sn.dest], dtype=np.int32)
        for name, arr in (("source_name", s.name), ("source_kind", s.kind)):
            v = ds.createVariable(name, str, ("source",))
            for j, x in enumerate(arr):
                v[j] = str(x)
        if with_discharge:
            q = ds.createVariable("discharge", "f4", ("time", "source"), zlib=complevel > 0,
                                  complevel=complevel)
            q.units = "m3 s-1"
            q.long_name = "source discharge as read from the dataset (all sources, kept or not)"
            q[:] = s.q
        if s.has_temp:
            st_ = ds.createVariable("source_temperature", "f4", ("time", "source"), zlib=complevel > 0,
                                    complevel=complevel, fill_value=np.float32(np.nan))
            st_.units = "degree_Celsius"
            st_.long_name = "source water temperature as assigned (NaN = unknown)"
            st_[:] = s.temp
        if with_operator:
            Wc = delivery.Wp.tocoo()
            ds.createDimension("link", Wc.nnz)
            for name, arr, typ in (("link_point", Wc.row, "i4"), ("link_source", Wc.col, "i4")):
                v = ds.createVariable(name, typ, ("link",), zlib=True)
                v[:] = arr
            w = ds.createVariable("link_weight", "f8", ("link",), zlib=True)
            w.units = "m-2"
            w.long_name = "runoff[point] = sum over links of link_weight * discharge[source]"
            w[:] = Wc.data
    return path


class SparseRunoff:
    """What read_sparse_netcdf returns."""

    def __init__(self, time, points, runoff, area, grid_shape, attrs, temperature=None):
        self.time, self.points, self.runoff, self.area = time, points, runoff, area
        self.grid_shape, self.attrs = tuple(int(x) for x in grid_shape), attrs
        self.temperature = temperature

    def to_dense(self, k, dtype=float):
        out = np.zeros(int(np.prod(self.grid_shape)), dtype=dtype)
        out[self.points] = self.runoff[k]
        return out.reshape(self.grid_shape)

    def temperature_dense(self, k, fill=0.0, dtype=float):
        """Dense temperature record (fill where NaN / no water), like the MITgcm file."""
        out = np.full(int(np.prod(self.grid_shape)), fill, dtype=dtype)
        t = np.asarray(self.temperature[k], dtype=float)
        out[self.points] = np.where(np.isfinite(t), t, fill)
        return out.reshape(self.grid_shape)

    def total_m3s(self, k):
        return float((np.asarray(self.runoff[k], float) * self.area).sum())


def read_sparse_netcdf(path):
    from netCDF4 import Dataset, num2date
    with Dataset(path) as ds:
        tv = ds["time"]
        t = num2date(tv[:], tv.units, getattr(tv, "calendar", "standard"),
                     only_use_cftime_datetimes=False, only_use_python_datetimes=True)
        time = np.array([np.datetime64(x) for x in t])
        temp = (np.asarray(np.ma.filled(ds["temperature"][:], np.nan))
                if "temperature" in ds.variables else None)
        return SparseRunoff(time, np.asarray(ds["point_index"][:]),
                            np.asarray(ds["runoff"][:]), np.asarray(ds["area"][:]),
                            np.asarray(ds.grid_shape).ravel(),
                            {k: ds.getncattr(k) for k in ds.ncattrs()}, temperature=temp)
