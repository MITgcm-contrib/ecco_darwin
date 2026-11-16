"""River water temperature -> SourceSet.temp (deg C).

Discharge products rarely carry water temperature (GloFAS has none), so temperature comes
from a separate gridded product and is SAMPLED at the river sources:

    clim = TemperatureClimatology.from_binary(sorted(glob("wattmp*.bin")), **WATTMP)   # or from_netcdf
    clim.save("river_temp_clim.npz")             # build once (it reads every file)
    src  = add_temperature(src, TemperatureClimatology.load("river_temp_clim.npz"))
    src, donors = fill_missing(src, max_km=150)  # sources with no value borrow a neighbour's

The climatology is what our production forcing used (make_glofas_river_temp_forcing.py):
  - every daily record of every file is placed on a CALENDAR slot (Jan 1 = 0 ... Dec 31 = 365,
    Feb 29 = 59, Mar 1 = 60 in every year), so 365- and 366-day files line up;
  - a value counts if it is finite, not a fill value, and 0-45 deg C after conversion;
  - a cell is usable if it has >= 30 good values and a non-zero mean; slots it lacks are filled
    with its overall mean;
  - each source takes the value of the NEAREST usable cell (great-circle distance) on the
    calendar slot of each record.
A climatology repeats the same seasonal cycle every year. TemperatureSeries samples a dated
(year-specific) product instead.

Our source was a JRA55-do-grid river temperature product (wattmp<year>.bin). Its format
broke every convention of the other files we used, and each mistake gave plausible-looking
values -- check a new product the same way (map the valid cells, compare known rivers):
    little-endian float32 (big-endian reads give ~85 % denormal garbage), (365|366, 720, 1440),
    row 0 = 89.875 N (north first), lon -179.875 .. 179.875 by 0.25, Kelvin, 1e20 over the sea.
"""
import datetime
import os

import numpy as np

from ..geometry import to_xyz, chord_to_km, km_to_chord

WATTMP = dict(nlon=1440, nlat=720, lon0=-179.875, lat_first=89.875, inc=0.25, dtype="<f4",
              north_first=True, kelvin=True)          # the wattmp<year>.bin layout (see above)
_CANON_LEAP = 2000


def calendar_slot(dates):
    """datetime64 dates -> calendar slot 0..365 (Jan 1 = 0, Feb 29 = 59, Mar 1 = 60, Dec 31 = 365)."""
    d = np.asarray(dates).astype("datetime64[D]")
    m = (d.astype("datetime64[M]").astype(int) % 12) + 1
    day = (d - d.astype("datetime64[M]")).astype(int) + 1
    base = datetime.date(_CANON_LEAP, 1, 1)
    return np.array([(datetime.date(_CANON_LEAP, int(a), int(b)) - base).days for a, b in zip(m, day)],
                    dtype=np.int64)


class TemperatureClimatology:
    """Daily calendar climatology of river water temperature on a regular lat-lon grid.

    lon (nx,), lat (ny,) ascending, clim (366, ny, nx) float32 deg C, valid (ny, nx) bool."""

    def __init__(self, lon, lat, clim, valid, attrs=None):
        self.lon = np.asarray(lon, float)
        self.lat = np.asarray(lat, float)
        self.clim = np.asarray(clim, np.float32)
        self.valid = np.asarray(valid, bool)
        self.attrs = dict(attrs or {})
        self._tree = None

    # ------------------------------------------------------------------ build
    @classmethod
    def _accumulate(cls, records, lon, lat, kelvin, valid_range, fill_values, min_good, verbose):
        ny, nx = lat.size, lon.size
        csum = np.zeros((366, ny, nx), dtype=np.float64)
        ccnt = np.zeros((366, ny, nx), dtype=np.int32)
        nrec = 0
        for slot, rec in records:
            rec = np.asarray(rec, dtype=np.float32)
            if kelvin:
                rec = rec - 273.15                    # float32 arithmetic, as in production
            bad = (~np.isfinite(rec)) | (rec < valid_range[0]) | (rec > valid_range[1]) \
                | np.isin(rec, fill_values)
            good = ~bad
            csum[slot] += np.where(good, rec, 0.0)
            ccnt[slot] += good
            nrec += 1
        clim = np.divide(csum, ccnt, out=np.zeros_like(csum), where=ccnt > 0)
        total = ccnt.sum(axis=0)
        cell_mean = np.divide(csum.sum(axis=0), total, out=np.zeros((ny, nx)), where=total > 0)
        nonzero = np.zeros((ny, nx), dtype=bool)
        for s in range(366):
            nonzero |= (ccnt[s] > 0) & (clim[s] != 0)
        valid = (total >= min_good) & nonzero
        for s in range(366):
            gap = valid & (ccnt[s] == 0)
            clim[s] = np.where(gap, cell_mean, clim[s])
        if verbose:
            mv = cell_mean[valid]
            print(f"  climatology from {nrec:,} records: {int(valid.sum()):,} usable cells of {valid.size:,}"
                  + (f"; cell means {mv.min():.2f} .. {mv.max():.2f} C" if mv.size else ""), flush=True)
        return clim.astype(np.float32), valid

    @classmethod
    def from_binary(cls, paths, nlon, nlat, lon0, lat_first, inc, dtype="<f4", north_first=True,
                    kelvin=True, valid_range=(0.0, 45.0),
                    fill_values=(-9999.0, -999.0, 1e20, 1e35), min_good=30, lat_range=None,
                    verbose=True):
        """From flat binary files of daily records (nrec, nlat, nlon), 365 or 366 records each
        (a 365-record file is placed on the non-leap calendar, 366 on the leap calendar).
        lat_range=(S, N) keeps only those rows (pad it: sources near the edge need neighbours)."""
        lat_rows = lat_first + (-inc if north_first else inc) * np.arange(nlat)   # file order
        lon = lon0 + inc * np.arange(nlon)
        sel = np.arange(nlat) if lat_range is None else \
            np.flatnonzero((lat_rows >= lat_range[0]) & (lat_rows <= lat_range[1]))
        order = sel[np.argsort(lat_rows[sel])]                                    # ascending lat
        lat = lat_rows[order]
        r0, r1 = int(sel.min()), int(sel.max()) + 1
        rec_bytes = nlon * nlat * np.dtype(dtype).itemsize

        def records():
            for p in paths:
                size = os.path.getsize(p)
                if size % rec_bytes:
                    raise ValueError(f"{p}: {size} bytes is not a whole number of records")
                n = size // rec_bytes
                if n not in (365, 366):
                    if verbose:
                        print(f"  {os.path.basename(p)}: {n} records, not 365/366 -- skipped")
                    continue
                d0 = datetime.date(_CANON_LEAP if n == 366 else _CANON_LEAP + 1, 1, 1)
                mm = np.memmap(p, dtype=dtype, mode="r", shape=(n, nlat, nlon))
                for k in range(n):
                    d = d0 + datetime.timedelta(days=k)
                    slot = (datetime.date(_CANON_LEAP, d.month, d.day) - datetime.date(_CANON_LEAP, 1, 1)).days
                    yield slot, np.asarray(mm[k, r0:r1, :])[order - r0]
                del mm
                if verbose:
                    print(f"  {os.path.basename(p)}: {n} records", flush=True)

        clim, valid = cls._accumulate(records(), lon, lat, kelvin, valid_range, fill_values,
                                      min_good, verbose)
        return cls(lon, lat, clim, valid, attrs=dict(source="binary", files=len(paths),
                                                     kelvin=int(kelvin)))

    @classmethod
    def from_netcdf(cls, paths, var, lon="lon", lat="lat", time="time", kelvin=None,
                    valid_range=(0.0, 45.0), fill_values=(-9999.0, -999.0, 1e20, 1e35),
                    min_good=30, lat_range=None, verbose=True):
        """From netCDF files of dated records (time, lat, lon). kelvin=None reads the units."""
        from netCDF4 import Dataset, num2date
        if isinstance(paths, str):
            paths = [paths]
        with Dataset(paths[0]) as ds:
            LON = np.asarray(ds[lon][:], float)
            LAT = np.asarray(ds[lat][:], float)
            u = getattr(ds[var], "units", "")
        if kelvin is None:
            kelvin = u.strip().lower() in ("k", "kelvin", "degk", "deg_k")
        sel = np.arange(LAT.size) if lat_range is None else \
            np.flatnonzero((LAT >= lat_range[0]) & (LAT <= lat_range[1]))
        order = sel[np.argsort(LAT[sel])]

        def records():
            for p in paths:
                with Dataset(p) as ds:
                    tv = ds[time]
                    tt = num2date(tv[:], tv.units, getattr(tv, "calendar", "standard"),
                                  only_use_cftime_datetimes=False, only_use_python_datetimes=True)
                    slots = calendar_slot(np.array([np.datetime64(x, "D") for x in tt]))
                    v = ds[var]
                    for k, s in enumerate(slots):
                        rec = np.ma.filled(v[k, :, :], np.nan)
                        yield s, np.asarray(rec)[order]
                if verbose:
                    print(f"  {os.path.basename(p)}: {len(slots)} records", flush=True)

        clim, valid = cls._accumulate(records(), LON, LAT[order], kelvin, valid_range, fill_values,
                                      min_good, verbose)
        return cls(LON, LAT[order], clim, valid, attrs=dict(source="netcdf", variable=var,
                                                             kelvin=int(kelvin)))

    # ------------------------------------------------------------------ io
    def save(self, path):
        np.savez(path, lon=self.lon, lat=self.lat, clim=self.clim, valid=self.valid,
                 attrs=np.array(repr(self.attrs)))
        return path

    @classmethod
    def load(cls, path):
        z = np.load(path, allow_pickle=False)
        return cls(z["lon"], z["lat"], z["clim"], z["valid"])

    # ------------------------------------------------------------------ sample
    def _lookup(self):
        if self._tree is None:
            from scipy.spatial import cKDTree
            rows, cols = np.nonzero(self.valid)
            LON = self.lon[cols]
            LAT = self.lat[rows]
            self._tree = (cKDTree(to_xyz(LON, LAT)), rows, cols)
        return self._tree

    def sample(self, lon, lat, time, max_km=None):
        """(temperature (nt, n) deg C, distance to the cell used (n,) km). NaN beyond max_km."""
        tree, rows, cols = self._lookup()
        d, k = tree.query(to_xyz(np.asarray(lon, float), np.asarray(lat, float)), k=1, workers=-1)
        dkm = chord_to_km(d)
        slots = calendar_slot(time)
        out = self.clim[slots[:, None], rows[k][None, :], cols[k][None, :]].astype(float)
        if max_km is not None:
            out[:, dkm > max_km] = np.nan
        return out, dkm


class TemperatureSeries:
    """Dated (year-specific) gridded river temperature, sampled at the nearest usable cell.
    Records are matched to the sources' dates exactly; missing dates give NaN."""

    def __init__(self, lon, lat, data, time, valid):
        self.lon, self.lat = np.asarray(lon, float), np.asarray(lat, float)
        self.data = data                    # (nt, ny, nx) deg C
        self.time = np.asarray(time).astype("datetime64[D]")
        self.valid = valid
        self._tree = None

    @classmethod
    def from_netcdf(cls, path, var, lon="lon", lat="lat", time="time", kelvin=None,
                    valid_range=(-2.0, 45.0), min_frac=0.5):
        from netCDF4 import Dataset, num2date
        with Dataset(path) as ds:
            LON = np.asarray(ds[lon][:], float)
            LAT = np.asarray(ds[lat][:], float)
            v = ds[var]
            u = getattr(v, "units", "")
            data = np.ma.filled(v[:], np.nan).astype(float)
            tv = ds[time]
            tt = num2date(tv[:], tv.units, getattr(tv, "calendar", "standard"),
                          only_use_cftime_datetimes=False, only_use_python_datetimes=True)
        if kelvin is None:
            kelvin = u.strip().lower() in ("k", "kelvin", "degk", "deg_k")
        if kelvin:
            data = data - 273.15
        data[(data < valid_range[0]) | (data > valid_range[1])] = np.nan
        valid = np.isfinite(data).mean(axis=0) >= min_frac
        if LON.ndim != 1 or LAT.ndim != 1:
            raise ValueError("TemperatureSeries expects 1-D lon/lat")
        return cls(LON, LAT, data, np.array([np.datetime64(x, "D") for x in tt]), valid)

    def sample(self, lon, lat, time, max_km=None):
        from scipy.spatial import cKDTree
        if self._tree is None:
            rows, cols = np.nonzero(self.valid)
            self._tree = (cKDTree(to_xyz(self.lon[cols], self.lat[rows])), rows, cols)
        tree, rows, cols = self._tree
        d, k = tree.query(to_xyz(np.asarray(lon, float), np.asarray(lat, float)), k=1, workers=-1)
        dkm = chord_to_km(d)
        pos = {t: i for i, t in enumerate(self.time)}
        idx = np.array([pos.get(t, -1) for t in np.asarray(time).astype("datetime64[D]")])
        out = np.full((idx.size, k.size), np.nan)
        ok = idx >= 0
        out[ok] = self.data[idx[ok][:, None], rows[k][None, :], cols[k][None, :]]
        if max_km is not None:
            out[:, dkm > max_km] = np.nan
        return out, dkm


def add_temperature(sources, sampler, max_km=None):
    """A copy of `sources` with temperature sampled from a TemperatureClimatology/Series."""
    t, dkm = sampler.sample(sources.lon, sources.lat, sources.time, max_km=max_km)
    s = sources.with_temp(t)
    s.attrs["temperature"] = (f"sampled at the nearest usable cell (median {np.median(dkm):.1f} km, "
                              f"max {dkm.max():.1f} km)")
    return s


def fill_missing(sources, max_km=150.0):
    """Sources whose temperature is unknown (NaN) take the daily temperature of the LARGEST
    (mean discharge) source with known temperature within max_km -- the rule the production
    forcing used for the rivers added by hand (their water was given the temperature of the
    largest real outlet within 1.5 degrees). Returns (new SourceSet, {i: donor or -1})."""
    from scipy.spatial import cKDTree
    if not sources.has_temp:
        raise ValueError("sources have no temperature to fill")
    T = sources.temp.copy()
    need = np.flatnonzero(~np.isfinite(T).all(axis=0))
    known = np.flatnonzero(np.isfinite(T).any(axis=0))
    donors = {}
    if need.size and known.size:
        mq = sources.mean_q()
        tree = cKDTree(to_xyz(sources.lon[known], sources.lat[known]))
        xyz = to_xyz(sources.lon[need], sources.lat[need])
        for i, cand in zip(need, tree.query_ball_point(xyz, km_to_chord(max_km))):
            cand = [known[c] for c in cand if known[c] != i]
            if not cand:
                donors[int(i)] = -1
                continue
            j = max(cand, key=lambda c: mq[c])
            donors[int(i)] = int(j)
            T[:, i] = np.where(np.isfinite(T[:, i]), T[:, i], T[:, j])
    else:
        donors = {int(i): -1 for i in need}
    s = sources.with_temp(T)
    left = int((np.array(list(donors.values())) < 0).sum()) if donors else 0
    s.attrs["temperature_fill"] = (f"{len(donors)} sources with missing temperature; "
                                   f"{len(donors) - left} filled from a neighbour within {max_km:g} km, "
                                   f"{left} left unknown")
    return s, donors
