"""SourceSet: the hand-over format between a dataset reader and the regridder.

Every reader in discharge2grid.readers, whatever the dataset (gridded routing model,
station table, coastal runoff field), produces the same thing: n point sources with a
location and a discharge time series in m3/s. The regridder never sees the dataset.

    lon, lat : (n,)      where the water ENTERS the ocean (or should be delivered)
    q        : (nt, n)   discharge, m3 s-1; NaN and negative values are set to 0
    time     : (nt,)     numpy datetime64 (start of each record; daily, monthly, ...)
    name     : (n,)      optional labels
    kind     : (n,)      optional provenance labels ('coastal', 'dead_end', 'station', ...)
    temp     : (nt, n)   optional river water temperature, deg C; NaN = unknown
                         (readers.temperature samples a gridded product at the sources)

Keep the set of sources FIXED across all records of a run (a source with 0 discharge on
some days is fine). That is what makes the sparse output's point list fixed in time, and
it is also the cure for outlets that drop in and out of a product from year to year.
"""
from dataclasses import dataclass, field

import numpy as np


@dataclass
class SourceSet:
    lon: np.ndarray
    lat: np.ndarray
    q: np.ndarray
    time: np.ndarray
    name: np.ndarray = None
    kind: np.ndarray = None
    attrs: dict = field(default_factory=dict)
    temp: np.ndarray = None

    def __post_init__(self):
        self.lon = np.asarray(self.lon, dtype=float).ravel()
        self.lat = np.asarray(self.lat, dtype=float).ravel()
        q = np.asarray(self.q, dtype=float)
        if q.ndim == 1:
            q = q[:, None] if self.lon.size == 1 else q[None, :]
        n = self.lon.size
        if self.lat.size != n or q.shape[1] != n:
            raise ValueError(f"lon {self.lon.size}, lat {self.lat.size}, q {q.shape}: "
                             f"q must be (ntime, nsource)")
        self.time = np.asarray(self.time)
        if self.time.size != q.shape[0]:
            raise ValueError(f"time has {self.time.size} records, q has {q.shape[0]}")
        bad = ~np.isfinite(q)
        neg = np.isfinite(q) & (q < 0)
        if bad.any() or neg.any():
            q = np.where(bad | neg, 0.0, q)
            self.attrs.setdefault("cleaned", f"{int(bad.sum())} non-finite and {int(neg.sum())} "
                                             f"negative values set to 0")
        self.q = q
        self.name = (np.asarray([f"src_{x:.3f}_{y:.3f}" for x, y in zip(self.lon, self.lat)])
                     if self.name is None else np.asarray(self.name).astype(str).ravel())
        self.kind = (np.full(n, "", dtype=object) if self.kind is None
                     else np.asarray(self.kind).astype(str).ravel())
        if self.name.size != n or self.kind.size != n:
            raise ValueError("name/kind must have one entry per source")
        if self.temp is not None:
            t = np.asarray(self.temp, dtype=float)
            if t.ndim == 1:
                t = t[:, None] if n == 1 else t[None, :]
            if t.shape != q.shape:
                raise ValueError(f"temp {t.shape} must have the shape of q {q.shape}")
            self.temp = t

    @property
    def n(self):
        return int(self.lon.size)

    @property
    def nt(self):
        return int(self.q.shape[0])

    @property
    def has_temp(self):
        return self.temp is not None

    def mean_q(self):
        return self.q.mean(axis=0) if self.nt else np.zeros(self.n)

    def _temp_or_nan(self):
        return self.temp if self.temp is not None else np.full(self.q.shape, np.nan)

    def subset(self, sel):
        """A new SourceSet with only the selected sources (bool mask or indices)."""
        sel = np.asarray(sel)
        if sel.dtype == bool:
            sel = np.flatnonzero(sel)
        return SourceSet(self.lon[sel], self.lat[sel], self.q[:, sel], self.time,
                         self.name[sel], self.kind[sel], dict(self.attrs),
                         temp=None if self.temp is None else self.temp[:, sel])

    def head(self, nrec):
        """The first nrec records (testing)."""
        return SourceSet(self.lon, self.lat, self.q[:nrec], self.time[:nrec], self.name,
                         self.kind, dict(self.attrs),
                         temp=None if self.temp is None else self.temp[:nrec])

    def with_temp(self, temp):
        """Same sources with a temperature array (nt, n), deg C, NaN = unknown."""
        return SourceSet(self.lon, self.lat, self.q, self.time, self.name, self.kind,
                         dict(self.attrs), temp=temp)

    def concat(self, other):
        """Append the sources of another SourceSet with the SAME time axis. If only one of
        the two has temperatures, the other's are unknown (NaN)."""
        if other.nt != self.nt or not np.array_equal(other.time, self.time):
            raise ValueError("concat needs identical time axes")
        temp = None
        if self.has_temp or other.has_temp:
            temp = np.concatenate([self._temp_or_nan(), other._temp_or_nan()], axis=1)
        return SourceSet(np.r_[self.lon, other.lon], np.r_[self.lat, other.lat],
                         np.concatenate([self.q, other.q], axis=1), self.time,
                         np.r_[self.name, other.name], np.r_[self.kind, other.kind],
                         {**self.attrs, **other.attrs}, temp=temp)

    def table(self):
        """List of dicts (name, kind, lon, lat, mean_q[, mean_temp]) -- for CSV/printing.
        mean_temp is discharge-weighted over the records with known temperature."""
        mq = self.mean_q()
        mt = None
        if self.has_temp:
            ok = np.isfinite(self.temp)
            w = np.where(ok, self.q, 0.0)
            num = (w * np.where(ok, self.temp, 0.0)).sum(axis=0)
            den = w.sum(axis=0)
            mt = np.where(den > 0, num / np.where(den > 0, den, 1.0), np.nan)
        rows = []
        for k in range(self.n):
            r = dict(name=self.name[k], kind=self.kind[k], lon=float(self.lon[k]),
                     lat=float(self.lat[k]), mean_q=float(mq[k]))
            if mt is not None:
                r["mean_temp"] = float(mt[k])
            rows.append(r)
        return rows

    def to_csv(self, path):
        import csv
        rows = self.table()
        fields = ["name", "kind", "lon", "lat", "mean_q"] + (["mean_temp"] if self.has_temp else [])
        with open(path, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=fields)
            w.writeheader()
            for r in rows:
                w.writerow(r)


def daily_axis(year, ndays=None):
    """datetime64[D] axis for a calendar year (or its first ndays days)."""
    t0 = np.datetime64(f"{int(year)}-01-01")
    n = int((np.datetime64(f"{int(year) + 1}-01-01") - t0).astype(int)) if ndays is None else ndays
    return t0 + np.arange(n).astype("timedelta64[D]")
