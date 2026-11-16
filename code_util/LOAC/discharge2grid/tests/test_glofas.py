"""GloFAS outlet finding on a synthetic GloFAS-like file. pytest, or `python tests/test_glofas.py`.

Synthetic world (0.05 deg, lat 60-62, lon 10-13, 60 days, chunked 29 x 7 x 9 + deflate):
  sea   : lon > 12.0, plus a strip of sea on the northern edge
  river A ends at the coast (cell next to the sea): a normal outlet, 500 m3/s
  river B runs ALONG the coast for 4 cells before its mouth: every one of those cells is
          coastal and carries ~300 m3/s -> along-coast duplicates
  river C ends inland 6 cells (~15 km) from the sea: a dead end, 400 m3/s
  background: 0.5 m3/s on every land cell (local runoff)
"""
import os
import sys
import tempfile

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from discharge2grid.readers import glofas  # noqa: E402


def make_file(path, nt=60):
    from netCDF4 import Dataset
    lat = np.round(np.arange(61.975, 59.99, -0.05), 3)          # north-to-south, like GloFAS
    lon = np.round(np.arange(10.025, 13.0, 0.05), 3)
    ny, nx = lat.size, lon.size
    base = np.full((ny, nx), 0.5)
    sea = (lon[None, :] > 12.0) | (lat[:, None] > 61.9)
    q = np.repeat(base[None], nt, axis=0)
    t = np.arange(nt)
    season = 1.0 + 0.5 * np.sin(2 * np.pi * t / nt)
    jA = int(np.argmin(np.abs(lat - 60.5)))
    iC = int(np.argmin(np.abs(lon - 12.0)))                     # last land column
    for i in range(iC - 20, iC + 1):                             # river A: west -> east
        q[:, jA, i] = 500 * season * (0.6 + 0.4 * (i - (iC - 20)) / 20)
    jB = int(np.argmin(np.abs(lat - 61.2)))
    for i in range(iC - 15, iC - 3):                             # river B inland part
        q[:, jB, i] = 300 * season * 0.95
    for k, j in enumerate(range(jB, jB + 4)):                    # then south along the coast
        q[:, j, iC] = 300 * season * (0.96 + 0.01 * k)
    jCc = int(np.argmin(np.abs(lat - 60.9)))
    for i in range(iC - 30, iC - 6):                             # river C: stops 6 cells short
        q[:, jCc, i] = 400 * season * (0.7 + 0.3 * (i - (iC - 30)) / 23)
    q[:, sea] = np.nan
    with Dataset(path, "w") as ds:
        ds.createDimension("valid_time", nt)
        ds.createDimension("latitude", ny)
        ds.createDimension("longitude", nx)
        tv = ds.createVariable("valid_time", "i8", ("valid_time",))
        tv.units = "days since 2015-01-01"
        tv.calendar = "proleptic_gregorian"
        tv[:] = t
        ds.createVariable("latitude", "f8", ("latitude",))[:] = lat
        ds.createVariable("longitude", "f8", ("longitude",))[:] = lon
        v = ds.createVariable("avg_dis", "f4", ("valid_time", "latitude", "longitude"),
                              zlib=True, complevel=1, chunksizes=(29, 7, 9), fill_value=np.nan)
        v.units = "m**3 s**-1"
        v[:] = q
    return lat, lon, q, dict(A=(jA, iC), B=(jB + 3, iC), C=(jCc, iC - 7), sea=sea)


def test_scan_outlets_deadends_and_points():
    with tempfile.TemporaryDirectory() as td:
        p = os.path.join(td, "glofas_original_2015.nc")
        lat, lon, q, key = make_file(p)
        scan = glofas.scan_raw(p, verbose=False)
        assert np.array_equal(scan.ocean, key["sea"])
        land = ~key["sea"]
        assert np.allclose(scan.qmean[land], q[:, land].mean(axis=0))
        out = glofas.coastal_outlets(scan)
        cells = set(zip(out["row"].tolist(), out["col"].tolist()))
        assert key["A"] in cells and key["B"] in cells
        assert key["C"] not in cells                  # the dead end is not coastal
        # the biggest coastal cell is river A
        assert (int(out["row"][0]), int(out["col"][0])) == key["A"]
        dup = glofas.flag_coastal_duplicates(out, scan)
        idx = {(int(r), int(c)): k for k, (r, c) in enumerate(zip(out["row"], out["col"]))}
        jB, iC = key["B"]
        along = [idx[(j, iC)] for j in range(jB - 3, jB + 1)]
        assert dup[along[-1]] == -1                   # the mouth (largest) is kept
        assert all(dup[k] >= 0 for k in along[:-1])   # the cells upstream along the coast
        assert dup[idx[key["A"]]] == -1
        dead = glofas.dead_end_candidates(scan, min_q=100, max_sea_km=50)
        assert any((d["row"], d["col"]) == key["C"] for d in dead)
        assert not any((d["row"], d["col"]) in (key["A"], key["B"]) for d in dead)
        # chunk-aware point reader == direct slicing
        pts_lon = [lon[key["A"][1]], lon[key["C"][1]], lon[3]]
        pts_lat = [lat[key["A"][0]], lat[key["C"][0]], lat[30]]
        qq, time = glofas.read_outlet_series(p, pts_lon, pts_lat)
        assert qq.shape == (60, 3) and time[0] == np.datetime64("2015-01-01")
        ref = np.nan_to_num(np.stack([q[:, key["A"][0], key["A"][1]], q[:, key["C"][0], key["C"][1]],
                                      q[:, 30, 3]], axis=1))
        assert np.allclose(qq, ref.astype(np.float32))
        s = glofas.sources_from_raw(os.path.join(td, "glofas_original_{year}.nc"), 2015,
                                    pts_lon[:2], pts_lat[:2], deliver_lon=[12.1, 12.1],
                                    deliver_lat=[60.5, 60.9])
        assert s.n == 2 and np.allclose(s.lon, 12.1)


if __name__ == "__main__":
    test_scan_outlets_deadends_and_points()
    print("ok test_scan_outlets_deadends_and_points")
