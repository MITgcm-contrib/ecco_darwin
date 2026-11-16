"""River temperature: mixing, heat conservation, outputs, filling, climatology building.
pytest, or `python tests/test_temperature.py`."""
import os
import sys
import tempfile

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from discharge2grid import (ModelGrid, Mapper, SourceSet, daily_axis, write_mitgcm_temperature,  # noqa
                            write_sparse_netcdf, read_sparse_netcdf)
from discharge2grid.geometry import regular_cell_area  # noqa
from discharge2grid.readers import temperature as rt, points  # noqa


def toy():
    lon = np.arange(0.025, 6.0, 0.05)
    lat = np.arange(60.025, 64.0, 0.05)
    LON, LAT = np.meshgrid(lon, lat)
    return ModelGrid(LON, LAT, regular_cell_area(lon, lat), LON > 2.0, dims=("j", "i"))


def sources(temp):
    nt = 6
    q = np.stack([np.linspace(500, 900, nt), np.linspace(100, 50, nt), np.full(nt, 40.0)], axis=1)
    return SourceSet([1.98, 1.98, 1.98], [61.0, 61.15, 63.0], q, daily_axis(2015, nt),
                     name=["A", "B", "C"], temp=temp)


def test_mixing_and_heat():
    g = toy()
    T = np.array([[5.0, 15.0, 8.0]] * 6)
    s = sources(T)
    d = Mapper(g, sigma_km=8, max_snap_km=20).map_sources(s)
    for k in range(s.nt):
        t, r = d.temperature_points(k, return_runoff=True)
        wet = r > 0
        assert np.all((t[wet] >= 5 - 1e-12) & (t[wet] <= 15 + 1e-12))
        # heat conservation: sum R T A == sum q T
        heat = float((r[wet] * t[wet] * d.point_area[wet]).sum())
        assert abs(heat / float((s.q[k] * s.temp[k]).sum()) - 1) < 1e-12
        # a cell reached only by C has C's temperature exactly
        onlyC = d.Wp[:, 2].toarray().ravel() > 0
        onlyC &= d.Wp[:, 0].toarray().ravel() == 0
        onlyC &= d.Wp[:, 1].toarray().ravel() == 0
        assert onlyC.any() and np.allclose(t[onlyC], 8.0, atol=1e-12)
    info = d.report()
    assert abs(info["heat_ratio"] - 1) < 1e-12 and info["cells_unknown_temp"] == 0
    with tempfile.TemporaryDirectory() as td:
        p = write_sparse_netcdf(os.path.join(td, "x.nc"), d)
        r = read_sparse_netcdf(p)
        b = os.path.join(td, "runoff_temp_x_2015")
        write_mitgcm_temperature(b, d)
        arr = np.fromfile(b, ">f4").reshape((s.nt,) + g.shape)
        for k in (0, 5):
            assert np.allclose(r.temperature[k], d.temperature_points(k), atol=1e-5)
            assert np.array_equal(arr[k], d.temperature_dense(k).astype(">f4"))
            assert np.allclose(r.temperature_dense(k), d.temperature_dense(k), atol=1e-5)
        assert arr[0][~g.wet].max() == 0.0


def test_unknown_temperature_and_fill():
    g = toy()
    T = np.array([[5.0, np.nan, 8.0]] * 6)        # B unknown
    s = sources(T)
    d = Mapper(g, sigma_km=8, max_snap_km=20).map_sources(s)
    t, r = d.temperature_points(0, return_runoff=True)
    onlyB = (d.Wp[:, 1].toarray().ravel() > 0) & (d.Wp[:, 0].toarray().ravel() == 0)
    if onlyB.any():                              # cells only B reaches: no known temperature
        assert np.all(np.isnan(t[onlyB]))
    assert d.report()["cells_unknown_temp"] >= 0
    f, donors = rt.fill_missing(s, max_km=50)
    assert donors == {1: 0}                      # B borrows A (largest within 50 km), not C (220 km)
    assert np.allclose(f.temp[:, 1], 5.0)
    f2, donors2 = rt.fill_missing(s, max_km=5)
    assert donors2 == {1: -1} and np.isnan(f2.temp[:, 1]).all()


def test_climatology_from_binary():
    # 8 x 12 grid, 1-degree cells, north-first rows, little-endian Kelvin, 1e20 over "sea"
    nlat, nlon = 8, 12
    lat_first, lon0, inc = 64.5, 0.5, 1.0
    lat_rows = lat_first - inc * np.arange(nlat)                 # 64.5 .. 57.5 (north first)
    sea = np.zeros((nlat, nlon), bool)
    sea[:, 8:] = True                                            # east is sea
    with tempfile.TemporaryDirectory() as td:
        files = []
        for year, n in ((2015, 365), (2016, 366)):
            base = 0.0 if year == 2015 else 1.0                  # years differ by 1 C
            t0 = np.datetime64(f"{year}-01-01")
            slots = rt.calendar_slot(t0 + np.arange(n).astype("timedelta64[D]"))
            rec = np.empty((n, nlat, nlon), dtype="<f4")
            for k in range(n):
                # value = 273.15 + 2 + base + slot/100 + row-latitude/1000  (all within 0..45 C)
                rec[k] = 273.15 + 2 + base + slots[k] / 100.0 + lat_rows[:, None] / 1000.0
                rec[k][sea] = 1e20
            p = os.path.join(td, f"tw{year}.bin")
            rec.tofile(p)
            files.append(p)
        c = rt.TemperatureClimatology.from_binary(files, nlon=nlon, nlat=nlat, lon0=lon0,
                                                  lat_first=lat_first, inc=inc, dtype="<f4",
                                                  north_first=True, kelvin=True, min_good=30,
                                                  verbose=False)
        assert np.all(np.diff(c.lat) > 0)                        # stored south to north
        assert c.valid[:, :8].all() and not c.valid[:, 8:].any()
        j = int(np.argmin(np.abs(c.lat - 60.5)))
        # Mar 1 (slot 60) is the mean of both years; Feb 29 (slot 59) only from 2016
        assert abs(c.clim[60, j, 3] - (2 + 0.5 + 0.60 + 60.5 / 1000)) < 1e-4
        assert abs(c.clim[59, j, 3] - (2 + 1.0 + 0.59 + 60.5 / 1000)) < 1e-4
        # sampling: a source at sea takes the nearest usable cell; dates map to slots
        t, dkm = c.sample([8.4], [60.6], np.array(["2015-03-01", "2016-02-29"], "datetime64[D]"))
        assert dkm[0] < 60 and abs(t[0, 0] - c.clim[60, j, 7]) < 1e-6 and abs(t[1, 0] - c.clim[59, j, 7]) < 1e-6
        t2, _ = c.sample([30.0], [60.5], np.array(["2015-03-01"], "datetime64[D]"), max_km=100)
        assert np.isnan(t2).all()                                # beyond max_km
        npz = c.save(os.path.join(td, "c.npz"))
        c2 = rt.TemperatureClimatology.load(npz)
        assert np.array_equal(c2.clim, c.clim) and np.array_equal(c2.valid, c.valid)


def test_table_with_temperature():
    with tempfile.TemporaryDirectory() as td:
        o, q, t = (os.path.join(td, f) for f in ("o.csv", "q.csv", "t.csv"))
        open(o, "w").write("name,lon,lat\nA,1.98,61.0\nB,1.98,62.0\n")
        open(q, "w").write("date,A,B\n2015-01-01,10,20\n2015-01-02,11,21\n")
        open(t, "w").write("date,A,B\n2015-01-01,1.5,\n2015-01-02,2.5,3.0\n")
        s = points.from_csv(o, q, temp_csv=t)
        assert s.has_temp and np.isnan(s.temp[0, 1]) and s.temp[1, 1] == 3.0


if __name__ == "__main__":
    for name, f in sorted(globals().items()):
        if name.startswith("test_") and callable(f):
            f()
            print("ok", name)
