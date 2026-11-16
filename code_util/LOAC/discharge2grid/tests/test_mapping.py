"""Mapping core on synthetic grids. Run with pytest, or `python tests/test_mapping.py`."""
import math
import os
import sys
import tempfile

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from discharge2grid import (ModelGrid, Mapper, SourceSet, daily_axis, write_mitgcm_runoff,  # noqa
                            write_sparse_netcdf, read_sparse_netcdf, KEPT, TOO_FAR, OUTSIDE)
from discharge2grid.geometry import regular_cell_area, gc_km  # noqa


def box_grid(lon0=0.0, lon1=6.0, lat0=60.0, lat1=64.0, d=0.05, coast_lon=2.0):
    """Regular grid; land west of coast_lon, sea east of it."""
    lon = np.arange(lon0 + d / 2, lon1, d)
    lat = np.arange(lat0 + d / 2, lat1, d)
    LON, LAT = np.meshgrid(lon, lat)
    wet = LON > coast_lon
    return ModelGrid(LON, LAT, regular_cell_area(lon, lat), wet, dims=("j", "i"), name="box")


def sources(lon, lat, nt=5, seed=0):
    rng = np.random.default_rng(seed)
    q = rng.uniform(10, 1000, size=(nt, len(lon)))
    q[1, 0] = 0.0                      # a zero day
    return SourceSet(lon, lat, q, daily_axis(2015, nt), name=[f"r{k}" for k in range(len(lon))])


def test_conservation_and_snap():
    g = box_grid()
    s = sources([1.98, 1.90, 1.0, 5.5], [61.0, 62.0, 62.0, 61.5])
    m = Mapper(g, sigma_km=10.0, max_snap_km=20.0)
    d = m.map_sources(s)
    st = d.snap.status
    assert st[0] == KEPT and st[1] == KEPT          # 1-6 km from water
    assert st[2] == TOO_FAR                          # ~55 km inland
    assert st[3] == KEPT                             # already at sea
    for k in range(s.nt):
        delivered = float((d.runoff_points(k) * d.point_area).sum())
        kept = float(s.q[k, st == KEPT].sum())
        assert abs(delivered / kept - 1.0) < 1e-12
    # dense and point outputs agree
    dense = d.runoff_dense(3)
    assert np.allclose(dense.ravel()[d.points], d.runoff_points(3))
    assert dense[~g.wet].max() == 0.0


def test_unclipped_footprint_area():
    g = box_grid(lon0=0, lon1=8, coast_lon=0.5)
    m = Mapper(g, sigma_km=12.0)
    kw = int(m.snap([4.0], [62.0]).dest[0])         # far from coast and edges
    a = m.footprint_area_m2(kw)
    ideal = 4 * math.pi * (12e3) ** 2
    assert abs(a / ideal - 1) < 0.03, a / ideal


def test_exclude_outside_domain():
    g = box_grid(lon0=0, lon1=6, lat0=60, lat1=64)
    lon, lat = [5.0], [64.3]                         # 33 km north of the grid edge, over "sea"
    on = Mapper(g, sigma_km=8, max_snap_km=60, exclude_outside=True).snap(lon, lat)
    off = Mapper(g, sigma_km=8, max_snap_km=60, exclude_outside=False).snap(lon, lat)
    assert on.status[0] == OUTSIDE and off.status[0] == KEPT
    # inland source inside the domain is kept with exclusion on (land cells count)
    inside = Mapper(g, sigma_km=8, max_snap_km=60, exclude_outside=True).snap([1.6], [62.0])
    assert inside.status[0] == KEPT


def test_nd_grid_and_sparse_roundtrip():
    # an LLC-like (tile, j, i) layout: three 2-D patches stacked
    tiles = [box_grid(lon0=x, lon1=x + 3, coast_lon=x + 1) for x in (0.0, 10.0, 20.0)]
    LON = np.stack([t.lon for t in tiles]); LAT = np.stack([t.lat for t in tiles])
    A = np.stack([t.area for t in tiles]); W = np.stack([t.wet for t in tiles])
    g = ModelGrid(LON, LAT, A, W, dims=("tile", "j", "i"), name="tiles")
    s = sources([0.98, 10.98, 20.98], [61.0, 62.0, 63.0], nt=4)
    d = Mapper(g, sigma_km=6, max_snap_km=10).map_sources(s)
    assert (d.snap.status == KEPT).all()
    tile_of_point = g.unravel(d.points)[0]
    assert set(tile_of_point.tolist()) == {0, 1, 2}
    with tempfile.TemporaryDirectory() as td:
        p = write_sparse_netcdf(os.path.join(td, "x.nc"), d)
        r = read_sparse_netcdf(p)
        assert r.grid_shape == g.shape
        for k in range(s.nt):
            dense = d.runoff_dense(k)
            assert np.allclose(r.to_dense(k), dense, rtol=1e-6, atol=0)
            assert abs(r.total_m3s(k) / s.q[k].sum() - 1) < 1e-6
        b = os.path.join(td, "runoff_x_2015")
        assert write_mitgcm_runoff(b, d) == s.nt
        arr = np.fromfile(b, dtype=">f4").reshape((s.nt,) + g.shape)
        assert np.array_equal(arr[2], d.runoff_dense(2).astype(">f4"))
        # the operator stored in the file reproduces the runoff
        from netCDF4 import Dataset
        with Dataset(p) as ds:
            lp, ls, lw = ds["link_point"][:], ds["link_source"][:], ds["link_weight"][:]
            q = ds["discharge"][:].astype(float)
        rr = np.zeros(d.n_point)
        np.add.at(rr, lp, lw * q[1, ls])
        assert np.allclose(rr, d.runoff_points(1), rtol=1e-6)


def test_regular_area_sums_to_sphere():
    lon = np.arange(-179.5, 180, 1.0)
    lat = np.arange(-89.5, 90, 1.0)
    tot = regular_cell_area(lon, lat).sum()
    assert abs(tot / (4 * math.pi * 6371e3 ** 2) - 1) < 1e-12


if __name__ == "__main__":
    for name, f in sorted(globals().items()):
        if name.startswith("test_") and callable(f):
            f()
            print("ok", name)
