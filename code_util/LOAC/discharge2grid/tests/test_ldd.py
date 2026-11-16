"""LDD-based outlets on a synthetic drainage map, north-up and south-up.
pytest, or `python tests/test_ldd.py`.

Land west of column 40, sea east of it (missing values). Filler land drains west (code 4),
away from the sea, so it never makes an outlet. Rivers, all flowing east (code 6):
  R1 ends in a PIT (5) on the last land column           -> 'mouth'
  R2's last cell points INTO the sea (6 -> sea cell)      -> 'mouth' (non-pit outlet)
  R3 ends in a pit 4 cells short of the coast (~10 km)    -> 'near_sea' (a GloFAS dead end)
  R4 ends in a pit 30 cells inland (~75 km)               -> 'inland' with max_sea_km = 50
"""
import os
import sys
import tempfile

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from discharge2grid.readers import glofas  # noqa: E402


def make_ldd(path, north_up=True):
    from netCDF4 import Dataset
    ny, nx, coast = 30, 60, 40
    lat = 61.0 - 0.05 * np.arange(ny)            # north-to-south
    lon = 10.025 + 0.05 * np.arange(nx)
    L = np.full((ny, nx), 4, dtype=np.int16)      # filler drains west
    L[:, coast:] = -1                             # sea (missing)
    rows = dict(R1=5, R2=10, R3=15, R4=20)
    L[rows["R1"], 0:coast] = 6
    L[rows["R1"], coast - 1] = 5                  # pit at the coast
    L[rows["R2"], 0:coast] = 6                    # last cell drains into the sea
    L[rows["R3"], 0:coast - 4] = 6
    L[rows["R3"], coast - 5] = 5                  # pit 4 cells short
    L[rows["R4"], 0:coast - 30] = 6
    L[rows["R4"], coast - 31] = 5                 # pit 30 cells short
    if not north_up:
        lat, L = lat[::-1], L[::-1]
    with Dataset(path, "w") as ds:
        ds.createDimension("lat", ny)
        ds.createDimension("lon", nx)
        ds.createVariable("lat", "f8", ("lat",))[:] = lat
        ds.createVariable("lon", "f8", ("lon",))[:] = lon
        v = ds.createVariable("Band1", "i2", ("lat", "lon"), fill_value=-1)
        v[:] = np.ma.masked_less(L, 0)
    return lat, lon, rows, coast


def check(north_up):
    with tempfile.TemporaryDirectory() as td:
        p = os.path.join(td, "ldd.nc")
        lat, lon, rows, coast = make_ldd(p, north_up)
        out = glofas.outlets_from_ldd(p, max_sea_km=50)
        got = {(round(float(y), 3), round(float(x), 3)): k
               for y, x, k in zip(out["lat"], out["lon"], out["kind"])}
        lat_ns = 61.0 - 0.05 * np.arange(30)
        want = {
            (round(lat_ns[rows["R1"]], 3), round(lon[coast - 1], 3)): "mouth",
            (round(lat_ns[rows["R2"]], 3), round(lon[coast - 1], 3)): "mouth",
            (round(lat_ns[rows["R3"]], 3), round(lon[coast - 5], 3)): "near_sea",
            (round(lat_ns[rows["R4"]], 3), round(lon[coast - 31], 3)): "inland",
        }
        assert got == want, (north_up, got, want)
        u = dict(zip(out["kind"], out["use"]))
        assert u["mouth"] == 1 and u["near_sea"] == 1 and u["inland"] == 0


def test_ldd_north_up():
    check(True)


def test_ldd_south_up():
    check(False)


if __name__ == "__main__":
    test_ldd_north_up()
    print("ok test_ldd_north_up")
    test_ldd_south_up()
    print("ok test_ldd_south_up")
