"""Table / gridded readers and the dropout scan. pytest, or `python tests/test_readers_checks.py`."""
import os
import sys
import tempfile

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from discharge2grid.readers import points, gridded  # noqa: E402
from discharge2grid.checks import scan_dropouts  # noqa: E402
from discharge2grid.geometry import regular_cell_area  # noqa: E402


def test_table_reader():
    with tempfile.TemporaryDirectory() as td:
        o, s = os.path.join(td, "o.csv"), os.path.join(td, "s.csv")
        open(o, "w").write("name,lon,lat\nAlpha,10.5,60.2\nBeta,-150.25,61.0\n")
        open(s, "w").write("date,Beta,Alpha\n2015-01,10,1\n2015-02,20,\n2015-03,30,3\n")
        src = points.from_csv(o, s)
        assert src.n == 2 and src.nt == 3
        assert list(src.name) == ["Alpha", "Beta"]
        assert np.allclose(src.q[:, 0], [1, 0, 3]) and np.allclose(src.q[:, 1], [10, 20, 30])
        assert src.time[1] == np.datetime64("2015-02-01")


def test_gridded_units():
    lon = np.arange(0.125, 2, 0.25)
    lat = np.arange(60.125, 61, 0.25)
    f = np.zeros((2, lat.size, lon.size))
    f[:, 1, 2] = [1.0, 2.0]                              # kg m-2 s-1
    s = gridded.from_field(lon, lat, f, np.array(["2015-01-01", "2015-01-02"], "datetime64[D]"),
                           "kg m-2 s-1")
    a = regular_cell_area(lon, lat)[1, 2]
    assert s.n == 1 and np.allclose(s.q[:, 0], np.array([1.0, 2.0]) * a / 1000.0)
    s2 = gridded.from_field(lon, lat, f, np.arange(2), "m3 s-1")
    assert np.allclose(s2.q[:, 0], [1.0, 2.0])


def test_dropout_scan():
    years = list(range(2000, 2010))
    tables = {}
    for y in years:
        lon, lat, q = [10.0, 20.0, 20.001, 30.0], [60.0, 61.0, 61.0, 62.0], [100.0, 50.0, 0.0, 80.0]
        if y in (2003, 2007):
            q[0] = 0.0                     # river 1 lost in two years
        if y == 2005:                      # river 2 moves 0.1 km (same river, other cell)
            q[1], q[2] = 0.0, 50.0
        if y == 2004:                      # river 3 moves 20 km: a gap that is NOT lost water
            lon[3], lat[3] = 30.35, 62.0
        tables[y] = (np.array(lon), np.array(lat), np.array(q))
    res = scan_dropouts(tables, merge_km=8, min_q=20)
    by = {round(r["lon"]): r for r in res}
    assert by[10]["gap_years"] == [2003, 2007] and by[10]["keep"]
    assert 20 not in by                        # merged within 8 km: no gap
    assert by[30]["gap_years"] == [2004] and not by[30]["keep"]


if __name__ == "__main__":
    for name, f in sorted(globals().items()):
        if name.startswith("test_") and callable(f):
            f()
            print("ok", name)
