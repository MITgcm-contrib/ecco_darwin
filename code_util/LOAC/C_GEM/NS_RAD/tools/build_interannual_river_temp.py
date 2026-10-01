"""
Derive a genuine multi-year (2005-2023) riverine water-temperature forcing from the
interannual air-temperature series, using the SAME regional relation build_river_temp.py
already validated for the single-2022-year forcing:

    T_river = max(0, 1.432 * (T_air_10day + 4.0))          [degrees C]

No new data source -- T_river has never been an independent observation in this project
(see build_river_temp.py's docstring: constrained fit against sparse USGS records,
applied regionally to all four rivers since only Kuparuk has a defensible site-specific
fit). This just applies that already-validated formula to
airtemp_interannual_2005-2023_degC.csv (tools/build_interannual_met.py) instead of the
single-year prda2h2022.txt.

ONE DIFFERENCE FROM THE SINGLE-YEAR VERSION: the 10-day smoothing there wraps
CIRCULARLY at the year boundary (defensible for a 365-day climatology that repeats).
Here the input is one continuous 19-year daily series, so the running mean is a
straight sequential mean with edge-padding only at the very start/end of the whole
record, not a per-year wrap -- there is no real physical discontinuity at each Dec
31 -> Jan 1 boundary in a genuine multi-year series, so wrapping there would be wrong.

Usage:  python tools/build_interannual_river_temp.py
"""
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent.parent
FORC = ROOT / "forcing"

SLOPE = 1.432
T0 = -4.0
SMOOTH_DAYS = 10

YEAR_START, YEAR_END = 2005, 2023


def smooth_edge_padded(x, window):
    """Running mean, edge-padded at the series boundaries (not circular -- see
    module docstring for why that differs from the single-year version)."""
    pad = np.r_[np.full(window, x[0]), x, np.full(window, x[-1])]
    return np.convolve(pad, np.ones(window) / window, mode="same")[window:-window]


def main():
    src = FORC / f"airtemp_interannual_{YEAR_START}-{YEAR_END}_degC.csv"
    air = np.genfromtxt(src, delimiter=",")
    n_expected = (YEAR_END - YEAR_START + 1) * 365
    assert air.size == n_expected, f"{src.name}: expected {n_expected}, got {air.size}"

    air_smooth = smooth_edge_padded(air, SMOOTH_DAYS)
    river_temp = np.maximum(0.0, SLOPE * (air_smooth + (-T0)))

    out = FORC / f"river_temp_interannual_{YEAR_START}-{YEAR_END}_degC.csv"
    np.savetxt(out, river_temp, fmt="%.4f")
    print(f"wrote {out.name} ({len(river_temp)} rows), "
          f"mean={river_temp.mean():.2f} degC, "
          f"days at floor (0 degC): {(river_temp == 0.0).sum()}")


if __name__ == "__main__":
    main()
