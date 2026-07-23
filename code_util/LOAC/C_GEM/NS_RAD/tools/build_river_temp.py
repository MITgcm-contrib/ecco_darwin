"""
Build river water-temperature (and air-temperature) forcings for NS-RAD.

WHY THIS EXISTS
---------------
The shipped `watertemp.csv` is NOT river temperature. It is NDBC buoy PRDA2
*sea-water* temperature at Prudhoe Bay -- verified r = 0.9922, mean |diff| 0.093 C
against prda2h2022.txt, and it reaches -2.00 C, the freezing point of seawater, which
fresh river water cannot do. Using it as a river forcing:

  * biased every temperature-dependent rate cold by 6-8 C, understating respiration
    1.84x (Fhet, Q10=2.75) and nitrification 2.63x (Fnit, Q10=5);
  * drove the ice gate off COASTAL SEA ICE, which clears weeks after river breakup,
    so the gate opened on day 173 while Colville and Kuparuk peak on day 153 --
    discarding 45% and 40% of annual water volume including the entire freshet.

This script replaces it with a river temperature modelled from air temperature.

THE RELATION
------------
    T_river = max(0, 1.432 * (T_air_10day + 4.0))          [degrees C]

where T_air_10day is a 10-day running mean of PRDA2 air temperature (ATMP), the
smoothing standing in for the river's thermal inertia.

Derivation and the choices behind it:

  * Fitted against USGS parameter 00010 (water temperature) daily values. Only two of
    the four gauges have any 2022 record, and both are summer-only: Sagavanirktok
    (15908000) 93 days, Kuparuk (15896000) 55 days. Colville and Canning have none.

  * The form is CONSTRAINED -- slope only, through a fixed zero-crossing T0 = -4 C --
    rather than a free linear fit. An unconstrained fit scores marginally better on
    the summer data but extrapolates catastrophically into winter: the Sagavanirktok
    fit crossed zero at T_air = -20.6 C, which would hold that river above freezing
    for 297 days a year. The constraint costs almost nothing in-sample (Kuparuk
    R2 0.707 -> 0.706) and makes winter behaviour physical.

  * The relation is REGIONAL, taken from Kuparuk and applied to all four rivers.
    Kuparuk is the only site whose temperature is well explained by air temperature
    (R2 = 0.706, RMSE 2.45 C, and improving monotonically with smoothing: 0.54 at
    1-day -> 0.71 at 10-day). The Sagavanirktok is NOT: constrained fits to its own
    data return NEGATIVE R2 at every T0 tested, i.e. worse than predicting its mean.
    It is mountain-fed with strong snowmelt and groundwater influence (and is the
    aufeis river), so it is decoupled from local air temperature. There is therefore
    no defensible site-specific fit for it, and no data at all for the other two.

VALIDATION (against observed river temperature)

    river            regional relation      sea-water forcing it replaces
    Kuparuk          RMSE 2.45  bias +0.05  RMSE 9.39  bias -8.30
    Sagavanirktok    RMSE 4.05  bias +2.61  RMSE 6.89  bias -6.00

So this is an improvement for both rivers, but it is NOT a good thermal model -- it
carries a +2.6 C warm bias on the Sagavanirktok. It trades a large cold bias for a
smaller warm one. A proper treatment needs an equilibrium-temperature or heat-budget
formulation with discharge dependence. Treat the absolute temperatures with caution;
treat the seasonal timing as much improved.

Usage:  python tools/build_river_temp.py
"""

import datetime as dt
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent.parent
FORCINGS = ROOT / "forcing"
NDBC = FORCINGS / "prda2h2022.txt"
YEAR = 2022

# Regional air->water relation. See module docstring for derivation.
SLOPE = 1.432
T0 = -4.0
SMOOTH_DAYS = 10

ATMP_COL = 13     # air temperature, degrees C
NDBC_MISSING = 99.0


def daily_air_temp():
    """365 daily-mean air temperatures from the NDBC record, gaps interpolated."""
    acc = {}
    with open(NDBC) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            col = line.split()
            if len(col) <= ATMP_COL:
                continue
            val = float(col[ATMP_COL])
            if val >= NDBC_MISSING:
                continue
            doy = dt.date(int(col[0]), int(col[1]), int(col[2])).timetuple().tm_yday - 1
            if doy < 365:
                acc.setdefault(doy, []).append(val)

    series = np.array([np.mean(acc[d]) if d in acc else np.nan for d in range(365)])
    ok = ~np.isnan(series)
    if not ok.any():
        raise SystemExit("no usable air temperature in NDBC file")
    series[~ok] = np.interp(np.flatnonzero(~ok), np.flatnonzero(ok), series[ok])
    return series


def smooth_circular(x, window):
    """Running mean that wraps at the year boundary -- the forcing is a climatology."""
    pad = np.r_[x[-window:], x, x[:window]]
    return np.convolve(pad, np.ones(window) / window, mode="same")[window:-window]


def write_csv(path, values):
    """Match the existing forcing format exactly: CRLF, 2 dp, no header, no final newline."""
    if len(values) != 365:
        raise SystemExit(f"{path.name}: expected 365 values, got {len(values)}")
    path.write_text("\r\n".join(f"{v:.2f}" for v in values), encoding="utf-8", newline="")


def main():
    air = daily_air_temp()
    air_s = smooth_circular(air, SMOOTH_DAYS)
    water = np.clip(SLOPE * (air_s - T0), 0.0, None)

    write_csv(FORCINGS / f"airtemp_{YEAR}_degC.csv", air)
    write_csv(FORCINGS / f"river_watertemp_{YEAR}_degC.csv", water)

    open_days = np.flatnonzero(water > 0)
    # Reproduce exfread's rolling ice-sum to report the resulting gate window.
    roll = np.convolve(water, np.ones(SMOOTH_DAYS), mode="same")
    gate = np.flatnonzero(roll > 0)

    print(f"air temperature   : min {air.min():6.1f}  max {air.max():5.1f}  mean {air.mean():6.2f} C")
    print(f"river temperature : min {water.min():6.1f}  max {water.max():5.1f}  "
          f"mean over open water {water[water > 0].mean():5.2f} C")
    print(f"  above 0 C on {len(open_days)} days (day {open_days.min()}-{open_days.max()})")
    print(f"  resulting ice gate open: day {gate.min()}-{gate.max()} ({len(gate)} days)")
    print(f"  previous sea-water gate: day 173-271 (99 days)")
    print()
    print(f"wrote {FORCINGS.name}/airtemp_{YEAR}_degC.csv")
    print(f"wrote {FORCINGS.name}/river_watertemp_{YEAR}_degC.csv")
    print()
    print("NOTE: watertemp.csv (sea water) is left in place but is no longer read by")
    print("      main.py. It is retained only as provenance for prda2h2022.txt.")


if __name__ == "__main__":
    main()
