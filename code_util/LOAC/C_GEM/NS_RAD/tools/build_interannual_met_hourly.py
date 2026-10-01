"""
Build genuine multi-year (2005-2023) HOURLY air-temperature AND wind-speed forcings
from NDBC PRDA2 -- resolves the diurnal cycle, unlike this project's original daily-mean
versions (build_interannual_met.py), which have no day/night structure at all. See
CLAUDE.md -> "Interannual forcing" -> "Sub-daily (diurnal) forcing" for why this exists
and file_module.py's generalized `row_interval_sec` support that makes the model able to
read it.

Wind is included here (not just air temperature) because there is no extra cost to it:
PRDA2's raw record already has both in the same file, at the same native hourly/sub-
hourly cadence -- the ORIGINAL daily-mean builder (build_interannual_met.py) averaged
both away, but nothing stops extracting wind hourly from the exact same cached files
air temperature already uses. (Humidity and storm surge are ALSO extended to hourly
this same way -- see build_interannual_humidity_hourly.py and
build_interannual_surge_hourly.py -- for the same reason: their raw sources are already
hourly too, so there was no real barrier, only an earlier, narrower scoping choice.)

Reuses forcing/_ndbc_cache/prda2h{year}.txt, already downloaded by
build_interannual_met.py -- no new download needed for either variable.

CALENDAR. One row per clock hour (any sub-hourly readings within the same hour are
averaged), 8760 hours/year (Feb 29's 24 hours dropped in leap years, same "every year
is exactly 365 days" convention as every other interannual forcing here, just at hourly
granularity: 19*8760 = 166,440 rows total).

GAPS. Short gaps (<= MAX_INTERP_GAP_HOURS) linearly interpolated; longer gaps filled
with the multi-year HOUR-OF-YEAR climatology (the mean of that same calendar hour
across all other years that have real data then) -- the hourly analogue of every other
builder's day-of-year climatology fill.

Usage:  python tools/build_interannual_met_hourly.py
"""
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent.parent
FORC = ROOT / "forcing"
CACHE = FORC / "_ndbc_cache"

YEAR_START = 2005
YEAR_END = 2023
N_YEARS = YEAR_END - YEAR_START + 1
HOURS_PER_YEAR = 365 * 24  # 8760, Feb 29 dropped
MAX_INTERP_GAP_HOURS = 6

COL_MONTH, COL_DAY, COL_HOUR = 1, 2, 3
COL_WSPD = 6
COL_ATMP = 13

_CUM_DAYS = [0, 31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334]


def _hour_of_year(month, day, hour):
    return (_CUM_DAYS[month - 1] + (day - 1)) * 24 + hour


def _parse_year_hourly(year):
    """Return (atmp_hourly, wspd_hourly) for one year."""
    path = CACHE / f"prda2h{year}.txt"
    is_leap = (year % 4 == 0 and year % 100 != 0) or (year % 400 == 0)
    nhours = 366 * 24 if is_leap else 365 * 24
    sums = {"atmp": np.zeros(nhours), "wspd": np.zeros(nhours)}
    counts = {"atmp": np.zeros(nhours), "wspd": np.zeros(nhours)}
    with open(path) as f:
        for line in f:
            if not line or not line[0].isdigit():
                continue
            parts = line.split()
            if len(parts) <= COL_ATMP:
                continue
            month, day, hour = int(parts[COL_MONTH]), int(parts[COL_DAY]), int(parts[COL_HOUR])
            idx = _hour_of_year(month, day, hour)
            atmp = float(parts[COL_ATMP])
            if atmp < 999.0:
                sums["atmp"][idx] += atmp
                counts["atmp"][idx] += 1
            wspd = float(parts[COL_WSPD])
            if wspd < 99.0:
                sums["wspd"][idx] += wspd
                counts["wspd"][idx] += 1
    with np.errstate(invalid="ignore", divide="ignore"):
        atmp_h = np.where(counts["atmp"] > 0, sums["atmp"] / counts["atmp"], np.nan)
        wspd_h = np.where(counts["wspd"] > 0, sums["wspd"] / counts["wspd"], np.nan)
    return atmp_h, wspd_h


def _drop_feb29_hours(arr, year):
    is_leap = (year % 4 == 0 and year % 100 != 0) or (year % 400 == 0)
    if is_leap and len(arr) == 366 * 24:
        feb29_start = 59 * 24
        return np.delete(arr, np.arange(feb29_start, feb29_start + 24))
    return arr[:HOURS_PER_YEAR]


def _interp_short_gaps(flat, max_gap):
    nan_mask = np.isnan(flat)
    if not nan_mask.any():
        return flat
    idx = np.arange(len(flat))
    d = np.diff(nan_mask.astype(int))
    starts = np.where(d == 1)[0] + 1
    ends = np.where(d == -1)[0] + 1
    if nan_mask[0]:
        starts = np.r_[0, starts]
    if nan_mask[-1]:
        ends = np.r_[ends, len(flat)]
    for s, e in zip(starts, ends):
        if (e - s) <= max_gap and s > 0 and e < len(flat):
            flat[s:e] = np.interp(idx[s:e], [s - 1, e], [flat[s - 1], flat[e]])
    return flat


def _climatology_fill(flat, n_years):
    grid = flat.reshape(n_years, HOURS_PER_YEAR)
    for h in range(HOURS_PER_YEAR):
        col = grid[:, h]
        nan_here = np.isnan(col)
        if nan_here.any() and not nan_here.all():
            col[nan_here] = np.nanmean(col[~nan_here])
    return grid.reshape(-1)


def _finish_and_write(per_year_list, label, out_name, unit_fmt="%.4f"):
    flat = np.concatenate(per_year_list)
    assert flat.size == N_YEARS * HOURS_PER_YEAR
    n_missing = int(np.isnan(flat).sum())
    flat = _interp_short_gaps(flat, MAX_INTERP_GAP_HOURS)
    flat = _climatology_fill(flat, N_YEARS)
    n_remaining = int(np.isnan(flat).sum())
    print(f"  {label}: {n_missing} missing hours total; {n_remaining} remain unfilled "
          f"(should be 0)")
    out = FORC / out_name
    np.savetxt(out, flat, fmt=unit_fmt)
    print(f"  wrote {out.name} ({len(flat)} rows), mean={flat.mean():.2f}")
    return flat


def main():
    print(f"Building hourly air temperature + wind speed {YEAR_START}-{YEAR_END} "
          f"from cached PRDA2 files...")
    atmp_years, wspd_years = [], []
    for year in range(YEAR_START, YEAR_END + 1):
        atmp, wspd = _parse_year_hourly(year)
        n_valid_atmp = int(np.sum(~np.isnan(atmp[:HOURS_PER_YEAR])))
        n_valid_wspd = int(np.sum(~np.isnan(wspd[:HOURS_PER_YEAR])))
        print(f"  {year}: {n_valid_atmp}/{HOURS_PER_YEAR} hours ATMP, "
              f"{n_valid_wspd}/{HOURS_PER_YEAR} hours WSPD")
        atmp_years.append(_drop_feb29_hours(atmp, year))
        wspd_years.append(_drop_feb29_hours(wspd, year))

    _finish_and_write(atmp_years, "airtemp",
                       f"airtemp_hourly_interannual_{YEAR_START}-{YEAR_END}_degC.csv")
    _finish_and_write(wspd_years, "windspeed",
                       f"windspeed_hourly_interannual_{YEAR_START}-{YEAR_END}_msec.csv")


if __name__ == "__main__":
    main()
