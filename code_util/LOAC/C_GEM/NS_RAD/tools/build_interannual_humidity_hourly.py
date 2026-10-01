"""
Build a genuine multi-year (2005-2023) HOURLY relative-humidity forcing from Deadhorse
Airport (NOAA ISD) -- resolves the diurnal cycle, unlike build_interannual_humidity.py's
daily-mean version. See CLAUDE.md -> "Interannual forcing" -> "Sub-daily (diurnal)
forcing".

Reuses forcing/_isd_cache/deadhorse_isd_{year}.csv, already downloaded by
build_interannual_humidity.py -- no new download needed (the ISD record is already
hourly; the daily version just averaged it away).

CALENDAR / GAPS. Same hourly convention as build_interannual_met_hourly.py: 8760
hours/year (Feb 29 dropped), short gaps (<= MAX_INTERP_GAP_HOURS) linearly
interpolated, longer gaps filled with the multi-year hour-of-year climatology.

Usage:  python tools/build_interannual_humidity_hourly.py
"""
import csv
import datetime as dt
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent.parent
FORC = ROOT / "forcing"
CACHE = FORC / "_isd_cache"

YEAR_START = 2005
YEAR_END = 2023
N_YEARS = YEAR_END - YEAR_START + 1
HOURS_PER_YEAR = 365 * 24
MAX_INTERP_GAP_HOURS = 6

_CUM_DAYS = [0, 31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334]


def _hour_of_year(month, day, hour):
    return (_CUM_DAYS[month - 1] + (day - 1)) * 24 + hour


def esat(T):
    return 6.112 * np.exp(17.67 * T / (T + 243.5))


def _parse_year_hourly(path, year):
    is_leap = (year % 4 == 0 and year % 100 != 0) or (year % 400 == 0)
    nhours = 366 * 24 if is_leap else 365 * 24
    sums = np.zeros(nhours)
    counts = np.zeros(nhours)
    with open(path) as f:
        for row in csv.DictReader(f):
            def parse(field):
                v = row[field].split(",")[0]
                return np.nan if v in ("+9999", "9999", "") else int(v) / 10.0
            T, Td = parse("TMP"), parse("DEW")
            if np.isnan(T) or np.isnan(Td) or T < -50:
                continue
            date = dt.datetime.strptime(row["DATE"], "%Y-%m-%dT%H:%M:%S")
            idx = _hour_of_year(date.month, date.day, date.hour)
            rh = np.clip(esat(Td) / esat(T), 0.1, 1.0)
            sums[idx] += rh
            counts[idx] += 1
    with np.errstate(invalid="ignore", divide="ignore"):
        return np.where(counts > 0, sums / counts, np.nan)


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


def main():
    print(f"Building hourly humidity {YEAR_START}-{YEAR_END} from cached Deadhorse ISD files...")
    per_year = []
    for year in range(YEAR_START, YEAR_END + 1):
        path = CACHE / f"deadhorse_isd_{year}.csv"
        rh = _parse_year_hourly(path, year)
        n_valid = int(np.sum(~np.isnan(rh[:HOURS_PER_YEAR])))
        print(f"  {year}: {n_valid}/{HOURS_PER_YEAR} hours with a valid RH mean")
        per_year.append(_drop_feb29_hours(rh, year))

    flat = np.concatenate(per_year)
    assert flat.size == N_YEARS * HOURS_PER_YEAR
    n_missing = int(np.isnan(flat).sum())
    flat = _interp_short_gaps(flat, MAX_INTERP_GAP_HOURS)
    flat = _climatology_fill(flat, N_YEARS)
    n_remaining = int(np.isnan(flat).sum())
    print(f"\n{n_missing} missing hours total; {n_remaining} remain unfilled (should be 0)")

    out = FORC / f"relhum_hourly_interannual_{YEAR_START}-{YEAR_END}_frac.csv"
    np.savetxt(out, flat, fmt="%.4f")
    print(f"wrote {out.name} ({len(flat)} rows), mean RH {100*flat.mean():.0f}%")


if __name__ == "__main__":
    main()
