"""
Build a genuine multi-year (2005-2023) daily relative-humidity forcing from Deadhorse
Airport (NOAA ISD), extending build_humidity.py's single-2022-year method across the
interannual common window (see CLAUDE.md -> "Interannual forcing"). Same station
(PASC, USAF 700637 / WBAN 27406), same Magnus-form RH calculation from hourly
temperature/dewpoint, just looped over years instead of one fixed file.

Deadhorse's ISD record goes back at least to 1980 (checked directly against
ncei.noaa.gov's per-year access files), so 2005 is comfortably inside it -- this
source is NOT the reason the common window starts at 2005 (that's NDBC PRDA2's own
wind/air-temp record, see build_interannual_met.py).

CALENDAR / GAPS. Same convention as the other interannual builders: daily mean
(skipping missing hourly obs), Feb 29 dropped, short gaps (<= MAX_INTERP_GAP days)
linearly interpolated, longer gaps filled with the multi-year day-of-year climatology.

Usage:  python tools/build_interannual_humidity.py
"""
import csv
import datetime as dt
import urllib.request
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent.parent
FORC = ROOT / "forcing"
CACHE = FORC / "_isd_cache"

YEAR_START = 2005
YEAR_END = 2023
N_YEARS = YEAR_END - YEAR_START + 1
MAX_INTERP_GAP = 5

ISD_URL_TMPL = "https://www.ncei.noaa.gov/data/global-hourly/access/{year}/70063727406.csv"

_CUM_DAYS = [0, 31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334]


def _day_of_year(month, day):
    return _CUM_DAYS[month - 1] + (day - 1)


def esat(T):
    """Saturation vapour pressure over water [hPa], Magnus form, T in degC."""
    return 6.112 * np.exp(17.67 * T / (T + 243.5))


def _fetch_year(year):
    CACHE.mkdir(exist_ok=True)
    cached = CACHE / f"deadhorse_isd_{year}.csv"
    if not cached.exists():
        print(f"  downloading {year} -> {cached.name}")
        urllib.request.urlretrieve(ISD_URL_TMPL.format(year=year), cached)
    return cached


def _parse_year(path, year):
    is_leap = (year % 4 == 0 and year % 100 != 0) or (year % 400 == 0)
    ndays = 366 if is_leap else 365
    sums = np.zeros(ndays)
    counts = np.zeros(ndays)
    with open(path) as f:
        for row in csv.DictReader(f):
            def parse(field):
                v = row[field].split(",")[0]
                return np.nan if v in ("+9999", "9999", "") else int(v) / 10.0
            T, Td = parse("TMP"), parse("DEW")
            if np.isnan(T) or np.isnan(Td) or T < -50:
                continue
            date = dt.datetime.strptime(row["DATE"], "%Y-%m-%dT%H:%M:%S")
            doy = _day_of_year(date.month, date.day)
            rh = np.clip(esat(Td) / esat(T), 0.1, 1.0)
            sums[doy] += rh
            counts[doy] += 1
    with np.errstate(invalid="ignore", divide="ignore"):
        return np.where(counts > 0, sums / counts, np.nan)


def _drop_feb29(arr, year):
    is_leap = (year % 4 == 0 and year % 100 != 0) or (year % 400 == 0)
    if is_leap and len(arr) == 366:
        return np.delete(arr, 59)
    return arr[:365]


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
    grid = flat.reshape(n_years, 365)
    for doy in range(365):
        col = grid[:, doy]
        nan_here = np.isnan(col)
        if nan_here.any() and not nan_here.all():
            col[nan_here] = np.nanmean(col[~nan_here])
    return grid.reshape(-1)


def main():
    print(f"Fetching Deadhorse ISD {YEAR_START}-{YEAR_END}...")
    per_year = []
    for year in range(YEAR_START, YEAR_END + 1):
        path = _fetch_year(year)
        rh = _parse_year(path, year)
        n_valid = int(np.sum(~np.isnan(rh[:365])))
        print(f"  {year}: {n_valid}/365 days with a valid RH mean")
        per_year.append(_drop_feb29(rh, year))

    flat = np.concatenate(per_year)
    n_missing = int(np.isnan(flat).sum())
    flat = _interp_short_gaps(flat, MAX_INTERP_GAP)
    flat = _climatology_fill(flat, N_YEARS)
    n_remaining = int(np.isnan(flat).sum())
    print(f"\n{n_missing} missing days total; {n_remaining} remain unfilled (should be 0)")

    out = FORC / f"relhum_interannual_{YEAR_START}-{YEAR_END}_frac.csv"
    np.savetxt(out, flat, fmt="%.4f")
    print(f"wrote {out.name} ({len(flat)} rows), mean RH {100*flat.mean():.0f}%")


if __name__ == "__main__":
    main()
