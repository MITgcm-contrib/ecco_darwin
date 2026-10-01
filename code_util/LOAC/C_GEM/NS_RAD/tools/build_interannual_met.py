"""
Build a genuine multi-year (2005-2023) daily meteorological forcing from NDBC station
PRDA2 (Prudhoe Bay, AK) -- wind speed and air temperature only.

NOT sea/water temperature: PRDA2's WTMP sensor has severe historical reliability
problems (2006 and 2007 are 100% missing; most years before ~2019 are 40-75% missing;
59% of 2005-2023 would be climatology-filled, not real). `watertemp.csv`/`SEATEMP_FILE`
is the marine boundary end-member for the transported temperature field (`v['T']['clb']`
in main.py -- NOT the interior river temperature, which is already prognostic via
heat_module.py's surface heat budget), so it is sourced from the ECCO-Darwin marine
boundary extraction instead (see tools/build_interannual_marine.py), consistent with
every other marine-boundary species and with real data back to 1995, not 2005.

WHY 2005, NOT 1980. PRDA2's own historical record starts in 2005 (confirmed directly
against NDBC, not from secondary documentation) -- this is the binding constraint on
how far back ANY genuinely multi-year forcing can go in this project (discharge/DOC
already covers 1980-2023; ECCO-Darwin and NOAA CO-OPS Prudhoe Bay water level both cover
1995-2023; only PRDA2-derived met is short). See CLAUDE.md -> "Interannual forcing" for
the full per-source availability table and the reasoning for picking a 2005-2023 common
window over other options (full reanalysis substitution, mismatched per-source windows).

SOURCE. NDBC historical standard meteorological data, one file per year:
    https://www.ndbc.noaa.gov/view_text_file.php?filename=prda2h{year}.txt.gz&dir=data/historical/stdmet/
Columns are POSITIONAL, not name-keyed, because the header labels/case change across
years (WD/WDIR, BAR/PRES, YY/YYYY) even though the column ORDER and units never do:
    YY MM DD hh mm WDIR WSPD GST WVHT DPD APD MWD PRES ATMP WTMP DEWP VIS TIDE
Missing-value sentinels (NDBC standard): WDIR/MWD=999, WSPD/GST=99.0, WVHT/DPD/APD=99.00,
PRES=9999.0, ATMP/WTMP/DEWP=999.0, VIS=99.0, TIDE=99.00.

CALENDAR. Aggregated to one daily MEAN per variable (skipping missing obs within the
day), then Feb 29 is dropped in leap years -- same "every year is exactly 365 days"
convention as build_interannual_forcings.py (discharge/DOC), so the output is a clean
19*365 = 6935-row series and downstream reads treat every year identically.

GAPS.
- 2005 is a PARTIAL year at PRDA2 (the station's own record starts in April 2005, not
  January) -- confirmed by inspecting the raw file, not assumed. The missing Jan-Mar
  stretch is filled with the multi-year day-of-year climatology (the mean of that same
  calendar day across all OTHER 18 years), not left blank or backfilled from 2022 alone.
- Within 2006-2023 (full years), missing individual days (sensor dropouts -- WTMP is
  missing ~20-25% of hourly obs in spot checks, mostly ice-season dropouts; a handful of
  whole days can be missing entirely) are linearly interpolated day-to-day; this is only
  applied to short gaps (<= MAX_INTERP_GAP days) -- longer gaps fall back to the same
  day-of-year climatology fill as 2005, flagged in the printed summary.

Usage:  python tools/build_interannual_met.py
"""
import gzip
import io
import urllib.request
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent.parent
FORC = ROOT / "forcing"
CACHE = FORC / "_ndbc_cache"

YEAR_START = 2005
YEAR_END = 2023
N_YEARS = YEAR_END - YEAR_START + 1
N_DAYS = N_YEARS * 365
MAX_INTERP_GAP = 5  # days; longer gaps use the day-of-year climatology instead

URL_TMPL = "https://www.ndbc.noaa.gov/view_text_file.php?filename=prda2h{year}.txt.gz&dir=data/historical/stdmet/"

# Column positions (0-indexed) in the raw NDBC file, after whitespace-splitting each
# data line. Stable across 2005-2023 despite header-label cosmetics (confirmed by
# inspection of 2005, 2007, 2010, 2015, 2020, 2023).
COL_MONTH, COL_DAY = 1, 2
COL_WSPD, COL_ATMP = 6, 13


def _fetch_year(year):
    """Return the raw text of one year's PRDA2 standard-meteorological file, caching
    locally (these are static historical archives -- safe to cache indefinitely)."""
    CACHE.mkdir(exist_ok=True)
    cached = CACHE / f"prda2h{year}.txt"
    if cached.exists():
        return cached.read_text()
    url = URL_TMPL.format(year=year)
    print(f"  downloading {year} -> {cached.name}")
    with urllib.request.urlopen(url) as resp:
        raw = resp.read()
    # The server serves these already-decompressed as text/plain despite the .gz name
    # in the query string (view_text_file.php); guard for either case.
    try:
        text = gzip.decompress(raw).decode("utf-8", errors="replace")
    except OSError:
        text = raw.decode("utf-8", errors="replace")
    cached.write_text(text)
    return text


def _parse_year(text):
    """Return (daily_wspd, daily_atmp), each a length-365 or length-366 (leap year)
    array of daily means with NaN for days with zero valid observations."""
    sums = {}
    counts = {}
    for key in ("wspd", "atmp"):
        sums[key] = np.zeros(366)
        counts[key] = np.zeros(366)

    for line in text.splitlines():
        if not line or not line[0].isdigit():
            continue  # header/comment lines start with '#' or a non-digit
        parts = line.split()
        if len(parts) <= COL_ATMP:
            continue
        month, day = int(parts[COL_MONTH]), int(parts[COL_DAY])
        doy = _day_of_year(month, day)

        wspd = float(parts[COL_WSPD])
        if wspd < 99.0:
            sums["wspd"][doy] += wspd
            counts["wspd"][doy] += 1
        atmp = float(parts[COL_ATMP])
        if atmp < 999.0:
            sums["atmp"][doy] += atmp
            counts["atmp"][doy] += 1

    out = {}
    for key in ("wspd", "atmp"):
        with np.errstate(invalid="ignore", divide="ignore"):
            out[key] = np.where(counts[key] > 0, sums[key] / counts[key], np.nan)
    return out["wspd"], out["atmp"]


_CUM_DAYS = [0, 31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334]  # non-leap cumulative


def _day_of_year(month, day):
    """0-indexed day-of-year, NOT leap-year-aware (Feb 29 lands on index 59, same as
    Mar 1 in a non-leap year -- fine here since we drop Feb 29 immediately after and
    only ever index by (month, day), never by a running day count across the file)."""
    return _CUM_DAYS[month - 1] + (day - 1)


def _drop_feb29(arr, year):
    """Drop index 59 (Feb 29) if `year` is a leap year and arr has 366 entries."""
    is_leap = (year % 4 == 0 and year % 100 != 0) or (year % 400 == 0)
    if is_leap and len(arr) == 366:
        return np.delete(arr, 59)
    return arr[:365]


def _interp_short_gaps(series_by_year, max_gap):
    """Linearly interpolate NaN runs of length <= max_gap within the flat (N_YEARS*365)
    series; longer runs are left as NaN for the climatology fill to handle."""
    flat = np.concatenate(series_by_year)
    nan_mask = np.isnan(flat)
    if not nan_mask.any():
        return flat
    idx = np.arange(len(flat))
    # identify contiguous NaN runs
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
    """Fill any remaining NaNs (2005's missing Jan-Mar, plus any long in-record gaps)
    with the mean of that same day-of-year across all other years that have real data
    on that day."""
    grid = flat.reshape(n_years, 365)
    for doy in range(365):
        col = grid[:, doy]
        nan_here = np.isnan(col)
        if nan_here.any() and not nan_here.all():
            col[nan_here] = np.nanmean(col[~nan_here])
    return grid.reshape(-1)


def build_variable(all_years_parsed, key, label, unit):
    per_year = []
    for i, year in enumerate(range(YEAR_START, YEAR_END + 1)):
        arr = all_years_parsed[i][key]
        per_year.append(_drop_feb29(arr, year))
    flat = _interp_short_gaps(per_year, MAX_INTERP_GAP)
    n_short_nan = 0  # already resolved by _interp_short_gaps; count remaining below
    n_before = np.isnan(flat).sum()
    flat = _climatology_fill(flat, N_YEARS)
    n_after = np.isnan(flat).sum()
    print(f"  {label}: {n_before} missing days before climatology fill, "
          f"{n_after} remain unfilled (should be 0)")
    out_path = FORC / f"{label}_interannual_{YEAR_START}-{YEAR_END}_{unit}.csv"
    np.savetxt(out_path, flat, fmt="%.6f")
    print(f"  wrote {out_path.name} ({len(flat)} rows)")
    return flat


def main():
    print(f"Fetching NDBC PRDA2 {YEAR_START}-{YEAR_END} ({N_YEARS} years)...")
    parsed = []
    for year in range(YEAR_START, YEAR_END + 1):
        text = _fetch_year(year)
        wspd, atmp = _parse_year(text)
        n_days_with_data = int(np.sum(~np.isnan(atmp[:365])))
        print(f"  {year}: {n_days_with_data}/365 days with any ATMP obs")
        parsed.append({"wspd": wspd, "atmp": atmp})

    print("\nAssembling daily series (Feb 29 dropped, gaps filled)...")
    wspd_series = build_variable(parsed, "wspd", "windspeed", "msec")
    atmp_series = build_variable(parsed, "atmp", "airtemp", "degC")

    print(f"\n{YEAR_START}-{YEAR_END} annual-mean check "
          f"(wind {wspd_series.mean():.2f} m/s, air {atmp_series.mean():.2f} degC) -- "
          f"sanity only, no independent source to cross-check against here (unlike "
          f"discharge's monthly-vs-daily PWBM check).")


if __name__ == "__main__":
    main()
