"""
Build a genuine multi-year (2005-2023) daily wind-driven storm-surge forcing from NOAA
CO-OPS station 9497645 (Prudhoe Bay) -- the same station/method as build_surge.py's
single 2022 year, extended across the whole interannual common window (see CLAUDE.md ->
"Interannual forcing"). Only Prudhoe Bay has a long verified water-level record on this
coast, so it remains the surge proxy for all three rivers, same as the single-year file.

    surge(t) = daily-mean observed water level - the WHOLE-RECORD mean (not a per-year
               mean, so any real multi-decadal trend is preserved rather than flattened
               into 19 separate zero-mean years)

SOURCE. `hourly_height` product (one request per year, no chunking needed -- unlike the
raw 6-min `water_level` product, which NOAA's API caps at 31 days per request). Verified
this station has NO real data before 1995 (both `hourly_height` and `water_level` return
"No data" for e.g. 1993 -- an earlier scratch check that seemed to show 1990 data via
`hourly_height` was a fluke/stale response, not reproduced on a clean re-check), so 2005
is comfortably inside its real record; 1995-2023 would also be available if the common
window were ever revisited (see CLAUDE.md).

CALENDAR / GAPS. Same convention as build_interannual_met.py: daily mean, Feb 29
dropped, short gaps (<= MAX_INTERP_GAP days) linearly interpolated, longer gaps (and any
day with zero valid hourly obs) filled with the multi-year day-of-year climatology.

Usage:  python tools/build_interannual_surge.py
"""
import json
import urllib.request
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent.parent
FORC = ROOT / "forcing"
CACHE = FORC / "_coops_cache"

YEAR_START = 2005
YEAR_END = 2023
N_YEARS = YEAR_END - YEAR_START + 1
STATION = "9497645"
MAX_INTERP_GAP = 5

URL_TMPL = ("https://api.tidesandcurrents.noaa.gov/api/prod/datagetter?"
            "product=hourly_height&datum=MSL&station=" + STATION +
            "&begin_date={year}0101&end_date={year}1231&time_zone=gmt&units=metric&format=json")

_CUM_DAYS = [0, 31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334]


def _day_of_year(month, day):
    return _CUM_DAYS[month - 1] + (day - 1)


def _fetch_year(year):
    CACHE.mkdir(exist_ok=True)
    cached = CACHE / f"prudhoe_hourly_height_{year}.json"
    if cached.exists():
        return json.loads(cached.read_text())
    url = URL_TMPL.format(year=year)
    print(f"  downloading {year}")
    with urllib.request.urlopen(url, timeout=90) as r:
        payload = json.load(r)
    cached.write_text(json.dumps(payload))
    return payload


def _daily_means(payload, year):
    """Return a length-365/366 array of daily-mean water level [m], NaN where no
    valid hourly observation exists that day."""
    is_leap = (year % 4 == 0 and year % 100 != 0) or (year % 400 == 0)
    ndays = 366 if is_leap else 365
    sums = np.zeros(ndays)
    counts = np.zeros(ndays)
    for rec in payload.get("data", []):
        v = rec["v"]
        if v in ("", "-"):
            continue
        date_str = rec["t"].split(" ")[0]  # 'YYYY-MM-DD'
        _, month, day = (int(x) for x in date_str.split("-"))
        doy = _day_of_year(month, day)
        sums[doy] += float(v)
        counts[doy] += 1
    with np.errstate(invalid="ignore", divide="ignore"):
        daily = np.where(counts > 0, sums / counts, np.nan)
    return daily


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
    print(f"Fetching NOAA CO-OPS {STATION} (Prudhoe Bay) {YEAR_START}-{YEAR_END}...")
    per_year = []
    for year in range(YEAR_START, YEAR_END + 1):
        payload = _fetch_year(year)
        daily = _daily_means(payload, year)
        n_valid = int(np.sum(~np.isnan(daily[:365])))
        print(f"  {year}: {n_valid}/365 days with a valid hourly-height mean")
        per_year.append(_drop_feb29(daily, year))

    flat = np.concatenate(per_year)
    n_missing = int(np.isnan(flat).sum())
    flat = _interp_short_gaps(flat, MAX_INTERP_GAP)
    flat = _climatology_fill(flat, N_YEARS)
    n_remaining = int(np.isnan(flat).sum())
    print(f"\n{n_missing} missing days total; {n_remaining} remain unfilled (should be 0)")

    surge = flat - np.mean(flat)  # single whole-record mean removed, not per-year
    out_path = FORC / f"surge_prudhoe_interannual_{YEAR_START}-{YEAR_END}_m.csv"
    np.savetxt(out_path, surge, fmt="%.4f")
    print(f"wrote {out_path.name} ({len(surge)} rows)")
    print(f"  surge range {surge.min():+.2f} to {surge.max():+.2f} m, "
          f"days>+0.4 m: {(surge > 0.4).sum()}")


if __name__ == "__main__":
    main()
