"""
Build a genuine multi-year (2005-2023) HOURLY solar-radiation forcing from NOAA GML's
Barrow (Utqiagvik) Atmospheric Baseline Observatory -- REAL 1-minute pyranometer
observations, not the NARR reanalysis build_interannual_solar_subdaily.py used.
Replaces that script as the solar source (kept for reference/comparison, not deleted).

WHY BARROW OVER NARR. Checked directly: NOAA GML's qcrad_v3 product
(gml.noaa.gov/aftp/data/radiation/baseline/qcrad_v3/brw/) is a clean, HEADER-LABELED,
one-file-per-day ASCII product with genuine 1-minute global shortwave (`GSW`, W/m^2)
-- real observations, not reanalysis, far finer native resolution than NARR's 3-hourly,
and MUCH smaller/more reliable to download (small per-day files, no multi-hundred-MB
single-file downloads with constant connection drops).

WHY NOT ALSO SWITCH WIND/AIR TEMP/HUMIDITY TO BARROW (this file's air temp/RH/wind
columns could supply those too). Checked directly: Barrow is ~330 km WEST of Prudhoe
Bay -- FARTHER from all three rivers than Prudhoe Bay already is, not closer. PRDA2 +
Deadhorse Airport (already built, tested, 91-98% real) remain the better regional
proxy for those three variables precisely because they are the closer station; Barrow
is used for solar ONLY because Prudhoe Bay has no solar-measuring instrument at all
(confirmed: PRDA2's raw file has no radiation column), and Barrow's BSRN status is
specifically why it has one. See CLAUDE.md -> "Interannual forcing" -> "Sub-daily
(diurnal) forcing" for the full reasoning trail.

SOURCE. One small file per day:
    https://gml.noaa.gov/aftp/data/radiation/baseline/qcrad_v3/brw/{year}/brw_{year}{mm}{dd}.qdat
Header-labeled whitespace columns; `GSW` is global (downwelling) shortwave, W/m^2,
1-minute resolution (1440 rows/day). -9999.9 is the missing-value sentinel.

CALENDAR. Aggregated 1-minute -> hourly mean (matching build_interannual_met_hourly.py's
air-temperature resolution, so every sub-daily forcing in this project shares one
common row interval and can be read by the same row_interval_sec=3600 code path).
8760 hours/year (Feb 29 dropped), 19*8760 = 166,440 rows total.

GAPS. Short gaps (<= MAX_INTERP_GAP_HOURS) linearly interpolated; longer gaps (or
whole missing days -- a day's file simply doesn't exist for some outage dates) filled
with the multi-year hour-of-year climatology, same convention as every other
sub-daily/daily builder here.

Usage:  python tools/build_interannual_solar_barrow.py
"""
import datetime as dt
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

import numpy as np
import urllib.request
import urllib.error

ROOT = Path(__file__).resolve().parent.parent
FORC = ROOT / "forcing"
CACHE = FORC / "_barrow_qcrad_cache"

YEAR_START = 2005
YEAR_END = 2023
N_YEARS = YEAR_END - YEAR_START + 1
HOURS_PER_YEAR = 365 * 24
MAX_INTERP_GAP_HOURS = 6
MISSING = -9999.9

URL_TMPL = ("https://gml.noaa.gov/aftp/data/radiation/baseline/qcrad_v3/brw/"
            "{year}/brw_{year}{month:02d}{day:02d}.qdat")


def _fetch_day(date):
    CACHE.mkdir(exist_ok=True)
    dest = CACHE / f"brw_{date:%Y%m%d}.qdat"
    if dest.exists():
        return dest
    url = URL_TMPL.format(year=date.year, month=date.month, day=date.day)
    try:
        with urllib.request.urlopen(url, timeout=30) as resp:
            data = resp.read()
        dest.write_bytes(data)
        return dest
    except urllib.error.HTTPError:
        return None  # day genuinely missing from the archive (station outage)


def _parse_day_gsw_hourly(path):
    """Return a length-24 array of hourly-mean GSW, NaN for hours with no valid
    1-minute reading."""
    with open(path) as f:
        header = f.readline().split()
        gsw_col = header.index("GSW")
        sums = np.zeros(24)
        counts = np.zeros(24)
        for line in f:
            parts = line.split()
            if len(parts) <= gsw_col:
                continue
            minute_of_day = int(parts[1])
            hour = minute_of_day // 60
            if hour >= 24:
                continue
            val = float(parts[gsw_col])
            if val > MISSING + 1e-6:
                sums[hour] += val
                counts[hour] += 1
    with np.errstate(invalid="ignore", divide="ignore"):
        return np.where(counts > 0, sums / counts, np.nan)


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
    all_dates = []
    for year in range(YEAR_START, YEAR_END + 1):
        d = dt.date(year, 1, 1)
        while d.year == year:
            if not (d.month == 2 and d.day == 29):  # Feb 29 dropped up front
                all_dates.append(d)
            d += dt.timedelta(days=1)
    print(f"Fetching {len(all_dates)} days of Barrow qcrad_v3 data "
          f"({YEAR_START}-{YEAR_END})...")

    with ThreadPoolExecutor(max_workers=10) as ex:
        futs = {ex.submit(_fetch_day, d): d for d in all_dates}
        paths = {}
        for n, fut in enumerate(as_completed(futs), 1):
            d = futs[fut]
            paths[d] = fut.result()
            if n % 500 == 0 or n == len(all_dates):
                print(f"  downloaded {n}/{len(all_dates)}", flush=True)

    n_missing_days = sum(1 for d in all_dates if paths[d] is None)
    print(f"{n_missing_days} days missing entirely from the archive")

    flat = np.full(len(all_dates) * 24, np.nan)
    for i, d in enumerate(all_dates):
        if paths[d] is not None:
            flat[i * 24:(i + 1) * 24] = _parse_day_gsw_hourly(paths[d])

    assert flat.size == N_YEARS * HOURS_PER_YEAR
    n_missing_hours = int(np.isnan(flat).sum())
    flat = _interp_short_gaps(flat, MAX_INTERP_GAP_HOURS)
    flat = _climatology_fill(flat, N_YEARS)
    n_remaining = int(np.isnan(flat).sum())
    print(f"{n_missing_hours} missing hours total; {n_remaining} remain unfilled "
          f"(should be 0)")

    flat = np.maximum(flat, 0.0)  # BSRN GSW can read slightly negative at night (instrument noise)
    out = FORC / f"solar_hourly_interannual_{YEAR_START}-{YEAR_END}_Wm2.csv"
    np.savetxt(out, flat, fmt="%.4f")
    print(f"wrote {out.name} ({len(flat)} rows), mean={flat.mean():.1f} W/m^2, "
          f"max={flat.max():.1f} W/m^2")


if __name__ == "__main__":
    main()
