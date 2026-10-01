"""
Build a genuine multi-year (2005-2023) 3-HOURLY solar-radiation forcing from NOAA
NARR -- resolves the diurnal cycle, unlike build_interannual_solar.py's daily-mean
version (and unlike every forcing in this project before this round). See CLAUDE.md ->
"Interannual forcing" -> "Sub-daily (diurnal) forcing".

WHY 3-HOURLY, NOT HOURLY. Checked directly (not assumed): NARR's `dswrf` has no truly
hourly product at all -- `Datasets/NARR/monolevel/dswrf.{year}.nc` (distinct from the
`Dailies/monolevel` daily-mean version already used) is 3-hourly (2920 records/year =
365*8, confirmed exact 3-hour spacing in the time coordinate). This is the finest
native resolution available from this no-auth source; true hourly would need a
different product (e.g. ERA5, which needs a Copernicus CDS credential this project has
otherwise avoided -- see CLAUDE.md's "Interannual forcing" ERA5 discussion).

SOURCE. Same download mechanics as build_interannual_solar.py (the connection drops
every ~10 MB regardless of which NARR product; resumed with the same retry loop) --
just a bigger file per year (~300 MB vs ~65 MB for the daily version, since sub-setting
one grid point still requires downloading the whole global 3-hourly field, same as the
daily case).

CALENDAR. 3-hour-resolution rows, 2920/year (365*8, Feb 29's 8 rows dropped in leap
years) -- 19*2920 = 55,480 rows total. NOT the same row count or spacing as any other
forcing file in this project; file_module.py's generalized `row_interval_sec` (10800 s
here) is what lets exfread() interpolate it correctly.

Usage:  python tools/build_interannual_solar_subdaily.py
"""
import subprocess
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent.parent
FORC = ROOT / "forcing"
CACHE = FORC / "_narr_hourly_cache"

YEAR_START = 2005
YEAR_END = 2023
N_YEARS = YEAR_END - YEAR_START + 1
STEPS_PER_DAY = 8  # 3-hourly
STEPS_PER_YEAR = 365 * STEPS_PER_DAY  # 2920, Feb 29 dropped

URL_TMPL = "https://psl.noaa.gov/thredds/fileServer/Datasets/NARR/monolevel/dswrf.{year}.nc"
TARGET_LAT, TARGET_LON = 70.4, -148.5
MAX_RESUME_ATTEMPTS = 150


def _is_valid_netcdf(path):
    """Byte count matching Content-Length is NOT sufficient (seen directly: a 2006
    download matched size exactly but was still HDF-corrupt, presumably a resume/
    proxy byte-alignment glitch) -- actually open the file."""
    try:
        import netCDF4 as nc
        ds = nc.Dataset(path)
        ds.variables["dswrf"]
        ds.close()
        return True
    except Exception:
        return False


def _download_with_resume(url, dest):
    for attempt in range(MAX_RESUME_ATTEMPTS):
        subprocess.run(
            ["curl", "-s", "--http1.1", "--max-time", "60", "-C", "-", "-o", str(dest), url],
            capture_output=True,
        )
        if dest.exists():
            head = subprocess.run(
                ["curl", "-sI", "--http1.1", "--max-time", "30", url],
                capture_output=True, text=True,
            )
            content_length = None
            for line in head.stdout.splitlines():
                if line.lower().startswith("content-length:"):
                    content_length = int(line.split(":", 1)[1].strip())
            if content_length and dest.stat().st_size >= content_length:
                if _is_valid_netcdf(dest):
                    return
                # Corrupt despite matching size -- start over, not just resume, since
                # a byte-alignment glitch could persist across further appends.
                dest.unlink()
    raise RuntimeError(f"failed to fully download {url} after {MAX_RESUME_ATTEMPTS} attempts")


def _fetch_year(year):
    CACHE.mkdir(exist_ok=True)
    dest = CACHE / f"dswrf3h.{year}.nc"
    # BUG FIXED (found live): this used to be `if not dest.exists()`, so a leftover
    # PARTIAL file from a previous run that was killed mid-download (not cleanly
    # exhausted its retry loop, so never got unlinked) was mistaken for a finished
    # download and skipped -- silently handing a corrupt file to the caller.
    if not dest.exists() or not _is_valid_netcdf(dest):
        print(f"  downloading {year} (resuming through connection drops, ~300 MB)...")
        _download_with_resume(URL_TMPL.format(year=year), dest)
    return dest


def _extract_point_series(path):
    import netCDF4 as nc
    ds = nc.Dataset(path)
    lat = ds.variables["lat"][:]
    lon = ds.variables["lon"][:]
    lon_adj = np.where(lon > 180, lon - 360, lon)
    dist = (lat - TARGET_LAT) ** 2 + (lon_adj - TARGET_LON) ** 2
    j, i = np.unravel_index(np.argmin(dist), dist.shape)
    return np.asarray(ds.variables["dswrf"][:, j, i], dtype=float)


def _drop_feb29_steps(arr, year):
    is_leap = (year % 4 == 0 and year % 100 != 0) or (year % 400 == 0)
    expected_full = 366 * STEPS_PER_DAY if is_leap else 365 * STEPS_PER_DAY
    assert len(arr) == expected_full, (len(arr), expected_full)
    if is_leap:
        feb29_start = 59 * STEPS_PER_DAY
        return np.delete(arr, np.arange(feb29_start, feb29_start + STEPS_PER_DAY))
    return arr


def main():
    print(f"Fetching NARR 3-hourly dswrf {YEAR_START}-{YEAR_END} at "
          f"({TARGET_LAT}N, {TARGET_LON}E)...")
    per_year = []
    for year in range(YEAR_START, YEAR_END + 1):
        path = _fetch_year(year)
        series = _extract_point_series(path)
        series = _drop_feb29_steps(series, year)
        print(f"  {year}: {len(series)} steps, mean={series.mean():.1f} W/m^2")
        per_year.append(series)

    flat = np.concatenate(per_year)
    assert flat.size == N_YEARS * STEPS_PER_YEAR, (flat.size, N_YEARS * STEPS_PER_YEAR)
    assert not np.isnan(flat).any(), "unexpected NaN in a reanalysis product"

    out = FORC / f"solar_3hourly_interannual_{YEAR_START}-{YEAR_END}_Wm2.csv"
    np.savetxt(out, flat, fmt="%.4f")
    print(f"\nwrote {out.name} ({len(flat)} rows), mean={flat.mean():.1f} W/m^2, "
          f"max={flat.max():.1f} W/m^2")


if __name__ == "__main__":
    main()
