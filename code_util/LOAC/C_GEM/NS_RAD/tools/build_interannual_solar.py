"""
Build a genuine multi-year (2005-2023) daily solar-radiation forcing from NOAA's North
American Regional Reanalysis (NARR), replacing solarradiation.csv (whose original
single-2022-year source/build script could not be identified or reproduced -- see
tools/build_interannual_met.py's docstring).

WHY NARR, NOT NSRDB. NSRDB's standard (PSM3) product's documented extent tops out
around 60N; Prudhoe Bay is 70.4N. NSRDB does have a dedicated Polar product above 60N,
but it only starts in 2013 -- 8 years short of this window -- AND its API host
(developer.nrel.gov) was unreachable (DNS NXDOMAIN) when this was built. NARR instead:
covers all of North America including northern Alaska (confirmed directly: grid point
(70.17N, -148.55E) is ~0.25 degrees from Prudhoe Bay, well inside the domain, which
reaches 85.3N), requires NO authentication, and its daily-mean downward shortwave
radiation flux (`dswrf`) at that point is physically sane (Jan mean 0.39 W/m^2 -- near
the polar-night floor, July mean 262 W/m^2 -- matching the existing single-year file's
summer peak of ~330-400 W/m^2 order of magnitude). Available 1979-present, so 2005 is
not this source's own constraint (PRDA2 wind/air-temp remains the binding one for the
overall common window -- see CLAUDE.md -> "Interannual forcing").

SOURCE. One NetCDF file per year, no auth:
    https://psl.noaa.gov/thredds/fileServer/Datasets/NARR/Dailies/monolevel/dswrf.{year}.nc
The HTTP/2 connection reliably drops partway through (confirmed: consistently around
9-11 MB into each ~65 MB file, both on HTTP/2 and forced HTTP/1.1 -- a server/path
issue, not a timeout, since re-requesting with Range picks up instantly) -- so this
downloads with an explicit resume-on-drop loop rather than a single request.

CALENDAR. NARR's own daily files are real calendar years (365 or 366 rows); Feb 29 is
dropped for leap years, same "every year is exactly 365 days" convention as every other
interannual forcing here.

Usage:  python tools/build_interannual_solar.py
"""
import http.client
import subprocess
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent.parent
FORC = ROOT / "forcing"
CACHE = FORC / "_narr_cache"

YEAR_START = 2005
YEAR_END = 2023
N_YEARS = YEAR_END - YEAR_START + 1

URL_TMPL = "https://psl.noaa.gov/thredds/fileServer/Datasets/NARR/Dailies/monolevel/dswrf.{year}.nc"
TARGET_LAT, TARGET_LON = 70.4, -148.5  # Prudhoe Bay, matches the other met sources' site
MAX_RESUME_ATTEMPTS = 15


def _download_with_resume(url, dest):
    """curl with -C - (resume) in a retry loop -- the server drops the connection
    every ~10 MB (confirmed on both HTTP/2 and HTTP/1.1), so a single request never
    completes a ~65 MB file; repeated resumes do."""
    for attempt in range(MAX_RESUME_ATTEMPTS):
        result = subprocess.run(
            ["curl", "-s", "--http1.1", "--max-time", "90", "-C", "-", "-o", str(dest), url],
            capture_output=True,
        )
        if dest.exists():
            # Compare against the server's declared size via a HEAD request once;
            # cheap enough to repeat, and avoids trusting curl's own exit code (a
            # dropped stream sometimes still exits 0 with a partial file).
            head = subprocess.run(
                ["curl", "-sI", "--http1.1", "--max-time", "30", url],
                capture_output=True, text=True,
            )
            content_length = None
            for line in head.stdout.splitlines():
                if line.lower().startswith("content-length:"):
                    content_length = int(line.split(":", 1)[1].strip())
            if content_length and dest.stat().st_size >= content_length:
                return
    raise RuntimeError(f"failed to fully download {url} after {MAX_RESUME_ATTEMPTS} attempts "
                       f"(got {dest.stat().st_size if dest.exists() else 0} bytes)")


def _fetch_year(year):
    CACHE.mkdir(exist_ok=True)
    dest = CACHE / f"dswrf.{year}.nc"
    if not dest.exists():
        print(f"  downloading {year} (resuming through connection drops)...")
        _download_with_resume(URL_TMPL.format(year=year), dest)
    return dest


def _extract_point_series(path, year):
    import netCDF4 as nc
    ds = nc.Dataset(path)
    lat = ds.variables["lat"][:]
    lon = ds.variables["lon"][:]
    lon_adj = np.where(lon > 180, lon - 360, lon)
    dist = (lat - TARGET_LAT) ** 2 + (lon_adj - TARGET_LON) ** 2
    j, i = np.unravel_index(np.argmin(dist), dist.shape)
    series = np.asarray(ds.variables["dswrf"][:, j, i], dtype=float)
    return series


def _drop_feb29(arr, year):
    is_leap = (year % 4 == 0 and year % 100 != 0) or (year % 400 == 0)
    if is_leap and len(arr) == 366:
        return np.delete(arr, 59)
    return arr[:365]


def main():
    print(f"Fetching NARR dswrf {YEAR_START}-{YEAR_END} at ({TARGET_LAT}N, {TARGET_LON}E)...")
    per_year = []
    for year in range(YEAR_START, YEAR_END + 1):
        path = _fetch_year(year)
        series = _extract_point_series(path, year)
        print(f"  {year}: {len(series)} days, mean={series.mean():.1f} W/m^2")
        per_year.append(_drop_feb29(series, year))

    flat = np.concatenate(per_year)
    assert flat.size == N_YEARS * 365, f"expected {N_YEARS*365}, got {flat.size}"
    assert not np.isnan(flat).any(), "unexpected NaN in a reanalysis product"

    out = FORC / f"solar_interannual_{YEAR_START}-{YEAR_END}_Wm2.csv"
    np.savetxt(out, flat, fmt="%.4f")
    print(f"\nwrote {out.name} ({len(flat)} rows), mean={flat.mean():.1f} W/m^2, "
          f"max={flat.max():.1f} W/m^2")


if __name__ == "__main__":
    main()
