"""
Build a genuine multi-year (2005-2023) marine boundary (clb) forcing -- salinity, sea
temperature, and carbonate/nutrient chemistry -- from ECCO-Darwin's time-extended output,
replacing the single climatological annual-mean `clb` each sites/<name>.py currently uses
(see that file's "MARINE (clb)" comment block).

SOURCE. data.nas.nasa.gov (the no-auth portal originally used for the climatology) was
unreachable when this was built (TCP timeout, confirmed transient/server-side, not a
network problem on this end -- general internet access was fine). Per the project
owner, this instead uses the JPL ECCO Drive's extended-output mirror, which needs HTTP
Basic Auth (credentials NOT stored here -- see `ECCO_JPL_USER`/`ECCO_JPL_PASS`
environment variables below):
    https://ecco.jpl.nasa.gov/drive/files/ECCO2/LLC270/ECCO-Darwin_extension/monthly/<VAR>/<VAR>.<step>.data
Confirmed to have the SAME grid/step numbering as the original climatology's source (the
first available step, 2232, matches exactly), just with more months appended -- 407
months total vs. the original's 272-293, i.e. a genuine time-extension of the same
underlying run, not a different product.

CALENDAR. `<step>` is a cumulative count of 1200 s (20-min) MITgcm timesteps; step/72
gives cumulative days elapsed. Anchored empirically (not assumed) by walking the leap-
year fingerprint in the gap pattern between consecutive files' cumulative-day count:
the very first gap is 29 days (a leap-year February), and Gregorian 1992 is the nearest
plausible leap year for an ECCO LLC270-family run (consistent with ECCO LLC270's own
documented 1992 start) -- so file index 0 = end of January 1992, file index m = end of
calendar month (1992-01 + m). This lines up an integer number of months later with zero
drift through all 407 files (every subsequent gap matches its calendar month's real
length, 28/29/30/31 days, not just the first one) -- i.e. it is a confirmed, not
guessed, mapping. Index 156 = Jan 2005, index 383 = Dec 2023 (228 months, the window
this script extracts).

SPECIES. SALTanom -> S (+35, MITgcm's standard salinity-anomaly convention), SST ->
T (marine boundary end-member for the transported temperature field, replacing
PRDA2/`watertemp.csv` -- see CLAUDE.md "Known defects" -> "distance" era discussion and
build_interannual_met.py's docstring for why), DIC/ALK/NO3/NH4/PO4/SiO2->dSi/O2/DOC->TOC
unconverted (already mmol/m^3, matching the model's own units directly, same as the
existing climatological `clb` values in sites/<name>.py). NOT sourced this way: DIA, pH,
RDOC/CH4/N2O/SPM -- same gaps as the climatological version, see sites/<name>.py.

GRID. Same nearest-wet-cell (j, i) LLC270 indices already validated for the
climatology (scratch/ecco_darwin/extract_all_rivers.py) -- reused here, not re-derived.

CALENDAR->DAILY. Each species/river gets a 228-point monthly series, then linearly
interpolated to daily (month value anchored at that month's center day), Feb 29 excluded
(the project's "every year is exactly 365 days" convention), giving a 19*365 = 6935-row
series matching every other interannual forcing file.

Usage:
    ECCO_JPL_USER=... ECCO_JPL_PASS=... python tools/build_interannual_marine.py
"""
import base64
import os
import time
import urllib.request
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor, as_completed

import numpy as np

ROOT = Path(__file__).resolve().parent.parent
FORC = ROOT / "forcing"
CACHE = FORC / "_ecco_jpl_cache"

BASE = "https://ecco.jpl.nasa.gov/drive/files/ECCO2/LLC270/ECCO-Darwin_extension/monthly"
NX, NY = 270, 3510
LEVEL_BYTES = NX * NY * 4
MISSING = -999.0

EPOCH_YEAR, EPOCH_MONTH = 1992, 1   # file index 0 == end of this calendar month
IDX_START = (2005 - EPOCH_YEAR) * 12 + (1 - EPOCH_MONTH)   # Jan 2005 -> 156
IDX_END = (2023 - EPOCH_YEAR) * 12 + (12 - EPOCH_MONTH)    # Dec 2023 -> 383 (inclusive)
N_MONTHS = IDX_END - IDX_START + 1  # 228
YEAR_START, YEAR_END = 2005, 2023
N_YEARS = YEAR_END - YEAR_START + 1

# (j, i) nearest wet LLC270 cell to each river mouth -- unchanged from
# scratch/ecco_darwin/extract_all_rivers.py.
POINTS = {
    "colville": (2460, 3),
    "kuparuk": (2475, 2),
    "sagavanirktok": (2481, 2),
}

# raw ECCO-Darwin variable -> (model species name, unit conversion)
SPECIES = {
    "SALTanom": ("S", lambda x: x + 35.0),
    "SST": ("T", lambda x: x),
    "DIC": ("DIC", lambda x: x),
    "ALK": ("ALK", lambda x: x),
    "NO3": ("NO3", lambda x: x),
    "NH4": ("NH4", lambda x: x),
    "PO4": ("PO4", lambda x: x),
    "SiO2": ("dSi", lambda x: x),
    "O2": ("O2", lambda x: x),
    "DOC": ("TOC", lambda x: x),
}


def _auth_header():
    """Preemptive HTTP Basic Auth header. NOTE: this server does not seem to issue a
    standard 401/WWW-Authenticate challenge on the first request, so urllib's normal
    HTTPBasicAuthHandler (which waits for that challenge before sending credentials)
    silently gets an unauthenticated response instead -- confirmed directly (it
    returned a 10 KB logged-out-looking fragment with zero file links, vs. 314 KB with
    407 real file links once the Authorization header is sent on the FIRST request,
    matching how `curl -u` behaves by default). Building the header manually instead."""
    user = os.environ["ECCO_JPL_USER"]
    pw = os.environ["ECCO_JPL_PASS"]
    token = base64.b64encode(f"{user}:{pw}".encode()).decode()
    return {"Authorization": f"Basic {token}"}


def _list_steps(var):
    """Return the sorted list of integer timesteps available for `var`."""
    url = f"{BASE}/{var}"
    req = urllib.request.Request(url, headers=_auth_header())
    with urllib.request.urlopen(req, timeout=60) as resp:
        html = resp.read().decode()
    import re
    steps = sorted(set(int(s) for s in re.findall(rf'{var}\.(\d+)\.data', html)))
    return steps


def _fetch_month(var, step):
    url = f"{BASE}/{var}/{var}.{step:010d}.data"
    headers = dict(_auth_header())
    headers["Range"] = f"bytes=0-{LEVEL_BYTES - 1}"
    req = urllib.request.Request(url, headers=headers)
    for attempt in range(5):
        try:
            with urllib.request.urlopen(req, timeout=90) as resp:
                data = resp.read(LEVEL_BYTES)
            if len(data) != LEVEL_BYTES:
                raise ValueError(f"got {len(data)} bytes, expected {LEVEL_BYTES}")
            arr = np.frombuffer(data, dtype=">f4").reshape(NY, NX)
            return step, {river: float(arr[j, i]) for river, (j, i) in POINTS.items()}
        except Exception:
            time.sleep(2 * (attempt + 1))
    return step, None


def _cache_path(var, step):
    CACHE.mkdir(exist_ok=True)
    return CACHE / f"{var}.{step}.npz"


def fetch_var_monthly(var, steps_window, max_workers=6):
    """Return {river: np.array([228 monthly values])} for one raw ECCO variable,
    fetching only the steps in `steps_window` (already sliced to 2005-2023), with a
    per-(var,step) on-disk cache."""
    per_river = {river: {} for river in POINTS}
    to_fetch = []
    for step in steps_window:
        cp = _cache_path(var, step)
        if cp.exists():
            d = np.load(cp)
            for river in POINTS:
                per_river[river][step] = float(d[river])
        else:
            to_fetch.append(step)

    if to_fetch:
        with ThreadPoolExecutor(max_workers=max_workers) as ex:
            futs = {ex.submit(_fetch_month, var, s): s for s in to_fetch}
            for n, fut in enumerate(as_completed(futs), 1):
                step, vals = fut.result()
                if vals is not None:
                    for river, v in vals.items():
                        per_river[river][step] = v
                    np.savez(_cache_path(var, step), **vals)
                if n % 20 == 0 or n == len(to_fetch):
                    print(f"    {var}: {n}/{len(to_fetch)} new months fetched", flush=True)

    out = {}
    for river in POINTS:
        vals = [per_river[river].get(s) for s in steps_window]
        missing = sum(v is None for v in vals)
        if missing:
            print(f"    WARNING {var}/{river}: {missing}/{len(vals)} months missing")
        arr = np.array([v if v is not None else np.nan for v in vals])
        out[river] = arr
    return out


_CUM_DAYS_365 = [0, 31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334, 365]


def _monthly_to_daily(monthly_vals, year_start, n_years):
    """Linearly interpolate a length-(n_years*12) monthly series to a length-
    (n_years*365) daily series, anchoring each month's value at that month's center
    day (day-of-year convention: no leap days, see module docstring)."""
    n_months = n_years * 12
    assert len(monthly_vals) == n_months
    month_center_day = []
    for k in range(n_months):
        y = k // 12
        m = k % 12
        start, end = _CUM_DAYS_365[m], _CUM_DAYS_365[m + 1]
        month_center_day.append(y * 365 + (start + end) / 2.0)
    daily_axis = np.arange(n_years * 365)
    return np.interp(daily_axis, month_center_day, monthly_vals)


def main():
    print(f"Window: index {IDX_START} (Jan {YEAR_START}) .. {IDX_END} (Dec {YEAR_END}), "
          f"{N_MONTHS} months")
    steps_all = _list_steps("DIC")  # any variable works; all confirmed same 407-length list
    assert len(steps_all) > IDX_END, f"only {len(steps_all)} months available remotely"
    steps_window = steps_all[IDX_START:IDX_END + 1]
    print(f"first/last step in window: {steps_window[0]} .. {steps_window[-1]}")

    for raw_var, (species, convert) in SPECIES.items():
        print(f"=== {raw_var} -> {species} ===")
        monthly = fetch_var_monthly(raw_var, steps_window)
        for river, vals in monthly.items():
            if np.isnan(vals).any():
                # short remote gaps (rare) -- interpolate before converting, same
                # spirit as the other interannual builders
                nan_mask = np.isnan(vals)
                idx = np.arange(len(vals))
                vals = vals.copy()
                vals[nan_mask] = np.interp(idx[nan_mask], idx[~nan_mask], vals[~nan_mask])
            converted = convert(vals)
            daily = _monthly_to_daily(converted, YEAR_START, N_YEARS)
            out_path = FORC / f"{river}_{species.lower()}_marine_interannual_{YEAR_START}-{YEAR_END}.csv"
            np.savetxt(out_path, daily, fmt="%.6f")
            print(f"  {river}: mean={converted.mean():.4f}, wrote {out_path.name} "
                  f"({len(daily)} rows)")


if __name__ == "__main__":
    main()
