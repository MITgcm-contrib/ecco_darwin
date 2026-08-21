"""
Fetch OBSERVED USGS daily discharge, full historical record (1980-2023), for the three
gauged North Slope rivers -- the ground truth to compare the interannual PWBM-derived
discharge forcing (tools/build_interannual_forcings.py) against in
docs/ns_rad_interannual.pdf.

This is a genuinely different comparison than tools/fetch_discharge.py, which fetches
ONE year (2022) to build the model's own discharge FORCING. This script fetches the
FULL multi-decade record from the same three gauges to VALIDATE that forcing's
interannual variability against reality, year by year, 1980-2023.

Canning has no gauge (see tools/fetch_discharge.py's RECONSTRUCTED note) so it is not
included here either -- there is nothing to validate against.

Output: docs/interannual_discharge_obs.json, ANNUAL statistics only (mean m3/s,
n_valid_days) per gauge per calendar year -- not the raw daily record, which would be
tens of thousands of rows duplicating USGS's own database for no benefit here (the
report only ever compares annual means, matching the model's own forcing granularity).
Years with fewer than 300 valid days are dropped (an incomplete year biases an annual
mean, especially since the missing days are disproportionately the ice-affected winter
months -- see tools/fetch_discharge.py's caveat on the 'A:e' estimated flag).

Usage:  python tools/fetch_interannual_discharge_obs.py
"""
import json
import sys
import urllib.request
from pathlib import Path

CFS_TO_CMS = 0.3048 ** 3   # exact: 1 ft = 0.3048 m by definition
OUT = Path(__file__).resolve().parent.parent / "docs" / "interannual_discharge_obs.json"
YEAR_START, YEAR_END = 1980, 2023
MIN_VALID_DAYS = 300

GAUGED = {
    "colville":      ("15875000", "Colville R at Umiat"),
    "kuparuk":       ("15896000", "Kuparuk R nr Deadhorse"),
    "sagavanirktok": ("15908000", "Sagavanirktok R nr Pump Sta 3"),
}

MISSING = {"", "Ice", "Eqp", "Ssn", "--", "Bkw", "Dis", "Rat", "Mnt", "ZFl"}


def _fetch(site, start, end):
    """Return {date_string: m3/s} for a USGS gauge over [start, end]. Same query/parse
    as tools/fetch_discharge.py, just a much longer date range."""
    url = ("https://waterservices.usgs.gov/nwis/dv/?format=rdb"
           f"&sites={site}&startDT={start}&endDT={end}"
           "&parameterCd=00060&statCd=00003")
    with urllib.request.urlopen(url, timeout=300) as fh:
        text = fh.read().decode("utf-8", "replace")

    out = {}
    for line in text.splitlines():
        if line.startswith("#"):
            continue
        col = line.split("\t")
        if len(col) < 4 or col[0] != "USGS":
            continue
        raw = col[3].strip()
        if raw in MISSING:
            continue
        out[col[2]] = float(raw) * CFS_TO_CMS
    return out


def main():
    result = {}
    for name, (gauge, label) in GAUGED.items():
        print(f"{name} <- USGS {gauge} ({label}), {YEAR_START}-{YEAR_END}")
        daily = _fetch(gauge, f"{YEAR_START}-01-01", f"{YEAR_END}-12-31")
        by_year = {}
        for date, val in daily.items():
            by_year.setdefault(date[:4], []).append(val)
        annual = {}
        for year, vals in sorted(by_year.items()):
            if len(vals) >= MIN_VALID_DAYS:
                annual[year] = {"mean_m3s": sum(vals) / len(vals), "n_valid_days": len(vals)}
        result[name] = {"gauge": gauge, "label": label, "annual": annual}
        print(f"  {len(annual)}/{YEAR_END - YEAR_START + 1} years with >= {MIN_VALID_DAYS} valid days "
              f"(of {len(by_year)} years with any data)")

    OUT.parent.mkdir(exist_ok=True)
    OUT.write_text(json.dumps(result, indent=1))
    print(f"wrote {OUT.relative_to(OUT.parent.parent)}")


if __name__ == "__main__":
    sys.exit(main())
