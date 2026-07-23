"""
Build C-GEM daily discharge forcings for the four major North Slope rivers.

Writes one <river>_river_discharge_2022_m3sec.csv per site into forcing/,
in the exact format C-GEM's file_module.exfread expects: 365 values, one per line,
2 decimal places, CRLF endings, no header, no trailing newline.

Provenance
----------
Colville, Kuparuk, Sagavanirktok  -- USGS NWIS daily values (parameterCd 00060,
statCd 00003), cfs -> m3/s via 0.3048**3. This reproduces the original
colville_river_discharge_2022_m3sec.csv byte-for-byte, which is the check that the
same recipe is being applied to all four sites.

Canning -- NO 2022 record. USGS 15955000 (Canning R ab Staines R) ran 2008-06-23 to
2012-09-30 only. Reconstructed from Hulahula R (15980000), the nearest gauge active in
2022, scaled by the ratio of means over their 731-day common record. See RECONSTRUCTED
below. This is a provisional forcing, not an observation -- see CLAUDE.md.

Gauge caveats that affect interpretation, not the mechanics:
  - Colville is gauged at Umiat, far upstream of the delta; Sagavanirktok at Pump
    Station 3, likewise well upstream. Both therefore UNDERSTATE discharge at the
    river mouth the model actually represents. Kuparuk near Deadhorse is close to
    tidewater and needs no such caveat.
  - ~250 of 365 days each year carry the USGS 'A:e' (approved, estimated) flag --
    these are the ice-affected winter months, modelled rather than measured. C-GEM's
    ice gate suppresses biogeochemistry across most of that period anyway.

Usage:  python tools/fetch_discharge.py
"""

import sys
import urllib.request
from pathlib import Path

# Exact: 1 ft = 0.3048 m by definition, so 1 ft^3 = 0.3048^3 m^3. Using the rounded
# 0.0283168 instead shifts ~10 days/yr by 0.01 m3/s and fails to reproduce the
# original Colville file byte-for-byte, so keep the exact form.
CFS_TO_CMS = 0.3048 ** 3
OUT_DIR = Path(__file__).resolve().parent.parent / "forcing"
YEAR = 2022

# site key -> (USGS gauge id, human label)
GAUGED = {
    "colville":      ("15875000", "Colville R at Umiat"),
    "kuparuk":       ("15896000", "Kuparuk R nr Deadhorse"),
    "sagavanirktok": ("15908000", "Sagavanirktok R nr Pump Sta 3"),
}

# site key -> (donor gauge id, donor label, target gauge id, overlap window)
RECONSTRUCTED = {
    "canning": ("15980000", "Hulahula R nr Kaktovik", "15955000",
                ("2010-10-01", "2012-09-30")),
}

MISSING = {"", "Ice", "Eqp", "Ssn", "--", "Bkw", "Dis", "Rat", "Mnt", "ZFl"}


def _fetch(site, start, end):
    """Return {date_string: m3/s} for a USGS gauge over [start, end]."""
    url = ("https://waterservices.usgs.gov/nwis/dv/?format=rdb"
           f"&sites={site}&startDT={start}&endDT={end}"
           "&parameterCd=00060&statCd=00003")
    with urllib.request.urlopen(url, timeout=120) as fh:
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


def _year_series(daily, year):
    """Order a {date: value} dict into a dense 365-value list, gap-filled."""
    import datetime as dt

    days = []
    d = dt.date(year, 1, 1)
    while d.year == year:
        if not (d.month == 2 and d.day == 29):  # C-GEM's axis is fixed at 365
            days.append(d.isoformat())
        d += dt.timedelta(days=1)

    vals, last = [], None
    for key in days:
        v = daily.get(key)
        if v is None:
            v = last if last is not None else 0.0
        vals.append(v)
        last = v
    return vals


def _write(name, values):
    if len(values) != 365:
        raise SystemExit(f"{name}: expected 365 values, got {len(values)}")
    path = OUT_DIR / f"{name}_river_discharge_{YEAR}_m3sec.csv"
    body = "\r\n".join(f"{v:.2f}" for v in values)  # CRLF, no trailing newline
    path.write_text(body, encoding="utf-8", newline="")
    print(f"  wrote {path.name:48s} "
          f"min={min(values):8.2f} max={max(values):9.2f} "
          f"mean={sum(values)/len(values):8.2f} m3/s")


def main():
    OUT_DIR.mkdir(exist_ok=True)

    for name, (gauge, label) in GAUGED.items():
        print(f"{name} <- USGS {gauge} ({label})")
        series = _year_series(_fetch(gauge, f"{YEAR}-01-01", f"{YEAR}-12-31"), YEAR)
        _write(name, series)

    for name, (donor, dlabel, target, (ov0, ov1)) in RECONSTRUCTED.items():
        print(f"{name} <- RECONSTRUCTED from USGS {donor} ({dlabel})")
        tgt = _fetch(target, ov0, ov1)
        don = _fetch(donor, ov0, ov1)
        both = sorted(set(tgt) & set(don))
        if not both:
            raise SystemExit(f"{name}: no overlap between {target} and {donor}")

        # Ratio of MEANS, not median of daily ratios: preserves annual water volume,
        # which is what the transport and carbonate budgets are sensitive to. The
        # median daily ratio (~3.7) is higher because low-flow days inflate it.
        tm = sum(tgt[k] for k in both) / len(both)
        dm = sum(don[k] for k in both) / len(both)
        scale = tm / dm
        print(f"  overlap={len(both)}d  target_mean={tm:.2f}  donor_mean={dm:.2f}  "
              f"scale={scale:.3f}")

        series = [v * scale for v in _year_series(_fetch(donor, f"{YEAR}-01-01",
                                                         f"{YEAR}-12-31"), YEAR)]
        _write(name, series)


if __name__ == "__main__":
    sys.exit(main())
