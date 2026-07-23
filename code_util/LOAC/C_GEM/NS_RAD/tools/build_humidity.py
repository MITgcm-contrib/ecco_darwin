"""
Build a daily relative-humidity forcing from Deadhorse Airport (NOAA ISD).

WHY NOT ERA5. The heat budget needs humidity at Prudhoe Bay to close its latent-heat
term. ERA5-Land would supply it but requires Copernicus CDS credentials. Deadhorse
Airport (PASC, USAF 700637 / WBAN 27406) sits at 70.19 N, -148.48 E -- essentially
colocated with the PRDA2 buoy that already supplies wind/air/water temperature -- and
its hourly dewpoint is in NOAA's freely accessible Integrated Surface Database. For a
single point that is a better humidity source than a reanalysis grid cell.

It replaces heat_module's assumed constant REL_HUMIDITY (0.85). Note the observed
open-water mean is ~0.87, so the constant was close; see docs/performance-of-heat... no,
see the heat-budget validation. Relative humidity from hourly T and Td via the Magnus
saturation ratio, daily-averaged, gaps filled by interpolation onto 365 days.

Usage:  python tools/build_humidity.py
"""
import csv
import datetime as dt
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent.parent
FORC = ROOT / "forcing"
ISD_URL = "https://www.ncei.noaa.gov/data/global-hourly/access/2022/70063727406.csv"


def esat(T):
    """Saturation vapour pressure over water [hPa], Magnus form, T in degC."""
    return 6.112 * np.exp(17.67 * T / (T + 243.5))


def main():
    import urllib.request
    local = FORC / "deadhorse_isd_2022.csv"
    if not local.exists():
        print(f"downloading Deadhorse ISD -> {local.name}")
        urllib.request.urlretrieve(ISD_URL, local)

    recs = {}
    with open(local) as f:
        for row in csv.DictReader(f):
            def parse(field):
                v = row[field].split(",")[0]
                return np.nan if v in ("+9999", "9999", "") else int(v) / 10.0
            T, Td = parse("TMP"), parse("DEW")
            if np.isnan(T) or np.isnan(Td) or T < -50:
                continue
            doy = dt.datetime.strptime(row["DATE"], "%Y-%m-%dT%H:%M:%S").timetuple().tm_yday - 1
            if doy < 365:
                recs.setdefault(doy, []).append(esat(Td) / esat(T))

    rh = np.array([np.clip(np.mean(recs[d]), 0.1, 1.0) if d in recs else np.nan
                   for d in range(365)])
    ok = ~np.isnan(rh)
    rh[~ok] = np.interp(np.flatnonzero(~ok), np.flatnonzero(ok), rh[ok])

    out = FORC / "relhum_2022_frac.csv"
    out.write_text("\r\n".join(f"{v:.4f}" for v in rh), encoding="utf-8", newline="")
    print(f"wrote {out.name}: RH {100*rh.min():.0f}-{100*rh.max():.0f}%, "
          f"mean {100*rh.mean():.0f}%, open-water mean {100*rh[140:304].mean():.0f}%")
    print(f"  (heat_module's assumed constant REL_HUMIDITY was 0.85)")


if __name__ == "__main__":
    main()
