"""
Build a wind-driven storm-surge forcing for the marine boundary from observed NOAA
CO-OPS water level, as the DAILY-MEAN sea level (the tide averages out over a day, so
what remains is the surge + seasonal sea-level signal). Written as a 365-value daily
series, mean-removed, the same format the other forcings use.

Only Prudhoe Bay (9497645) has continuous 2022 water level on this coast, so it is the
surge record for Kuparuk (its own tidal station) and the regional proxy for the others.
The Beaufort coast is microtidal (~0.3 m tide) but wind surges reach 1-3 m in big storm
years; 2022 was moderate (peak daily-mean +0.75 m). This is the physical saltwater-
intrusion driver the harmonic tide alone omits.

    surge(t) = daily-mean observed water level - annual mean

Output: forcing/surge_<proxy>_2022_m.csv  (365 daily values)
Wire in via sites/<name>.py: SURGE_FILE = "surge_prudhoe_2022_m.csv".
"""
import json
import urllib.request
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent.parent
OUT = ROOT / "forcing" / "surge_prudhoe_2022_m.csv"
STATION = "9497645"   # Prudhoe Bay


def main():
    url = ("https://api.tidesandcurrents.noaa.gov/api/prod/datagetter?"
           "product=hourly_height&datum=MSL&station=" + STATION +
           "&begin_date=20220101&end_date=20221231&time_zone=gmt&units=metric&format=json")
    with urllib.request.urlopen(url, timeout=90) as r:
        data = json.load(r)["data"]
    # hourly water level; keep gaps as NaN then fill
    wl = np.array([float(x["v"]) if x["v"] not in ("", "-") else np.nan for x in data])
    # pad/truncate to whole days
    n = (len(wl) // 24) * 24
    daily = np.nanmean(wl[:n].reshape(-1, 24), axis=1)      # tide averages out -> surge+seasonal
    # to exactly 365 days
    if len(daily) >= 365:
        daily = daily[:365]
    else:
        daily = np.pad(daily, (0, 365 - len(daily)), mode="edge")
    daily = np.where(np.isfinite(daily), daily, np.nanmean(daily))
    surge = daily - np.mean(daily)                          # mean-removed fluctuation
    np.savetxt(OUT, surge, fmt="%.4f")
    print(f"wrote {OUT.relative_to(ROOT)}  (365 daily values)")
    print(f"  surge range {surge.min():+.2f} to {surge.max():+.2f} m,  "
          f"days>+0.4 m: {(surge > 0.4).sum()}")


if __name__ == "__main__":
    main()
