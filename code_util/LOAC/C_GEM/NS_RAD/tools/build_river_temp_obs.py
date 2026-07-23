"""
Build observation-blended upstream river-temperature forcings for the gauged sites.

For Sagavanirktok and Kuparuk, USGS parameter 00010 gives OBSERVED daily river
temperature in 2022 (93 and 55 days, open-water only). This uses that observation as
the upstream boundary where it exists, and falls back to the modelled air-temp
regression (river_watertemp_2022_degC.csv) elsewhere -- winter, and any gaps outside
the observed span. A short linear taper at each end of the observed span avoids a
discontinuity where the two sources meet.

Why: the shared regression is +2.6 C warm on the Sagavanirktok (a mountain/
groundwater-fed river it was not fitted to); the observed record removes that bias
directly for the days it covers. Kuparuk, which the regression WAS fitted to, barely
changes -- a useful consistency check.

Colville and Canning have no 00010 record and keep the pure regression.

This series drives BOTH the upstream T boundary and the ice gate (via exfread's
rolling sum). The observed data is summer-only, so winter stays the regression and
the ice-gate timing -- set by the shoulder-season zero crossings -- is unchanged.
Verified in the smoke test below.

Usage:  python tools/build_river_temp_obs.py
"""
import datetime as dt
import urllib.request
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent.parent
FORC = ROOT / "forcing"
REGRESSION = FORC / "river_watertemp_2022_degC.csv"
GAUGE = {"sagavanirktok": "15908000", "kuparuk": "15896000"}
TAPER = 5   # days over which to blend observed -> regression at each span edge


def usgs_daily(gauge):
    u = (f"https://waterservices.usgs.gov/nwis/dv/?format=rdb&sites={gauge}"
         "&startDT=2022-01-01&endDT=2022-12-31&parameterCd=00010&statCd=00003")
    obs = {}
    for l in urllib.request.urlopen(u, timeout=60).read().decode().splitlines():
        c = l.split("\t")
        if len(c) < 4 or c[0] != "USGS":
            continue
        x = c[3].strip()
        if x in ("", "Ice", "Ssn", "Eqp"):
            continue
        d = dt.date(*map(int, c[2].split("-"))).timetuple().tm_yday - 1
        if d < 365:
            obs[d] = float(x)
    return obs


def write_csv(path, values):
    if len(values) != 365:
        raise SystemExit(f"{path.name}: {len(values)} values")
    path.write_text("\r\n".join(f"{v:.2f}" for v in values),
                    encoding="utf-8", newline="")


def main():
    reg = np.genfromtxt(open(REGRESSION, encoding="utf-8-sig"), delimiter=",")
    for site, gauge in GAUGE.items():
        obs = usgs_daily(gauge)
        days = np.array(sorted(obs))
        lo, hi = days.min(), days.max()
        # interpolate observed across its internal gaps, within [lo, hi] only
        obs_interp = np.interp(np.arange(lo, hi + 1), days,
                               [obs[d] for d in days])

        out = reg.copy()
        # blend weight: 1 inside the observed span, linear taper to 0 over TAPER
        # days beyond each edge (clamped to the span so we never invent obs)
        w = np.zeros(365)
        w[lo:hi + 1] = 1.0
        for k in range(1, TAPER + 1):
            f = 1.0 - k / (TAPER + 1)
            if lo - k >= 0:
                w[lo - k] = f
            if hi + k < 365:
                w[hi + k] = f
        # observed value to use in the taper zone = nearest span-edge observed value
        obs_full = out.copy()
        obs_full[lo:hi + 1] = obs_interp
        obs_full[:lo] = obs_interp[0]
        obs_full[hi + 1:] = obs_interp[-1]

        blended = w * obs_full + (1.0 - w) * out
        fname = FORC / f"{site}_watertemp_obs_2022_degC.csv"
        write_csv(fname, blended)

        # report the effect over the observed days
        bias_reg = np.mean([reg[d] - obs[d] for d in days])
        print(f"{site:15s} obs {len(obs)} d (day {lo}-{hi})  "
              f"regression bias was {bias_reg:+.2f} C -> now observed on those days  "
              f"-> {fname.name}")


if __name__ == "__main__":
    main()
