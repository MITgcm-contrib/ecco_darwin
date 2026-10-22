"""
Builds CLIMATOLOGICAL forcing files (day-of-year, or hour-of-year, mean across all
19 real years, 2005-2023) for a "FORTE campaign year" site variant -- a single
representative repeating year standing in for whichever real year the actual FORTE
Expedition field campaign (2027, 2028, ...) occurs in, built FROM the already-
validated 19-year full-forcing interannual record (tools/build_interannual_*.py)
rather than from any one specific historical year.

WHY A CLIMATOLOGY, NOT ONE ANALOG YEAR (project owner's explicit choice,
2026-10-19): runs/osse_nature (2007 low-flow / 2008 average / 2009 high-flow) and
runs/interannual_full (the full 19-year daily record) both stay exactly as they are,
for their own existing purposes (the OSSE sampling/ranking report, the publication
figures). This is a SEPARATE, NEW framework purpose-built for the FORTE campaign's
real monitored window (day-of-year 152-273: Phase 1a's Jun 1 start through the
mooring-recovery end of September -- see tools/build_osse_observations.py's PHASES
and MOORING_WINDOW) -- a smooth, multi-year-averaged "typical campaign season"
forcing is the right baseline for that: reusable across every future campaign year
without re-choosing an analog year each time, and its whole purpose is to be the
thing REAL field-campaign data (especially time-varying upstream `cub` chemistry,
which does not exist as real data for ANY historical year, synthetic or otherwise)
gets compared against, or eventually substituted into, once collected.

OUTPUT: one 365-ROW (daily forcings) or 8760-ROW (hourly forcings) CSV per
variable, named <same stem>_climatology_<same unit suffix> -- feeds
sites/<river>_forte_campaign.py, a thin site variant pointing at these instead of
the 19-year series. A 365-day/8760-hour file is exactly what every ORIGINAL
single-year site (colville.py etc.) already expects -- repeatYear wraps it
automatically, the same mechanism every non-interannual site already uses, so no
new model-code changes are needed, only new forcing files + a new site variant.

Usage:  python tools/build_forte_campaign_forcing.py
"""
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent.parent
FORC = ROOT / "forcing"
SITES = ["colville", "kuparuk", "sagavanirktok"]
MARINE_SPECIES = ["s", "dic", "alk", "no3", "nh4", "po4", "dsi", "o2", "toc"]


def _climatology_daily(vals):
    """Day-of-year mean (365 values) across an N*365-row daily series (no leap
    days, matching every interannual forcing file's own convention -- verified
    directly: every *_interannual_2005-2023_*.csv is exactly 19*365 rows)."""
    n_years = vals.size // 365
    return vals[:n_years * 365].reshape(n_years, 365).mean(axis=0)


def _climatology_hourly(vals):
    """Hour-of-year mean (8760 values) across an N*8760-row hourly series
    (verified: every *_hourly_interannual_2005-2023_*.csv is exactly 19*8760)."""
    n_years = vals.size // 8760
    return vals[:n_years * 8760].reshape(n_years, 8760).mean(axis=0)


def _build_one(src_name, dst_name, hourly):
    src = FORC / src_name
    if not src.exists():
        print(f"  SKIP (not found): {src_name}")
        return
    vals = np.genfromtxt(open(src, encoding="utf-8-sig"), delimiter=",")
    clim = _climatology_hourly(vals) if hourly else _climatology_daily(vals)
    np.savetxt(FORC / dst_name, clim, delimiter=",", fmt="%.6f")
    print(f"  wrote forcing/{dst_name} ({clim.size} rows, "
          f"{'hourly' if hourly else 'daily'})")


def main():
    print("Shared (met/tide) climatologies:")
    _build_one("windspeed_hourly_interannual_2005-2023_msec.csv",
                "windspeed_hourly_climatology_msec.csv", hourly=True)
    _build_one("airtemp_hourly_interannual_2005-2023_degC.csv",
                "airtemp_hourly_climatology_degC.csv", hourly=True)
    _build_one("relhum_hourly_interannual_2005-2023_frac.csv",
                "relhum_hourly_climatology_frac.csv", hourly=True)
    _build_one("solar_hourly_interannual_2005-2023_Wm2.csv",
                "solar_hourly_climatology_Wm2.csv", hourly=True)
    _build_one("surge_prudhoe_interannual_2005-2023_m.csv",
                "surge_prudhoe_climatology_m.csv", hourly=False)
    _build_one("river_temp_interannual_2005-2023_degC.csv",
                "river_temp_climatology_degC.csv", hourly=False)

    for site in SITES:
        print(f"{site}:")
        _build_one(f"{site}_river_discharge_interannual_2005-2023_m3sec.csv",
                    f"{site}_river_discharge_climatology_m3sec.csv", hourly=False)
        _build_one(f"{site}_toc_interannual_2005-2023_mmolC_m3.csv",
                    f"{site}_toc_climatology_mmolC_m3.csv", hourly=False)
        _build_one(f"{site}_t_marine_interannual_2005-2023.csv",
                    f"{site}_t_marine_climatology.csv", hourly=False)
        for sp in MARINE_SPECIES:
            _build_one(f"{site}_{sp}_marine_interannual_2005-2023.csv",
                        f"{site}_{sp}_marine_climatology.csv", hourly=False)


if __name__ == "__main__":
    main()
