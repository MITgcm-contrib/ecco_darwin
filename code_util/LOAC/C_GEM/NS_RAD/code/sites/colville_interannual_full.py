"""
Colville, FULL-FORCING interannual variant: EVERY forcing category is genuinely
multi-year for the same 2005-2023 window (19 years), not just discharge/DOC like
`colville_interannual.py`. `import *` inherits colville.py's geometry/parameters
unchanged; only the forcing FILE paths below differ. `colville.py`,
`colville_interannual.py`, runs/definitive, runs/interannual, and every existing report/
validation tool are completely unaffected by this file's existence.

WHY 2005-2023, NOT 1980-2023. This is the SHORTEST COMMON WINDOW across every forcing
source that has a real, non-repeating multi-year record -- see CLAUDE.md -> "Interannual
forcing" -> "Full-forcing interannual (2005-2023)" for the full per-source availability
table and the reasoning (chosen explicitly over per-source mismatched windows or
switching meteorology to a reanalysis product). Discharge/DOC (otherwise 1980-2023) are
trimmed to match, not left at their own longer native range, so every variable in this
run is real for the same 19 years.

Per-category source (all built by tools/build_interannual_*.py, 2026-09-30):
  - Discharge + riverine TOC: PWBM, sliced to 2005-2023 from the same 1980-2023 series
    colville_interannual.py uses.
  - Wind + air temperature: NDBC PRDA2 (shared across all rivers, same station as
    colville.py's single-year files), HOURLY -- resolves the diurnal cycle, unlike
    every forcing in this project's history before this round. Both built from the
    SAME already-downloaded PRDA2 raw archive in one pass (no new source needed) by
    tools/build_interannual_met_hourly.py: wind was extended to hourly alongside air
    temperature because its raw source is exactly as available, at zero extra
    download cost -- the daily-only choice for these two was never a hard
    limitation, just an initial, narrower scoping choice. ~89-93% real observations.
  - Humidity: Deadhorse Airport NOAA ISD (shared, same station as colville.py's
    single-year file), HOURLY -- same reasoning as wind: the raw ISD record is
    already hourly, so build_interannual_humidity_hourly.py extracts it directly
    from the already-cached files build_interannual_humidity.py downloaded. ~95%
    real observations.
  - Storm surge: NOAA CO-OPS 9497645, Prudhoe Bay (shared regional proxy, same
    station as colville.py's single-year surge file), DAILY -- deliberately NOT
    extended to hourly. NOAA's raw `hourly_height` product is tide-INCLUDED water
    level; the daily mean is specifically what cancels the ~12 h tidal oscillation
    out, leaving the true non-tidal surge residual. Reading it hourly without first
    subtracting the harmonic tide would double-count the tide against the model's
    own harmonic `Tide()` formula -- a real methodological step, not a free
    upgrade, and storm surges themselves evolve over many hours to days anyway
    (not a diurnal process), so this was left at daily resolution. ~99.4% real
    observations.
  - Solar radiation: NOAA GML Barrow (Utqiagvik) Atmospheric Baseline Observatory
    (shared), HOURLY real 1-minute pyranometer observations aggregated up -- NOT
    the same source as colville.py's single-year `solarradiation.csv` (whose
    original derivation could not be identified or reproduced -- see
    tools/build_interannual_met.py's docstring), and NOT a reanalysis product
    (NSRDB doesn't reach this latitude in its standard product and its polar
    product only starts in 2013; NOAA NARR reanalysis was tried first and works
    but is a MODELED estimate -- see tools/build_interannual_solar_subdaily.py,
    superseded, kept for reference). Barrow is ~330 km from Prudhoe Bay -- farther
    from all three rivers than Prudhoe Bay itself -- but is used anyway because it
    is the only nearby site that actually measures solar radiation at all (PRDA2
    has no radiation instrument); wind/air-temp/humidity stay sourced from the
    closer Prudhoe Bay stations. ~98% real observations. See
    tools/build_interannual_solar_barrow.py and CLAUDE.md -> "Interannual forcing"
    -> "Sub-daily (diurnal) forcing" for the full reasoning trail, including why
    hourly (not true 1-minute) was chosen as this project's diurnal resolution.
  - Marine boundary (S, sea temperature, DIC, ALK, NO3, NH4, PO4, dSi, O2, TOC) AND sea
    temperature (replacing PRDA2's `watertemp.csv` -- see
    tools/build_interannual_met.py's docstring for why): ECCO-Darwin v5 time-extended
    output, this river's own nearest-wet-cell point, monthly native resolution
    interpolated to daily. 100% real months, zero gaps, at all three rivers.
  - Riverine temperature: no independent source (never an independent observation
    even in the single-year version -- see build_river_temp.py); the same regional
    air->water formula applied to the new multi-year air temperature instead of the
    single 2022 year (tools/build_interannual_river_temp.py).
  - NOT extended: ice-model forcing has no separate source (it consumes air temp/wind,
    both already covered, plus discharge, also already covered) -- nothing left to do
    there. Nothing else in this model reads a forcing file at all.

    CGEM_SITE=colville_interannual_full CGEM_MAXT_DAYS=6935 CGEM_WARMUP_DAYS=365 \\
        PYTHONPATH=code python code/main.py

or tools/run_interannual_full.sh colville. Writes to runs/interannual_full/colville/ by
convention (not enforced here -- see the run wrapper).
"""
from .colville import *  # noqa: F401,F403

DISCHARGE_FILE = "colville_river_discharge_interannual_2005-2023_m3sec.csv"
WIND_FILE = "windspeed_hourly_interannual_2005-2023_msec.csv"
WIND_FREQ_SEC = 3600
AIRTEMP_FILE = "airtemp_hourly_interannual_2005-2023_degC.csv"
AIRTEMP_FREQ_SEC = 3600
RELHUM_FILE = "relhum_hourly_interannual_2005-2023_frac.csv"
RELHUM_FREQ_SEC = 3600
SOLAR_FILE = "solar_hourly_interannual_2005-2023_Wm2.csv"
SOLAR_FREQ_SEC = 3600
SURGE_FILE = "surge_prudhoe_interannual_2005-2023_m.csv"
SEATEMP_FILE = "colville_t_marine_interannual_2005-2023.csv"
WATERTEMP_FILE = "river_temp_interannual_2005-2023_degC.csv"

BOUNDARY_FORCING = {
    "TOC": {"cub": "colville_toc_interannual_2005-2023_mmolC_m3.csv",
            "clb": "colville_toc_marine_interannual_2005-2023.csv"},
    "S":   {"clb": "colville_s_marine_interannual_2005-2023.csv"},
    "DIC": {"clb": "colville_dic_marine_interannual_2005-2023.csv"},
    "ALK": {"clb": "colville_alk_marine_interannual_2005-2023.csv"},
    "NO3": {"clb": "colville_no3_marine_interannual_2005-2023.csv"},
    "NH4": {"clb": "colville_nh4_marine_interannual_2005-2023.csv"},
    "PO4": {"clb": "colville_po4_marine_interannual_2005-2023.csv"},
    "dSi": {"clb": "colville_dsi_marine_interannual_2005-2023.csv"},
    "O2":  {"clb": "colville_o2_marine_interannual_2005-2023.csv"},
}

LABEL = "Colville-interannual-full"
