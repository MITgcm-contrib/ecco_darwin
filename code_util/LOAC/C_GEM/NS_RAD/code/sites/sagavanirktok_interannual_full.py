"""
Sagavanirktok, FULL-FORCING interannual variant: EVERY forcing category is genuinely
multi-year for the same 2005-2023 window (19 years), not just discharge/DOC like
`sagavanirktok_interannual.py`. See `colville_interannual_full.py` for the full
per-category source documentation (shared met/surge/solar/humidity sources, this
river's own ECCO-Darwin marine-boundary point) -- not repeated here to avoid drift
between the three rivers' variants. `sagavanirktok.py`, `sagavanirktok_interannual.py`,
runs/definitive, runs/interannual, and every existing report/validation tool are
completely unaffected by this file's existence.

    CGEM_SITE=sagavanirktok_interannual_full CGEM_MAXT_DAYS=6935 CGEM_WARMUP_DAYS=365 \\
        PYTHONPATH=code python code/main.py

or tools/run_interannual_full.sh sagavanirktok. Writes to
runs/interannual_full/sagavanirktok/ by convention (not enforced here -- see the run
wrapper).
"""
from .sagavanirktok import *  # noqa: F401,F403

DISCHARGE_FILE = "sagavanirktok_river_discharge_interannual_2005-2023_m3sec.csv"
WIND_FILE = "windspeed_hourly_interannual_2005-2023_msec.csv"
WIND_FREQ_SEC = 3600
AIRTEMP_FILE = "airtemp_hourly_interannual_2005-2023_degC.csv"
AIRTEMP_FREQ_SEC = 3600
RELHUM_FILE = "relhum_hourly_interannual_2005-2023_frac.csv"
RELHUM_FREQ_SEC = 3600
SOLAR_FILE = "solar_hourly_interannual_2005-2023_Wm2.csv"
SOLAR_FREQ_SEC = 3600
SURGE_FILE = "surge_prudhoe_interannual_2005-2023_m.csv"
SEATEMP_FILE = "sagavanirktok_t_marine_interannual_2005-2023.csv"
WATERTEMP_FILE = "river_temp_interannual_2005-2023_degC.csv"

BOUNDARY_FORCING = {
    "TOC": {"cub": "sagavanirktok_toc_interannual_2005-2023_mmolC_m3.csv",
            "clb": "sagavanirktok_toc_marine_interannual_2005-2023.csv"},
    "S":   {"clb": "sagavanirktok_s_marine_interannual_2005-2023.csv"},
    "DIC": {"clb": "sagavanirktok_dic_marine_interannual_2005-2023.csv"},
    "ALK": {"clb": "sagavanirktok_alk_marine_interannual_2005-2023.csv"},
    "NO3": {"clb": "sagavanirktok_no3_marine_interannual_2005-2023.csv"},
    "NH4": {"clb": "sagavanirktok_nh4_marine_interannual_2005-2023.csv"},
    "PO4": {"clb": "sagavanirktok_po4_marine_interannual_2005-2023.csv"},
    "dSi": {"clb": "sagavanirktok_dsi_marine_interannual_2005-2023.csv"},
    "O2":  {"clb": "sagavanirktok_o2_marine_interannual_2005-2023.csv"},
}

LABEL = "Sagavanirktok-interannual-full"
