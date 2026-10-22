"""
FORTE campaign year: a thin site variant pointing at the CLIMATOLOGICAL forcing
files (tools/build_forte_campaign_forcing.py -- day-of-year/hour-of-year mean
across all 19 real years, 2005-2023) instead of kuparuk_interannual_full.py's
19-year real series. `import *` inherits geometry/parameters/FREQ_SEC settings
unchanged; only the forcing FILE paths below differ.

WHY A CLIMATOLOGY, NOT ONE ANALOG YEAR: this framework exists specifically for the
FORTE Expedition field campaign's real monitored window (day-of-year 152-273,
Phase 1a's Jun 1 start through the mooring-recovery end of September -- see
tools/build_osse_observations.py's PHASES/MOORING_WINDOW), as a reusable baseline
for every future campaign year (2027, 2028, ...) -- and, eventually, the thing REAL
field-collected upstream (`cub`) chemistry gets compared against or substituted
into, once collected (no historical year has real cub data to be an "analog" for
in the first place). See CLAUDE.md -> "Interannual forcing" -> "FORTE campaign
framework" and tools/build_forte_campaign_forcing.py's docstring for the full
reasoning. runs/osse_nature/ (2007/2008/2009 analog years) and
runs/interannual_full/ (the real 19-year record) are both completely unaffected by
this file's existence -- kept exactly as-is for their own existing purposes.

    CGEM_SITE=kuparuk_forte_campaign CGEM_MAXT_DAYS=638 CGEM_WARMUP_DAYS=365 \
        CGEM_TS=12 PYTHONPATH=code python code/main.py

or tools/run_forte_campaign.sh kuparuk. 638 = 365 (warmup, one climatological year)
+ 273 (day-of-year the mooring recovery/campaign season ends, late September) of a
second climatological year. Writes to runs/forte_campaign/kuparuk/ by convention
(not enforced here -- see the run wrapper).
"""
from .kuparuk_interannual_full import *  # noqa: F401,F403

DISCHARGE_FILE = "kuparuk_river_discharge_climatology_m3sec.csv"
WIND_FILE = "windspeed_hourly_climatology_msec.csv"
AIRTEMP_FILE = "airtemp_hourly_climatology_degC.csv"
RELHUM_FILE = "relhum_hourly_climatology_frac.csv"
SOLAR_FILE = "solar_hourly_climatology_Wm2.csv"
SURGE_FILE = "surge_prudhoe_climatology_m.csv"
SEATEMP_FILE = "kuparuk_t_marine_climatology.csv"
WATERTEMP_FILE = "river_temp_climatology_degC.csv"

BOUNDARY_FORCING = {
    "TOC": {"cub": "kuparuk_toc_climatology_mmolC_m3.csv",
            "clb": "kuparuk_toc_marine_climatology.csv"},
    "S":   {"clb": "kuparuk_s_marine_climatology.csv"},
    "DIC":   {"clb": "kuparuk_dic_marine_climatology.csv"},
    "ALK":   {"clb": "kuparuk_alk_marine_climatology.csv"},
    "NO3":   {"clb": "kuparuk_no3_marine_climatology.csv"},
    "NH4":   {"clb": "kuparuk_nh4_marine_climatology.csv"},
    "PO4":   {"clb": "kuparuk_po4_marine_climatology.csv"},
    "dSi":   {"clb": "kuparuk_dsi_marine_climatology.csv"},
    "O2":   {"clb": "kuparuk_o2_marine_climatology.csv"},
}

LABEL = "Kuparuk-forte-campaign"
