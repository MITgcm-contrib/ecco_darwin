"""
IDEALIZED North Slope river -- a verification fixture, NOT a real river.

Purpose
-------
A generic, fully analytic estuary for exercising the model end to end with a KNOWN,
reproducible input over a full seasonal cycle:

  * round-number geometry (30 km, 2 m deep, a simple 400->100 m flare);
  * analytic, time-varying meteorology and discharge (a spring freshet large enough
    to trigger hydraulic ice break-up, an air temperature that crosses 0 C twice);
  * TIME-VARYING BOUNDARY CHEMISTRY -- riverine DOC flushed high on the freshet,
    alkalinity/DIC diluted by snowmelt, nitrate drawn down in summer, and a marine
    salinity that freshens under summer sea-ice melt.

The last of these is the point of the fixture: it drives BOTH boundary ends (upstream
`cub` and downstream `clb`) from forcing files via config.BOUNDARY_FORCING, so the
file-driven time-varying-boundary machinery is tested, not just the constant-boundary
path the four real rivers use.

All series are built by tools/build_idealized_forcings.py from a single PARAMS block,
so the fixture is reproducible. See docs/idealized_verification.md and run it with
tools/verify_idealized.py.

Nothing here is measured. IS_IDEALIZED = True makes config.py emit a note so idealized
output is never mistaken for a river result.
"""
from ._baseline import AMPL, pfun  # noqa: F401  (single-sinusoid tidal fallback)

LABEL = "Idealized"
IS_IDEALIZED = True
GEOMETRY_IS_PLACEHOLDER = False   # the idealization is deliberate, not a stand-in
ARCTIC_BGC = True                 # run the Arctic biogeochemistry extension (this fixture tests it)

# --- geometry: clean round numbers -------------------------------------------
EL = 30000        # 30 km domain  -> M = EL/DELXI+1 = 151 -> forced even = 150
DEPTH_lb = 2.0    # constant 2 m depth, both ends
DEPTH_ub = 2.0
B_lb = 400        # seaward width [m]
B_ub = 100        # prismatic width upstream of the flare [m]
L_FLARE = 8000    # width converges over 0..8 km, prismatic beyond
distance = 1

# --- tides: use the single-sinusoid fallback (AMPL/pfun from _baseline) -------
# The idealized site is deliberately NOT in forcing/tidal_constituents.json,
# so fun_module.Tide falls back to the microtidal single sinusoid. Keeps the fixture
# self-contained (no dependency on the harmonic-reconstruction JSON).

# --- analytic forcings (built by tools/build_idealized_forcings.py) ----------
DISCHARGE_FILE = "idealized_discharge_m3sec.csv"
WATERTEMP_FILE = "idealized_river_watertemp_degC.csv"   # upstream T boundary
SEATEMP_FILE = "idealized_seatemp_degC.csv"             # downstream T boundary
WIND_FILE = "idealized_wind_msec.csv"
SOLAR_FILE = "idealized_solar_Wm2.csv"
AIRTEMP_FILE = "idealized_airtemp_degC.csv"
RELHUM_FILE = "idealized_relhum_frac.csv"
PCO2_FILE = "idealized_pCO2_uatm.csv"

# --- boundary chemistry ------------------------------------------------------
# BOUNDARIES seeds the t=0 profile and supplies the CONSTANT value for any end not
# listed in BOUNDARY_FORCING below. Species/ends that ARE time-varying are set here
# to their day-0 value (main.py overwrites them every timestep anyway). Units follow
# the model: mmol m^-3 for solutes, g L^-1 for SPM, PSU for S, degC for T.
#                (clb = downstream/marine, cub = upstream/riverine)
BOUNDARIES = {
    "S":   (30.0, 0.0),        # clb varies (sea-ice melt); river is fresh
    "DIA": (1.0, 0.5),
    "dSi": (5.0, 60.0),
    "NO3": (0.5, 6.9),         # cub varies (winter-high)
    "NH4": (0.2, 1.0),
    "PO4": (0.6, 0.3),
    "O2":  (350.0, 400.0),
    "TOC": (60.0, 250.0),      # labile DOC; cub varies (freshet DOC flush)
    "RDOC": (40.0, 450.0),     # refractory + chromophoric DOC; cub varies (freshet flush)
    "CH4": (0.005, 0.6),       # dissolved methane [mmol m^-3]; river supersaturated
    "N2O": (0.012, 0.03),      # dissolved N2O [mmol m^-3]
    "SPM": (0.05, 1.0),
    "DIC": (2050.0, 880.0),    # cub varies (snowmelt dilution)
    "ALK": (2200.0, 900.0),    # cub varies (snowmelt dilution)
    "pH":  (8.1, 7.6),         # seed only; diagnosed from S/DIC/ALK
    "T":   (-1.8, 0.0),        # seed only; main.py drives both ends each timestep
}

# --- time-varying boundaries (the feature under test) ------------------------
# Each entry makes that species' clb and/or cub follow a 365-day forcing series,
# refreshed every timestep in main.py. Everything else stays at its BOUNDARIES value.
BOUNDARY_FORCING = {
    "S":   {"clb": "idealized_S_clb_marine.csv"},   # marine end (downstream)
    "TOC": {"cub": "idealized_TOC_cub_river.csv"},  # riverine end (upstream)
    "RDOC": {"cub": "idealized_RDOC_cub_river.csv"},  # refractory DOC, new-species BC test
    "NO3": {"cub": "idealized_NO3_cub_river.csv"},
    "ALK": {"cub": "idealized_ALK_cub_river.csv"},
    "DIC": {"cub": "idealized_DIC_cub_river.csv"},
}

# --- distributed lateral loading (the Arctic-extension feature under test) ----
# Tundra/thermokarst inputs entering ALONG the channel, not just at the upstream
# boundary. Total inflow spread over the domain; concentrations are the lateral
# water's chemistry (DOC/DIC/CH4-rich, as thermokarst drainage is). Injected as a
# mixing source in main.py. The four real rivers leave this off (unconstrained).
LATERAL_INFLOW = 8.0           # total lateral inflow [m^3 s^-1] (~15% of mean discharge)
LATERAL_CONC = {               # chemistry of the lateral (tundra/thermokarst) water
    "RDOC": 600.0, "TOC": 300.0, "DIC": 1100.0, "NO3": 5.0,
    "CH4": 2.0, "N2O": 0.04, "O2": 300.0,
}
