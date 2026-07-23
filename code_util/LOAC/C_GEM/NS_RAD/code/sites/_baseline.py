"""
Shared starting values for all four North Slope sites.

IMPORTANT -- read before trusting any run.

The geometry below is NOT a per-river measurement. It is the single hybrid
configuration the model shipped with, which mixed sources: EL / B_lb / B_ub came
from SWORD delta segments, DEPTH_lb / DEPTH_ub were cited to the Sagavanirktok
gauge (USGS 15908000), and the discharge forcing was Colville. No single river is
described by this combination.

It is reproduced here so that all four sites are immediately runnable and the
plumbing (forcings, output separation, ice gating) can be exercised end to end.
Each site module is expected to OVERRIDE these with values for its own river
before any run is treated as scientifically meaningful. Sites that have not yet
been given real geometry declare GEOMETRY_IS_PLACEHOLDER = True, and config.py
emits a warning naming the site when one is used.

BOUNDARIES maps species -> (clb, cub):
    clb = downstream / marine boundary concentration
    cub = upstream   / riverine boundary concentration
Units follow the model: mmol m^-3 for solutes, g L^-1 for SPM, PSU for S,
pH units for pH. These are likewise the original single-site values and want
per-river replacement -- riverine DOC, alkalinity and nutrients differ a lot
across the North Slope.
"""

# --- channel geometry (PLACEHOLDER -- see module docstring) -------------------
EL = 27175        # Estuarine length [m]
DEPTH_lb = 15     # Depth at downstream boundary [m]
DEPTH_ub = 15     # Depth at upstream boundary [m]
B_lb = 1215       # Width at downstream boundary [m]
B_ub = 859        # Width at upstream boundary [m]
distance = 1      # Grid points in the saline zone

# --- water temperature -------------------------------------------------------
# Modelled river temperature, shared across all four sites. The relation is
# REGIONAL (fitted on Kuparuk, the only gauge whose temperature air temperature
# actually explains) and its only input -- PRDA2 air temperature -- is regional, so
# there is currently no basis for a per-river series. Overridable per site if a
# site-specific record ever becomes available. Built by tools/build_river_temp.py.
#
# This replaces the shipped watertemp.csv, which was PRDA2 SEA-water temperature and
# biased every rate cold by 6-8 C while driving the ice gate off coastal sea ice.
WATERTEMP_FILE = "river_watertemp_2022_degC.csv"

# --- shared meteorological forcings (per-site OVERRIDABLE) --------------------
# These are the regional records every real North Slope site uses in common (one
# wind/solar/air/humidity/pCO2 series for the whole coast, plus the PRDA2 sea-water
# temperature that forms the DOWNSTREAM temperature boundary). They live here as
# defaults so a site can override any of them -- the idealized verification site
# (sites/idealized.py) points all of them at its own analytic series. config.py
# reads each with getattr(_site, ..., <this default>), so the four real rivers,
# which do not set them, are completely unaffected. Resolved to absolute paths in
# main.py against forcing/.
WIND_FILE = "windspeed.csv"
SOLAR_FILE = "solarradiation.csv"
AIRTEMP_FILE = "airtemp_2022_degC.csv"
RELHUM_FILE = "relhum_2022_frac.csv"
PCO2_FILE = "pCO2_Barrow_2022.csv"
SEATEMP_FILE = "watertemp.csv"    # PRDA2 sea-water temp = downstream T boundary

# --- time-varying boundary conditions (per-site, OPTIONAL) -------------------
# Empty here: the four real rivers hold every solute boundary CONSTANT in time
# (only temperature, handled separately in main.py, is seasonal). A site may declare
#
#     BOUNDARY_FORCING = {"TOC": {"cub": "some_river_series.csv"},
#                         "S":   {"clb": "some_marine_series.csv"}, ...}
#
# to make any species' downstream (clb) and/or upstream (cub) boundary follow a
# 365-day forcing series, refreshed every timestep exactly as temperature already is
# (openbound reads v[s]["clb"/"cub"] fresh each call, so updating them per step is all
# that is needed). The BOUNDARIES table above still seeds t=0 and supplies the
# constant value for any end NOT listed here. Used by the idealized site.
BOUNDARY_FORCING = {}

# --- tides -------------------------------------------------------------------
# Tides are now a per-river multi-constituent harmonic reconstruction from the
# nearest NOAA CO-OPS station (forcing/tidal_constituents.json, built by
# tools/build_tides.py; loaded in config, summed in fun_module.Tide). AMPL/pfun below
# are only the FALLBACK single-sinusoid, used if no constituent data is found for the
# site. The Beaufort coast is microtidal -- observed short-period range ~0.3-0.4 m,
# M2 amplitude ~0.06-0.07 m -- so the shipped 0.2 m was a reasonable amplitude guess;
# the constituent data mainly adds spring-neap and diurnal-inequality structure.
AMPL = 0.2        # fallback tidal amplitude at mouth [m]
pfun = 0.0787     # fallback tidal frequency [cycle/hr]

# --- boundary chemistry (PLACEHOLDER -- see module docstring) ----------------
BOUNDARIES = {
    "S":   (27.1, 0.0),
    "DIA": (1.6, 1.1),
    "dSi": (2.35, 174.6),
    "NO3": (0.03, 7.68),
    "NH4": (0.01, 3.9),
    "PO4": (0.56, 0.01),
    "O2":  (392.0, 367.0),
    "TOC": (69.0, 1582.0),      # labile / bioavailable DOC (existing pool)
    # Arctic extension pools/gases (placeholder Arctic values; refine per-river with data):
    "RDOC": (40.0, 800.0),      # refractory + chromophoric DOC -- dominant in Arctic rivers
    "CH4": (0.003, 0.5),        # dissolved methane [mmol m^-3]; rivers are supersaturated
    "N2O": (0.012, 0.02),       # dissolved N2O [mmol m^-3]; rivers near-to-super-saturated
    "SPM": (0.15, 2.0),
    "DIC": (1869.0, 1717.0),
    "ALK": (1957.0, 1596.0),
    "pH":  (8.1, 7.5),
    # Water temperature [degC]. These two values seed the initial profile ONLY --
    # main.py overwrites v['T']['clb'] and ['cub'] every timestep from the sea and
    # river temperature forcings, because unlike salinity the thermal boundaries are
    # strongly seasonal (-2 to +18 degC). The values here are representative of
    # t = 0 (1 January): sea at its freezing point, river at 0.
    "T":   (-1.8, 0.0),
}
