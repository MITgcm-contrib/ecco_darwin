"""
Central configuration and single source of all model parameters (translated from the
original C-GEM headers, then heavily extended for the North Slope multi-site setup).

Structure:
  * SITE SELECTION -- picks a river from CGEM_SITE and imports sites/<name>.py, from which
    only the per-site quantities (geometry, tides, discharge/temperature filenames, the
    BOUNDARIES table, optional surge) are pulled; everything below is SHARED across all
    four rivers, so a change here applies to every site at once.
  * Then, in order: run controls (MAXT/WARMUP/DELTI/TS/output format, all env-overridable),
    grid (M from EL/DELXI), geometry & width/dispersion model, hydrodynamic & sediment
    parameters, the biogeochemical reaction-network constants, the ice-model block, and
    the tidal-constituent load.

Env-var overrides (defaults in parentheses): CGEM_SITE (colville), CGEM_MAXT_DAYS (730),
CGEM_WARMUP_DAYS (365), CGEM_TS (12), CGEM_OUTPUT (nc), CGEM_ICE (on). See CLAUDE.md.
"""
import math
import os
import importlib
import warnings

# SITE SELECTION
# ---------------------------------------------------------------------------
# Which river this run represents. Set via the CGEM_SITE environment variable:
#     CGEM_SITE=kuparuk python main.py
# Valid keys are listed in sites/__init__.py. Defaults to colville, which is the
# configuration the model originally shipped with.
#
# Everything below the site block is SHARED across all four rivers -- edit it once
# and every site sees the change. Only geometry, tides, the discharge filename and
# the boundary concentrations are per-site; those live in sites/<name>.py.
SITE = os.environ.get("CGEM_SITE", "colville")

try:
    _site = importlib.import_module(f"sites.{SITE}")
except ModuleNotFoundError as exc:
    from sites import SITES
    raise SystemExit(
        f"Unknown CGEM_SITE={SITE!r}. Valid sites: {', '.join(SITES)}"
    ) from exc

SITE_LABEL = _site.LABEL
DISCHARGE_FILE = _site.DISCHARGE_FILE
WATERTEMP_FILE = _site.WATERTEMP_FILE
BOUNDARIES = _site.BOUNDARIES

# Meteorological forcing filenames. Shared/regional for the four real rivers (the
# defaults below match what main.py used to hardcode, so those runs are unchanged);
# a site -- e.g. the idealized verification site -- may override any of them in
# sites/<name>.py. getattr keeps sites that don't set them working as before.
WIND_FILE = getattr(_site, "WIND_FILE", "windspeed.csv")
SOLAR_FILE = getattr(_site, "SOLAR_FILE", "solarradiation.csv")
AIRTEMP_FILE = getattr(_site, "AIRTEMP_FILE", "airtemp_2022_degC.csv")
RELHUM_FILE = getattr(_site, "RELHUM_FILE", "relhum_2022_frac.csv")
PCO2_FILE = getattr(_site, "PCO2_FILE", "pCO2_Barrow_2022.csv")
SEATEMP_FILE = getattr(_site, "SEATEMP_FILE", "watertemp.csv")
# Optional wind-driven storm-surge forcing at the marine boundary (daily sea-level
# residual, m; tools/build_surge.py). None = harmonic tide only. Added to Tide() in
# fun_module each timestep. See sites/<name>.py to enable per river.
SURGE_FILE = getattr(_site, "SURGE_FILE", None)

# Optional per-species time-varying boundaries, {species: {"clb"|"cub": filename}}.
# Empty for the real rivers (constant solute boundaries); main.py refreshes any
# entries here each timestep. See sites/_baseline.py for the contract.
BOUNDARY_FORCING = getattr(_site, "BOUNDARY_FORCING", {})

# Informational flag for the idealized verification fixture, so its NetCDF output is
# never mistaken for a real river. Not a placeholder warning -- the idealization is
# deliberate. config emits a one-line note below.
IS_IDEALIZED = getattr(_site, "IS_IDEALIZED", False)

if getattr(_site, "GEOMETRY_IS_PLACEHOLDER", False):
    warnings.warn(
        f"[{SITE_LABEL}] channel geometry and boundary chemistry are still the "
        "shipped placeholder values, not values for this river -- see "
        "sites/_baseline.py. Plumbing is valid; results are not.",
        stacklevel=2,
    )
if getattr(_site, "DISCHARGE_IS_RECONSTRUCTED", False):
    warnings.warn(
        f"[{SITE_LABEL}] discharge is RECONSTRUCTED from a donor gauge, not "
        f"observed (no 2022 record exists) -- see sites/{SITE}.py.",
        stacklevel=2,
    )
if IS_IDEALIZED:
    warnings.warn(
        f"[{SITE_LABEL}] this is the IDEALIZED verification fixture -- analytic, "
        "time-varying forcings and boundaries, not a real river. Output is for "
        "testing the model machinery, not for interpretation. See "
        "sites/idealized.py and docs/idealized_verification.md.",
        stacklevel=2,
    )

# GEOMETRICAL AND PHYSICAL PARAMETERS  (per-site, from sites/<name>.py)
EL = _site.EL              # Estuarine length [m]
DEPTH_lb = _site.DEPTH_lb  # Depth at downstream boundary [m]
DEPTH_ub = _site.DEPTH_ub  # Depth at upstream boundary [m]
B_lb = _site.B_lb          # Width at seaward boundary [m]
B_ub = _site.B_ub          # PRISMATIC width upstream of the flare [m]
L_FLARE = _site.L_FLARE    # Flare length [m] -- see WIDTH_MODEL below
distance = int(os.environ.get("CGEM_DISTANCE", _site.distance))  # Grid points in saline
# zone (downstream plateau of the Chezy/sediment ramps). Env-overridable for sensitivity
# sweeps -- see the distance-vs-salinity test in docs. Default per site (1 for all four).

# WIDTH MODEL
# -----------
# 'flare'  : width converges exponentially over 0..L_FLARE, then is PRISMATIC.
# 'expo'   : the original C-GEM law, exponential convergence over the whole domain.
#
# SWORD v17c node widths say 'flare' is the right shape for these rivers. Formal
# model comparison on 1 km binned log-width over the 27 km domain (AIC, lower better):
#
#     river           exponential      constant      flare+prismatic   best
#     colville        R2 0.06 / -28.3     -29.0      R2 0.39 / -35.4   flare
#     kuparuk         R2 0.23 / -74.2     -68.9      R2 0.43 / -80.7   flare
#     canning         R2 0.15 / -65.6     -63.0      R2 0.35 / -70.9   flare
#     sagavanirktok   R2 0.46 / -26.7     -23.2      (9 bins, gap)     expo*
#
# * the Sagavanirktok has a 5-20 km data gap and only 9 usable bins, so its
#   apparent preference for the exponential is not trustworthy; it uses borrowed
#   flare parameters instead (see sites/sagavanirktok.py).
#
# Physically: these are rivers with small deltas, not tide-dominated funnel
# estuaries (observed microtidal range ~0.3-0.4 m, distance = 1, negligible saline
# intrusion). The exponential convergence law comes from Savenije's alluvial-estuary
# theory and is being applied outside its regime. Set 'expo' to recover the original.
WIDTH_MODEL = "flare"

RS = 1.0  # Storage width ratio
rho_w = 1000.0  # Density of pure water [kg/m^3]
G = 9.81  # Gravity acceleration [m/s^2]

# HYDRODYNAMIC AND SEDIMENT PARAMETERS
Chezy_lb = 60  # Chezy coefficient (downstream) [m^-1/2 s^-1] Panchenko and Alabyan, 2022 https://doi.org/10.1016/j.mex.2022.101669
Chezy_ub = 40  # Chezy coefficient (upstream) [m^-1/2 s^-1] Panchenko and Alabyan, 2022 https://doi.org/10.1016/j.mex.2022.101669
#ws = 1.0e-3  # Settling velocity for sediment [m/s]
kISS = 51.0e-3  # Inorganic sinking velocity half saturation concentration [g L^-1] - Clark et al 2022 https://doi.org/10.1029/2022JG007139
wMAX = 2.0 / 86400  #	Maximum inorganic suspended sediment sinking velocity [m s^−1] - Clark et al 2022 https://doi.org/10.1029/2022JG007139
tau_ero_lb = 0.005  # Erosion shear stress (downstream) [N/m^2] - Clark et al 2022 https://doi.org/10.1029/2022JG007139
tau_dep_lb = 0.4  # Deposition shear stress (downstream) [N/m^2]
tau_ero_ub = 0.005  # Erosion shear stress (upstream) [N/m^2] - Clark et al 2022 https://doi.org/10.1029/2022JG007139
tau_dep_ub = 1.0  # Deposition shear stress (upstream) [N/m^2]
Mero_lb = 5.787037037037037e-05  # Erosion coefficient (downstream) [mg/m^2 s] - Clark et al 2020 https://doi.org/10.1029/2019JG005442
Mero_ub = 5.787037037037037e-05  # Erosion coefficient (upstream) [mg/m^2 s] - Clark et al 2020 https://doi.org/10.1029/2019JG005442

# BIOGEOCHEMICAL PARAMETERS
Pbmax = 1.3888888888888888e-05  # Max photosynthetic rate [s^-1] - Le Fouest 2013 www.biogeosciences.net/10/4785/2013/
Chla2CMIN = 0.0125  # Minimum Chla : C ratio for phytoplankton [g/g] - Le Fouest 2013 www.biogeosciences.net/10/4785/2013/
KE = 8.0  # Photoacclimation parameter [E/m^2/d^1] - Le Fouest 2013 www.biogeosciences.net/10/4785/2013/
alpha = 2.0  # Photosynthetic efficiency [mg C/(mg Chl)/(E/m^2/d^1)] - Le Fouest 2013 www.biogeosciences.net/10/4785/2013/
KdSi = 1.07  # Silica constant [muM Si]
KPO4 = 0.2  # Phosphate constant [muM P]
KNH4 = 0.5  # Ammonium constant [muM N] - Le Fouest 2013 www.biogeosciences.net/10/4785/2013/
KNO3 = 1  # Nitrate constant [muM N] - Le Fouest 2013 www.biogeosciences.net/10/4785/2013/
KTOC = 186.25  # Organic matter constant [muM C]
KO2_ox = 31.0  # Oxygen constant [muM O2]
KO2_nit = 51.25  # Oxygen constant [muM O2]
KN = 1.13  # Dissolved nitrogen constant [muM N]
KinO2 = 33.0  # Denitrification inhibition term [muM O2]
redsi = 16.0 / 80.0  # Redfield ratio for silica
redn = 16.0 / 106.0  # Redfield ratio for nitrogen
redp = 1.0 / 106.0  # Redfield ratio for phosphorus
kmaint = 4.6e-7  # Maintenance rate constant [s^-1]
kmort = 1.1574074074074074e-07  # Phytoplankton mortality [s^-1] - Le Fouest 2013 www.biogeosciences.net/10/4785/2013/
kexcr = 5e-2  # Excretion constant [-]
kgrowth = 2.9e-1  # Growth constant [-]
KD1 = 1.3  # Background extinction coefficient [m^-1]
KD2 = 0.06  # Suspended matter attenuation [mg^-1 m^-1]
kox = 6.08e-4  # Aerobic degradation rate [muM C s^-1]
kdenit = 5.05e-4  # Denitrification rate [muM C s^-1]
knit = 3.472222222222222e-06  # Nitrification rate [muM C s^-1] - Le Fouest 2013 www.biogeosciences.net/10/4785/2013/

# ARCTIC BIOGEOCHEMISTRY EXTENSION  (docs/arctic_biogeochemistry.md)
# ------------------------------------------------------------------
# Four process groups added on top of the shipped C-GEM network, motivated by what
# controls Arctic land-to-ocean CO2/CH4/N2O: (1) a refractory + chromophoric DOC pool
# with slow oxidation AND photomineralisation; (2) CH4 and N2O cycling + gas exchange;
# (3) a benthic (sediment) DIC/alkalinity/CH4 efflux; (4) distributed lateral loading.
# All are ADDITIVE: the existing TOC pool keeps its role as labile/bioavailable DOC,
# and every new term is gated to open water + the ice/`react` branch like the rest.

# MASTER SWITCH. The extension REACTIONS (RDOC oxidation/photolysis, CH4/N2O cycling,
# benthic efflux) run only where a site opts in with ARCTIC_BGC = True. The three new
# tracers are still transported everywhere (so they exist for every site), but with the
# switch off they are inert passive tracers and the DIC/O2/ALK/NO3 updates reduce EXACTLY
# to the shipped network -- so the four real rivers' carbon results are bit-identical
# until their RDOC/CH4/N2O boundary chemistry is constrained (they carry placeholder
# values today). The idealized verification fixture sets it True. Lateral loading has its
# own switch (LATERAL_INFLOW), likewise off for the real rivers.
ARCTIC_BGC = getattr(_site, "ARCTIC_BGC", False)

# (1) refractory/chromophoric DOC (RDOC) + photomineralisation
krefr = 1.0e-5     # RDOC aerobic oxidation MAX rate [mmol C m^-3 s^-1] (~60x slower than kox);
                   # refractory DOC turns over on months-years (Hansell 2013, Ann Rev Mar Sci)
K_RDOC = 100.0     # RDOC half-saturation for oxidation/photolysis [mmol C m^-3]
PHOTO_EFF = 5.0e-8 # apparent photomineralisation efficiency [mmol C m^-2 s^-1 per
                   # (muE m^-2 s^-1) absorbed]; tuned so summer surface photo-oxidation is a
                   # sizeable fraction of respiration, per Cory et al. 2014 (Science 345:925)
                   # where photochemistry rivals microbial DIC production in Arctic surface water
f_photo_lab = 0.2  # fraction of photoproduct returned as labile DOC (rest -> DIC directly);
                   # Cory & Kling 2018 (photo-oxidation yields both CO2 and partially-oxidised DOC)

# (2) methane and nitrous oxide
k_ch4ox = 1.2e-6   # methanotrophy MAX rate [mmol m^-3 s^-1] (turnover days-weeks; Bussmann 2013)
K_ch4 = 0.5        # CH4 half-saturation for oxidation [mmol m^-3]
K_ch4O2 = 5.0      # O2 half-saturation for methanotrophy [mmol m^-3]
pCH4_atm = 1.9e-6  # atmospheric CH4 mole fraction [atm] (~1.9 ppm, 2022; NOAA GML)
y_nit = 0.003      # N2O-N yield per N nitrified [mol/mol] (~0.3%; Beaulieu 2011 PNAS 108:214)
y_denit = 0.005    # N2O-N yield per N denitrified [mol/mol] (~0.5%; incomplete denitrification)
pN2O_atm = 0.335e-6  # atmospheric N2O mole fraction [atm] (~335 ppb, 2022; NOAA GML)

# (3) benthic (sediment) efflux -- a first-order sediment-oxygen-demand closure (no OM pool yet)
k_sod = 2.3e-4     # sediment O2 demand MAX rate [mmol O2 m^-2 s^-1] (~20 mmol m^-2 d^-1;
                   # river/estuary SOD 5-50, Cai & Sayles 1996)
K_sod = 10.0       # O2 half-saturation for SOD [mmol m^-3]
k_bdenit = 5.0e-5  # benthic denitrification MAX rate [mmol N m^-2 s^-1]
k_methano = 2.0e-8 # benthic methanogenesis MAX rate [mmol CH4 m^-2 s^-1], only where O2 is low

# (4) distributed lateral loading (per-site; OFF unless the site sets LATERAL_INFLOW).
# Tundra/thermokarst/tributary inputs entering ALONG the channel rather than only at the
# upstream boundary. Total lateral inflow is spread uniformly over the domain and injected
# as a mixing source in main.py; the four real rivers leave it off (unconstrained), so
# their results are unchanged. See sites/idealized.py for the fixture that exercises it.
LATERAL_INFLOW = getattr(_site, "LATERAL_INFLOW", 0.0)     # total lateral inflow [m^3 s^-1], 0 = off
LATERAL_CONC = getattr(_site, "LATERAL_CONC", {})          # {species: concentration} of the lateral water

# EXTERNAL FORCINGS
# Qr = -177.0  # River discharge [m^3/s]
AMPL = _site.AMPL    # Tidal amplitude at mouth [m]  -- fallback single-sinusoid only
pfun = _site.pfun    # Tidal frequency [cycle/hr]    -- fallback single-sinusoid only

# TIDES
# -----
# Sea-surface elevation at the mouth is a sum of harmonic constituents from the
# nearest NOAA CO-OPS station (tools/build_tides.py), reconstructed in
# fun_module.Tide as  eta(t) = sum_i A_i cos(speed_i*(t/3600) - G_i).
# TIDE_AMP [m], TIDE_SPEED [deg/hr], TIDE_PHASE [deg] are parallel arrays.
# If no constituent data exists for the site, these stay empty and Tide falls back
# to the single-sinusoid (AMPL, pfun) above. The Beaufort coast is microtidal
# (~0.4 m range), so this is a realistic-structure refinement, not a big amplitude
# change from the shipped 0.2 m idealisation.
import json as _json
_tide_file = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                          "..", "forcing", "tidal_constituents.json")
TIDE_AMP, TIDE_SPEED, TIDE_PHASE = [], [], []
TIDE_STATION = ""
if os.path.exists(_tide_file):
    with open(_tide_file) as _tf:
        _td = _json.load(_tf).get(SITE)
    if _td:
        TIDE_STATION = _td["station_name"]
        for _c in _td["constituents"]:
            TIDE_AMP.append(_c["amp_m"])
            TIDE_SPEED.append(_c["speed_deg_hr"])
            TIDE_PHASE.append(_c["phase_deg"])
#Uw_sal = 8  # Wind speed in saline estuary [m/s]
#Uw_tid = 8  # Wind speed in tidal river [m/s]
#water_temp = 12  # Water temperature [C]
#pCO2 = 380 * 1e-6  # CO2 partial pressure in the atmosphere [atm]
repeatYear = 1  # equal to 1 if 365-day climatology (repeat forcings)
nbday_ice = 10  # number of day of negative water temperature that will stop all reactions (transport, biogeo, sed)


# SURFACE HEAT BUDGET
# -------------------
# Drives the transported temperature field v['T'] (heat_module.py). Without this,
# T is only a mixing tracer between the sea and river boundary values and nothing
# warms the channel interior.
#
#     rho_w cp_w H dT/dt = Q_sw + Q_lw + Q_sens + Q_lat
#
# THE HUMIDITY ASSUMPTION. Q_lat needs vapour pressure, and the PRDA2 buoy's DEWP
# field is entirely empty for 2022 -- there is no humidity observation for this
# coast in the data we hold. REL_HUMIDITY below is therefore an ASSUMED constant.
# 0.85 is representative of an Arctic coastal summer, but it is a parameter, not a
# measurement, and latent heat is not a small term in this budget. Treat modelled
# temperature as uncertain in proportion to how wrong this is; source ERA5-Land
# humidity if that uncertainty matters. Set HEAT_BUDGET = False to disable the whole
# term set and fall back to pure advection.
HEAT_BUDGET = True

REL_HUMIDITY = 0.85     # ASSUMED -- no observation available. See above.
albedo_water = 0.06     # open-water shortwave albedo [-]
emiss_water = 0.97      # water longwave emissivity [-]
sigma_SB = 5.67e-8      # Stefan-Boltzmann [W m^-2 K^-4]
rho_air = 1.3           # air density [kg m^-3]
cp_air = 1005.0         # air specific heat [J kg^-1 K^-1]
cp_water = 4186.0       # water specific heat [J kg^-1 K^-1]
L_vap = 2.5e6           # latent heat of vaporisation [J kg^-1]
C_sens = 1.3e-3         # bulk sensible transfer coefficient [-]
C_lat = 1.3e-3          # bulk latent transfer coefficient [-]
P_atm = 1013.0          # reference surface pressure [hPa]
T_FREEZE = 0.0          # freshwater freezing point [degC]
H_MIN_HEAT = 0.05       # floor on depth in the heat budget [m]; without it a cell
                        # approaching bottom-fast conditions gives dT/dt -> infinity

# RIVER ICE MODEL  (item 2 / docs/ice_model_plan.md step c)
# ---------------------------------------------------------
# Prognostic ice thickness per cell, grown from the heat-budget energy deficit and by
# conduction through the ice slab (Stefan-type), melted by surface heat at breakup.
# Ice insulates the water (heat budget skips ice-covered cells), stops gas exchange,
# attenuates under-ice light, and where it reaches the bed (BOTTOM-FAST) closes the
# channel. Breakup is HYDRAULIC: the freshet discharge surge mechanically clears ice.
# Set ICE_MODEL = False to fall back to the crude previousdays gate.
ICE_MODEL = os.environ.get("CGEM_ICE", "on") != "off"
rho_ice = 917.0         # ice density [kg m^-3]
L_fusion = 3.34e5       # latent heat of fusion [J kg^-1]
k_ice = 2.2             # thermal conductivity of ice [W m^-1 K^-1]
albedo_ice = 0.6        # ice/snow shortwave albedo [-]
k_ice_PAR = 1.5         # PAR attenuation coefficient in ice [m^-1] (under-ice light)
H_MIN_ICE = 0.01        # floor on ice thickness in the conduction term [m]
ICE_FORM_THRESH = 0.005 # ice thickness above which a cell counts as ice-covered [m]
# Hydraulic breakup: clear ice where discharge exceeds this multiple of the ANNUAL
# MEAN discharge (the freshet surge). North Slope breakup is mechanical, not thermal.
# Referenced to the mean, not the winter baseflow, because the ice-affected winter
# percentile is ~0 on these gauges; freshets run 6-30x the mean, so 3x fires only then.
BREAKUP_Q_FACTOR = 3.0

# OTHER PARAMETERS
Euler = 0.5772156649  # Euler's constant
PI = math.pi  # Pi value
pH_ite = 50  # number of iterations to converge to pH following Follows et al., 2006
mass_mol_B = 10.8110  # molar mass of Boron g/mol

# CARBONATE UNITS -- now solved consistently in mol/kg (see biogeo_module /
# fun_module._h_solve_kg, ported from C-GEM v2). The state (DIC/ALK, mmol/m^3) is
# converted to mol/kg via the local in-situ density before speciation, so it shares the
# mol/kg basis of K0/K1/K2/KB and boron is a plain 0.1336*S mg/kg -> mol/kg with no
# scale factor. This REPLACES the earlier `CARB_BT_SCALE` patch, which only rescaled the
# borate term of the mmol/m^3-vs-mol/kg-mixed legacy solve. Because boron scales with
# SALINITY and these rivers run ~fresh (S~0) in the domain, the correction moves FCO2
# only slightly (see docs / CLAUDE.md). No toggle: the unit-correct solve is the path.

# NUMERICAL INTEGRATION
# MAXT/WARMUP may be overridden in days via the environment, which is how to run a
# smoke test without editing (and later forgetting to revert) this file:
#     CGEM_MAXT_DAYS=2 CGEM_WARMUP_DAYS=1 CGEM_SITE=kuparuk python main.py
# Biogeochemistry and sediment are skipped until t > WARMUP, so a run shorter than
# WARMUP exercises only hydrodynamics and transport.
MAXT = int(os.environ.get("CGEM_MAXT_DAYS", 365 * 2)) * 24 * 60 * 60  # Max time [s]
WARMUP = int(os.environ.get("CGEM_WARMUP_DAYS", 365)) * 24 * 60 * 60  # Warmup [s]
DELTI = 75  # Delta t [s]
# Save every TS timesteps. Default 12 -> every 900 s (~15 min), which produces ~4.6 GB
# of .dat output per 2-year site. Override via CGEM_TS to thin it: CGEM_TS=48 gives
# hourly output at ~1/4 the disk (plenty for diagnostics; still resolves the ~12 h
# tide). Output cadence only -- does not affect the integration.
TS = int(os.environ.get("CGEM_TS", 12))

# OUTPUT FORMAT
# -------------
# 'nc'   -> compact NetCDF (one output.nc per site; time x distance; ~5-10x smaller
#           than .dat, self-describing, no trailing-tab-NaN quirk).  [default]
# 'dat'  -> the legacy tab-separated .dat files (one per field).
# 'both' -> write both. Override via CGEM_OUTPUT.
OUTPUT_FORMAT = os.environ.get("CGEM_OUTPUT", "nc")
DELXI = 200  # Delta x [m]
TOL = 1e-10  # Convergence criterion
M = int(EL / DELXI) + 1  # Max even grid points
if M % 2 == 0:
    pass
else:
    M = M - 1
M1 = M - 1  # Max odd grid points
M2 = M - 2  # Last even grid point
M3 = M - 3  # Last odd grid point
MAXV = 13  # Max species in the registry (DIA,dSi,NO3,NH4,PO4,O2,TOC,S,SPM,DIC,ALK,pH,T).
           # Informational only -- the species live in the `variables.v` dict, keyed by
           # name; nothing sizes an array from MAXV. Was 12 before T was added.

# DISPERSION
# ----------
# 'seo'      : Seo & Cheong (1998) longitudinal dispersion for RIVERS,
#              K = 5.915 (W/H)^0.620 (U/u*)^1.428 H u*,  u* = sqrt(g) U / Chezy.
# 'savenije' : the original C-GEM Van der Burgh / Savenije estuary formulation.
#
# Why the change. The Savenije form collapses under the observation-based geometry:
# with realistic shallow depths, beta = -(K LC Qr)/(D0 B_lb DEPTH_lb) explodes and
# dispersion[i] is clamped to zero after only a few grid points --
#
#     geometry                          beta     dispersion non-zero over
#     shipped (B=1215/859, D=15 m)       0.5     all 137 pts (27.4 km)
#     colville (1074/400, D=2.25)        8.2     16 pts (3.2 km)
#     kuparuk  (102/48,  D=1.34)        69.5      3 pts (0.6 km)
#     sagavanirktok (98/54, D=0.98)    138.8      2 pts (0.4 km)
#
# i.e. transport became purely advective. That is the same root cause as the width
# misfit: an estuary parameterisation applied to rivers with negligible tides.
#
# Seo & Cheong rather than Fischer (1979): Fischer assumes W/H ~ 10-100, but Colville
# runs W/H = 477, where it returns 21661 m2/s -- above the numerical bound below.
# Seo & Cheong was regressed across a wider W/H range and gives 342-474 m2/s here.
DISPERSION_MODEL = "seo"

# Numerical ceiling on dispersion. schemes_module.disp_sch builds
#     r[i] = 1 - (c1 + c2) * DELTI / (8 DELXI^2)
# and r < 0 makes the Crank-Nicolson right-hand side oscillatory. With c1 ~ c2 ~ K
# that requires K <= 4 DELXI^2 / DELTI. Exceeding it is a numerical failure, not a
# physical statement, so the value is capped and the cap is reported once per run.
DISP_MAX = 4.0 * DELXI ** 2 / DELTI