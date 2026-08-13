"""Colville River -- largest North Slope drainage (~60,000 km2), Beaufort coast."""
from ._baseline import WATERTEMP_FILE, AMPL, pfun, BOUNDARIES as _BASE_BOUNDARIES  # noqa: F401

# --- riverine boundary chemistry from sufficiently-upstream WQP grabs ----------
# Upstream (cub) concentrations from the "Colville R nr Nuiqsut" mainstem station
# (WQP, ~32 km inland at the delta head -- i.e. essentially the model's upstream
# boundary location, and freshwater not brackish). Open-water medians, converted to
# model units (mmol/m^3): alkalinity mg/L CaCO3 x 19.98; DOC mg/L x 83.3; N species
# mg/L-as-N x 71.4; PO4 mg/L-as-PO4 x 10.5. Colville drains carbonate-rich terrain,
# so alkalinity is genuinely high; the large correction was TOC (358 vs placeholder
# 1582, ~5x lower). The marine clb end is unchanged.
#
# CARBONATE: DIC is solved to reproduce the observed PAIRED-sample pCO2. The 74
# co-located ALK+pH grabs at USGS-15880000 give a median water pCO2 of 658 uatm
# (SUPERSATURATED -> the river outgasses), so DIC is back-solved from the median ALK
# (1319) and that pCO2 -> 1337. This replaces the earlier median-ALK + median-pH
# (separate) solve, which gave DIC 1302 and a spurious pCO2 ~459 (near atmospheric):
# taking pH and ALK medians independently understates the flux-relevant pCO2, since
# pH and ALK covary and the CO2-rich samples are what drive the air-sea flux.
# tools/usgs_carbonate_boundary.py siteid USGS-15880000.
BOUNDARIES = dict(_BASE_BOUNDARIES)
for _sp, _cub in [("NO3", 5.0), ("NH4", 4.0), ("PO4", 0.5), ("TOC", 358.0),
                  ("pH", 7.82), ("ALK", 1318.9), ("DIC", 1336.6)]:
    BOUNDARIES[_sp] = (_BASE_BOUNDARIES[_sp][0], _cub)

# --- O2 and SPM (cub): USGS Water Quality Portal, same gauge as above (siteid USGS-15880000)
# ---------------------------------------------------------------------------------------
# Open-water (Jun-Sep) medians. O2: n=63 discrete DO samples, 10.4 mg/L -> 325.0 mmol/m^3
# (32 g/mol). SPM: n=74 Suspended Sediment Concentration samples, 16.5 mg/L -> 0.0165 g/L
# -- two orders of magnitude below the shared placeholder (2.0 g/L), consistent with these
# being clear, low-sediment-load tundra rivers, not the turbid estuary C-GEM shipped with.
for _sp, _cub in [("O2", 325.00), ("SPM", 0.0165)]:
    BOUNDARIES[_sp] = (_BASE_BOUNDARIES[_sp][0], _cub)

# --- DIA (cub): epilithic->water-column proxy, NOT a phytoplankton measurement ----------
# North Slope rivers have essentially no water-column phytoplankton monitoring (checked
# WQP + literature) -- what IS measured is EPILITHIC chlorophyll (algae attached to the
# streambed), a different compartment than DIA (this model's transported, advected tracer).
# Kuparuk reference (unfertilized) reach epilithic chlorophyll = 2.8 mg Chl/m2 on rock (85%
# of bed cover), Slavik et al. 2004 (Ecology 85:939-954) Table 2, July 1998 -- the only
# quantitative North Slope epilithic value found, applied to all four rivers for lack of a
# site-specific one. Converted to a rough water-column-equivalent by DIVIDING BY THIS SITE'S
# DEPTH (as if the entire bed stock were resuspended and evenly mixed through the water
# column -- NOT how these compartments actually behave; a proxy, not a measurement), then to
# carbon via C:Chl = 75 gC/gChla (user-specified; close to the model's own Chla2CMIN-implied
# maximum of 1/0.0125 = 80 gC/gChla) and to mmol via /12.011 g/mol:
#   DIA = (2.8 / 2.25) * 75 / 12.011 = 7.771 mmol C/m^3   (this site's DEPTH_ub = 2.25 m,
#   set below in the geometry section -- computed here as a literal since DEPTH_ub isn't
#   defined yet at this point in the file)
BOUNDARIES["DIA"] = (_BASE_BOUNDARIES["DIA"][0], 7.771)

# --- MARINE (clb): ECCO-Darwin v5, nearest wet LLC270 cell to the delta mouth ---------
# Cell (70.568N, -150.789E), ~14.6 km offshore of the SWORD mouth point (70.44N, -150.7W).
# Climatological annual mean over the full archived record (272-292 months depending on
# variable, native LLC270 monthly output, public NASA data.nas.nasa.gov/ecco portal, no
# auth). S is SALTanom + 35 (MITgcm's standard salinity-anomaly convention). DIC/ALK/NO3/
# NH4/PO4/dSi/O2/TOC units match the model's mmol m^-3 directly, no conversion. Built by
# scratch/ecco_darwin/extract_all_rivers.py (not tracked -- rerun to reproduce or extend).
# NOT sourced this way: DIA (ECCO-Darwin gives Chl1-5, not carbon biomass -- would need an
# uncertain C:Chl ratio), pH (not in the archived output; could instead be diagnosed from
# this DIC/ALK/S pair via the model's own carbonate solve), RDOC/CH4/N2O/SPM (no analog).
for _sp, _clb in [("S", 30.6588), ("DIC", 2042.7683), ("ALK", 2143.9555), ("NO3", 5.1300),
                   ("NH4", 0.0338), ("PO4", 0.7973), ("dSi", 8.1513), ("O2", 372.9380),
                   ("TOC", 39.3797)]:
    BOUNDARIES[_sp] = (_clb, BOUNDARIES[_sp][1])   # override marine clb, keep river cub

BOUNDARY_CHEM_SOURCE = ("WQP Colville R nr Nuiqsut 15880000 (paired-pCO2 DIC); "
                        "marine (clb): ECCO-Darwin v5 nearest-cell climatology")

LABEL = "Colville"
DISCHARGE_FILE = "colville_river_discharge_2022_m3sec.csv"
# Wind-driven storm surge at the marine boundary. Prudhoe Bay (9497645) is the only
# continuous 2022 water-level record on this coast, used as the REGIONAL proxy for all
# four rivers (see tools/build_surge.py); added to the harmonic tide in fun_module.Tide.
SURGE_FILE = "surge_prudhoe_2022_m.csv"

# USGS 15875000 (Colville R at Umiat). Gauged FAR upstream of the delta -- Umiat
# drains roughly half the basin -- so this understates discharge at the mouth.
GAUGE = "15875000"
DISCHARGE_IS_UPSTREAM_PROXY = True

# --- geometry: SWORD v17b nodes + USGS channel surveys -----------------------
# Width from SWORD v17b node width / n_chan_mod (per conveyance channel; raw SWORD
# width is the sum across braids and this river runs 3-4 braided channels upstream,
# which made it appear to WIDEN inland before the correction). Lakes excluded,
# restricted to the model domain, 104 nodes.
EL = 4150         # observed estuary length [m] (4.15 km)
# Grid refined for the shorter EL: the shared default DELXI=200/DELTI=75 collapses this
# domain to M=20 points. Colville's own 200-day/freshet diagnostic run didn't show as
# unambiguous a blow-up as Kuparuk's (whose upstream boundary spiked to 4x nominal depth
# under the shared default -- see kuparuk.py), but the same few-points-at-a-short-domain
# structural issue applies, so the same precautionary refinement is used here too. DELXI
# halved (M=42) with DELTI cut to keep dispersion's Crank-Nicolson stability margin
# (DISP_MAX = 4 DELXI^2/DELTI) well above the ~350-650 m^2/s Seo & Cheong range. See
# config.py's "GRID SPACING & CFL" note.
DELXI = 100       # Delta x [m] (was 200 shared default) -> M=42
DELTI = 30        # Delta t [s] (was 75 shared default)
# DELTA-MOUTH WIDTH = raw SWORD v17b distributary SUM at the Harrison Bay delta
# (1310 + 240 m across the two seaward-most channels, tools/extract_sword.py). The
# fresh SWORD extraction shows Colville's mouth is a wide two-channel delta, not the
# single 1207 m channel the earlier one-off extraction gave -- so it is widened here too.
B_lb_perchan = 1207     # earlier single-channel seaward width [m]
SWORD_MOUTH_SUM = 1550  # raw SWORD distributary sum [m]
B_lb = SWORD_MOUTH_SUM
B_ub = 423         # PRISMATIC width upstream of the flare [m]
L_FLARE = 5500    # flare length [m]; width converges over 0..L_FLARE,
                   # then is prismatic. Fit R2 = 0.392 (vs exponential-over-27km).

# Multi-channel geometry, used ONLY when CGEM_MULTICHANNEL=on (see config.MULTICHANNEL).
# N_CHAN_LB = the 2 seaward distributaries SWORD resolves at Harrison Bay (whose widths
# sum to B_lb). B_UB_TOTAL = the median RAW SWORD width in the prismatic reach (>5.5 km,
# main stem, 102 nodes) -- i.e. the braided TOTAL conveyance, measured directly rather
# than reconstructed from B_ub x a braid count. config derives the thread count from it
# (1052/423 = 2.49); note that is NOT the median n_chan_mod of 3.0, because B_ub is a
# median of per-node ratios and the median of a ratio is not the ratio of medians.
#
# The consequence: the total prismatic conveyance is 1052 m against a 1550 m mouth, so
# roughly HALF of Colville's apparent "flare" is the B_lb(sum)/B_ub(per-channel) definition
# change rather than real convergence -- and the pre-adoption 3x velocity acceleration just
# inside the mouth (0.156 -> 0.467 m/s over 5.4 km) was an artifact of it.
# tools/extract_sword.py.
N_CHAN_LB = 2.0
B_UB_TOTAL = 1052.0     # SWORD v17b raw prismatic median [m]; IQR 792-1354

# Depth from at-a-station hydraulic geometry D = 0.360*Q^0.297 (208 USGS ADCP
# surveys, 1953-2026) evaluated at the 2022 open-water mean discharge, 484.5 m3/s.
# CAVEAT: the gauge is far upstream of the delta, so this likely UNDERSTATES depth
# in the estuarine channel the model represents. It replaces a shipped 15 m whose
# cited source (the Sagavanirktok gauge) shows a 0.88 m median.
DEPTH_lb = 2.25
DEPTH_ub = 2.25

distance = 1      # grid points in saline zone -- almost no saline intrusion
GEOMETRY_IS_PLACEHOLDER = False
