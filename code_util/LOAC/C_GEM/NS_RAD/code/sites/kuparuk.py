"""Kuparuk River -- ~8,100 km2, discharges near Prudhoe Bay."""
from ._baseline import AMPL, pfun, BOUNDARIES as _BASE_BOUNDARIES  # noqa: F401

# Observed USGS 00010 (blended to the regression outside the 55-day
# record). The regression was already unbiased here (+0.05 C); this is a
# consistency check. See tools/build_river_temp_obs.py.
WATERTEMP_FILE = "kuparuk_watertemp_obs_2022_degC.csv"

# --- riverine boundary chemistry from Arctic LTER (headwater proxy) ----------
# The upstream (cub) concentrations for these species come from the Arctic LTER
# Streams Chemistry record (EDI knb-lter-arc.10303), NATURAL "Reference" reach only
# -- the fertilized reaches are a decades-long P-enrichment EXPERIMENT and must not
# be used. Open-water medians, 1978-2019. Units: uM = mmol/m^3 = model units.
#
# BIG CAVEAT, do not lose it: the LTER site is at the Dalton Highway crossing
# (68.6 N, 720 m), ~190 km / ~266 river-km upstream of the mouth -- ~163 km ABOVE
# this model's own upstream boundary. It is the same river but its HEADWATERS.
# Downstream the river integrates tributaries, mineral weathering and in-stream
# processing, so these are a first-order proxy for the delta boundary, not a
# measurement of it. Alkalinity especially (274 vs the old placeholder 1596) is
# soft tundra headwater and is likely a LOWER bound for the delta. Treat as
# "far better than the shipped placeholder, still not delta-observed."
#
BOUNDARIES = dict(_BASE_BOUNDARIES)

# --- NUTRIENTS: Arctic LTER headwater reference reach (see caveat above) ------
# NO3/NH4/PO4/TOC still come from the LTER headwater record -- no delta-proximal
# nutrient series is available. Species not covered (dSi, DIA, S) stay placeholder.
# See tools/lter_boundary.py. (As of 2026-08, the EDI/PASTA data endpoint this tool reads
# returns 403 for anonymous access -- new bot protection, not a dataset change; re-running
# it may need a session-authenticated download instead of a bare urllib GET.)
for _sp, _cub in [("NO3", 3.46), ("NH4", 0.37), ("PO4", 0.05), ("TOC", 295.6)]:
    BOUNDARIES[_sp] = (_BASE_BOUNDARIES[_sp][0], _cub)   # keep marine clb, override river cub

# --- O2 and SPM (cub): USGS Water Quality Portal, delta gauge (siteid USGS-15896000) -----
# Open-water (Jun-Sep) medians. O2: n=32 discrete DO samples, 12.15 mg/L -> 379.7 mmol/m^3
# (32 g/mol). SPM: n=52 Suspended Sediment Concentration samples, 7.5 mg/L -> 0.0075 g/L --
# two orders of magnitude below the shared placeholder (2.0 g/L); Kuparuk is the clearest
# of the four (lowest SPM), consistent with it being closest to tidewater/best-constrained.
for _sp, _cub in [("O2", 379.69), ("SPM", 0.0075)]:
    BOUNDARIES[_sp] = (_BASE_BOUNDARIES[_sp][0], _cub)

# --- DIA (cub): epilithic->water-column proxy, NOT a phytoplankton measurement ----------
# North Slope rivers have essentially no water-column phytoplankton monitoring (checked
# WQP + literature) -- what IS measured is EPILITHIC chlorophyll (algae attached to the
# streambed), a different compartment than DIA (this model's transported, advected tracer).
# Kuparuk's OWN reference (unfertilized) reach epilithic chlorophyll = 2.8 mg Chl/m2 on
# rock (85% of bed cover), Slavik et al. 2004 (Ecology 85:939-954) Table 2, July 1998 --
# this is the site the LTER fertilization study was actually conducted on, so unlike the
# other three rivers this value is not borrowed. Converted to a rough water-column-
# equivalent by DIVIDING BY THIS SITE'S DEPTH (as if the entire bed stock were resuspended
# and evenly mixed through the water column -- NOT how these compartments actually behave;
# a proxy, not a measurement), then to carbon via C:Chl = 75 gC/gChla (user-specified;
# close to the model's own Chla2CMIN-implied maximum of 1/0.0125 = 80 gC/gChla) and to
# mmol via /12.011 g/mol:
#   DIA = (2.8 / 1.34) * 75 / 12.011 = 13.048 mmol C/m^3   (this site's DEPTH_ub = 1.34 m,
#   set below in the geometry section -- computed here as a literal since DEPTH_ub isn't
#   defined yet at this point in the file)
BOUNDARIES["DIA"] = (_BASE_BOUNDARIES["DIA"][0], 13.048)

# --- CARBONATE: the river's OWN near-tidewater gauge (DELTA-PROXIMAL) ---------
# ALK and pH are open-water (Jun-Sep) samples at gauge 15896000 (nr Deadhorse,
# ~tidewater) -- NOT the 163-km-upstream LTER headwater. DIC is solved to reproduce
# the observed PAIRED-sample pCO2: 48 co-located ALK+pH grabs give median water pCO2
# 796 uatm (SUPERSATURATED -> outgas), and DIC is back-solved from the median ALK
# (1099) and that pCO2 -> 1141 (tools/usgs_carbonate_boundary.py siteid USGS-15896000).
# The delta alkalinity (1099) is ~4x the soft-water headwater value (274) -- downstream
# mineral weathering, exactly the "LTER is a lower bound for the delta" caveat above --
# so the river outgasses as Arctic rivers physically do. The retired headwater proxy
# (pH 7.32, ALK 274.4, DIC 284.8) sat at pCO2 ~ 421 uatm (== atmosphere), holding the
# model spuriously near equilibrium.
for _sp, _cub in [("pH", 7.64), ("ALK", 1099.1), ("DIC", 1141.1)]:
    BOUNDARIES[_sp] = (_BASE_BOUNDARIES[_sp][0], _cub)   # keep marine clb, override river cub

# --- MARINE (clb): ECCO-Darwin v5, nearest wet LLC270 cell to the delta mouth ---------
# Cell (70.516N, -148.819E), ~11 km offshore of the SWORD mouth point (70.42N, -148.9W).
# Climatological annual mean over the full archived record (272-292 months depending on
# variable, native LLC270 monthly output, public NASA data.nas.nasa.gov/ecco portal, no
# auth). S is SALTanom + 35 (MITgcm's standard salinity-anomaly convention -- the +35
# reconstruction gives a physically sane Beaufort-shelf value, cross-checked against
# O2 ~374 mmol/m3 landing close to the old placeholder 392). DIC/ALK/NO3/NH4/PO4/dSi/O2/TOC
# units match the model's mmol m^-3 directly, no conversion. Built by
# scratch/ecco_darwin/extract_point.py (not tracked -- rerun to reproduce or extend).
# NOT sourced this way: DIA (ECCO-Darwin gives Chl1-5, not carbon biomass -- would need an
# uncertain C:Chl ratio), pH (not in the archived output; could instead be diagnosed from
# this DIC/ALK/S pair via the model's own carbonate solve), RDOC/CH4/N2O/SPM (no analog).
for _sp, _clb in [("S", 30.2522), ("DIC", 2025.1029), ("ALK", 2123.3752), ("NO3", 5.6323),
                   ("NH4", 0.0358), ("PO4", 0.8545), ("dSi", 8.1139), ("O2", 373.7506),
                   ("TOC", 32.7940)]:
    BOUNDARIES[_sp] = (_clb, BOUNDARIES[_sp][1])   # override marine clb, keep river cub

BOUNDARY_CHEM_SOURCE = ("carbonate: USGS 15896000 delta gauge (paired-pCO2 DIC); "
                        "nutrients: Arctic LTER reference reach, knb-lter-arc.10303; "
                        "marine (clb): ECCO-Darwin v5 nearest-cell climatology")

LABEL = "Kuparuk"
DISCHARGE_FILE = "kuparuk_river_discharge_2022_m3sec.csv"

# Wind-driven storm-surge forcing at the marine boundary, from OBSERVED 2022 water level
# at Prudhoe Bay (NOAA CO-OPS 9497645 -- Kuparuk's own tidal station), daily-mean
# residual (tools/build_surge.py). The Beaufort coast is microtidal (~0.3 m tide), so
# this non-tidal surge (up to +0.66 m in 2022) is the real saltwater-intrusion driver.
# Only Kuparuk has an observed water-level station on this coast.
SURGE_FILE = "surge_prudhoe_2022_m.csv"

# USGS 15896000 (Kuparuk R nr Deadhorse). Close to tidewater, so unlike Colville and
# Sagavanirktok this needs no upstream-proxy caveat -- best constrained of the four.
GAUGE = "15896000"
DISCHARGE_IS_UPSTREAM_PROXY = False

# --- geometry: SWORD v17b nodes + USGS channel surveys -----------------------
# Width from SWORD per-channel node width, 149 nodes. This is the ONLY one of the
# four whose width profile actually fits C-GEM's exponential convergence with any
# skill (R2 = 0.23 per-channel, 0.63 on raw width over the model domain).
EL = 7156         # observed estuary length [m] (7.156 km)
# Grid refined for the shorter EL: the shared default DELXI=200/DELTI=75 collapses this
# domain to M=36 points and, empirically (a 200-day run through spring freshet), is
# numerically unstable -- the upstream boundary spiked to U=3.97 m/s, depth=6.3 m against
# a nominal 1.34 m, a transient blow-up, not physical. DELXI halved (M=72) with DELTI cut
# to keep dispersion's Crank-Nicolson stability margin (DISP_MAX = 4 DELXI^2/DELTI) well
# above the ~350-650 m^2/s Seo & Cheong range. See config.py's "GRID SPACING & CFL" note.
DELXI = 100       # Delta x [m] (was 200 shared default) -> M=72
DELTI = 30        # Delta t [s] (was 75 shared default)
# DELTA-MOUTH WIDTH = distributary SUM (see CLAUDE.md "braided vs deltaic"). At the
# delta the distributaries diverge to separate sea outlets and never rejoin, so they
# convey the total discharge to the coast in PARALLEL -- the estuarine/salt-exchange
# width is the SUM, not the per-channel width. Using per-channel width with total Q
# over-estimated the mouth velocity ~n_chan-fold and flushed the salt intrusion to ~0.
# INTERIM: per-channel x N_CHAN_MOUTH, a morphological distributary count for Gwydyr Bay
# (delta ~5 km wide). To be replaced by the exact raw SWORD node-width sum (SWORD_MOUTH_SUM
# below, None until the v17 extraction lands). Converges to the per-channel width upstream.
B_lb_perchan = 123       # SWORD per-channel width at the seaward boundary [m]
N_CHAN_MOUTH = 4         # (interim distributary count, superseded by SWORD below)
# Raw SWORD v17b delta conveyance = SUM of the seaward-most distinct distributary reaches
# in Gwydyr Bay (335 + 182 m), tools/extract_sword.py. Replaces the interim 4x estimate
# (492 m) -- which happened to be within 5% of the recovered value.
SWORD_MOUTH_SUM = 516    # raw SWORD distributary sum [m]
B_lb = SWORD_MOUTH_SUM or (B_lb_perchan * N_CHAN_MOUTH)
B_ub = 58         # PRISMATIC width upstream of the flare [m] (single channel)
L_FLARE = 7000    # flare length [m]; width converges over 0..L_FLARE,
                   # then is prismatic. Fit R2 = 0.431 (vs exponential-over-27km).
# Depth from D = 0.291*Q^0.309 (323 surveys, 1970-2026) at open-water mean 141.7 m3/s.
# Gauge sits near tidewater, so this depth is the most representative of the four.
#
# KEPT at the mean-discharge value deliberately, after the EL=7156 m shortening (below)
# was found to interact badly with this depth during freshet: the domain now reaches
# just past L_FLARE, where the channel is fully narrowed to its ~60 m prismatic width
# (confirmed correct -- see width_at() and the actual model "width" output field, both
# checked directly), and Kuparuk's freshet discharge peaks at ~23.7x its mean (day 152,
# 2022). At fixed 1.34 m depth that implies a physically-impossible ~19 m/s velocity
# through that cross-section, so the model instead raises the water level (observed:
# depth up to ~5x nominal, velocity up to ~4 m/s) to actually pass the flow -- a
# self-consistent hydrodynamic response, not a numerical instability (grid refinement to
# DELXI=100/DELTI=30, below, was tested and does not change this behavior). There is no
# spatially-resolved survey data to derive a better depth specifically at this narrow
# boundary (D=c*Q^f is a single-gauge, at-a-station relation, not a spatial profile), and
# no floodplain/overbank storage in this 1-D fixed-bank model -- so an extreme freshet
# response here is accepted as a real consequence of a narrow channel carrying a very
# large flood, not something to paper over with an untethered depth value.
DEPTH_lb = 1.34
DEPTH_ub = 1.34

# Multi-channel geometry for CGEM_MULTICHANNEL=on (config.MULTICHANNEL). B_UB_TOTAL is
# the OBSERVED raw braided total beyond the flare; config derives the thread count from it. Two seaward
# distributaries at the delta (SWORD, widths 335+182 = B_lb); the prismatic reach is
# SINGLE-THREAD (see the per-channel note above), so N_CHAN_UP = 1 and the flare here
# is a REAL convergence, not a definition change. This site is the control.
N_CHAN_LB = 2.0
B_UB_TOTAL = 60.0      # SWORD v17b raw prismatic median [m]; IQR 42-78, 127 nodes.
                        # 60 vs B_ub 58 -> thread count 1.03: SINGLE-THREAD, confirmed by
                        # median n_chan_mod = 1.0. This site's flare is REAL convergence,
                        # not a definition change -- the control for the experiment.

distance = 1
GEOMETRY_IS_PLACEHOLDER = False
