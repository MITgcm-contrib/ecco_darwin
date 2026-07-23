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
# DIC is NOT in the LTER record. Rather than guess a ratio, it is SOLVED from the
# two carbonate variables that ARE observed -- pH 7.32 and ALK 274.4 -- via the
# model's own carbonate system at S=0, T~10 C, giving DIC = 284.8 (which reproduces
# pH 7.32 exactly, confirmed). So DIC is derived, not independently measured.
# Only Kuparuk gets this -- the other three rivers have no LTER data and keep the
# shared placeholder. Species not covered here (dSi, O2, DIA, S, SPM) also stay
# placeholder. See tools/lter_boundary.py.
BOUNDARIES = dict(_BASE_BOUNDARIES)
for _sp, _cub in [("NO3", 3.46), ("NH4", 0.37), ("PO4", 0.05), ("TOC", 295.6),
                  ("pH", 7.32), ("ALK", 274.4), ("DIC", 284.8)]:
    BOUNDARIES[_sp] = (_BASE_BOUNDARIES[_sp][0], _cub)   # keep marine clb, override river cub
BOUNDARY_CHEM_SOURCE = "Arctic LTER headwater (reference reach), knb-lter-arc.10303"

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

# --- geometry: SWORD v17c nodes + USGS channel surveys -----------------------
# Width from SWORD per-channel node width, 149 nodes. This is the ONLY one of the
# four whose width profile actually fits C-GEM's exponential convergence with any
# skill (R2 = 0.23 per-channel, 0.63 on raw width over the model domain).
EL = 27175
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
DEPTH_lb = 1.34
DEPTH_ub = 1.34

distance = 1
GEOMETRY_IS_PLACEHOLDER = False
