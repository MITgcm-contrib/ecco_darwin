"""Sagavanirktok River -- ~14,900 km2, enters the Beaufort east of Prudhoe Bay."""
from ._baseline import AMPL, pfun, BOUNDARIES as _BASE_BOUNDARIES  # noqa: F401

# --- riverine boundary chemistry: upstream WQP grabs + literature TOC ----------
# Upstream (cub) from WQP mainstem stations Sagwon / Franklin Bluffs (~55-137 km
# inland, unambiguously freshwater). The Sagavanirktok drains Brooks Range limestone,
# so it is a HARD-water river: ALK 2078 (n=23) and DIC 2017 (solved from ALK + the
# observed pH 8.00, n=24) are genuinely high, above the placeholder. NO3 11 (n=27) is
# well sampled. TOC = 250 is LITERATURE-augmented (Rember & Trefry 2004 report Sag
# DOC 167-742 uM over the season; the WQP upstream DOC has n=1 and is unreliable) --
# flagged as such. NH4/PO4 are near zero upstream. Marine clb unchanged.
# See tools/build_boundary_chem.py.
BOUNDARIES = dict(_BASE_BOUNDARIES)
for _sp, _cub in [("NO3", 11.0), ("NH4", 1.0), ("PO4", 0.1), ("TOC", 250.0),
                  ("pH", 8.00), ("ALK", 2078.0), ("DIC", 2017.0)]:
    BOUNDARIES[_sp] = (_BASE_BOUNDARIES[_sp][0], _cub)
BOUNDARY_CHEM_SOURCE = "WQP Sagwon/Franklin Bluffs (upstream) + Rember 2004 TOC (lit)"

# Observed USGS 00010 (blended to the regression outside the 93-day
# open-water record). Removes the +2.6 C regression warm bias directly. See
# tools/build_river_temp_obs.py.
WATERTEMP_FILE = "sagavanirktok_watertemp_obs_2022_degC.csv"

LABEL = "Sagavanirktok"
DISCHARGE_FILE = "sagavanirktok_river_discharge_2022_m3sec.csv"
# Wind-driven storm surge at the marine boundary. Prudhoe Bay (9497645) is the only
# continuous 2022 water-level record on this coast, used as the REGIONAL proxy for all
# four rivers (see tools/build_surge.py); added to the harmonic tide in fun_module.Tide.
SURGE_FILE = "surge_prudhoe_2022_m.csv"

# USGS 15908000 (nr Pump Sta 3), ~130 km inland, under half the basin -- 2022 mean
# only 47.6 m3/s, which reflects gauge placement rather than the delta.
GAUGE = "15908000"
DISCHARGE_IS_UPSTREAM_PROXY = True

# --- geometry: SWORD v17c nodes + USGS surveys, WITH A BORROWED RATIO --------
# B_lb = 98 m is observed (SWORD per-channel near the mouth). B_ub is NOT observed.
# SWORD shows this river WIDENING upstream even after the per-channel correction
# (98 -> 119 m, i.e. a negative convergence length), from only 46 nodes with a data
# gap between 5 and 20 km. A negative LC would make the model funnel backwards, so
# the convergence RATIO 0.55 is borrowed from the Canning -- the nearest
# single-channel Brooks Range river -- and applied to the observed B_lb.
# Same donor-transfer caveat as the Canning discharge: the shape is assumed.
EL = 27175
# DELTA-MOUTH WIDTH = distributary SUM (see kuparuk.py / CLAUDE.md "braided vs deltaic").
# The Sagavanirktok has a heavily braided/deltaic mouth near Prudhoe Bay; the per-channel
# width with total Q over-flushes salt. INTERIM: per-channel x N_CHAN_MOUTH, pending the
# exact raw SWORD summed width. Converges to the per-channel width upstream.
B_lb_perchan = 105       # SWORD per-channel width at the seaward boundary [m]
N_CHAN_MOUTH = 4         # (interim distributary count, superseded by SWORD below)
# Raw SWORD v17b delta conveyance = SUM of the seaward-most distinct distributary reaches
# near Prudhoe Bay (320 + 233 m), tools/extract_sword.py.
SWORD_MOUTH_SUM = 553    # raw SWORD distributary sum [m]
B_lb = SWORD_MOUTH_SUM or (B_lb_perchan * N_CHAN_MOUTH)
B_ub = 102        # PRISMATIC width [m] = SWORD per-channel median in the prismatic reach
                   # (>7 km, main stem, braided total / n_chan_mod). Replaces the value
                   # borrowed from Canning -- the fresh SWORD v17b per-channel data is a
                   # direct observation of the Sagavanirktok's single-channel width.
# Only L_FLARE (the flare LENGTH) is still borrowed from Canning: SWORD's Sag flare fit
# scores R2 = -0.65 (worse than a constant) over a mid-reach data gap, so the convergence
# LENGTH is not resolved. Both endpoints are now SWORD-observed: B_lb (delta sum) and
# B_ub (per-channel median) above.
L_FLARE = 7000    # flare length [m]; width converges over 0..L_FLARE, then prismatic.
WIDTH_RATIO_IS_BORROWED = True

# Depth from D = 0.280*Q^0.259 (189 surveys, 1991-2026) at open-water mean 127.6 m3/s.
# These same 189 surveys are what contradict the shipped DEPTH = 15 m, which was
# cited to THIS gauge but whose record spans 0.5-1.9 m, median 0.88 m.
DEPTH_lb = 0.98
DEPTH_ub = 0.98

distance = 1
GEOMETRY_IS_PLACEHOLDER = False
