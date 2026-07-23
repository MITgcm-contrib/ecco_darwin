"""Canning River -- ~6,200 km2, western boundary of ANWR."""
from ._baseline import WATERTEMP_FILE, AMPL, pfun, BOUNDARIES as _BASE_BOUNDARIES  # noqa: F401

# --- riverine boundary chemistry: CARBONATE only (weakest-constrained site) ----
# The Canning itself has NO discrete carbonate samples. Its region (eastern North
# Slope / ANWR) is pooled instead: 40 open-water ALK grabs (28 co-located with pH)
# across 20 stations near Kaktovik (Hulahula/Sadlerochit/Marsh-Creek drainages, the
# Canning's neighbours and the same terrain its discharge is proxied from). Median
# ALK 1797 mmol/m^3 and median PAIRED water pCO2 613 uatm (SUPERSATURATED -> outgas);
# DIC back-solved -> 1804. All other species keep the shared placeholder -- consistent
# with this being the site with reconstructed discharge and borrowed temperature too.
# tools/usgs_carbonate_boundary.py bbox -146.5,69.0,-143.5,70.6.
BOUNDARIES = dict(_BASE_BOUNDARIES)
for _sp, _cub in [("pH", 7.97), ("ALK", 1796.6), ("DIC", 1803.7)]:
    BOUNDARIES[_sp] = (_BASE_BOUNDARIES[_sp][0], _cub)
BOUNDARY_CHEM_SOURCE = "WQP eastern-ANWR regional pool (paired-pCO2 DIC); other species placeholder"

LABEL = "Canning"
DISCHARGE_FILE = "canning_river_discharge_2022_m3sec.csv"
# Wind-driven storm surge at the marine boundary. Prudhoe Bay (9497645) is the only
# continuous 2022 water-level record on this coast, used as the REGIONAL proxy for all
# four rivers (see tools/build_surge.py); added to the harmonic tide in fun_module.Tide.
SURGE_FILE = "surge_prudhoe_2022_m.csv"

# NO 2022 DISCHARGE OBSERVATIONS. USGS 15955000 ran only 2008-06-23 to 2012-09-30.
# Reconstructed from Hulahula R (15980000) scaled by 2.971, the ratio of means over
# their 731-day common record (log-space r = 0.87). Rebuild: tools/fetch_discharge.py.
GAUGE = None
DISCHARGE_IS_RECONSTRUCTED = True
DISCHARGE_DONOR_GAUGE = "15980000"
DISCHARGE_IS_UPSTREAM_PROXY = False

# --- geometry: SWORD v17c nodes + USGS channel surveys -----------------------
# Width from SWORD per-channel node width, 136 nodes. Single-channel river, so the
# braiding correction barely applies here. Its 0.55 convergence ratio is also the
# donor for the Sagavanirktok.
EL = 27175
# DELTA-MOUTH WIDTH = raw SWORD v17b distributary SUM at the Camden Bay delta
# (562 + 236 m across the two seaward-most channels, tools/extract_sword.py). The fresh
# SWORD extraction shows a wide two-channel delta mouth, not the 272 m the earlier
# one-off extraction gave (which had traced only the main channel), so it is widened.
B_lb_perchan = 272      # earlier single-channel seaward width [m]
SWORD_MOUTH_SUM = 797   # raw SWORD distributary sum [m]
B_lb = SWORD_MOUTH_SUM
B_ub = 64          # PRISMATIC width [m] = SWORD per-channel median in the prismatic reach
                    # (>7 km, main stem, braided total / n_chan_mod), SWORD v17b. Replaces
                    # the earlier 132 m, which sat above the observed per-channel width.
L_FLARE = 7000    # flare length [m]; width converges over 0..L_FLARE,
                   # then is prismatic. Fit R2 = 0.347 (vs exponential-over-27km).
# Depth from D = 0.224*Q^0.325 at open-water mean 138.3 m3/s. Note this rests on only
# 28 surveys (2008-2012) -- by far the thinnest survey record of the four -- though
# the fit itself is the best (R2 = 0.89). So Canning has observed GEOMETRY but
# reconstructed DISCHARGE and borrowed TEMPERATURE: the weakest-constrained site.
DEPTH_lb = 1.11
DEPTH_ub = 1.11

distance = 1
GEOMETRY_IS_PLACEHOLDER = False
