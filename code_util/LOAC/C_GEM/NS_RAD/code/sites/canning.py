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

# --- O2 and SPM (cub): USGS Water Quality Portal, same bbox as above -------------------
# Open-water (Jun-Sep) medians. O2: n=24 discrete DO samples, 11.6 mg/L -> 362.5 mmol/m^3
# (32 g/mol). SPM: n=9 Suspended Sediment Concentration samples, 4.0 mg/L -> 0.0040 g/L --
# two orders of magnitude below the shared placeholder (2.0 g/L); the lowest sediment load
# of the four, consistent with a mountain-fed, single-thread, less-braided system.
for _sp, _cub in [("O2", 362.50), ("SPM", 0.0040)]:
    BOUNDARIES[_sp] = (_BASE_BOUNDARIES[_sp][0], _cub)

# --- DIA (cub): epilithic->water-column proxy, NOT a phytoplankton measurement ----------
# North Slope rivers have essentially no water-column phytoplankton monitoring (checked
# WQP + literature) -- what IS measured is EPILITHIC chlorophyll (algae attached to the
# streambed), a different compartment than DIA (this model's transported, advected tracer).
# Kuparuk reference (unfertilized) reach epilithic chlorophyll = 2.8 mg Chl/m2 on rock (85%
# of bed cover), Slavik et al. 2004 (Ecology 85:939-954) Table 2, July 1998 -- the only
# quantitative North Slope epilithic value found, applied here for lack of a site-specific
# one (consistent with Canning already being the weakest-constrained site). Converted to a
# rough water-column-equivalent by DIVIDING BY THIS SITE'S DEPTH (as if the entire bed
# stock were resuspended and evenly mixed through the water column -- NOT how these
# compartments actually behave; a proxy, not a measurement), then to carbon via C:Chl = 75
# gC/gChla (user-specified; close to the model's own Chla2CMIN-implied maximum of
# 1/0.0125 = 80 gC/gChla) and to mmol via /12.011 g/mol:
#   DIA = (2.8 / 1.11) * 75 / 12.011 = 15.751 mmol C/m^3   (this site's DEPTH_ub = 1.11 m,
#   set below in the geometry section -- computed here as a literal since DEPTH_ub isn't
#   defined yet at this point in the file)
BOUNDARIES["DIA"] = (_BASE_BOUNDARIES["DIA"][0], 15.751)

# --- MARINE (clb): ECCO-Darwin v5, nearest wet LLC270 cell to the delta mouth ---------
# Cell (70.203N, -145.916E), ~8.5 km offshore of the SWORD mouth point (70.13N, -145.85W).
# Climatological annual mean over the full archived record (272-292 months depending on
# variable, native LLC270 monthly output, public NASA data.nas.nasa.gov/ecco portal, no
# auth). S is SALTanom + 35 (MITgcm's standard salinity-anomaly convention). DIC/ALK/NO3/
# NH4/PO4/dSi/O2/TOC units match the model's mmol m^-3 directly, no conversion. Built by
# scratch/ecco_darwin/extract_all_rivers.py (not tracked -- rerun to reproduce or extend).
# Unlike the other three sites, NO3/NH4/PO4 here still carry the SHARED PLACEHOLDER cub
# (Canning has no discrete nutrient record of its own -- see the module docstring), so
# this is the first time those two species get any site-specific value for Canning at all.
# NOT sourced this way: DIA (ECCO-Darwin gives Chl1-5, not carbon biomass -- would need an
# uncertain C:Chl ratio), pH (not in the archived output; could instead be diagnosed from
# this DIC/ALK/S pair via the model's own carbonate solve), RDOC/CH4/N2O/SPM (no analog).
for _sp, _clb in [("S", 29.8775), ("DIC", 2009.8823), ("ALK", 2107.3902), ("NO3", 5.8853),
                   ("NH4", 0.0274), ("PO4", 0.8892), ("dSi", 7.8052), ("O2", 375.6487),
                   ("TOC", 27.2495)]:
    BOUNDARIES[_sp] = (_clb, BOUNDARIES[_sp][1])   # override marine clb, keep river cub

BOUNDARY_CHEM_SOURCE = ("WQP eastern-ANWR regional pool (paired-pCO2 DIC); other species placeholder; "
                        "marine (clb): ECCO-Darwin v5 nearest-cell climatology")

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

# --- geometry: SWORD v17b nodes + USGS channel surveys -----------------------
# Width from SWORD per-channel node width, 136 nodes. Single-channel river, so the
# braiding correction barely applies here. Its 0.55 convergence ratio is also the
# donor for the Sagavanirktok.
EL = 0            # PLACEHOLDER pending an observed value -- see module docstring caveat.
                  # M = int(EL/DELXI)+1, forced even (config.py) -> M=0 for this site: the
                  # grid degenerates and Canning will NOT run until this is set to a real length.
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

# Multi-channel geometry for CGEM_MULTICHANNEL=on (config.MULTICHANNEL). B_UB_TOTAL is
# the OBSERVED raw braided total beyond the flare; config derives the thread count from it. Two seaward
# distributaries (SWORD, 562+236 = B_lb); single-channel river upstream, so N_CHAN_UP = 1
# and this flare is a real convergence.
N_CHAN_LB = 2.0
B_UB_TOTAL = 63.0      # SWORD v17b raw prismatic median [m]; IQR 53-78, 68 nodes.
                        # 63 vs B_ub 64 -> thread count 0.98, i.e. single-thread (median
                        # n_chan_mod = 1.0). Clamped to 1.0 in config; real flare.

distance = 1
GEOMETRY_IS_PLACEHOLDER = False
