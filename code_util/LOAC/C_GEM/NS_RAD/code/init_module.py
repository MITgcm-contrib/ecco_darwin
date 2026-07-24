"""
Initialization module (translated from init.c)
"""

import math
from file_module import exfread
from config import (M, PI, EL, DELXI, G, Chezy_lb, Chezy_ub, DEPTH_lb, DEPTH_ub,
                    B_ub, B_lb, distance, BOUNDARIES, WIDTH_MODEL, L_FLARE, LATERAL_INFLOW,
                    MULTICHANNEL, N_CHAN_LB, N_CHAN_UP)
from variables import v, D, Dold, DEPTH, E, Y, U, TU, B, B_thread, n_chan, ZZ, Z, fl, dispersion, GPP, NPP, NPP_NO3, NPP_NH4, phy_death, aer_deg, denit, nit, O2_ex, Chezy, C, FCO2, Hplus, q_lat
import numpy as np

# Prismatic-end TOTAL conveyance width. With MULTICHANNEL off this is the site's
# B_ub unchanged; with it on, B_ub is the PER-CHANNEL median so the total across
# N_CHAN_UP braids is the conveyance the river actually offers. See config.MULTICHANNEL.
B_ub_total = float(B_ub) * N_CHAN_UP


def width_at(s):
    """
    TOTAL conveyance/surface width [m] at distance s [m] upstream of the seaward
    boundary -- summed across parallel threads where the channel braids or splits.

    WIDTH_MODEL == 'flare' (default): width converges exponentially over the first
    L_FLARE metres from B_lb down to the prismatic width B_ub_total, and is constant
    thereafter. SWORD v17b node widths support this shape over the original
    whole-domain exponential for 3 of the 4 rivers by AIC -- these are rivers with
    short delta flares, not tide-dominated funnel estuaries. See config.WIDTH_MODEL.

    WIDTH_MODEL == 'expo': the original C-GEM law, exponential over the whole domain.
    """
    if WIDTH_MODEL == "expo":
        LC = 1.0 / -(1.0 / float(EL) * math.log(B_ub_total / float(B_lb)))
        return B_lb * math.exp(-s / LC) if s > 0 else B_lb

    if s <= 0.0:
        return float(B_lb)
    if s >= L_FLARE:
        return B_ub_total
    # exponential across the flare only: B(L_FLARE) == B_ub_total exactly
    LC_flare = -float(L_FLARE) / math.log(B_ub_total / float(B_lb))
    return B_lb * math.exp(-s / LC_flare)


def n_chan_at(s):
    """
    Number of parallel threads at distance s [m] upstream of the seaward boundary:
    delta distributaries at the mouth, braids in the prismatic reach.

    Blended GEOMETRICALLY over the same L_FLARE the width uses, so that the derived
    per-thread width B/n_chan is itself a smooth exponential between its two
    endpoints (B_lb/N_CHAN_LB at the mouth, exactly B_ub upstream) rather than
    picking up a kink from mixing an exponential with a linear ramp.

    Returns 1.0 everywhere unless config.MULTICHANNEL is on.
    """
    if not MULTICHANNEL:
        return 1.0
    # blend over whatever length the width itself converges over
    L = float(EL) if WIDTH_MODEL == "expo" else float(L_FLARE)
    if s <= 0.0:
        return N_CHAN_LB
    if s >= L:
        return N_CHAN_UP
    return N_CHAN_LB * (N_CHAN_UP / N_CHAN_LB) ** (s / L)


def init():
    """Initialize model arrays and set up boundary conditions."""
    include_constantDEPTH = 0

    # (The Van der Burgh coefficient and the commented-out Savenije dispersion block
    # that used to sit here have been removed. Dispersion is now computed per
    # timestep from local width/depth/velocity by fun_module.river_dispersion --
    # see config.DISPERSION_MODEL. The K defined here was already dead code.)

    for i in range(M + 1):
        # Hydrodynamics: Initialize arrays
        E[i] = 0.0
        Y[i] = -5.0
        # Chezy ramps linearly from Chezy_lb at the seaward end (i = distance) to Chezy_ub
        # at the head (i = M). The reference was a hardcoded 50 -- a leftover from the
        # original 160 km site where `distance` was ~50 -- which, with distance=1 here, put
        # Chezy ABOVE Chezy_lb for i<50 (too little near-mouth friction). Now keyed off
        # `distance`, consistent with sed_module's tau ramps.
        Chezy[i] = Chezy_lb - (Chezy_lb - Chezy_ub) * (i - distance) / (M - distance) if i >= distance else Chezy_lb

        U[i] = 0.0
        TU[i] = 0.0
        _s = i * float(DELXI)
        B[i] = width_at(_s)
        # Per-thread width for the shear-dispersion closure only; == B[i] unless
        # MULTICHANNEL is on. Everything else (cross-section, depth, gas exchange,
        # surface area) uses the TOTAL width B[i]. See config.MULTICHANNEL.
        n_chan[i] = n_chan_at(_s)
        B_thread[i] = B[i] / n_chan[i]
        ZZ[i] = B[i] * (DEPTH_lb - (DEPTH_lb - DEPTH_ub) / float(EL) * _s)
        Z[i] = 0.0

        if include_constantDEPTH == 1:
            ZZ[i] = B[i] * DEPTH_lb
            Chezy[i] = Chezy_ub if i >= distance else Chezy_lb

        for t in range(5):
            C[i][t] = 0.0

        H = B[i] * 0.0
        TH = B[i] * 0.0
        D[i] = H + ZZ[i]
        Dold[i] = H + ZZ[i]
        DEPTH[i] = D[i] / B[i]

        # Transport: Initialize arrays
        fl[i] = 0.0
        dispersion[i] = 0.0

        # Biogeochemistry: Initialize arrays
        GPP[i] = NPP[i] = NPP_NO3[i] = NPP_NH4[i] = phy_death[i] = 0.0
        aer_deg[i] = denit[i] = nit[i] = O2_ex[i] = FCO2[i] = 0.0


    # Boundary conditions for substances, from the active site's BOUNDARIES table
    # (sites/<name>.py). clb = downstream/marine, cub = upstream/riverine.
    for name, (clb, cub) in BOUNDARIES.items():
        initialize_substance(name, f"{name}.dat", clb, cub)

    # Distributed lateral inflow [m^3 s^-1 per cell], spread uniformly over the M
    # interior cells. Zero unless the active site sets LATERAL_INFLOW (real rivers
    # leave it off). Consumed by lateral_module each timestep.
    q_lat[:] = 0.0
    if LATERAL_INFLOW > 0.0:
        for i in range(1, M + 1):
            q_lat[i] = LATERAL_INFLOW / float(M)


def initialize_substance(name, filename, clb, cub):
    """Initialize the given substance with boundary conditions."""
    v[name]["name"] = filename
    v[name]["env"] = 1
    v[name]["clb"] = clb
    v[name]["cub"] = cub

    if name == "S" or name == "SPM":
        for i in range(0, M + 1):
            v[name]["c"][i] = clb + (cub - clb) / float(M) * i
            if v[name]["c"][i] < 0:
                v[name]["c"][i] = 0
            v[name]["avg"][i] = v[name]["concflux"][i] = 0.0
            v[name]["advflux"][i] = v[name]["disflux"][i] = 0.0
    else:
        for i in range(1, M + 1):
            v[name]["c"][i] = clb + (cub - clb) / float(M) * i
            v[name]["avg"][i] = v[name]["concflux"][i] = 0.0
            v[name]["advflux"][i] = v[name]["disflux"][i] = 0.0

    # pH is not affected by advection/diffusion
    # it changes only based on S, DIC, ALK which are advected/diffused
    if name == "pH":
        v[name]["env"] = 0
