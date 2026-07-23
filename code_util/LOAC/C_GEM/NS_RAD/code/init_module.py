"""
Initialization module (translated from init.c)
"""

import math
from file_module import exfread
from config import (M, PI, EL, DELXI, G, Chezy_lb, Chezy_ub, DEPTH_lb, DEPTH_ub,
                    B_ub, B_lb, distance, BOUNDARIES, WIDTH_MODEL, L_FLARE, LATERAL_INFLOW)
from variables import v, D, Dold, DEPTH, E, Y, U, TU, B, ZZ, Z, fl, dispersion, GPP, NPP, NPP_NO3, NPP_NH4, phy_death, aer_deg, denit, nit, O2_ex, Chezy, C, FCO2, Hplus, q_lat
import numpy as np

def width_at(s):
    """
    Channel width [m] at distance s [m] upstream of the seaward boundary.

    WIDTH_MODEL == 'flare' (default): width converges exponentially over the first
    L_FLARE metres from B_lb down to the prismatic width B_ub, and is constant at
    B_ub thereafter. SWORD v17c node widths support this shape over the original
    whole-domain exponential for 3 of the 4 rivers by AIC -- these are rivers with
    short delta flares, not tide-dominated funnel estuaries. See config.WIDTH_MODEL.

    WIDTH_MODEL == 'expo': the original C-GEM law, exponential over the whole domain.
    """
    if WIDTH_MODEL == "expo":
        LC = 1.0 / -(1.0 / float(EL) * math.log(float(B_ub) / float(B_lb)))
        return B_lb * math.exp(-s / LC) if s > 0 else B_lb

    if s <= 0.0:
        return float(B_lb)
    if s >= L_FLARE:
        return float(B_ub)
    # exponential across the flare only: B(L_FLARE) == B_ub exactly
    LC_flare = -float(L_FLARE) / math.log(float(B_ub) / float(B_lb))
    return B_lb * math.exp(-s / LC_flare)


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
        B[i] = width_at(i * float(DELXI))
        ZZ[i] = B[i] * (DEPTH_lb - (DEPTH_lb - DEPTH_ub) / float(EL) * i * float(DELXI))
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
