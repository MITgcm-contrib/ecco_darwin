"""
Initialization module (translated from init.c)
"""

import math
from config import M, PI, EL, DELXI, G, Chezy_lb, Chezy_ub, DEPTH_lb, DEPTH_ub, B_ub, B_lb, distance
from variables import v, D, Dold, DEPTH, E, Y, U, TU, B, ZZ, Z, fl, dispersion, GPP, NPP, NPP_NO3, NPP_NH4, phy_death, aer_deg, denit, nit, O2_ex, Chezy, C, FCO2, Hplus
import numpy as np
from forcings_module import get_discharge, initialize_forcings
from config import SIM_START_DATETIME, Qr, DELTI, MAXT


def init():
    """Initialize model arrays and set up boundary conditions."""
    include_constantDEPTH = 0

    initialize_forcings(SIM_START_DATETIME, DELTI, MAXT)

    # Convergence length [m]
    LC = 1.0 / -(1.0 / float(EL) * math.log(float(B_ub) / float(B_lb)))

    # Van der Burgh's coefficient [-]
    K = 4.38 * DEPTH_lb**0.36 * B_lb**-0.21 * LC**-0.14

    Qr0 = -get_discharge(0, SIM_START_DATETIME)
    #print(f"[DEBUG] Qr0={Qr0}")

    # Dispersion coefficient [m^2/s]
    for i in range(M + 1):
        N = -PI * Qr0 / (DEPTH_lb * B_lb)
        D0 = 26 * math.sqrt(N * G) * DEPTH_lb**1.5  # use abs(N) to prevent domain error *maybe*
        beta = -(K * LC * Qr0) / (D0 * B_lb * DEPTH_lb)
        d = D0 * (1.0 - beta * (math.exp((i * DELXI) / LC) - 1.0))
        dispersion[i] = max(d, 0.0)

    for i in range(M + 1):
        # Hydrodynamics: Initialize arrays
        E[i] = 0.0
        Y[i] = -5.0
        Chezy[i] = Chezy_lb - (Chezy_lb - Chezy_ub) * (i - 50) / (M - 50) if i >= distance else Chezy_lb

        U[i] = 0.0
        TU[i] = 0.0
        B[i] = B_lb * math.exp(-i * float(DELXI) / float(LC)) if i > 0 else B_lb
        ZZ[i] = B[i] * (DEPTH_lb - (DEPTH_lb - DEPTH_ub) / float(EL) * i * float(DELXI))
        Z[i] = 0.0

        if include_constantDEPTH == 1:
            ZZ[i] = B[i] * DEPTH_lb
            Chezy[i] = Chezy_ub if i >= distance else Chezy_lb

        for t in range(5):
            C[i][t] = 0.0

        # FIX: Set proper initial water depths instead of zero
        H = B[i] * 2.0  # 2 meters initial water depth
        TH = B[i] * 2.0
        D[i] = H + ZZ[i]
        Dold[i] = H + ZZ[i]
        DEPTH[i] = D[i] / B[i]

        # Transport: Initialize arrays
        fl[i] = 0.0
        dispersion[i] = 0.0

        # Biogeochemistry: Initialize arrays
        GPP[i] = NPP[i] = NPP_NO3[i] = NPP_NH4[i] = phy_death[i] = 0.0
        aer_deg[i] = denit[i] = nit[i] = O2_ex[i] = FCO2[i] = 0.0



    # Boundary conditions for substances
    initialize_substance("S", "S.dat", 34.0, 0.0)
    initialize_substance("DIA", "DIA.dat", 1.0, 10.0)
    initialize_substance("dSi", "dSi.dat", 9.0, 87.0)
    initialize_substance("NO3", "NO3.dat", 132.0, 24.0)
    initialize_substance("NH4", "NH4.dat", 12.0, 4.0)
    initialize_substance("PO4", "PO4.dat", 30.0, 2.0)
    initialize_substance("O2", "O2.dat", 280.0, 280.0)
    initialize_substance("TOC", "TOC.dat", 0.0, 545.0)
    initialize_substance("SPM", "SPM.dat", 0.0, 0.01)
    initialize_substance("DIC", "DIC.dat", 2000, 1837)
    initialize_substance("ALK", "ALK.dat", 2223, 1749)
    initialize_substance("pH", "pH.dat", 8.2, 7.67)

    #print(f"[DEBUG] After init(): D[M]={D[M]}, B[M]={B[M]}, DEPTH[M]={DEPTH[M]}")
    #print(f"[DEBUG] After init(): Qr0={Qr0}, U[M]={U[M]}")

    # init_sub(substance name, filename, concentration lower bound, concentration upper bound)
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
