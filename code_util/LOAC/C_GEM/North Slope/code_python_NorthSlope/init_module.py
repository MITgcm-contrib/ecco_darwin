"""
Initialization module (translated from init.c)
"""

import math
from file_module import exfread
from config import M, PI, EL, DELXI, G, Chezy_lb, Chezy_ub, DEPTH_lb, DEPTH_ub, B_ub, B_lb, distance
from variables import v, D, Dold, DEPTH, E, Y, U, TU, B, ZZ, Z, fl, dispersion, GPP, NPP, NPP_NO3, NPP_NH4, phy_death, aer_deg, denit, nit, O2_ex, Chezy, C, FCO2, Hplus
import numpy as np

def init():
    """Initialize model arrays and set up boundary conditions."""
    include_constantDEPTH = 0

    # Convergence length [m]
    LC = 1.0 / -(1.0 / float(EL) * math.log(float(B_ub) / float(B_lb)))

    # Van der Burgh's coefficient [-]
    K = 4.38 * DEPTH_lb**0.36 * B_lb**-0.21 * LC**-0.14

    # Dispersion coefficient [m^2/s]
    #for i in range(M + 1):
        #N = -PI * Qr / (DEPTH_lb * B_lb)
        #D0 = 26 * math.sqrt(N * G) * DEPTH_lb**1.5
        #beta = -(K * LC * Qr) / (D0 * B_lb * DEPTH_lb)
        #d = D0 * (1.0 - beta * (math.exp((i * DELXI) / LC) - 1.0))
        #dispersion[i] = max(d, 0.0)

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


    # Boundary conditions for substances
    initialize_substance("S", "S.dat", 27.1, 0.0)
    initialize_substance("DIA", "DIA.dat", 1.6, 1.1)
    initialize_substance("dSi", "dSi.dat", 2.35, 174.6)
    initialize_substance("NO3", "NO3.dat", 0.03, 7.68)
    initialize_substance("NH4", "NH4.dat", 0.01, 3.9)#/(7.68/0.03)
    initialize_substance("PO4", "PO4.dat", 0.56, 0.01)
    initialize_substance("O2", "O2.dat", 392.0, 367.0)
    initialize_substance("TOC", "TOC.dat", 69.0, 1582.0)
    initialize_substance("SPM", "SPM.dat", 0.15, 2.0)
    initialize_substance("DIC", "DIC.dat", 1869, 1717)
    initialize_substance("ALK", "ALK.dat", 1957, 1596)
    initialize_substance("pH", "pH.dat", 8.1, 7.5)


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
