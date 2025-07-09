"""
Initialization module (translated from init.c)
"""

import math
from config import M, PI, EL, DELXI, G, Qr, Chezy_lb, Chezy_ub, DEPTH_lb, DEPTH_ub, B_ub, B_lb, distance
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
    for i in range(M + 1):
        N = -PI * Qr / (DEPTH_lb * B_lb)
        D0 = 26 * math.sqrt(N * G) * DEPTH_lb**1.5
        beta = -(K * LC * Qr) / (D0 * B_lb * DEPTH_lb)
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
    initialize_substance("S", "S.dat", 27.55, 6.05) # Twilley et al
    initialize_substance("DIA", "DIA.dat", 43.0, 130.0) #Cifuentes et al PC:Chl a ratio=104 (coastal value for limited influence of OC from terrestrial plants (mangroves)) for converting chla from Twilley et al -> 520 mg C m-3 -> 43 [mmol C m-3] 
    initialize_substance("dSi", "dSi.dat", 83.0, 179.5) # Twilley et al
    initialize_substance("NO3", "NO3.dat", 5.25, 17.6) # Twilley et al [mmol m^-3]
    initialize_substance("NH4", "NH4.dat", 1.8, 3.8) # Twilley et al 1 µM = 1 µmol L-1 = 0.001 mmol L-1 = 1 mmol m-3 [mmol m^-3]
    initialize_substance("PO4", "PO4.dat", 1.75, 2.8) # Twilley et al
    initialize_substance("O2", "O2.dat", 203.125, 167.4) # Twilley et al [g/l] 3.7 ml L-1 / 22400 ml per mol (molar volume of an ideal gas)-> mol L-1 * 1000  * 1000 -> [mmol m^-3]
    initialize_substance("TOC", "TOC.dat", 92.95, 278.75) #Cifuentes et al PC:Chl a ratio=104 (coastal value for limited influence of OC from terrestrial plants (mangroves)) for converting chla from Twilley et al and then substracted from chla conversion to PC with ratio for Guayas river (327) [mmol C m-3] 
    initialize_substance("SPM", "SPM.dat", 0.214, 0.475) # Twilley et al [g/l]
    initialize_substance("DIC", "DIC.dat", 1806.91, 1278.0) # BElliard et al [mmol m^-3]
    initialize_substance("ALK", "ALK.dat", 1798.0, 1100.25) # BElliard et al [mmol m^-3]
    initialize_substance("pH", "pH.dat", 8.0, 7.88) # BElliard et al


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