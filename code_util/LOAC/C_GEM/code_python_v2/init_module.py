"""
Initialization module (translated from init.c)
"""
import math
from config import M, PI, EL, DELXI, G, Qr, Chezy_lb, Chezy_ub, DEPTH_lb, DEPTH_ub, distance, B_ub, B_lb
from variables import v, Dold, D, DEPTH, TH, E, Y, U, TU, B, ZZ, Z, fl, dispersion, GPP, NPP, NPP_NO3, NPP_NH4, phy_death, aer_deg, denit, nit, O2_ex, Chezy, C
import numpy as np

def calculate_convergence_length(EL, B_ub, B_lb):
    """
    Calculate the convergence length LC [m].

    Parameters:
    - EL: length scale (float)
    - B_ub: width at upper boundary (float)
    - B_lb: width at lower boundary (float)

    Returns:
    - LC: convergence length (float)
    """
    return 1.0 / (-(1.0 / float(EL)) * math.log(float(B_ub) / float(B_lb)))

def calculate_dispersion():
    """Calculate dispersion coefficient array (m^2/s) with vectorized operations."""
    LC = calculate_convergence_length(EL, B_ub, B_lb)
    # Van der Burgh's coefficient [-]
    K = 4.38 * DEPTH_lb**0.36 * B_lb**-0.21 * LC**-0.14
    # Compute constants
    N = -PI * Qr / (DEPTH_lb * B_lb)
    D0 = 26.0 * math.sqrt(N * G) * DEPTH_lb**1.5
    beta = -(K * LC * Qr) / (D0 * B_lb * DEPTH_lb)

    # Create index array
    i_arr = np.arange(M + 1)

    # Compute dispersion with exponential term
    exp_term = np.exp((i_arr * DELXI) / LC) - 1.0
    d = D0 * (1.0 - beta * exp_term)

    # Ensure non-negative dispersion values
    dispersion[:] = np.maximum(d, 0.0)
    print(dispersion)

def init():
    """Initialize model arrays and set up boundary conditions."""
    include_constantDEPTH = 0

    # Convergence length [m]
    LC = calculate_convergence_length(EL, B_ub, B_lb)
    # Dispersion coefficient [m^2/s]
    calculate_dispersion()

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
        aer_deg[i] = denit[i] = nit[i] = O2_ex[i] = 0.0

    # Boundary conditions for substances
    initialize_substance("S", "S.dat", 30.3, 0.0)
    initialize_substance("DIA", "DIA.dat", 10.0, 50.0)
    initialize_substance("dSi", "dSi.dat", 10.0, 250.0)
    initialize_substance("NO3", "NO3.dat", 50.0, 198.0)
    initialize_substance("NH4", "NH4.dat", 0.1, 520.0)
    initialize_substance("PO4", "PO4.dat", 0.1, 17.0)
    initialize_substance("O2", "O2.dat", 250.0, 106.0)
    initialize_substance("TOC", "TOC.dat", 0.1, 393.0)
    initialize_substance("SPM", "SPM.dat", 0.03, 0.07)

def initialize_substance(name, filename, clb, cub):
    """Initialize the given substance with boundary conditions."""
    v[name]["name"] = filename
    v[name]["env"] = 1
    v[name]["clb"] = clb
    v[name]["cub"] = cub

    indices = np.arange(M + 1)
    profile = clb + (cub - clb) * indices / float(M)

    if name in ("S", "SPM"):
        # Ensure no negative values
        v[name]["c"][:] = np.maximum(profile, 0.0)
    else:
        # Set values only from 1 to M
        v[name]["c"][1:] = profile[1:]

    # Zero-initialize other fields
    v[name]["avg"][:] = 0.0
    v[name]["concflux"][:] = 0.0
    v[name]["advflux"][:] = 0.0
    v[name]["disflux"][:] = 0.0

