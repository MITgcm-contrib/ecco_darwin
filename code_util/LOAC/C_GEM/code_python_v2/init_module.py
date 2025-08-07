import math
import numpy as np
from config import M, PI, EL, DELXI, G, Qr, Chezy_lb, Chezy_ub, DEPTH_lb, DEPTH_ub, distance, B_ub, B_lb
from variables import (
    v, Dold, D, DEPTH, TH, E, Y, U, TU, B, ZZ, Z, fl, dispersion,
    GPP, NPP, NPP_NO3, NPP_NH4, phy_death, aer_deg, denit, nit, O2_ex,
    Chezy, C
)

def calculate_convergence_length(EL, B_ub, B_lb):
    return 1.0 / (-(1.0 / EL) * math.log(B_ub / B_lb))

def calculate_dispersion():
    LC = calculate_convergence_length(EL, B_ub, B_lb)
    K = 4.32 * (DEPTH_lb**0.36 / (B_lb**0.21 * LC**0.14))
    N = -PI * Qr / (DEPTH_lb * B_lb)
    D0 = 26.0 * math.sqrt(N * G) * DEPTH_lb**1.5
    beta = -(K * LC * Qr) / (D0 * B_lb * DEPTH_lb)

    i_arr = np.arange(M + 1)
    exp_term = np.exp((i_arr * DELXI) / LC) - 1.0
    d = D0 * (1.0 - beta * exp_term)

    dispersion[:] = np.maximum(d, 0.0)

def init():
    include_constantDEPTH = 0
    LC = calculate_convergence_length(EL, B_ub, B_lb)
    calculate_dispersion()

    i = np.arange(M + 1)

    # Hydrodynamics
    E[:] = 0.0
    Y[:] = -5.0
    Chezy[:] = Chezy_lb
    Chezy[i >= distance] = Chezy_lb - (Chezy_lb - Chezy_ub) * (i[i >= distance] - 50) / (M - 50)

    U[:] = 0.0
    TU[:] = 0.0
    B[:] = B_lb * np.exp(-i * DELXI / LC)
    ZZ[:] = B * (DEPTH_lb - (DEPTH_lb - DEPTH_ub) * i * DELXI / EL)
    Z[:] = 0.0

    if include_constantDEPTH == 1:
        ZZ[:] = B * DEPTH_lb
        Chezy[:] = Chezy_lb
        Chezy[i >= distance] = Chezy_ub

    C[:] = 0.0

    H = np.zeros_like(B)
    TH[:] = B * 0.0
    D[:] = H + ZZ
    Dold[:] = D.copy()
    DEPTH[:] = D / B

    fl[:] = 0.0

    # Biogeochemistry
    GPP[:] = 0.0
    NPP[:] = 0.0
    NPP_NO3[:] = 0.0
    NPP_NH4[:] = 0.0
    phy_death[:] = 0.0
    aer_deg[:] = 0.0
    denit[:] = 0.0
    nit[:] = 0.0
    O2_ex[:] = 0.0

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
    v[name]["name"] = filename
    v[name]["env"] = 1
    v[name]["clb"] = clb
    v[name]["cub"] = cub

    i = np.arange(M + 1)
    profile = clb + (cub - clb) * i / M

    if name in ("S", "SPM"):
        v[name]["c"][:] = np.maximum(profile, 0.0)
    else:
        v[name]["c"][1:] = profile[1:]

    v[name]["avg"][:] = 0.0
    v[name]["concflux"][:] = 0.0
    v[name]["advflux"][:] = 0.0
    v[name]["disflux"][:] = 0.0