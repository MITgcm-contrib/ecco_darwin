import numpy as np
from numba import njit
from config import (
    M, PI, EL, DELXI, G, Qr, Chezy_lb, Chezy_ub,
    DEPTH_lb, DEPTH_ub, distance, B_ub, B_lb, USE_CO2_FLUX
)
from variables import (
    v, Dold, D, DEPTH, TH, E, Y, U, TU, B, ZZ, Z, fl, dispersion,
    GPP, NPP, NPP_NO3, NPP_NH4, phy_death, aer_deg, denit, nit, O2_ex,
    Chezy, C
)


@njit
def calculate_convergence_length(EL, B_ub, B_lb):
    return 1.0 / (-(1.0 / EL) * np.log(B_ub / B_lb))


@njit
def calculate_dispersion(EL, B_ub, B_lb, DELXI, Qr, G, DEPTH_lb, B_lb_val, dispersion_out):
    LC = calculate_convergence_length(EL, B_ub, B_lb)
    K = 4.32 * (DEPTH_lb**0.36 / (B_lb_val**0.21 * LC**0.14))
    N = -PI * Qr / (DEPTH_lb * B_lb_val)
    D0 = 26.0 * np.sqrt(N * G) * DEPTH_lb**1.5
    beta = -(K * LC * Qr) / (D0 * B_lb_val * DEPTH_lb)

    for i in range(M + 1):
        exp_term = np.exp((i * DELXI) / LC) - 1.0
        d = D0 * (1.0 - beta * exp_term)
        dispersion_out[i] = max(d, 0.0)

@njit
def init_hydro(DELXI, EL, LC, distance, include_constantDEPTH,
                      E, Y, Chezy, U, TU, B, ZZ, Z, C, TH, D, Dold, DEPTH, fl):
    i = np.arange(M + 1)

    E[:] = 0.0
    Y[:] = -5.0
    Chezy[:] = Chezy_lb
    for idx in range(M + 1):
        if idx >= distance:
            Chezy[idx] = Chezy_lb - (Chezy_lb - Chezy_ub) * (idx - 50) / (M - 50)

    U[:] = 0.0
    TU[:] = 0.0

    for idx in range(M + 1):
        B[idx] = B_lb * np.exp(-idx * DELXI / LC)
        ZZ[idx] = B[idx] * (DEPTH_lb - (DEPTH_lb - DEPTH_ub) * idx * DELXI / EL)

    if include_constantDEPTH == 1:
        for idx in range(M + 1):
            ZZ[idx] = B[idx] * DEPTH_lb
            Chezy[idx] = Chezy_ub if idx >= distance else Chezy_lb

    Z[:] = 0.0
    C[:] = 0.0
    H = np.zeros_like(B)
    TH[:] = B * 0.0
    D[:] = H + ZZ
    Dold[:] = D
    for idx in range(M + 1):
        if B[idx] != 0:
            DEPTH[idx] = D[idx] / B[idx]
        else:
            DEPTH[idx] = 0.0
    fl[:] = 0.0

def init():
    include_constantDEPTH = 0

    # Compute convergence length
    LC = calculate_convergence_length(EL, B_ub, B_lb)

    # Dispersion array
    calculate_dispersion(EL, B_ub, B_lb, DELXI, Qr, G, DEPTH_lb, B_lb, dispersion)

    # Initialize hydrodynamic fields
    init_hydro(DELXI, EL, LC, distance, include_constantDEPTH,
                      E, Y, Chezy, U, TU, B, ZZ, Z, C, TH, D, Dold, DEPTH, fl)

    # Zero out biogeochemical arrays
    for arr in [GPP, NPP, NPP_NO3, NPP_NH4, phy_death, aer_deg, denit, nit, O2_ex]:
        arr[:] = 0.0

    # Init all tracers
    initialize_substance("S", "S.dat", 34.0, 0.0)
    initialize_substance("DIA", "DIA.dat", 1.0, 10.0)
    initialize_substance("dSi", "dSi.dat", 9.0, 87.0)
    initialize_substance("NO3", "NO3.dat", 5.0, 72.0)
    initialize_substance("NH4", "NH4.dat", 1.0, 18.0)
    initialize_substance("PO4", "PO4.dat", 1.0, 3.0)
    initialize_substance("O2", "O2.dat", 280.0, 280.0)
    initialize_substance("TOC", "TOC.dat", 0.0, 545.0)
    initialize_substance("SPM", "SPM.dat", 0.0, 0.01)

    if USE_CO2_FLUX:
        initialize_substance("DIC", "DIC.dat", 2000, 1837)
        initialize_substance("ALK", "ALK.dat", 2223, 1749)
        initialize_substance("pH", "pH.dat", 8.2, 7.67)

def initialize_substance(name, filename, clb, cub):
    v[name]["name"] = filename
    v[name]["env"] = 1
    v[name]["clb"] = clb
    v[name]["cub"] = cub

    i = np.arange(M + 1)
    j = np.arange(1,M + 1,1)

    if name in ("S", "SPM"):
        v[name]["c"][:] = clb + (cub - clb) * i / M
        v[name]["avg"][:] = 0.0
        v[name]["concflux"][:] = 0.0
        v[name]["advflux"][:] = 0.0
        v[name]["disflux"][:] = 0.0
    else:
        v[name]["c"][1:] = clb + (cub - clb) * j / M
        v[name]["avg"][1:] = 0.0
        v[name]["concflux"][1:] = 0.0
        v[name]["advflux"][1:] = 0.0
        v[name]["disflux"][1:] = 0.0

    v[name]["avg"][:] = 0.0
    v[name]["concflux"][:] = 0.0
    v[name]["advflux"][:] = 0.0
    v[name]["disflux"][:] = 0.0