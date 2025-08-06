"""
Update hydrodynamics module (modernized from uphyd.c)
"""
import numpy as np
from config import M, M1, M2, M3, Qr
from variables import B, H, TH, D, ZZ, DEPTH, U, TU, Dold
from fun_module import Tide


def new_bc(t):
    """Set new boundary conditions at the start of time step."""
    H[1] = B[1] * Tide(t)
    TH[1] = H[1]
    D[1] = H[1] + ZZ[1]
    DEPTH[1] = D[1] / max(B[1], 1e-6)

    # Outlet velocity and transport (assume flow is imposed via Qr)
    D[M] = max(D[M], 1e-6)
    U[M] = Qr / D[M]
    TU[M] = U[M]


def update():
    """Vectorized update of TH, D, DEPTH, and TU arrays after solving system."""
    # Even indices: i = 2, 4, ..., M2 (update TH, D, DEPTH at U points)
    even = np.arange(2, M2 + 1, 2)
    TH[even] = 0.5 * (TH[even - 1] + TH[even + 1])

    tmp_left = TH[even - 1] + ZZ[even - 1]
    tmp_right = TH[even + 1] + ZZ[even + 1]
    D[even] = 0.5 * (tmp_left + tmp_right)

    B_safe_left = np.maximum(B[even - 1], 1e-6)
    B_safe_right = np.maximum(B[even + 1], 1e-6)
    DEPTH[even] = 0.5 * ((tmp_left / B_safe_left) + (tmp_right / B_safe_right))

    # Odd indices: i = 3, 5, ..., M1 (update D, DEPTH, TU at H points)
    odd = np.arange(3, M1 + 1, 2)
    D[odd] = TH[odd] + ZZ[odd]
    DEPTH[odd] = D[odd] / np.maximum(B[odd], 1e-6)
    TU[odd] = 0.5 * (TU[odd - 1] + TU[odd + 1])

    # Outlet: i = M
    TH_M1 = TH[M1]
    TH_M3 = TH[M3]
    ZZ_M1 = ZZ[M1]
    ZZ_M3 = ZZ[M3]

    D[M] = 0.5 * (3.0 * (TH_M1 + ZZ_M1) - (TH_M3 + ZZ_M3))
    DEPTH[M] = D[M] / max(B[M], 1e-6)
    TH[M] = 0.5 * (3.0 * TH_M1 - TH_M3)

    # Inlet extrapolation
    TU[1] = TU[2]

def new_uh():
    """Finalize velocity (U) and surface height (H) fields."""
    U[1:M + 1] = np.where(TU[1:M + 1] == 0.0, 1e-4, TU[1:M + 1])
    H[1:M + 1] = TH[1:M + 1]
    Dold[1:M + 1] = D[1:M + 1]
