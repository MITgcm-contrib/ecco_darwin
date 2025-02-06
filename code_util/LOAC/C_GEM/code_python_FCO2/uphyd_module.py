"""
Update hydrodynamics module (translated from uphyd.c)
"""

from config import M, M1, M2, M3, Qr
from variables import B, H, TH, D, ZZ, DEPTH, U, TU, Dold
from fun_module import Tide

def new_bc(t):
    """Set new boundary conditions."""
    H[1] = B[1] * Tide(t)
    TH[1] = H[1]
    D[1] = H[1] + ZZ[1]
    DEPTH[1] = D[1] / B[1]
    U[M] = Qr / D[M]
    TU[M] = U[M]

def update():
    """Update TH, D, DEPTH, and TU during iteration."""
    tmp = [0.0] * (M + 1)

    # Update even-indexed points
    for i in range(2, M2 + 1, 2):
        tmp[i - 1] = TH[i - 1] + ZZ[i - 1]
        tmp[i + 1] = TH[i + 1] + ZZ[i + 1]
        TH[i] = (TH[i - 1] + TH[i + 1]) / 2.0
        D[i] = (tmp[i - 1] + tmp[i + 1]) / 2.0
        DEPTH[i] = ((tmp[i - 1] / B[i - 1]) + (tmp[i + 1] / B[i + 1])) / 2.0

    # Update odd-indexed points
    for i in range(3, M1 + 1, 2):
        D[i] = TH[i] + ZZ[i]
        DEPTH[i] = D[i] / B[i]
        TU[i] = (TU[i + 1] + TU[i - 1]) / 2.0

    # Update the last point
    D[M] = (3.0 * (TH[M1] + ZZ[M1]) - (TH[M3] + ZZ[M3])) / 2.0
    DEPTH[M] = D[M] / B[M]
    TH[M] = (3.0 * TH[M1] - TH[M3]) / 2.0
    TU[1] = TU[2]

def new_uh():
    """Update H and U arrays."""
    for i in range(1, M + 1):
        U[i] = TU[i]
        if U[i] == 0:
            U[i] = 0.0001  # Avoid zero velocity
        H[i] = TH[i]
        Dold[i] = D[i]
