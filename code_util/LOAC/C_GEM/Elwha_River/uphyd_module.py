"""
Update hydrodynamics module (translated from uphyd.c)
"""

from config import M, M1, M2, M3, Qr, SIM_START_DATETIME
from variables import B, H, TH, D, ZZ, DEPTH, U, TU, Dold
from fun_module import Tide
from forcings_module import get_discharge

def new_bc(t):
    """Set new boundary conditions."""
    dynamic_Qr = -get_discharge(t, SIM_START_DATETIME)
    #print(f"[DEBUG] new_bc: t={t}, dynamic_Qr={dynamic_Qr}")

    H[1] = B[1] * Tide(t)
    TH[1] = H[1]
    D[1] = H[1] + ZZ[1]
    DEPTH[1] = D[1] / B[1]
    U[M] = dynamic_Qr / D[M]
    TU[M] = U[M]
    #print("DEPTH id after new_bc():", id(DEPTH), "first 5 values:", DEPTH[:5])

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

    # Debug: Check for zero or minimum depth after update
    #for idx, d in enumerate(DEPTH):
    #    if d == 0:
    #        print(f"Zero depth detected at index {idx} in update(), B={B[idx]}, D={D[idx]}, TH={TH[idx]}, ZZ={ZZ[idx]}")
    #print(f"Minimum depth after update: {min(DEPTH)}")
    #print("DEPTH id after update():", id(DEPTH), "first 5 values:", DEPTH[:5])

def new_uh():
    """Update H and U arrays."""
    for i in range(1, M + 1):
        U[i] = TU[i]
        if U[i] == 0:
            U[i] = 0.0001  # Avoid zero velocity
        H[i] = TH[i]
        Dold[i] = D[i]

#print("DEPTH values at start of coeff_a:", DEPTH)
