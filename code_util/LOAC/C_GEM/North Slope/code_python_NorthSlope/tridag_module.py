"""
Tridiagonal matrix solver module (translated from tridag.c)
"""

import math
from config import M, M1, M2, M3, DELTI, DELXI, G, RS
from variables import D, Z, C, U, B, H, TH, TU, Chezy, DEPTH

def coeff_a(t, Qr):
    """Calculate coefficients for tridiagonal matrix."""
    i = 0
    for j in range(3, M3 + 1, 2):
        i = j + 1
        Z[j] = RS * H[j] / DELTI
        C[j][1] = -D[j - 1] / (2.0 * DELXI)
        C[j][2] = RS / DELTI
        C[j][3] = D[j + 1] / (2.0 * DELXI)
        C[j][4] = 0.0

        Z[i] = 1.0 / (G * DELTI) * U[i]
        C[i][1] = -1.0 / (2.0 * DELXI * B[i - 1])
        C[i][2] = (
            1.0 / (G * DELTI)
            + (1.0 / (Chezy[i]**2)) * abs(U[i]) / DEPTH[i]
            + (U[i + 2] - U[i - 2]) / (4.0 * G * DELXI)
        )
        C[i][3] = 1.0 / (2.0 * DELXI * B[i + 1])
        C[i][4] = 0.0

    Z[2] = 1.0 / (2.0 * DELXI) * (H[1] / B[1]) + 1.0 / (G * DELTI) * U[2]
    Z[M1] = RS * TH[M1] / DELTI - Qr / (2.0 * DELXI)
    C[2][1] = 0.0
    C[2][2] = (
        1.0 / (G * DELTI)
        + (1.0 / (Chezy[2]**2)) * abs(U[2]) / DEPTH[2]
        + (U[4] - U[2]) / (2.0 * G * DELXI)
    )
    C[2][3] = 1.0 / (2.0 * DELXI * B[3])
    C[2][4] = 0.0
    C[M1][1] = -D[M2] / (2.0 * DELXI)
    C[M1][2] = RS / DELTI
    C[M1][3] = 0.0
    C[M1][4] = 0.0

def tridag():
    """Solve the tridiagonal matrix system."""
    gam = [0.0] * (M + 1)
    var = [0.0] * (M + 1)

    bet = C[2][2]
    var[2] = Z[2] / bet

    for j in range(3, M1 + 1):
        gam[j] = C[j - 1][3] / bet
        bet = C[j][2] - C[j][1] * gam[j]
        var[j] = (Z[j] - C[j][1] * var[j - 1]) / bet

    for j in range(M2, 1, -1):
        var[j] -= gam[j + 1] * var[j + 1]

    for j in range(2, M2 + 1, 2):
        TU[j] = var[j]
        TH[j + 1] = var[j + 1]

def conv(s, e, toler, xarray, yarray):
    """Test for convergence between two arrays."""
    t = 0.0

    for i in range(s, e + 1, 2):
        diff = xarray[i] - yarray[i]
        t = max(t, abs(diff))

    result = 0.0 if abs(t) >= toler else 1.0

    for i in range(s, e + 1, 2):
        yarray[i] = xarray[i]

    return result
