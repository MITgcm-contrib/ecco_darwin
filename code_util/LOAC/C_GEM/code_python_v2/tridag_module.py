"""
Tridiagonal matrix solver module (translated from tridag.c)
"""
import math
import numpy as np
from config import M, M1, M2, M3, DELTI, DELXI, G, RS, Qr, TH_ABS_FLOOR, TU_ABS_FLOOR, TH_REL, TU_REL, LOOSE_CAP
from variables import D, Z, C, U, B, H, TH, TU, Chezy, DEPTH
from scipy.linalg import solve_banded
from scipy.linalg.lapack import get_lapack_funcs
from numba import njit

def build_tridiag_system():
    """
    Build the banded matrix (ab) and RHS for scipy.linalg.solve_banded.
    Only includes rows [2:M1] which are the active equations.
    """
    i_start = 2
    i_end = M1
    n = i_end - i_start + 1  # Number of active equations

    lower = np.zeros(n)
    diag = np.zeros(n)
    upper = np.zeros(n)
    rhs = np.zeros(n)

    for idx, i in enumerate(range(i_start, i_end + 1)):
        if i % 2 == 1:
            rhs[idx] = RS * H[i] / DELTI
            lower[idx] = -D[i - 1] / (2.0 * DELXI)
            diag[idx] = RS / DELTI
            upper[idx] = D[i + 1] / (2.0 * DELXI)
        else:
            u = U[i]
            chezy_sq = max(Chezy[i], 1e-6)**2
            depth_i = max(DEPTH[i], 1e-6)
            B_left = max(B[i - 1], 1e-6)
            B_right = max(B[i + 1], 1e-6)

            rhs[idx] = u / (G * DELTI)
            lower[idx] = -1.0 / (2.0 * DELXI * B_left)
            diag[idx] = (
                1.0 / (G * DELTI)
                + abs(u) / (chezy_sq * depth_i)
                + (U[i + 2] - U[i - 2]) / (4.0 * G * DELXI)
            )
            upper[idx] = 1.0 / (2.0 * DELXI * B_right)
    # after the loop that fills lower/diag/upper
    lower[0] = 0.0  # matches C[2][1] = 0
    upper[-1] = 0.0  # matches C[M1][3] = 0
    # Boundary overrides
    diag[0] = (
        1.0 / (G * DELTI)
        + abs(U[2]) / (max(Chezy[2], 1e-6)**2 * max(DEPTH[2], 1e-6))
        + (U[4] - U[2]) / (2.0 * G * DELXI)
    )
    rhs[0] = H[1] / (2.0 * DELXI * B[1]) + U[2] / (G * DELTI)
    upper[0] = 1.0 / (2.0 * DELXI * max(B[3], 1e-6))

    diag[-1] = RS / DELTI
    rhs[-1] = RS * TH[M1] / DELTI - Qr / (2.0 * DELXI)
    lower[-1] = -D[M2] / (2.0 * DELXI)

    # Assemble banded matrix for solve_banded
    ab = np.zeros((3, n))
    ab[0, 1:] = upper[:-1]  # upper diagonal
    ab[1, :] = diag         # main diagonal
    ab[2, :-1] = lower[1:]  # lower diagonal

    return ab, rhs, i_start

def conv(s, e, x, y):
    """
    Max |x - y| over s..e step 2; then copy x->y on those indices.
    Returns the residual (float).
    """
    sl = slice(s, e + 1, 2)
    r = np.max(np.abs(x[sl] - y[sl]))
    np.copyto(y[sl], x[sl], casting='no')
    return float(r)

def adaptive_tols():
    # Magnitudes to scale relative tolerances
    th_mag = np.max(np.abs(TH[3:M1+1:2])) or 1.0   # avoid 0
    tu_mag = np.max(np.abs(TU[2:M2+1:2])) or 1.0

    tol_th = max(TH_ABS_FLOOR, TH_REL * th_mag)
    tol_tu = max(TU_ABS_FLOOR, TU_REL * tu_mag)

    # Roughness indicator: CFL close to 1 or very shallow cells
    # CFL ≈ |u| * dt / dx = |U| * DELTI / (1/DELXI) = |U| * DELTI * DELXI
    cfl = (np.max(np.abs(U[2:M2+1])) * DELTI * (1.0/ (1.0/DELXI)))  # = |U|*dt/dx
    # Simpler: cfl = np.max(np.abs(U[2:M2+1])) * DELTI * DELXI
    cfl = np.max(np.abs(U[2:M2+1])) * DELTI * DELXI

    shallow = np.any(DEPTH[2:M2+1] < 1e-3)  # or a domain-specific threshold

    if cfl > 0.7 or shallow:
        # allow a bit looser stopping to avoid thrashing;
        # still cap so we don't reintroduce transport drift
        tol_th = min(max(tol_th, 1e-9), LOOSE_CAP)
        tol_tu = min(max(tol_tu, 1e-9), LOOSE_CAP)

    return tol_th, tol_tu

def tridag():
    """Solve tridiagonal system and assign results to TU and TH."""
    ab, rhs, i_start = build_tridiag_system()
    var = solve_banded((1, 1), ab, rhs)

    for idx, i in enumerate(range(i_start, i_start + len(var))):
        if i % 2 == 0:
            TU[i] = var[idx]
        else:
            TH[i] = var[idx]

