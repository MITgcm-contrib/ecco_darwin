"""
Advection and dispersion schemes module (translated from schemes.c)
"""
import math
from variables import U, v, D, Dold2, fl, dispersion
import numpy as np
from numba import njit, prange
from config import M, M1, M2, M3, DELTI, DELXI

# -------------------------------
# TVD Advection with Minmod Limiter
# -------------------------------

@njit
def tvd_advection_twin(co, U, DELXI, DELTI):
    N = co.shape[0]
    q_new = co.copy()
    flux = np.zeros(N)
    b = 2.0

    # Serial loop instead of parallel prange
    for j in range(1, N - 3):
        vx = U[j + 1]
        cfl = abs(vx) * DELTI / (2 * DELXI)

        if vx > 0.0:
            f = co[j + 2] - co[j]
            if abs(f) > 1e-35 and j > 1:
                rg = (co[j] - co[j - 2]) / f
                phi = (2.0 - cfl + rg * (1.0 + cfl)) / 3.0
                philen = max(0.0, min(b, min(b * rg, phi)))
            else:
                philen = 0.0
            q_int = co[j] + 0.5 * (1.0 - cfl) * philen * f
        else:
            f = co[j] - co[j + 2]
            if abs(f) > 1e-35 and j < N - 5:
                rg = (co[j + 2] - co[j + 4]) / f
                phi = (2.0 - cfl + rg * (1.0 + cfl)) / 3.0
                philen = max(0.0, min(b, min(b * rg, phi)))
            else:
                philen = 0.0
            q_int = co[j + 2] + 0.5 * (1.0 - cfl) * philen * f

        flux[j + 1] = vx * q_int * D[j + 1]

    # Serial loop with step 2 (still not parallel)
    for j in range(3, N - 3, 2):
        if D[j] > 1e-12:  # avoid division by zero or near zero
            q_new[j] = co[j] - (DELTI / (2.0 * DELXI)) * (flux[j + 1] - flux[j - 1]) / D[j]
        else:
            q_new[j] = co[j]

        # Prevent negative or unphysical concentrations
        if q_new[j] < 0.0:
            q_new[j] = 0.0

    return q_new

# -------------------------------
# Crank–Nicolson Dispersion Solver
# -------------------------------

from scipy.linalg import solve_banded

def crank_nicholson_dispersion(co, D, dispersion, DELXI, DELTI):
    N = len(co)
    alpha = DELTI / (8.0 * DELXI**2)

    a = np.zeros(N)
    b = np.ones(N)
    c = np.zeros(N)
    rhs = np.zeros(N)

    # Interior points vectorized (except edges)
    c1 = dispersion[:-2] * D[:-2] / D[1:-1]
    c2 = dispersion[2:] * D[2:] / D[1:-1]

    a[1:-1] = -c1 * alpha
    b[1:-1] = 1.0 + (c1 + c2) * alpha
    c[1:-1] = -c2 * alpha

    rhs[1:-1] = (1.0 - (c1 + c2) * alpha) * co[1:-1]
    rhs[1:-1] += c1 * alpha * co[:-2] + c2 * alpha * co[2:]

    # Boundary conditions (Dirichlet)
    b[0] = b[-1] = 1.0
    rhs[0] = co[0]
    rhs[-1] = co[-1]

    # Assemble banded matrix for solve_banded
    ab = np.zeros((3, N))
    ab[0, 1:] = c[0:-1]  # upper diagonal
    ab[1, :] = b         # main diagonal
    ab[2, :-1] = a[1:]   # lower diagonal

    from scipy.linalg import solve_banded
    q_new = solve_banded((1,1), ab, rhs)

    return q_new


# -------------------------------
# Open Boundary Conditions
# -------------------------------

def openbound(co, s, v, U, DELXI, DELTI):
    """
    Modernized open boundary condition update.
    clb, cub: external boundary conditions
    """
    if U[2] >= 0.0:
        co[1] -= (co[1] - v[s]["clb"]) * U[2] * DELTI / DELXI
    else:
        co[1] -= (co[3] - co[1]) * U[2] * DELTI / DELXI

    if U[-2] >= 0.0:
        co[-1] -= (co[-1] - co[-3]) * U[-2] * DELTI / DELXI
    else:
        co[-1] -= (v[s]["cub"] - co[-1]) * U[-2] * DELTI / DELXI

    co[-2] = co[-1]  # enforce ghost cell update
