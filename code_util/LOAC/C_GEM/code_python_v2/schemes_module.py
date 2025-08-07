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
def minmod(a, b):
    if a * b <= 0.0:
        return 0.0
    else:
        return math.copysign(min(abs(a), abs(b)), a)

@njit
def tvd_advection_twin(co, U, DELXI, DELTI):
    N = co.shape[0]
    q_new = co.copy()
    flux = np.zeros(N)
    b = 2.0
    max_cfl = np.max(np.abs(U) * DELTI / DELXI)
    if max_cfl > 1.0:
        print(f"Warning: CFL condition violated: {max_cfl}")
    # Serial loop instead of parallel prange
    #for j in range(1, N - 3):
    for j in range(1, N - 2, 2):
        vx = U[j + 1]
        cfl = abs(vx) * DELTI / (2 * DELXI)

        if vx > 0.0:
            f = co[j + 2] - co[j]
            if abs(f) > 1e-35 and j > 1:
                #rg = (co[j] - co[j - 2]) / f
                #phi = (2.0 - cfl + rg * (1.0 + cfl)) / 3.0
                #philen = max(0.0, min(b, min(b * rg, phi)))
                a = (co[j] - co[j - 2]) / 2.0
                b = f / 2.0
                philen = minmod(a, b)
            else:
                philen = 0.0
            #q_int = co[j] + 0.5 * (1.0 - cfl) * philen * f
            q_new[j+1] = co[j] + (1.0 - cfl) * philen
        else:
            f = co[j] - co[j + 2]
            if abs(f) > 1e-35 and j < N - 2:
                #rg = (co[j + 2] - co[j + 4]) / f
                #phi = (2.0 - cfl + rg * (1.0 + cfl)) / 3.0
                #philen = max(0.0, min(b, min(b * rg, phi)))
                a = (co[j + 2] - co[j + 4]) / 2.0
                b = f / 2.0
                philen = minmod(a, b)
            else:
                philen = 0.0
            #q_int = co[j + 2] + 0.5 * (1.0 - cfl) * philen * f
            q_new[j+1] = co[j + 2] + (1.0 - cfl) * philen

        flux[j + 1] = vx * q_new[j+1] * D[j + 1]

    # Serial loop with step 2 (still not parallel)
    #for j in range(3, N - 3, 2):
    for j in range(3, N - 2, 2):
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

    # Solve only up to N-3; keep last two values fixed
    c1 = dispersion[0:N-3] * D[0:N-3] / D[1:N-2]
    c2 = dispersion[2:N-1] * D[2:N-1] / D[1:N-2]

    a[1:N-2] = -c1 * alpha
    b[1:N-2] = 1.0 + (c1 + c2) * alpha
    c[1:N-2] = -c2 * alpha

    rhs[1:N-2] = (1.0 - (c1 + c2) * alpha) * co[1:N-2]
    rhs[1:N-2] += c1 * alpha * co[0:N-3] + c2 * alpha * co[2:N-1]

    # BC at j=0: Dirichlet (or Neumann if desired)
    b[0] = 1.0
    rhs[0] = co[0]

    # j = N-2 and N-1: no update
    b[N-2:] = 1.0
    a[N-2:] = 0.0
    c[N-2:] = 0.0
    rhs[N-2:] = co[N-2:]

    # Assemble banded matrix
    ab = np.zeros((3, N))
    ab[0, 1:] = c[:-1]
    ab[1, :] = b
    ab[2, :-1] = a[1:]

    q_new = solve_banded((1,1), ab, rhs)

    return q_new


# -------------------------------
# Open Boundary Conditions
# -------------------------------

@njit
def openbound(co, clb, cub, U, DELXI, DELTI):
    """
    Boundary condition update similar to original, robust under inflow/outflow.
    M: last valid index
    M1: M - 1
    M2: M - 2
    """
    # Lower boundary (j=1)
    if U[2] >= 0.0:
        co[1] -= (co[1] - clb) * U[2] * DELTI / DELXI
    else:
        co[1] -= (co[3] - co[1]) * U[2] * DELTI / DELXI

    # Upper boundary (j=M)
    if U[M1] >= 0.0:
        co[M] -= (co[M] - co[M2]) * U[M1] * DELTI / DELXI
    else:
        co[M] -= (cub - co[M]) * U[M1] * DELTI / DELXI

    # Ghost cell update
    co[M1] = co[M]