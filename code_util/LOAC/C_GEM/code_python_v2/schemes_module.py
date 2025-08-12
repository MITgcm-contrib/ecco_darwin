"""
Advection and dispersion schemes module (translated from schemes.c)
"""

import math
import numpy as np
from numba import njit
from config import M, M1, M2, M3, DELTI, DELXI
from variables import U, v, D, Dold2, fl, dispersion

@njit
def openbound(co, clb, cub):
    """Set boundary conditions with ghost cells."""
    if U[2] >= 0.0:
        co[1] = co[1] - (co[1] - clb) * U[2] * (float(DELTI) / float(DELXI))
    else:
        co[1] = co[1] - (co[3] - co[1]) * U[2] * (float(DELTI) / float(DELXI))

    if U[M1] >= 0.0:
        co[M] = co[M] - (co[M] - co[M2]) * U[M1] * (float(DELTI) / float(DELXI))
    else:
        co[M] = co[M] - (cub - co[M]) * U[M1] * (float(DELTI) / float(DELXI))
    co[M1] = co[M]

def openbound_carbonate(co, clb, cub, U, M):
    """
    Enforce DIC/ALK BCs at both ends with sign-aware inflow/outflow.
    - co: tracer array indexed 0..M (and optionally M+1 if you keep a right ghost)
    - clb, cub: scalar or array-like; if array-like we use clb[1], cub[M]
    - U: velocity array (use your conventional indices; here we take U[1] and U[M])
    """

    # helpers to accept scalar or array-like bcs
    def _left_val(bc):
        try:
            return float(bc[1])   # array-like with 1..M interior
        except Exception:
            return float(bc)      # scalar

    def _right_val(bc):
        try:
            return float(bc[M])
        except Exception:
            return float(bc)

    lbc = _left_val(clb)
    rbc = _right_val(cub)

    # face velocities (adjust if your U is staggered differently)
    uL = float(U[1])   # near left face
    uR = float(U[M])   # near right face

    # ---- Left boundary (upstream) ----
    if uL > 0.0:                 # inflow from left
        co[1] = lbc              # Dirichlet into first interior cell
        # left ghost (if present) consistent with Dirichlet
        if co.shape[0] >= 2:     # co[0] exists for 0..M
            co[0] = 2.0*lbc - co[1]
    else:                        # outflow at left → zero-gradient
        if co.shape[0] >= 2:
            co[0] = co[1]

    # ---- Right boundary (downstream) ----
    if uR < 0.0:                 # inflow from right
        co[M] = rbc              # Dirichlet into last interior cell
        # right ghost only if it exists (0..M+1 layout)
        if co.shape[0] > M+1:
            co[M+1] = 2.0*rbc - co[M]
    else:                        # outflow at right → zero-gradient ghost
        if co.shape[0] > M+1:
            co[M+1] = co[M]

# -------------------------------
# TVD Advection with Minmod Limiter
# -------------------------------
@njit
def tvd(co, cold, fl, U, D, Dold2, DELTI, DELXI, M, M2):
    b = 2.0
    dt_dx2 = DELTI / (2.0 * DELXI)

    for j in range(1, M2 + 1, 2):
        vx = U[j + 1]
        cfl = abs(vx) * dt_dx2

        if vx > 0.0:
            f = cold[j + 2] - cold[j]
            if abs(f) > 1e-35 and j != 1:
                rg = (cold[j] - cold[j - 2]) / f
                phi = (2.0 - cfl + rg * (1.0 + cfl)) / 3.0
                philen = max(0.0, min(b, min(b * rg, phi)))
            else:
                philen = 0.0

            co[j + 1] = cold[j] + 0.5 * (1.0 - cfl) * philen * f
            fl[j + 1] = vx * D[j + 1] * co[j + 1]

        else:
            f = cold[j] - cold[j + 2]
            if abs(f) > 1e-35 and j != M2:
                if j + 4 > M:
                    rg = (cold[j + 2] - cold[j + 3]) / f
                else:
                    rg = (cold[j + 2] - cold[j + 4]) / f
                phi = (2.0 - cfl + rg * (1.0 + cfl)) / 3.0
                philen = max(0.0, min(b, min(b * rg, phi)))
            else:
                philen = 0.0

            co[j + 1] = cold[j + 2] + 0.5 * (1.0 - cfl) * philen * f
            fl[j + 1] = vx * D[j + 1] * co[j + 1]

    # Flux divergence and diffusion
    for j in range(3, M2 + 1, 2):
        rat = Dold2[j] / D[j]
        rat1 = dt_dx2 / D[j]
        cold[j] = rat * cold[j] - rat1 * (fl[j + 1] - fl[j - 1])
        co[j] = cold[j]
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
