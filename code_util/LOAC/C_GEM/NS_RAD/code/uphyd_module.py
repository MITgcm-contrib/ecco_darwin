"""
Update hydrodynamics module (translated from uphyd.c)

`update` runs once per hydrodynamic solver iteration (~44 per timestep), so it is
@njit-compiled alongside the tridag kernels -- see tridag_module for why the arrays
are passed in explicitly rather than read from module globals. Arithmetic and
ordering are unchanged, so results are bit-identical to the pure-Python version.
"""

import numpy as np
from numba import njit

from config import M, M1, M2, M3
from variables import B, H, TH, D, ZZ, DEPTH, U, TU, Dold
from fun_module import Tide


def new_bc(t, Qr):
    """Set new boundary conditions."""
    H[1] = B[1] * Tide(t)
    TH[1] = H[1]
    D[1] = H[1] + ZZ[1]
    DEPTH[1] = D[1] / B[1]
    U[M] = Qr / D[M]
    TU[M] = U[M]


# After the tridiagonal solve, the new water level TH and velocity TU are known only on
# their native staggered nodes. This fills in the intermediate nodes and the derived
# fields (ZZ is the fixed cross-section at reference level, so D = TH + ZZ is the total
# section, and DEPTH = D / width):
#   * even nodes 2..M2: TH, D and DEPTH by averaging the two odd (elevation) neighbours
#   * odd nodes 3..M1:  D and DEPTH from TH+ZZ; velocity TU by averaging even neighbours
#   * last node M: TH, D, DEPTH by 2nd-order linear EXTRAPOLATION from the two interior odd
#     nodes, (3*x[M1] - x[M3]) / 2 -- there is no node beyond the boundary to average.
@njit(cache=True)
def _update(TH, ZZ, D, DEPTH, B, TU, M, M1, M2, M3):
    """Jitted post-solve update of the derived hydrodynamic geometry: recompute the total
    cross-section D and depth DEPTH from the new free cross-section TH and width B (and
    the reference cross-section ZZ), on the staggered grid. Mutates in place."""
    tmp = np.zeros(M + 1)

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


def update():
    """Update TH, D, DEPTH, and TU during iteration."""
    _update(TH, ZZ, D, DEPTH, B, TU, M, M1, M2, M3)


@njit(cache=True)
def _new_uh(U, TU, H, TH, Dold, D, M):
    """Jitted commit of the converged iteration: copy the temporary velocity TU and free
    cross-section TH into the live U and H, and roll the cross-section history (Dold <- D)
    for the next timestep. Mutates U, H, Dold in place."""
    for i in range(1, M + 1):
        U[i] = TU[i]
        if U[i] == 0:
            U[i] = 0.0001  # Avoid zero velocity
        H[i] = TH[i]
        Dold[i] = D[i]


def new_uh():
    """Update H and U arrays."""
    _new_uh(U, TU, H, TH, Dold, D, M)
