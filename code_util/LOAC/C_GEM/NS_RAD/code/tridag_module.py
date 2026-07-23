"""
Tridiagonal matrix solver module (translated from tridag.c)

The three routines here are the hydrodynamic inner loop: `hyd` iterates
coeff_a -> tridag -> conv -> update until convergence, which under the shallow
observation-based geometry takes ~45-65 passes per timestep (site-dependent: ~64 on
shallow Kuparuk, D=1.34 m; it was ~21 with the old 15 m depth -- shallower channels
converge more slowly). Profiling put this loop at ~60% of total runtime, so the numeric
kernels are @njit-compiled (the whole convergence loop is now fused into one njit driver,
hyd_module._hyd_iterate, which calls these cores natively).

The public functions keep their original signatures and still mutate the shared
arrays in variables.py in place. The jitted kernels take those arrays as explicit
arguments because numba treats module globals as compile-time constants, which
would freeze the array contents. Arithmetic and its ordering are unchanged, so
results are bit-identical to the pure-Python version (verified by output hash).
"""

import numpy as np
from numba import njit

from config import M, M1, M2, M3, DELTI, DELXI, G, RS
from variables import D, Z, C, U, B, H, TH, TU, Chezy, DEPTH


# Assemble the tridiagonal system for one iteration of C-GEM's implicit 1-D hydrodynamics
# (coupled continuity + momentum, the Saint-Venant equations; Volta et al. 2014). The grid
# is STAGGERED: odd nodes carry water level / cross-section (H), even nodes carry velocity
# (U). Each interior pair (j odd, i = j+1 even) writes two matrix rows into Z (right-hand
# side) and C[.,1..4] (sub / diagonal / super / spare diagonals):
#   * odd row j  -> CONTINUITY  dA/dt + dQ/dx = 0:  C[j,1]=-D[j-1]/2dx, C[j,2]=RS/dt,
#     C[j,3]=D[j+1]/2dx, Z[j]=RS*H[j]/dt   (RS = storage-width ratio, config.RS)
#   * even row i -> MOMENTUM  dU/dt + U dU/dx + g dH/dx + friction = 0, where C[i,2] carries
#     bed friction (1/Chezy^2)|U|/DEPTH and advection (U[i+2]-U[i-2])/(4 g dx)
# The four hand-written rows Z[2]/C[2,*] (downstream, the tidal mouth) and Z[M1]/C[M1,*]
# (upstream, where river discharge enters as -Qr/2dx) are the boundary closures. Coefficients
# are transcribed verbatim from the C original -- see Volta et al. 2014 for the derivation.
@njit(cache=True)
def _coeff_a(Z, C, D, U, B, H, TH, Chezy, DEPTH,
             M1, M2, M3, DELTI, DELXI, G, RS, Qr):
    """Jitted assembly of the tridiagonal coefficient matrix `C` and RHS `Z` for the
    hydrodynamic solve: discretises the 1-D shallow-water continuity + momentum equations
    (with Chezy friction) on the staggered grid for this iteration. Mutates C, Z in place."""
    for j in range(3, M3 + 1, 2):
        i = j + 1
        Z[j] = RS * H[j] / DELTI
        C[j, 1] = -D[j - 1] / (2.0 * DELXI)
        C[j, 2] = RS / DELTI
        C[j, 3] = D[j + 1] / (2.0 * DELXI)
        C[j, 4] = 0.0

        Z[i] = 1.0 / (G * DELTI) * U[i]
        C[i, 1] = -1.0 / (2.0 * DELXI * B[i - 1])
        C[i, 2] = (
            1.0 / (G * DELTI)
            + (1.0 / (Chezy[i] ** 2)) * abs(U[i]) / DEPTH[i]
            + (U[i + 2] - U[i - 2]) / (4.0 * G * DELXI)
        )
        C[i, 3] = 1.0 / (2.0 * DELXI * B[i + 1])
        C[i, 4] = 0.0

    Z[2] = 1.0 / (2.0 * DELXI) * (H[1] / B[1]) + 1.0 / (G * DELTI) * U[2]
    Z[M1] = RS * TH[M1] / DELTI - Qr / (2.0 * DELXI)
    C[2, 1] = 0.0
    C[2, 2] = (
        1.0 / (G * DELTI)
        + (1.0 / (Chezy[2] ** 2)) * abs(U[2]) / DEPTH[2]
        + (U[4] - U[2]) / (2.0 * G * DELXI)
    )
    C[2, 3] = 1.0 / (2.0 * DELXI * B[3])
    C[2, 4] = 0.0
    C[M1, 1] = -D[M2] / (2.0 * DELXI)
    C[M1, 2] = RS / DELTI
    C[M1, 3] = 0.0
    C[M1, 4] = 0.0


def coeff_a(t, Qr):
    """Calculate coefficients for tridiagonal matrix."""
    _coeff_a(Z, C, D, U, B, H, TH, Chezy, DEPTH,
             M1, M2, M3, DELTI, DELXI, G, RS, Qr)


@njit(cache=True)
def _tridag(C, Z, TU, TH, M, M1, M2):
    """Jitted Thomas-algorithm solve of the tridiagonal system `C·x = Z` from _coeff_a,
    writing the new velocities (TU, even indices) and free cross-sections / levels (TH,
    odd indices) of the staggered grid in place."""
    gam = np.zeros(M + 1)
    var = np.zeros(M + 1)

    bet = C[2, 2]
    var[2] = Z[2] / bet

    for j in range(3, M1 + 1):
        gam[j] = C[j - 1, 3] / bet
        bet = C[j, 2] - C[j, 1] * gam[j]
        var[j] = (Z[j] - C[j, 1] * var[j - 1]) / bet

    for j in range(M2, 1, -1):
        var[j] -= gam[j + 1] * var[j + 1]

    for j in range(2, M2 + 1, 2):
        TU[j] = var[j]
        TH[j + 1] = var[j + 1]


def tridag():
    """Solve the tridiagonal matrix system."""
    _tridag(C, Z, TU, TH, M, M1, M2)


@njit(cache=True)
def _conv(s, e, toler, xarray, yarray):
    """Jitted convergence test over staggered indices s..e (step 2): returns 1.0 if the
    max |xarray - yarray| is within `toler`, else 0.0. The hyd loop iterates until both
    the TH and TU tests pass (sum == 2.0)."""
    t = 0.0
    for i in range(s, e + 1, 2):
        diff = abs(xarray[i] - yarray[i])
        if diff > t:
            t = diff
    result = 0.0 if t >= toler else 1.0
    for i in range(s, e + 1, 2):
        yarray[i] = xarray[i]
    return result


def conv(s, e, toler, xarray, yarray):
    """Test for convergence between two arrays."""
    return _conv(s, e, toler, xarray, yarray)
