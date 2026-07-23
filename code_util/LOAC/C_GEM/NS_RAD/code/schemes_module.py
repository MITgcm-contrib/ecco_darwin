"""
Advection and dispersion schemes module (translated from schemes.c)

tvd and disp_sch are the per-species-per-timestep advection/dispersion hot loops.
They are @njit-compiled following the same kernel pattern as tridag_module /
uphyd_module: an inner _kernel taking arrays + scalars explicitly (numba freezes
module globals as compile-time constants, which would freeze array contents), wrapped
by a thin same-named public function that reads the globals and forwards them, keeping
the call signatures tvd(co, s, t) / disp_sch(co) unchanged. openbound stays pure
Python -- it reads the v[s] boundary dict, which does not jit, and only touches the
two boundary cells (cheap). Results are bit-identical (verified by output hash).
"""

import numpy as np
from numba import njit

from config import M, M1, M2, M3, DELTI, DELXI
from variables import U, v, D, Dold2, fl, dispersion

def openbound(co, s):
    """Set boundary conditions."""

    if U[2] >= 0.0:
        co[1] = co[1] - (co[1] - v[s]["clb"]) * U[2] * (float(DELTI) / float(DELXI))
    else:
        co[1] = co[1] - (co[3] - co[1]) * U[2] * (float(DELTI) / float(DELXI))

    if U[M1] >= 0.0:
        co[M] = co[M] - (co[M] - co[M2]) * U[M1] * (float(DELTI) / float(DELXI))
    else:
        co[M] = co[M] - (v[s]["cub"] - co[M]) * U[M1] * (float(DELTI) / float(DELXI))
    co[M1] = co[M]

# Total Variation Diminishing (TVD) advection on the staggered grid (C-GEM; Volta et al.
# 2014). Loop 1 builds a flux-limited estimate of each cell-face concentration:
#   cfl        Courant number |u| dt / (2 dx)
#   rg         ratio of consecutive concentration gradients (upwind / local) -- the
#              smoothness monitor: ~1 in smooth flow, far from 1 near a front
#   philen     the flux limiter, max(0, min(b, min(b*rg, phi))) with b = 2 -- it caps the
#              anti-diffusive correction so no new maxima/minima are created (the TVD
#              property), degrading toward first-order upwind (philen->0) at sharp fronts
# co[j+1] is then the limited face value and fl[j+1] = u * D * co the advective flux; the
# two velocity signs are handled symmetrically. Loop 2 is the conservative update: cell j
# advances from the flux divergence (fl[j+1]-fl[j-1]), scaled by Dold2/D for the
# time-varying cross-section.
@njit(cache=True)
def _tvd(co, U, D, Dold2, fl, M, M2, DELTI, DELXI):
    """Jitted TVD advection of one scalar field `co` on the staggered grid: a flux-limited
    (superbee-type) upwind scheme that is monotone (no spurious over/undershoots). Writes
    the advected field back into `co` in place; `fl` is scratch for the advective flux."""
    b = 2.0
    # Snapshot of the field at the start of the step. The original aliased (`cold = co`)
    # and shadowed the module-level `cold`; a real copy is bit-identical here (the
    # access patterns are disjoint -- loop 1 writes only even slots & reads only odd,
    # loop 2 touches only its own odd index) and is required for the jitted kernel.
    cold = co.copy()

    for j in range(1, M2 + 1, 2):
        vx = U[j + 1]
        cfl = abs(vx) * (float(DELTI) / (2.0 * float(DELXI)))

        # For positive velocity, vx > 0
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

        # For negative velocity, vx <= 0
        if vx <= 0.0:
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

    for j in range(3, M2 + 1, 2):
        rat = Dold2[j] / D[j]
        rat1 = (float(DELTI) / (2.0 * float(DELXI))) / D[j]
        cold[j] = rat * cold[j] - rat1 * (fl[j + 1] - fl[j - 1])
        co[j] = cold[j]


def tvd(co, s, t):
    """Total Variation Diminishing (TVD) advection scheme. `s`/`t` unused by the
    numerics -- they exist only for the transport call signature."""
    _tvd(co, U, D, Dold2, fl, M, M2, DELTI, DELXI)


# Longitudinal dispersion  d/dx( K A dc/dx )  advanced implicitly by Crank-Nicolson and
# solved with the Thomas algorithm (C-GEM; Volta et al. 2014). Per interior node: a and c
# are the sub- and super-diagonals, b the diagonal, r the mirrored right-hand-side
# coefficients, and di the explicit (old-time) RHS vector. The factor dt/(8 dx^2) is the
# CN half-weight on the staggered 2*dx spacing; config.DISP_MAX bounds K <= 4 dx^2/dt so
# the RHS coefficient r stays non-negative (non-oscillatory). Ends are held fixed
# (a[1]=b[1]=1, di[1]=co[1], likewise at M). The final two loops are the Thomas forward
# elimination and back substitution, writing the updated concentrations back into co in place.
@njit(cache=True)
def _disp_sch(co, dispersion, D, M, M1, M3, DELTI, DELXI):
    """Jitted longitudinal dispersion of one scalar field `co`: a Crank-Nicolson (implicit)
    solve of the diffusion operator with the per-cell dispersion coefficient, via a
    tridiagonal system. Updates `co` in place. See config.DISP_MAX for the stability cap."""
    a = np.zeros(M + 1)
    b = np.zeros(M + 1)
    c = np.zeros(M + 1)
    r = np.zeros(M + 1)
    di = np.zeros(M + 1)
    gam = np.zeros(M + 1)
    a[1] = b[1] = 1.0
    di[1] = co[1]
    a[M] = b[M] = 1.0
    di[M] = co[M]

    for i in range(2, M1 + 1):
        c1 = dispersion[i - 1] * D[i - 1] / D[i]
        c2 = dispersion[i + 1] * D[i + 1] / D[i]
        a[i] = -c1 * (float(DELTI) / (8.0 * float(DELXI)**2))
        c[i] = -c2 * (float(DELTI) / (8.0 * float(DELXI)**2))
        b[i] = 1.0 + c2 * (float(DELTI) / (8.0 * float(DELXI)**2)) + c1 * (float(DELTI) / (8.0 * float(DELXI)**2))
        r[i] = 1.0 - c2 * (float(DELTI) / (8.0 * float(DELXI)**2)) - c1 * (float(DELTI) / (8.0 * float(DELXI)**2))

    for i in range(3, M1 + 1, 2):
        if i + 2 > M:
            di[i] = -c[i] * co[i + 1] + r[i] * co[i] - a[i] * co[i - 2]
        else:
            di[i] = -c[i] * co[i + 2] + r[i] * co[i] - a[i] * co[i - 2]

    if b[1] == 0.0:
        raise ValueError("Error: b[1] cannot be 0")

    bet = b[1]
    co[1] = di[1] / bet

    for j in range(3, M1 + 1, 2):
        gam[j] = c[j - 2] / bet
        bet = b[j] - a[j] * gam[j]
        co[j] = (di[j] - a[j] * co[j - 2]) / bet

    for j in range(M3, 0, -2):
        co[j] = co[j] - gam[j + 2] * co[j + 2]


def disp_sch(co):
    """Dispersion scheme (Thomas solve)."""
    _disp_sch(co, dispersion, D, M, M1, M3, DELTI, DELXI)
