"""
Hydrodynamics module (translated from hyd.c)
"""

from config import M, M1, M2, M3, TOL, TS, DELTI, DELXI, G, RS
from variables import D, Dold2, TH, E, TU, Y, Z, C, U, B, H, Chezy, DEPTH, ZZ
from uphyd_module import new_bc, new_uh, _update
from tridag_module import _coeff_a, _tridag, _conv
from file_module import hydwrite
from numba import njit

# Cap on the tridiagonal-solve iterations per timestep. conv() returns exactly 0.0 or
# 1.0, so `rsum == 2.0` is a safe test WHEN it converges -- but if the solution
# oscillates (possible at extreme discharge/geometry) the loop would spin forever with
# no cap. The cap makes that a bounded, warned event instead of a silent hang; the
# model proceeds with the last iterate. ~64 iters/step is typical, so 500 is generous.
HYD_MAXITER = 500
_hyd_cap_warned = [False]  # warn once per run, not per timestep


# The convergence loop, fused into a single @njit driver. It runs ~64 passes per
# timestep, and each pass used to cross the Python<->numba boundary five times --
# coeff_a, tridag, conv, conv and update were each a thin Python wrapper around a jitted
# kernel, i.e. ~320 dispatches per timestep (~270 million over a full run) on kernels
# that are tiny ~136-point loops. Running the whole while-loop inside numba, calling the
# already-jitted cores natively, removes all of that boundary crossing. The arithmetic,
# its ordering and the exact convergence test (rsum != 2.0) are unchanged, so results
# are bit-identical -- verified by the output-hash check in docs/performance.md.
#
# Plain @njit (no cache=True): the driver calls cache=True kernels imported from other
# modules and is entered once per timestep, so the one-time ~1 s compile is negligible
# and there is no reason to take on a cross-module cached-call cache dependency.
@njit
def _hyd_iterate(Z, C, D, U, B, H, TH, Chezy, DEPTH, TU, E, Y, ZZ,
                 M, M1, M2, M3, DELTI, DELXI, G, RS, TOL, Qr, maxiter):
    """Jitted hydrodynamic solve for one timestep: iterate coeff_a -> tridiagonal solve
    -> convergence test -> update to convergence (both TH and TU within TOL), capped at
    `maxiter` to prevent an infinite loop. Returns the iteration count; mutates the
    staggered-grid arrays (TH, TU, D, U, H, ...) in place."""
    rsum = 0.0
    it = 0
    while rsum != 2.0 and it < maxiter:
        it += 1
        _coeff_a(Z, C, D, U, B, H, TH, Chezy, DEPTH,
                 M1, M2, M3, DELTI, DELXI, G, RS, Qr)
        _tridag(C, Z, TU, TH, M, M1, M2)
        rsum = _conv(3, M1, TOL, TH, E) + _conv(2, M2, TOL, TU, Y)
        _update(TH, ZZ, D, DEPTH, B, TU, M, M1, M2, M3)
    return it, rsum


def hyd(t, Qr):
    """Main hydrodynamics routine."""
    # Store previous depths (copy D[1..M] -> Dold2[1..M]; vectorised, same values)
    Dold2[1:] = D[1:]

    # Set new boundary conditions
    new_bc(t, Qr)

    # Iterative solve to convergence, fully inside numba (see _hyd_iterate). Returns the
    # iteration count and the final rsum so the once-per-run cap warning stays in Python.
    i, rsum = _hyd_iterate(Z, C, D, U, B, H, TH, Chezy, DEPTH, TU, E, Y, ZZ,
                           M, M1, M2, M3, DELTI, DELXI, G, RS, TOL, Qr, HYD_MAXITER)

    if i >= HYD_MAXITER and rsum != 2. and not _hyd_cap_warned[0]:
        print(f"  [hyd] convergence not reached in {HYD_MAXITER} iters at "
              f"t={t/86400:.2f} d (rsum={rsum}); proceeding with last iterate. "
              f"(warned once)")
        _hyd_cap_warned[0] = True

    # Update U and H values
    new_uh()

    # Write hydrodynamic data to files at specified intervals
    if (float(t) / float(TS * DELTI)) % 1 == 0:
        hydwrite(t)
