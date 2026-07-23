"""
Sediment module (translated from sed.c)
"""

from config import (
    M, rho_w, G, Chezy_lb, Chezy_ub, Mero_lb, Mero_ub,
    tau_ero_lb, tau_ero_ub, tau_dep_lb, tau_dep_ub, distance, DELTI,
    wMAX, kISS,
)
from config import ICE_MODEL
from variables import (U, DEPTH, v, tau_b, Mero, tau_ero, tau_dep, erosion, deposition,
                       Chezy, include_constantDEPTH, ice_frac)
from numba import njit


# The per-grid-point loop is @njit-compiled. It was ~0.31 s of pure-Python time in the
# phase-2 profile (the second-largest un-jitted item after transport), and is a plain
# numeric loop over ~136 points -- the same thin-wrapper / arrays-as-arguments pattern
# as biogeo_module and the hydro kernels. The transcription is verbatim, so results are
# bit-identical (verified by the output-hash check in docs/performance.md).
#
# @njit WITHOUT cache=True, exactly as biogeo_module: the kernel reads its config
# constants (rho_w, G, Mero_ub/lb, tau_*_lb/ub, distance, DELTI, wMAX, kISS) as module
# globals, and numba freezes globals as compile-time constants -- caching keys on
# bytecode, not on the captured values, so an edit to config.py would silently keep the
# previously compiled numbers. Recompiling each run costs ~1 s against a ~20 min run.
# The runtime-varying input (previousdays) and the run flags (ICE_MODEL,
# include_constantDEPTH) are passed as arguments instead.
#
# wISS(t, i, SPM) = wMAX * SPM/(SPM + kISS) is inlined here; it was used only by sed.
@njit
def _sed_loop(c_SPM, U, DEPTH, Chezy, tau_b, Mero, tau_ero, tau_dep,
              erosion, deposition, ice_frac, M, previousdays, ice_on, const_depth):
    """Jitted per-cell SPM (suspended sediment) update: bed shear stress -> erosion
    (scaled by open-water fraction under ice) and deposition (settling), then the SPM
    concentration change. Conserves SPM year-round under the ice model. Mutates in place."""
    for i in range(1, M + 1):
        # Settling velocity (wISS inlined)
        ws = wMAX * (c_SPM[i] / (c_SPM[i] + kISS))
        # Erosion and deposition rates for SPM [mg m^-2 s^-1]
        tau_b[i] = rho_w * G * U[i]**2 / Chezy[i]**2
        Mero[i] = Mero_ub if i >= distance else Mero_lb

        tau_ero[i] = (
            tau_ero_lb + (tau_ero_ub - tau_ero_lb) * (i - distance) / (M - distance)
            if i >= distance else tau_ero_lb
        )

        tau_dep[i] = (
            tau_dep_lb + (tau_dep_ub - tau_dep_lb) * (i - distance) / (M - distance)
            if i >= distance else tau_dep_lb
        )

        erosion[i] = (
            0.0 if tau_ero[i] >= tau_b[i]
            else Mero[i] * (tau_b[i] / tau_ero[i] - 1.0)
        )
        # An ice cover armours the bed against resuspension (no wind-wave stress, and
        # bottom-fast ice seals it entirely): scale erosion by the open-water fraction.
        if ice_on:
            erosion[i] *= (1.0 - ice_frac[i])
        deposition[i] = (
            ws * (1.0 - tau_b[i] / tau_dep[i]) * c_SPM[i]
            if tau_dep[i] >= tau_b[i] else 0.0
        )

        if const_depth == 1:
            tau_ero[i] = tau_ero_ub if i >= distance else tau_ero_lb
            tau_dep[i] = tau_dep_ub if i >= distance else tau_dep_lb

        # Conserve SPM year-round under the ice model (erosion already gated above);
        # fall back to the crude winter-zeroing gate when the ice model is off.
        if ice_on or previousdays > 0:
            # Update SPM concentration [g/l]
            c_SPM[i] = c_SPM[i] + (1.0 / DEPTH[i]) * (erosion[i] - deposition[i]) * DELTI
        else:
            c_SPM[i] = 0.0


def sed(t, previousdays):
    """Calculate sediment erosion and deposition rates."""
    _sed_loop(v['SPM']['c'], U, DEPTH, Chezy, tau_b, Mero, tau_ero, tau_dep,
              erosion, deposition, ice_frac, M, previousdays,
              ICE_MODEL, include_constantDEPTH)
