"""
Sediment module (translated from sed.c)
"""

from config import (
    M, rho_w, G, Chezy_lb, Chezy_ub, Mero_lb, Mero_ub,
    tau_ero_lb, tau_ero_ub, distance, DELTI,
    wMAX, kISS,
)
from config import ICE_MODEL
from variables import (U, DEPTH, v, tau_b, Mero, tau_ero, erosion, deposition,
                       Chezy, include_constantDEPTH, ice_frac)
from numba import njit


# The per-grid-point loop is @njit-compiled. It was ~0.31 s of pure-Python time in the
# phase-2 profile (the second-largest un-jitted item after transport), and is a plain
# numeric loop over ~136 points -- the same thin-wrapper / arrays-as-arguments pattern
# as biogeo_module and the hydro kernels.
#
# @njit WITHOUT cache=True, exactly as biogeo_module: the kernel reads its config
# constants (rho_w, G, Mero_ub/lb, tau_ero_lb/ub, distance, DELTI, wMAX, kISS) as module
# globals, and numba freezes globals as compile-time constants -- caching keys on
# bytecode, not on the captured values, so an edit to config.py would silently keep the
# previously compiled numbers. Recompiling each run costs ~1 s against a ~20 min run.
# The runtime-varying input (previousdays) and the run flags (ICE_MODEL,
# include_constantDEPTH) are passed as arguments instead.
#
# wISS(t, i, SPM) = wMAX * SPM/(SPM + kISS) is inlined here; it was used only by sed.
#
# TWO STACKED DEFECTS, BOTH FIXED -- see CLAUDE.md -> "Known defects" (search
# "tau_dep" and "MG_TO_G") for the full citation trail and verification numbers.
#
# (1) DEPOSITION WAS SHEAR-STRESS-GATED BY AN UNCITED `tau_dep` THRESHOLD (removed).
# Erosion here (`Mero*(tau_b/tau_ero - 1)`) is Clark et al. 2022's actual Yukon-delta
# formulation (https://doi.org/10.1029/2022JG007139, Eq. 4: `M_tau*(tau_b - tau_crit)`,
# algebraically identical once M_tau = Mero/tau_ero) -- confirmed against the paper's
# own Table 1 (tau_crit = 0.005 Pa, M_tau = 1.0e-5 g m^-2 s^-1 Pa^-1). But that paper's
# OWN deposition is governed purely by the concentration-dependent settling velocity
# below (`ws`), with NO separate shear-stress threshold -- the `tau_dep`-gated
# `(1 - tau_b/tau_dep)` factor previously here had no citation in config.py and no
# counterpart in the cited paper, and it suppressed deposition during exactly the
# higher-flow conditions when erosion was most active, breaking the negative feedback
# that should keep the two in balance. `tau_dep`/`tau_dep_lb`/`tau_dep_ub` are removed
# entirely (not just disabled).
#
# (2) EROSION WAS MISSING AN mg->g UNIT CONVERSION (fixed, `MG_TO_G` below). `Mero` is
# documented AND cited in mg m^-2 s^-1 (matching Clark et al. 2022's Table 1 once
# converted), but `c_SPM` -- the state variable erosion feeds directly into -- is in
# g/L (see variables.py, config.py's own unit comments, and every sites/<name>.py
# BOUNDARIES value). `deposition[i] = ws*c_SPM[i]` is correctly scaled (ws is true
# SI m/s, c_SPM true g/L), but erosion[i] was combined with it UNCONVERTED -- a
# thousandfold mismatch between the two terms of the same subtraction. Caught only
# after fixing (1): with the tau_dep gate removed, SPM should have settled at the
# 0.0004-0.0053 g/L an analytic equilibrium check predicted (matching the real USGS
# river values, 0.0075-0.0165 g/L) -- instead a definitive rerun still spiked to
# 232 g/L at Kuparuk's mouth during its documented extreme-freshet event (day ~518,
# see "Geometry" -> "Kuparuk's extreme freshet response"). Reproducing the exact
# per-cell arithmetic at that event's bed shear stress converged to a true
# equilibrium of ~5440 g/L as literally coded (unconverted) vs. ~5.5 g/L with the
# mg->g factor applied -- a clean 1000x, confirming the unit mismatch rather than a
# further physics gap. (5.5 g/L at this one documented-extreme event is still
# elevated vs. the ~0.01 g/L baseline, which is expected: erosion scales with
# velocity squared, and this event's velocity is already flagged elsewhere as a
# physically extreme consequence of Kuparuk's geometry, not a numerics problem.)
MG_TO_G = 1.0e-3


@njit
def _sed_loop(c_SPM, U, DEPTH, Chezy, tau_b, Mero, tau_ero,
              erosion, deposition, ice_frac, M, previousdays, ice_on, const_depth):
    """Jitted per-cell SPM (suspended sediment) update: bed shear stress -> erosion
    (scaled by open-water fraction under ice) and deposition (settling), then the SPM
    concentration change. Conserves SPM year-round under the ice model. Mutates in place."""
    for i in range(1, M + 1):
        # Settling velocity (wISS inlined)
        ws = wMAX * (c_SPM[i] / (c_SPM[i] + kISS))
        # Bed shear stress [N/m^2]; Mero (erosion coefficient) stays in its documented/
        # cited mg m^-2 s^-1 for anyone reading it back -- the mg->g conversion is
        # applied only where erosion[i] is computed, right before it meets c_SPM (g/L).
        tau_b[i] = rho_w * G * U[i]**2 / Chezy[i]**2
        Mero[i] = Mero_ub if i >= distance else Mero_lb

        tau_ero[i] = (
            tau_ero_lb + (tau_ero_ub - tau_ero_lb) * (i - distance) / (M - distance)
            if i >= distance else tau_ero_lb
        )

        # erosion[i] is in g m^-2 s^-1 (MG_TO_G-converted), matching deposition[i]
        # below and the g/L state c_SPM -- see (2) above.
        erosion[i] = (
            0.0 if tau_ero[i] >= tau_b[i]
            else Mero[i] * MG_TO_G * (tau_b[i] / tau_ero[i] - 1.0)
        )
        # An ice cover armours the bed against resuspension (no wind-wave stress, and
        # bottom-fast ice seals it entirely): scale erosion by the open-water fraction.
        if ice_on:
            erosion[i] *= (1.0 - ice_frac[i])
        # Concentration-dependent settling only -- see (1) above for why this is NOT
        # also shear-stress-gated.
        deposition[i] = ws * c_SPM[i]

        if const_depth == 1:
            tau_ero[i] = tau_ero_ub if i >= distance else tau_ero_lb

        # Conserve SPM year-round under the ice model (erosion already gated above);
        # fall back to the crude winter-zeroing gate when the ice model is off.
        if ice_on or previousdays > 0:
            # Update SPM concentration [g/l]
            c_SPM[i] = c_SPM[i] + (1.0 / DEPTH[i]) * (erosion[i] - deposition[i]) * DELTI
        else:
            c_SPM[i] = 0.0


def sed(t, previousdays):
    """Calculate sediment erosion and deposition rates."""
    _sed_loop(v['SPM']['c'], U, DEPTH, Chezy, tau_b, Mero, tau_ero,
              erosion, deposition, ice_frac, M, previousdays,
              ICE_MODEL, include_constantDEPTH)
