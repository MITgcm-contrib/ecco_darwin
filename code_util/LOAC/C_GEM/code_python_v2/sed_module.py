import numpy as np
from config import (
    M, rho_w, G, Mero_lb, Mero_ub, tau_ero_lb, tau_ero_ub, tau_dep_lb, tau_dep_ub,
    distance, ws, DELTI, TS
)
from variables import U, DEPTH, v, tau_b, Mero, tau_ero, tau_dep, erosion, deposition, Chezy, include_constantDEPTH

def sed(t, io):
    sl = slice(1, M + 1)
    i  = np.arange(1, M + 1)

    Uloc   = U[sl]
    Cz     = Chezy[sl]
    depth  = DEPTH[sl]
    spm    = v['SPM']['c'][sl]

    eps = 1e-30

    # bed shear τb [Pa]
    tau_b_arr = rho_w * G * Uloc**2 / np.maximum(Cz**2, eps)

    # spatial fields
    Mero_arr      = np.where(i >= distance, Mero_ub, Mero_lb)
    ramp          = (i - distance) / (M - distance)
    tau_ero_ramp  = np.where(i >= distance, tau_ero_lb + (tau_ero_ub - tau_ero_lb) * ramp, tau_ero_lb)
    tau_dep_ramp  = np.where(i >= distance, tau_dep_lb + (tau_dep_ub - tau_dep_lb) * ramp, tau_dep_lb)

    # sign-preserving denominators
    den_ero = np.where(tau_ero_ramp >= 0.0, np.maximum(tau_ero_ramp,  eps), np.minimum(tau_ero_ramp, -eps))
    den_dep = np.where(tau_dep_ramp >= 0.0, np.maximum(tau_dep_ramp,  eps), np.minimum(tau_dep_ramp, -eps))

    # erosion [kg m^-2 s^-1]
    erosion_arr = np.where(
        tau_ero_ramp >= tau_b_arr, 0.0, Mero_arr * (tau_b_arr / den_ero - 1.0)
    )

    # smooth deposition switch
    alpha = np.clip(1.0 - tau_b_arr / den_dep, 0.0, 1.0)     # ∈ [0,1]

    # legacy deposition diagnostic (for output): ws * alpha * SPM
    deposition_arr = ws * alpha * spm

    # semi-implicit update for SPM (stable for stiff sink)
    coef   = (DELTI / np.maximum(depth, 1e-12))
    numer  = spm + coef * erosion_arr
    denom  = 1.0 + coef * ws * alpha
    spm_new = numer / denom

    # nonnegative and tiny floor
    spm_new = np.maximum(spm_new, 0.0)

    # commit
    v['SPM']['c'][sl] = spm_new

    # publish diagnostics / fields
    tau_b[sl]   = tau_b_arr
    Mero[sl]    = Mero_arr
    tau_ero[sl] = tau_ero_ramp
    tau_dep[sl] = tau_dep_ramp
    erosion[sl] = erosion_arr
    deposition[sl] = deposition_arr

    if include_constantDEPTH == 1:
        tau_ero[sl] = np.where(i >= distance, tau_ero_ub, tau_ero_lb)
        tau_dep[sl] = np.where(i >= distance, tau_dep_ub, tau_dep_lb)

    # integer write cadence
    write_every = max(1, int(TS * DELTI))
    if (t % write_every) == 0:
        io.write_rates(erosion,    "erosion",    t)
        io.write_rates(deposition, "deposition", t)
