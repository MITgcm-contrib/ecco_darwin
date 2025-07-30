"""
Sediment module (translated from sed.c)
"""

from config import (
    M, rho_w, G, Chezy_lb, Chezy_ub, Mero_lb, Mero_ub,
    tau_ero_lb, tau_ero_ub, tau_dep_lb, tau_dep_ub, distance, ws, DELTI, TS, WARMUP
)
from variables import U, DEPTH, v, tau_b, Mero, tau_ero, tau_dep, erosion, deposition, Chezy, include_constantDEPTH
from forcings_module import get_sediment
from file_module import Rates

def sed(t):
    """Calculate sediment erosion and deposition rates."""
    for i in range(1, M + 1):
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
        deposition[i] = (
            ws * (1.0 - tau_b[i] / tau_dep[i]) * v['SPM']['c'][i]
            if tau_dep[i] >= tau_b[i] else 0.0
        )

        if include_constantDEPTH == 1:
            tau_ero[i] = tau_ero_ub if i >= distance else tau_ero_lb
            tau_dep[i] = tau_dep_ub if i >= distance else tau_dep_lb

        # Update SPM concentration [g/l]
        #v['SPM']['c'][i] = v['SPM']['c'][i] + (1.0 / DEPTH[i]) * (erosion[i] - deposition[i]) * DELTI
        #v['SPM']['c'][i] = get_sediment(t) + (1.0 / DEPTH[i]) * (erosion[i] - deposition[i]) * DELTI
        if t <= WARMUP:
            # pure physical model during warmup - let system equilibrate
            v['SPM']['c'][i] = v['SPM']['c'][i] + (1.0 / DEPTH[i]) * (erosion[i] - deposition[i]) * DELTI
        else:
            # hybrid model after warmup - time series + physical processes
            v['SPM']['c'][i] = get_sediment(t) + (1.0 / DEPTH[i]) * (erosion[i] - deposition[i]) * DELTI

    # Write erosion/deposition process rates [mg m^-2 s^-1]
    if (float(t) / float(TS * DELTI)) % 1 == 0:
        Rates(erosion, "erosion.dat", t)
        Rates(deposition, "deposition.dat", t)




