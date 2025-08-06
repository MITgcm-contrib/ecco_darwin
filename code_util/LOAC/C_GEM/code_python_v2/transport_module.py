"""
Transport module (translated from transport.c)
"""
import numpy as np
from config import MAXV, M, TS, DELTI, DELXI
from variables import v, U, C, fl, Z, D , dispersion
from file_module import transwrite
from schemes_module import tvd_advection_twin, crank_nicholson_dispersion, openbound

def transport(t):
    """Main transport routine."""
    CFL = max(abs(U)) * DELTI / DELXI
    assert CFL < 1.0, f"CFL = {CFL:.3f} → too large! Reduce DELTI or increase DELXI."
    names = list(v.keys())
    for name in names:
        if v[name]["env"] == 1:
            co = v[name]["c"]
            # Apply open boundary conditions
            openbound(co, name, v, U, DELXI, DELTI)
            # Advection step
            co = tvd_advection_twin(co, U, DELXI, DELTI)
            # Dispersion step
            co = crank_nicholson_dispersion(co, D, dispersion, DELXI, DELTI)
            # Update concentration
            v[name]["c"][:] = co
            # Accumulate average concentrations
            v[name]["avg"] += v[name]["c"]

            # Write output at intervals
            if (float(t) / float(TS * DELTI)) % 1 == 0:
                transwrite(v[name]["c"], v[name]["name"], t)