"""
Transport module (translated from transport.c)
"""

from config import MAXV, M, TS, DELTI, ICE_MODEL
from variables import v, U, C, fl, Z
from schemes_module import openbound, tvd, disp_sch
from file_module import transwrite

def transport(t, previousdays):
    """Main transport routine."""
    names = list(v.keys())
    for names in names:
        if v[names]["env"] == 1:
            # With the prognostic ice model on, transport runs YEAR-ROUND: water (and
            # its tracers) keeps flowing under a floating ice cover, so the fields are
            # CONSERVED and advected/dispersed by whatever tidal + baseflow velocity
            # exists rather than being zeroed each winter. Zeroing was the crude gate's
            # stand-in for "no exchange under ice" and is what left DIC/ALK at ~0 at
            # ice-out (the carbonate-solve blow-up). With ICE_MODEL off, fall back to it.
            if ICE_MODEL or previousdays > 0:
                # Apply boundary conditions, advection, and dispersion
                openbound(v[names]["c"], names)
                tvd(v[names]["c"], names, t)
                disp_sch(v[names]["c"])
            else:
                # In-place, not a rebind. `c` is a numpy array shared by reference
                # with biogeo/sed/schemes; assigning a new list here would detach
                # every one of those references (see variables.py).
                v[names]["c"][:] = 0.0

            # Accumulate average concentrations (vectorised over indices 1..M, so
            # arithmetically identical to the per-element loop it replaces)
            v[names]["avg"][1:] += v[names]["c"][1:]

            # Write concentration output files at specific intervals
            if (float(t) / float(TS * DELTI)) % 1 == 0:
                transwrite(v[names]["c"], v[names]["name"], t)

        else:  # write pH output without being affected by transport
            # pH is not affected by advection/diffusion
            # it changes only based on S, DIC, ALK which are advected/diffused

            # Accumulate average concentrations (vectorised over indices 1..M, so
            # arithmetically identical to the per-element loop it replaces)
            v[names]["avg"][1:] += v[names]["c"][1:]

            # Write concentration output files at specific intervals
            if (float(t) / float(TS * DELTI)) % 1 == 0:
                transwrite(v[names]["c"], v[names]["name"], t)



