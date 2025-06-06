"""
Transport module (translated from transport.c)
"""

from config import MAXV, M, TS, DELTI
from variables import v, U, C, fl, Z
from schemes_module import openbound, tvd, disp_sch
from file_module import transwrite

def transport(t, previousdays):
    """Main transport routine."""
    names = list(v.keys())
    for names in names:
        if v[names]["env"] == 1:
            if previousdays > 0:
                # Apply boundary conditions, advection, and dispersion
                openbound(v[names]["c"], names)
                tvd(v[names]["c"], names, t)
                disp_sch(v[names]["c"])
            else:
                v[names]["c"] = [0.0] * (M + 1)

            # Accumulate average concentrations
            for i in range(1, M + 1):
                v[names]["avg"][i] += v[names]["c"][i]

            # Write concentration output files at specific intervals
            if (float(t) / float(TS * DELTI)) % 1 == 0:
                transwrite(v[names]["c"], v[names]["name"], t)

        else:  # write pH output without being affected by transport
            # pH is not affected by advection/diffusion
            # it changes only based on S, DIC, ALK which are advected/diffused

            # Accumulate average concentrations
            for i in range(1, M + 1):
                v[names]["avg"][i] += v[names]["c"][i]

            # Write concentration output files at specific intervals
            if (float(t) / float(TS * DELTI)) % 1 == 0:
                transwrite(v[names]["c"], v[names]["name"], t)



