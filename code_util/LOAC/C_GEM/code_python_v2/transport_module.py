"""
Transport module (translated from transport.c)
"""

from config import MAXV, M, TS, DELTI, DELXI, M2
from variables import v, U, C, fl, Z, dispersion, D, U, Dold2
from schemes_module import openbound, tvd, crank_nicholson_dispersion
from file_module import transwrite

def transport(t):
    for name, tracer in v.items():
        if tracer["env"] != 1:
            continue

        co = tracer["c"]
        clb = tracer["clb"]
        cub = tracer["cub"]

        # Apply boundary conditions
        openbound(co, clb, cub)

        # Advection step (copy for cold state)
        cold = co.copy()
        co_advected = co.copy()
        tvd(co_advected, cold, fl, U, D, Dold2, DELTI, DELXI, M, M2)

        # Dispersion step
        co_new = crank_nicholson_dispersion(co_advected, D, dispersion, DELXI, DELTI)
        tracer["c"][:] = co_new

        # Accumulate average with optional Kahan
        tracer["avg"][1:M + 1] += co_new[1:M + 1]

        if t % (TS * DELTI) == 0:
            transwrite(co_new, tracer["name"], t)