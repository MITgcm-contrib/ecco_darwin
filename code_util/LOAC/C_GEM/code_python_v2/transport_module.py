"""
Transport module (translated from transport.c)
"""

from config import M, TS, DELTI
from variables import v
from schemes_module import openbound, disp_sch, tvd

def transport(t, io):
    """Main transport routine (optimized old-style)."""
    # Bind globals to locals
    vloc = v
    Mloc = M
    openb = openbound
    adv   = tvd
    disp  = disp_sch

    step = TS * DELTI
    write_now = (t % step == 0.0)
    if not write_now:
        r = t / step
        write_now = abs(r - round(r)) < 1e-12

    for name, tr in vloc.items():
        co = tr["c"]

        if tr["env"] == 1:
            # BCs, advection, dispersion for transported tracers
            openb(co, name)
            adv(co, name, t)
            disp(co)
            # accumulate averages
            tr["avg"][1:Mloc+1] += co[1:Mloc+1]

            if write_now:
                io.write_tracer(co, tr["name"], t)

        else:
            # Non-transported (pH): just average + optional write
            tr["avg"][1:Mloc+1] += co[1:Mloc+1]
            if write_now:
                io.write_tracer(co, tr["name"], t)
