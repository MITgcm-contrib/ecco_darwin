"""
Hydrodynamics module (translated from hyd.c)
"""

from config import M, DELTI, TS, M1, M2, ITE
from variables import D, Dold2, TH, TU, E, Y
from uphyd_module import new_bc, new_uh, update
from tridag_module import tridag, conv, adaptive_tols

def hyd(t, io):
    Dold2[1:M+1] = D[1:M+1]
    new_bc(t)

    for it in range(ITE):  # safety cap
        tridag()
        tol_th, tol_tu = adaptive_tols()  # uses TU, not U

        r_th = conv(3, M1, TH, E)   # returns residual and copies TH->E
        r_tu = conv(2, M2, TU, Y)   # returns residual and copies TU->Y

        update()                    # keep before convergence check (legacy behavior)

        if (r_th < tol_th) and (r_tu < tol_tu):
            break

    else:
        # didn’t converge: optional fallback (halve dt, extra damped pass, log warning)
        raise RuntimeError("hyd: could not converge")
        pass

    new_uh()

    if t % (TS * DELTI) == 0:
        io.write_hyd(t)
