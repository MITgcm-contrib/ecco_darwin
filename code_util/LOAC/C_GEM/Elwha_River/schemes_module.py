"""
Advection and dispersion schemes module (translated from schemes.c)
"""

import math
from config import M, M1, M2, M3, DELTI, DELXI
from variables import U, v, D, Dold2, fl, dispersion, names, cold

def openbound(co, s):
    """Set boundary conditions."""

    if U[2] >= 0.0:
        co[1] = co[1] - (co[1] - v[s]["clb"]) * U[2] * (float(DELTI) / float(DELXI))
    else:
        co[1] = co[1] - (co[3] - co[1]) * U[2] * (float(DELTI) / float(DELXI))

    if U[M1] >= 0.0:
        co[M] = co[M] - (co[M] - co[M2]) * U[M1] * (float(DELTI) / float(DELXI))
    else:
        co[M] = co[M] - (v[s]["cub"] - co[M]) * U[M1] * (float(DELTI) / float(DELXI))
    co[M1] = co[M]

def tvd(co, s, t):
    """Total Variation Diminishing (TVD) advection scheme."""
    b = 2.0
    cold = co

    for j in range(1, M2 + 1, 2):
        vx = U[j + 1]
        cfl = abs(vx) * (float(DELTI) / (2.0 * float(DELXI)))

        # For positive velocity, vx > 0
        if vx > 0.0:
            f = cold[j + 2] - cold[j]
            if abs(f) > 1e-35 and j != 1:
                rg = (cold[j] - cold[j - 2]) / f
                phi = (2.0 - cfl + rg * (1.0 + cfl)) / 3.0
                philen = max(0.0, min(b, min(b * rg, phi)))
            else:
                philen = 0.0

            co[j + 1] = cold[j] + 0.5 * (1.0 - cfl) * philen * f
            fl[j + 1] = vx * D[j + 1] * co[j + 1]

        # For negative velocity, vx <= 0
        if vx <= 0.0:
            f = cold[j] - cold[j + 2]
            if abs(f) > 1e-35 and j != M2:
                if j + 4 > M:
                    rg = (cold[j + 2] - cold[j + 3]) / f
                else:
                    rg = (cold[j + 2] - cold[j + 4]) / f
                phi = (2.0 - cfl + rg * (1.0 + cfl)) / 3.0
                philen = max(0.0, min(b, min(b * rg, phi)))
            else:
                philen = 0.0

            co[j + 1] = cold[j + 2] + 0.5 * (1.0 - cfl) * philen * f
            fl[j + 1] = vx * D[j + 1] * co[j + 1]

    for j in range(3, M2 + 1, 2):
        rat = Dold2[j] / D[j]
        rat1 = (float(DELTI) / (2.0 * float(DELXI))) / D[j]
        cold[j] = rat * cold[j] - rat1 * (fl[j + 1] - fl[j - 1])
        co[j] = cold[j]

def disp_sch(co):
    """Dispersion scheme."""
    a, b, c, r, di, gam = [[0.0] * (M + 1) for _ in range(6)]
    a[1] = b[1] = 1.0
    di[1] = co[1]
    a[M] = b[M] = 1.0
    di[M] = co[M]

    for i in range(2, M1 + 1):
        c1 = dispersion[i - 1] * D[i - 1] / D[i]
        c2 = dispersion[i + 1] * D[i + 1] / D[i]
        a[i] = -c1 * (float(DELTI) / (8.0 * float(DELXI)**2))
        c[i] = -c2 * (float(DELTI) / (8.0 * float(DELXI)**2))
        b[i] = 1.0 + c2 * (float(DELTI) / (8.0 * float(DELXI)**2)) + c1 * (float(DELTI) / (8.0 * float(DELXI)**2))
        r[i] = 1.0 - c2 * (float(DELTI) / (8.0 * float(DELXI)**2)) - c1 * (float(DELTI) / (8.0 * float(DELXI)**2))

    for i in range(3, M1 + 1, 2):
        if i + 2 > M:
            di[i] = -c[i] * co[i + 1] + r[i] * co[i] - a[i] * co[i - 2]
        else:
            di[i] = -c[i] * co[i + 2] + r[i] * co[i] - a[i] * co[i - 2]

    if b[1] == 0.0:
        raise ValueError("Error: b[1] cannot be 0")

    bet = b[1]
    co[1] = di[1] / bet

    for j in range(3, M1 + 1, 2):
        gam[j] = c[j - 2] / bet
        bet = b[j] - a[j] * gam[j]
        co[j] = (di[j] - a[j] * co[j - 2]) / bet

    for j in range(M3, 0, -2):
        co[j] = co[j] - gam[j + 2] * co[j + 2]
