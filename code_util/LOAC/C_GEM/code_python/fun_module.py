"""
Utility functions module (translated from fun.c)
"""

import math
from config import PI, AMPL, pfun, Uw_sal, Uw_tid, water_temp, WARMUP, distance, M
from variables import v, U, DEPTH, kflow, kwind, vp

def D_O2(t):
    """Molecular diffusion coefficient for O2 [m^2/s]."""
    abstemp = Tabs(t)
    d_o2 = (6.35 * abstemp - 1664.0) * 1.0e-11
    return d_o2

def Sc(t, i):
    """Schmidt number [-]."""
    temp = Tabs(t) - 273.15
    sc0 = 1800.6 - 120.1 * temp + 3.7818 * (temp**2) - 0.047608 * (temp**3)
    sc = sc0 * (1 + 3.14e-3 * v['S']['c'][i])
    return sc

def piston_velocity(t):
    """Calculate piston velocity for O2 exchange across air-water interface."""
    for i in range(M + 1):
        kflow[i] = math.sqrt(abs(U[i]) * D_O2(t) / DEPTH[i])
        if i <= distance:
            kwind[i] = (1.0 / 3.6e5) * 0.31 * (Uw_sal**2) * (Sc(t, i) / 660)**-0.5
        else:
            kwind[i] = (1.0 / 3.6e5) * 0.31 * (Uw_tid**2) * (Sc(t, i) / 660)**-0.5
        vp[i] = kflow[i] + kwind[i]

def Tide(t):
    """Calculate tidal elevation [m]."""
    omega = 2.0 * PI / 3600.0
    v0 = 0.5 * AMPL * math.sin(pfun * omega * t)
    return v0

def Tabs(t):
    """Calculate absolute temperature [K]."""
    abstemp = 273.15 + water_temp
    return abstemp

def I0(t):
    """Calculate light intensity [muE m^-2 s^-1]."""
    ncloud = 0.6
    nday = math.ceil(t / (24 * 60 * 60)) - 1 - float(WARMUP) / (24 * 60 * 60)
    h = math.ceil(t / (60 * 60)) - 1 - float(WARMUP) / (60 * 60) - nday * 24
    r = 16.5
    light = (PI / 2.0) * math.sin(PI * ((t / (60 * 60)) - (12.0 - r / 2.0)) / r)
    irr = max(1000 * ncloud * light, 0.0)
    return irr

def Fhet(t):
    """Calculate temperature dependence for heterotrophic reactions."""
    a = (Tabs(t) - 278.0) / 10.0
    f = 2.75**a
    return f

def Fnit(t):
    """Calculate temperature dependence for nitrification."""
    a = (Tabs(t) - 278.0) / 10.0
    f = 5.0**a
    return f

def O2sat(t, i):
    """Calculate O2 saturation [mmol O2 m^-3]."""
    abstemp = Tabs(t)
    lno2sat = (-1.3529996 * 100 + 157228.8 / abstemp - 66371490 / (abstemp**2)
               + 12436780000 / (abstemp**3) - 8621061 * 100000 / (abstemp**4))
    f = -0.020573 + 12.142 / abstemp - 2363.1 / (abstemp**2)
    o2 = math.exp(lno2sat + f * v['S']['c'][i])
    return o2
