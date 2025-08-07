import numpy as np
from numba import njit, prange
from config import PI, AMPL, pfun, Uw_sal, Uw_tid, water_temp, WARMUP, distance, M
from variables import v, U, DEPTH, kflow, kwind, vp

@njit(cache=True)
def Tabs(t):
    return 273.15 + water_temp

@njit(cache=True)
def D_O2(t):
    abstemp = Tabs(t)
    return (6.35 * abstemp - 1664.0) * 1.0e-11

@njit(cache=True)
def Tide(t):
    omega = 2.0 * PI / 3600.0
    return 0.5 * AMPL * np.sin(pfun * omega * t)

@njit(cache=True)
def Fhet(t):
    a = (Tabs(t) - 278.0) / 10.0
    return 2.75 ** a

@njit(cache=True)
def Fnit(t):
    a = (Tabs(t) - 278.0) / 10.0
    return 5.0 ** a

@njit(cache=True)
def I0(t):
    ncloud = 0.6
    nday = np.ceil(t / (24 * 3600.0)) - 1 - float(WARMUP) / (24 * 3600.0)
    h = np.ceil(t / 3600.0) - 1 - float(WARMUP) / 3600.0 - nday * 24
    r = 16.5
    light = (PI / 2.0) * np.sin(PI * ((t / 3600.0) - (12.0 - r / 2.0)) / r)
    return max(1000 * ncloud * light, 0.0)

@njit(parallel=True, cache=True)
def Sc(t,S):
    temp = Tabs(t) - 273.15
    sc0 = 1800.6 - 120.1 * temp + 3.7818 * temp**2 - 0.047608 * temp**3
    sc = np.empty_like(S)
    for i in prange(S.shape[0]):
        sc[i] = sc0 * (1 + 3.14e-3 * S[i])
    return sc

@njit(parallel=True, cache=True)
def O2sat(t,S):
    abstemp = Tabs(t)
    lno2sat = (-1.3529996 * 100 + 157228.8 / abstemp - 66371490 / (abstemp**2)
               + 12436780000 / (abstemp**3) - 8621061 * 100000 / (abstemp**4))
    f = -0.020573 + 12.142 / abstemp - 2363.1 / (abstemp**2)
    out = np.empty_like(S)
    for i in prange(S.shape[0]):
        out[i] = np.exp(lno2sat + f * S[i])
    return out

@njit(parallel=True, cache=True)
def piston_velocity(t,S,kflow,kwind,vp):
    d_o2 = D_O2(t)
    sc = Sc(t,S)
    for i in prange(M + 1):
        kflow[i] = np.sqrt(np.abs(U[i]) * d_o2 / DEPTH[i])
        wind_speed = Uw_sal if i <= distance else Uw_tid
        kwind[i] = (1.0 / 3.6e5) * 0.31 * wind_speed**2 * (sc[i] / 660.0)**-0.5
        vp[i] = kflow[i] + kwind[i]
