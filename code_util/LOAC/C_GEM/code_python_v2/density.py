# --- coding: utf-8 ---
"""
Seawater density (UNESCO 1983) — Numba-optimized.

API (unchanged for scalars):
    dens(S, T, P_dbar)  ->  rho [kg/m^3]

Extras (vector kernels):
    dens_all(S, T, P_dbar, out_rho)
    svan_all(...), sigma_all(...)

All vector kernels expect 1D arrays of the same length and write into 'out_' arrays.
"""

import math
import numpy as np
from numba import njit, prange

# --------------------- scalar core (njit) ---------------------

@njit(cache=True)
def _dens0_scalar(S, T):
    """Density at zero pressure."""
    # pure water
    a0 = 999.842594
    a1 = 6.793952e-2
    a2 = -9.095290e-3
    a3 = 1.001685e-4
    a4 = -1.120083e-6
    a5 = 6.536332e-9
    SMOW = a0 + (a1 + (a2 + (a3 + (a4 + a5*T)*T)*T)*T)*T

    # salinity polynomials
    b0 = 8.24493e-1
    b1 = -4.0899e-3
    b2 = 7.6438e-5
    b3 = -8.2467e-7
    b4 = 5.3875e-9
    RB = b0 + (b1 + (b2 + (b3 + b4*T)*T)*T)*T

    c0 = -5.72466e-3
    c1 = 1.0227e-4
    c2 = -1.6546e-6
    RC = c0 + (c1 + c2*T)*T

    d0 = 4.8314e-4

    return SMOW + RB*S + RC*(S**1.5) + d0*S*S

@njit(cache=True)
def _seck_scalar(S, T, P_bar):
    """Secant bulk modulus K(S,T,P[bar]). P must be in bar."""
    # pure water terms
    h0 = 3.239908;  h1 = 1.43713e-3; h2 = 1.16092e-4; h3 = -5.77905e-7
    AW = h0 + (h1 + (h2 + h3*T)*T)*T

    k0 = 8.50935e-5; k1 = -6.12293e-6; k2 = 5.2787e-8
    BW = k0 + (k1 + k2*T)*T

    e0 = 19652.21; e1 = 148.4206; e2 = -2.327105; e3 = 1.360477e-2; e4 = -5.155288e-5
    KW = e0 + (e1 + (e2 + (e3 + e4*T)*T)*T)*T

    # seawater, P = 0
    SR = math.sqrt(S)
    i0 = 2.2838e-3; i1 = -1.0981e-5; i2 = -1.6078e-6; j0 = 1.91075e-4
    A = AW + (i0 + (i1 + i2*T)*T + j0*SR)*S

    f0 = 54.6746; f1 = -0.603459; f2 = 1.09987e-2; f3 = -6.1670e-5
    g0 = 7.944e-2; g1 = 1.6483e-2; g2 = -5.3009e-4
    K0 = KW + (f0 + (f1 + (f2 + f3*T)*T)*T + (g0 + (g1 + g2*T)*T)*SR)*S

    # general expression
    m0 = -9.9348e-7; m1 = 2.0816e-8; m2 = 9.1697e-10
    B = BW + (m0 + (m1 + m2*T)*T)*S

    return K0 + (A + B*P_bar)*P_bar

@njit(cache=True)
def dens_scalar(S, T, P_dbar):
    """
    UNESCO density with pressure P in dbar.
    Returns rho [kg/m^3].
    """
    P_bar = 0.1 * P_dbar  # convert dbar -> bar
    K = _seck_scalar(S, T, P_bar)
    rho0 = _dens0_scalar(S, T)
    denom = 1.0 - P_bar / K
    return rho0 / denom

# --------------------- vector kernels (njit parallel) ---------------------

@njit(parallel=True, cache=True)
def dens_all(S, T, P_dbar, out_rho):
    n = out_rho.shape[0]
    for i in prange(n):
        out_rho[i] = dens_scalar(S[i], T[i], P_dbar[i])

@njit(parallel=True, cache=True)
def sigma_all(S, T, P_dbar, out_sigma):
    n = out_sigma.shape[0]
    for i in prange(n):
        out_sigma[i] = dens_scalar(S[i], T[i], P_dbar[i]) - 1000.0

@njit(parallel=True, cache=True)
def svan_all(S, T, P_dbar, out_svan):
    # specific volume anomaly = 1/rho(S,T,P) - 1/rho(35,0,P)
    n = out_svan.shape[0]
    for i in prange(n):
        rho = dens_scalar(S[i], T[i], P_dbar[i])
        rho_std = dens_scalar(35.0, 0.0, P_dbar[i])
        out_svan[i] = 1.0 / rho - 1.0 / rho_std

# --------------------- Python convenience wrappers ---------------------

def dens(S, T, P_dbar=0.0):
    """
    Python-friendly wrapper.
    Accepts scalars or 1D arrays. For arrays, returns a new array.
    """
    if np.isscalar(S):
        return float(dens_scalar(float(S), float(T), float(P_dbar)))
    S = np.asarray(S, dtype=np.float64)
    T = np.asarray(T, dtype=np.float64)
    P = np.asarray(P_dbar, dtype=np.float64)
    out = np.empty_like(S)
    dens_all(S, T, P, out)
    return out

def sigma(S, T, P_dbar=0.0):
    if np.isscalar(S):
        return dens(S, T, P_dbar) - 1000.0
    S = np.asarray(S, dtype=np.float64)
    T = np.asarray(T, dtype=np.float64)
    P = np.asarray(P_dbar, dtype=np.float64)
    out = np.empty_like(S)
    sigma_all(S, T, P, out)
    return out

def svan(S, T, P_dbar=0.0):
    if np.isscalar(S):
        return 1.0/dens(S, T, P_dbar) - 1.0/dens(35.0, 0.0, P_dbar)
    S = np.asarray(S, dtype=np.float64)
    T = np.asarray(T, dtype=np.float64)
    P = np.asarray(P_dbar, dtype=np.float64)
    out = np.empty_like(S)
    svan_all(S, T, P, out)
    return out
