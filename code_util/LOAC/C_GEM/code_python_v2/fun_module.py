import numpy as np
import math
from numba import njit, prange
from config import PI, AMPL, pfun, Uw_sal, Uw_tid, water_temp, WARMUP, distance, M
from variables import v, U, DEPTH, kflow, kwind, vp
from config import G, MOLAR_MASS_B, PH_ITERS, CO2_PISTON_FROM_O2, pCO2, M, DELTI
from density import dens_scalar  # <- njit-friendly core, NOT dens()

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

# ---------- pressure at mid-water (bar) ----------
@njit(parallel=True, cache=True)
def p_bar_all(t, S, depth, out_p):
    """
    Hydrostatic pressure at mid-depth in *bar*.
    Inputs:
      S[i]      : salinity
      depth[i]  : water depth [m]
    """
    T_C = Tabs(t) - 273.15
    n = depth.shape[0]
    for i in prange(1, n):
        h = 0.5 * depth[i]
        # first guess: rho at P=0
        rho0 = dens_scalar(S[i], T_C, 0.0)
        P_dbar = rho0 * 9.81 * h / 1.0e4
        # one refinement with pressure-dependent rho
        rho = dens_scalar(S[i], T_C, P_dbar)
        P_dbar = rho * 9.81 * h / 1.0e4
        out_p[i] = 0.1 * P_dbar  # convert dbar -> bar

# ---------- Henry's K0 (CO2) ----------
@njit(parallel=True, cache=True)
def K0_CO2(t, S, out_K0):
    """
    Henry constant K0 for CO2 (Regnier-like parametrization)
    Units consistent with your flux formula (mmol m^-3 atm^-1 or mol kg^-1 atm^-1),
    depending on how you treat DIC. This mirrors your scalar K0_CO2.
    """
    tk = Tabs(t)
    tk100 = tk / 100.0
    tk1002 = tk100 * tk100
    ln_part = 93.4517 / tk100 - 60.2409 + 23.3585 * math.log(tk100)
    for i in prange(1, S.shape[0]):
        out_K0[i] = math.exp(ln_part + S[i] * (0.023517 - 0.023656 * tk100 + 0.0047036 * tk1002))

# ---------- K1, K2 with pressure correction ----------
@njit(parallel=True, cache=True)
def K1_CO2(t, S, pbar, out_K1):
    """
    First dissociation constant (Miller 1995 / Darwin notes).
    Uses simple pressure correction with water_temp.
    """
    tk = Tabs(t)
    invtk = 1.0 / tk
    dlogtk = math.log(tk)
    tC = tk - 273.15
    for i in prange(1, S.shape[0]):
        s  = S[i]; s2 = s * s
        K1_p0 = 10.0 ** (-1.0 * (3670.7 * invtk - 62.008 + 9.7944 * dlogtk - 0.0118 * s + 0.000116 * s2))
        # pressure correction (bar-based gas constant 83.143)
        out_K1[i] = K1_p0 * math.exp((24.2 - 0.085 * tC) * pbar[i] / (83.143 * tk))

@njit(parallel=True, cache=True)
def K2_CO2(t, S, pbar, out_K2):
    """
    Second dissociation constant (Miller 1995 / Darwin notes).
    """
    tk = Tabs(t)
    invtk = 1.0 / tk
    tC = tk - 273.15
    for i in prange(1, S.shape[0]):
        s  = S[i]; s2 = s * s
        K2_p0 = 10.0 ** (-1.0 * (1394.7 * invtk + 4.777 - 0.0184 * s + 0.000118 * s2))
        out_K2[i] = K2_p0 * math.exp((16.4 - 0.040 * tC) * pbar[i] / (83.143 * tk))

# ---------- KB (boric acid) ----------
@njit(parallel=True, cache=True)
def KB(t, S, out_KB):
    tk = Tabs(t)
    invtk = 1.0 / tk
    dlogtk = math.log(tk)
    for i in prange(1, S.shape[0]):
        s = S[i]
        sqrts = math.sqrt(s)
        s2 = s * s
        s15 = s * sqrts
        out_KB[i] = math.exp(
            (-8966.90 - 2890.53 * sqrts - 77.942 * s + 1.728 * s15 - 0.0996 * s2) * invtk
            + (148.0248 + 137.1942 * sqrts + 1.62142 * s)
            + (-24.4344 - 25.085 * sqrts - 0.2474 * s) * dlogtk
            + 0.053105 * sqrts * tk
        )

# ---------- pH solver (fixed small number of iterations) ----------
@njit(cache=True)
def _pH_one(sal, H_init, dic, alk, KB_val, K1, K2):
    """
    One-cell pH fixed-point loop (Follows 2006-like).
    Returns hydrogen ion concentration H+ (same units as inputs, mol/kg or mmol m^-3).
    """
    # total boron from salinity: bt [mol/kg]
    # bt_g (mg/kg) = 0.1336 * S  -> mol/kg = bt_g * 1e-3 / MOLAR_MASS_B
    bt = (0.1336 * sal) * 1e-3 / MOLAR_MASS_B
    hg = H_init
    for _ in range(PH_ITERS):
        # borate from boric acid
        bohg = bt * KB_val / (hg + KB_val)  # mol/kg
        # carbonate alkalinity
        cag = alk - bohg
        # gamma = DIC / CA
        gamm = dic / max(cag, 1e-30)
        # quadratic for carbonate system
        K1g = K1
        K2g = K2
        A = (1.0 - gamm) * (1.0 - gamm) * K1g * K1g
        B = 4.0 * K1g * K2g * (1.0 - 2.0 * gamm)
        disc = A - B
        if disc < 0.0:
            disc = 0.0
        hg = 0.5 * ((gamm - 1.0) * K1g + math.sqrt(disc))
        if hg <= 0.0:
            hg = 1e-12
    return hg

@njit(parallel=True, cache=True)
def pH(S, DIC, ALK, KB_v, K1_v, K2_v, H_guess, out_H):
    for i in prange(1, S.shape[0]):
        out_H[i] = _pH_one(S[i], H_guess[i], DIC[i], ALK[i], KB_v[i], K1_v[i], K2_v[i])

# ---------- complete CO2 exchange step (vectorized & in-place) ----------
@njit(parallel=True, cache=True)
def co2_flux_step(t,
                  S, DIC, ALK, pH_arr, depth, vp,
                  NPP_NO3, NPP_NH4, aer_deg, denit, nit,
                  out_FCO2, pCO2_air):
    """
    Computes air–sea CO2 exchange and updates DIC, ALK, pH in-place.
    Writes volumetric flux to out_FCO2 (mmol C m^-3 s^-1).
    All arrays expected length >= M+1 with index 0 unused, just like the rest of the model.
    """
    n = depth.shape[0]
    # work arrays
    pbar = np.zeros(n)
    K0   = np.zeros(n)
    K1   = np.zeros(n)
    K2   = np.zeros(n)
    KBv  = np.zeros(n)
    H    = np.zeros(n)

    # 1) pressure and constants
    p_bar_all(t, S, depth, pbar)
    K0_CO2(t, S, K0)
    K1_CO2(t, S, pbar, K1)
    K2_CO2(t, S, pbar, K2)
    KB(t, S, KBv)

    # 2) pH (solve for H+) using current pH as initial guess
    for i in prange(1, n):
        H[i] = 10.0 ** (-pH_arr[i])
    pH(S, DIC, ALK, KBv, K1, K2, H, H)

    # 3) CO2* and flux
    for i in prange(1, n):
        denom = 1.0 + (K1[i] / H[i]) + (K1[i] * K2[i] / (H[i] * H[i]))
        co2s  = DIC[i] / max(denom, 1e-30)
        vpCO2 = CO2_PISTON_FROM_O2 * vp[i]
        RCO2  = - vpCO2 * (co2s - K0[i] * pCO2_air)   # mmol C m^-2 s^-1
        out_FCO2[i] = RCO2 / max(depth[i], 1e-6)

        # 4) Update tracers (mmol C m^-3 s^-1 * Δt)
        DIC[i] += (out_FCO2[i] - NPP_NO3[i] - NPP_NH4[i] + aer_deg[i] + denit[i]) * DELTI
        ALK[i] += (((15.0/106.0) * aer_deg[i]) + ((93.4/106.0) * denit[i]) - (2.0 * nit[i])
                   - ((15.0/106.0) * NPP_NH4[i]) + ((17.0/106.0) * NPP_NO3[i])) * DELTI
        # 5) store pH
        if not (H[i] > 0.0) or not math.isfinite(H[i]):
            # fallback to previous pH (or a default) and small H+
            H[i] = 1e-12
        pH_arr[i] = -math.log10(max(H[i], 1e-12))