import numpy as np
import math
from numba import njit, prange
from config import PI, AMPL, pfun, Uw_sal, Uw_tid, water_temp, WARMUP, distance, M, USE_CO2_FLUX
from config import G, MOLAR_MASS_B, PH_ITERS, CO2_PISTON_FROM_O2, pCO2, M, DELTI
from density import dens_scalar  # <- njit-friendly core, NOT dens()
if USE_CO2_FLUX:
    from variables import Hplus

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
def piston_velocity(t,S,U_in,DEPTH_in,kflow,kwind,vp):
    d_o2 = D_O2(t)
    sc = Sc(t,S)
    for i in prange(M + 1):
        kflow[i] = np.sqrt(np.abs(U_in[i]) * d_o2 / DEPTH_in[i])
        wind_speed = Uw_sal if i <= distance else Uw_tid
        kwind[i] = (1.0 / 3.6e5) * 0.31 * wind_speed**2 * (sc[i] / 660.0)**-0.5
        vp[i] = kflow[i] + kwind[i]

# ---------- pressure at mid-water (bar) ----------
@njit(parallel=True, cache=True)
def p_bar_all(t, S, depth, out_p):
    """
    Hydrostatic pressure at mid-depth (bar), using a 2-step pressure/density refinement.
    Inputs:
      S[i]     : salinity (PSS-78)
      depth[i] : water depth [m]
    Writes:
      out_p[i] : mid-column pressure [bar]
    """
    T_C = Tabs(t) - 273.15
    n = depth.shape[0]

    # i=0 often unused in your grid; keep it deterministic
    if n > 0:
        out_p[0] = 0.0

    for i in prange(1, n):
        # mid-depth (avoid zero depth)
        h = 0.5 * depth[i]
        if h <= 0.0:
            out_p[i] = 0.0
            continue

        # first guess: density at P=0 dbar
        rho = dens_scalar(S[i], T_C, 0.0)
        P_dbar = rho * G * h / 1.0e4

        # one refinement with pressure-dependent density
        rho = dens_scalar(S[i], T_C, P_dbar)
        P_dbar = rho * G * h / 1.0e4

        # convert dbar -> bar
        out_p[i] = 0.1 * P_dbar

# ---------- Henry's K0 (CO2) ----------
@njit(cache=True)
def _K0_CO2_scalar(tk, s):
    tk100  = tk / 100.0
    tk1002 = tk100 * tk100
    return math.exp(93.4517 / tk100 - 60.2409 + 23.3585 * math.log(tk100)
                    + s * (0.023517 - 0.023656 * tk100 + 0.0047036 * tk1002))

@njit(cache=True)
def _K1_CO2_scalar(tk, tC, s, pbar):
    invtk = 1.0 / tk
    dlogk = math.log(tk)
    s2    = s * s
    K1_p0 = 10.0 ** (-1.0 * (3670.7*invtk - 62.008 + 9.7944*dlogk - 0.0118*s + 0.000116*s2))
    return K1_p0 * math.exp((24.2 - 0.085 * tC) * pbar / (83.143 * tk))

@njit(cache=True)
def _K2_CO2_scalar(tk, tC, s, pbar):
    invtk = 1.0 / tk
    s2    = s * s
    K2_p0 = 10.0 ** (-1.0 * (1394.7*invtk + 4.777 - 0.0184*s + 0.000118*s2))
    return K2_p0 * math.exp((16.4 - 0.040 * tC) * pbar / (83.143 * tk))

@njit(cache=True)
def _KB_scalar(tk, s):
    invtk = 1.0 / tk
    dlogk = math.log(tk)
    sqs   = math.sqrt(s)
    s2    = s * s
    s15   = s * sqs
    return math.exp(
        (-8966.90 - 2890.53*sqs - 77.942*s + 1.728*s15 - 0.0996*s2) * invtk
        + (148.0248 + 137.1942*sqs + 1.62142*s)
        + (-24.4344 - 25.085*sqs - 0.2474*s) * dlogk
        + 0.053105 * sqs * tk
    )

# ---- safe H+ (mol/kg) solve (Follows-style) ----
@njit(cache=True)
def _H_solve(sal, H_init, dic_kg, alk_kg, KB, K1, K2):
    # total boron: 0.1336*S mg/kg -> mol/kg
    bt = (0.1336 * sal) * 1e-3 / MOLAR_MASS_B
    hg = H_init if H_init > 0.0 else 1e-12
    for _ in range(PH_ITERS):
        bohg = bt * KB / (hg + KB)
        cag  = alk_kg - bohg
        if cag < 1e-30: cag = 1e-30
        gamma = dic_kg / cag
        disc = (1.0 - gamma)*(1.0 - gamma)*K1*K1 - 4.0*K1*K2*(1.0 - 2.0*gamma)
        if disc < 0.0: disc = 0.0
        hg_new = 0.5 * ((gamma - 1.0) * K1 + math.sqrt(disc))
        hg = hg_new if hg_new > 0.0 else 1e-12
    return hg

# ---- SEMI-IMPLICIT carbonate step (like SPM fix) ----
@njit(parallel=True, cache=True)
def co2_step_semi_implicit(
    t,
    S, depth, vp,
    DIC_v, ALK_v, pH_arr,             # state in mmol/m^3; pH dimensionless
    NPP_NO3, NPP_NH4, aer_deg, denit, nit,
    FCO2_out                           # output mmol/m^3/s
):
    """
    Semi-implicit DIC update with frozen linearization of FCO2:
      FCO2 ≈ -(vp/depth)/denom * DIC_v + (vp/depth)*(K0_v*pCO2_air)
    ALK explicit. pH recomputed after update. Uses p_bar_all for K1/K2 and rho.
    """
    n  = DIC_v.shape[0]
    tk = Tabs(t)
    tC = tk - 273.15

    # mid-column pressure (bar) for all cells
    pbar = np.zeros(n)
    p_bar_all(t, S, depth, pbar)

    rho   = np.zeros(n)     # kg/m^3
    H     = np.zeros(n)     # mol/kg
    denom = np.zeros(n)     # dimensionless buffer factor
    k     = np.zeros(n)     # s^-1
    b     = np.zeros(n)     # mmol/m^3/s

    # 1) coefficients & density & H from OLD state
    for i in prange(1, n):
        s        = S[i]
        pbar_i   = pbar[i]               # bar
        p_dbar   = 10.0 * pbar_i         # dbar
        # density at mid-depth pressure
        rho_i    = dens_scalar(s, tC, p_dbar)
        rho[i]   = rho_i

        # equilibrium constants in mol/kg basis
        K0 = _K0_CO2_scalar(tk, s)
        K1 = _K1_CO2_scalar(tk, tC, s, pbar_i)
        K2 = _K2_CO2_scalar(tk, tC, s, pbar_i)
        KB = _KB_scalar(tk, s)

        # convert state to mol/kg for speciation
        dic_kg = (DIC_v[i] * 1e-3) / rho_i
        alk_kg = (ALK_v[i] * 1e-3) / rho_i

        # H+ from old state
        H0 = 10.0 ** (-pH_arr[i]) if pH_arr[i] > 0.0 else 1e-8
        Hi = _H_solve(s, H0, dic_kg, alk_kg, KB, K1, K2)
        H[i] = Hi

        # buffer factor
        d = 1.0 + (K1 / Hi) + (K1 * K2 / (Hi * Hi))
        if d < 1e-12: d = 1e-12
        denom[i] = d

        # linearization slope & intercept (per volume)
        vpCO2     = CO2_PISTON_FROM_O2 * vp[i]
        inv_depth = 1.0 / (depth[i] if depth[i] > 1.0e-6 else 1.0e-6)
        k[i]      = vpCO2 * inv_depth / d
        K0_v      = K0 * rho_i * 1e3                    # mol/kg → mmol/m^3
        b[i]      = vpCO2 * inv_depth * (K0_v * pCO2)

    # 2) semi-implicit DIC; explicit ALK (per-volume, mmol/m^3)
    for i in prange(1, n):
        rhs = DIC_v[i] + DELTI * (-NPP_NO3[i] - NPP_NH4[i] + aer_deg[i] + denit[i] + b[i])
        DIC_v[i] = rhs / (1.0 + DELTI * k[i])
        ALK_v[i] += DELTI * ((15.0/106.0)*aer_deg[i] + (93.4/106.0)*denit[i]
                             - 2.0*nit[i] - (15.0/106.0)*NPP_NH4[i] + (17.0/106.0)*NPP_NO3[i])

    # 3) recompute pH from NEW state (mol/kg), reuse pbar for consistency
    for i in prange(1, n):
        s        = S[i]
        pbar_i   = pbar[i]
        p_dbar   = 10.0 * pbar_i
        rho_i    = rho[i]  # already at mid-depth pressure

        dic_kg = (DIC_v[i] * 1e-3) / rho_i
        alk_kg = (ALK_v[i] * 1e-3) / rho_i

        H0   = 10.0 ** (-pH_arr[i]) if pH_arr[i] > 0.0 else 1e-8
        H[i] = _H_solve(
            s, H0, dic_kg, alk_kg,
            _KB_scalar(tk, s),
            _K1_CO2_scalar(tk, tC, s, pbar_i),
            _K2_CO2_scalar(tk, tC, s, pbar_i)
        )
        pH_arr[i] = -math.log10(H[i] if H[i] > 1e-12 else 1e-12)

    # 4) diagnostic FCO2 from NEW DIC and frozen (k,b)
    for i in prange(1, n):
        FCO2_out[i] = -k[i] * DIC_v[i] + b[i]