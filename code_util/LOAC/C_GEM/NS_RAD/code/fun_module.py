"""
Utility functions module (translated from fun.c)
"""

import math
from config import (PI, AMPL, pfun, WARMUP, distance, M, pH_ite, G, mass_mol_B,
                    kISS, wMAX, DISPERSION_MODEL, DISP_MAX,
                    TIDE_AMP, TIDE_SPEED, TIDE_PHASE)
from variables import v, U, DEPTH, kflow, kwind, vp, Hplus, B, B_thread, Chezy, dispersion
from density import dens
from numba import njit
import random

def D_O2(t, water_temp):
    """Molecular diffusion coefficient for O2 [m^2/s]."""
    abstemp = Tabs(t, water_temp)
    d_o2 = (6.35 * abstemp - 1664.0) * 1.0e-11
    return d_o2

def Sc(t, i, water_temp):
    """Schmidt number [-]."""
    sc0 = 1800.6 - 120.1 * water_temp + 3.7818 * (water_temp**2) - 0.047608 * (water_temp**3)
    sc = sc0 * (1 + 3.14e-3 * v['S']['c'][i])
    return sc

def piston_velocity(t, Uw_sal, Uw_tid, T):
    """
    Calculate piston velocity for O2 exchange across air-water interface.

    `T` is now the transported temperature FIELD (variables.v['T']['c']), indexed per
    grid point, rather than the single scalar that used to be applied uniformly.

    The per-cell loop is @njit-compiled (_piston_velocity_loop, with the scalar D_O2
    and Sc helpers jitted as _d_o2 / _sc); it ran every timestep from biogeo and was
    ~0.16 s of pure-Python time in the phase-2 profile. The Python D_O2/Sc above are
    kept unchanged; they were used only here.
    """
    _piston_velocity_loop(U, DEPTH, T, v['S']['c'], kflow, kwind, vp,
                          Uw_sal, Uw_tid, M, distance)

# Wind-driven storm-surge offset [m] added to the mouth elevation. Refreshed each
# timestep from main.py (the daily surge series, tools/build_surge.py) when the active
# site sets SURGE_FILE; otherwise stays 0 and the mouth sees the harmonic tide only.
# The Beaufort coast is microtidal (~0.3 m), so this non-tidal residual -- surges reach
# 1-3 m in big storm years, +0.66 m in 2022 -- is the dominant saltwater-intrusion
# driver, which the astronomical tide alone omits.
SURGE = 0.0


def set_surge(eta_surge):
    """Set the module-level wind-surge offset [m] that Tide() adds to the harmonic elevation (called each timestep from main.py when the site has a SURGE_FILE)."""
    global SURGE
    SURGE = float(eta_surge)


def Tide(t):
    """
    Sea-surface elevation at the mouth [m] = harmonic tide + wind-driven surge.

    Multi-constituent harmonic tide from the nearest NOAA CO-OPS station (per-site
    constituents in config; built by tools/build_tides.py):

        eta(t) = sum_i  A_i * cos( speed_i * (t/3600) - G_i )  +  SURGE

    A_i [m], speed_i [deg/hr], G_i [deg]. The phase reference is the station's GMT
    equilibrium argument; for this idealised climatology run the absolute calendar
    phase is not pinned, but the amplitudes, frequencies and RELATIVE phases (hence
    spring-neap and diurnal inequality) are real. SURGE is the observed daily-mean
    sea-level residual (0 unless the site provides a surge forcing).

    Falls back to the shipped single sinusoid if no constituent data is loaded.
    """
    if TIDE_AMP:
        th = t / 3600.0                         # hours
        eta = 0.0
        for a, sp, g in zip(TIDE_AMP, TIDE_SPEED, TIDE_PHASE):
            eta += a * math.cos(math.radians(sp * th - g))
        return eta + SURGE
    # fallback: single-constituent idealisation
    omega = 2.0 * PI / 3600.0
    return 0.5 * AMPL * math.sin(pfun * omega * t) + SURGE

def Tabs(t, water_temp):
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

def Fhet(t, water_temp):
    """Calculate temperature dependence for heterotrophic reactions."""
    a = (Tabs(t, water_temp) - 278.0) / 10.0
    f = 2.75**a
    return f

def Fnit(t, water_temp):
    """Calculate temperature dependence for nitrification."""
    a = (Tabs(t, water_temp) - 278.0) / 10.0
    f = 5.0**a
    return f

def O2sat(t, i, water_temp):
    """Calculate O2 saturation [mmol O2 m^-3]."""
    abstemp = Tabs(t, water_temp)
    lno2sat = (-1.3529996 * 100 + 157228.8 / abstemp - 66371490 / (abstemp**2)
               + 12436780000 / (abstemp**3) - 8621061 * 100000 / (abstemp**4))
    f = -0.020573 + 12.142 / abstemp - 2363.1 / (abstemp**2)
    o2 = math.exp(lno2sat + f * v['S']['c'][i])
    return o2

def p_bar(t, i, water_temp):
    """Mid-column pressure [bar] at cell `i` (non-jitted form, reads DEPTH/salinity globally); the jitted biogeo path uses _pbar_rho instead."""
    # water density approx [kg/m**3] based on T, S, D/P in the middle of the water column
    density = dens(v['S']['c'][i], water_temp, DEPTH[i] / 2)  # (assuming depth=p_bar)
    p = DEPTH[i] / 2 * density * G  # N/m^2 (10e-4 dbar)
    p = p * 0.1  # bar
    return p

def K0_CO2(t, i, water_temp):
    """ Calculate K0 Henry's constant for CO2 exchanges as a function of temperature and salinity mmol m^−3 atm^−1 or
    mol kg-1 atm-1 (Regnier et al, 2013)"""
    tk = Tabs(t, water_temp)
    tk100 = tk / 100
    tk1002 = tk100 * tk100
    KO_TS = math.exp(93.4517 / tk100 - 60.2409 + 23.3585 * math.log(tk100) + v['S']['c'][i] *
                     (0.023517 - 0.023656 * tk100 + 0.0047036 * tk1002))
    return KO_TS

def K1_CO2(t, i, water_temp):
    """Estimate first dissociation constant of carbonic acid K1
! in mol/kg-SW on the SWS pH-scale (Miller 1995) from solvesaphe (Munhoven, 2013)"""
    # terms used more than once for:
    # temperature
    tk = Tabs(t, water_temp)
    tk100 = tk / 100
    invtk = 1.0 / tk
    dlogtk = math.log(tk)
    # salinity
    s = v['S']['c'][i]
    s2 = s * s
    # at constant pressure (zero)
    K1_p0 = 10**(-1 * (3670.7 * invtk - 62.008 + 9.7944 * dlogtk - 0.0118 * s + 0.000116 * s2))
    # Pressure correction
    K1 = K1_p0 * math.exp((24.2 - 0.085 * water_temp) * p_bar(t, i, water_temp) / (83.143 * tk))
    # # at constant pressure (zero)
    # K1_p0 = 2.18867 - 2275.0360 / Tabs(t) - 1.468591 * math.log(Tabs(t)) +\
    #         (-0.138681 - 9.33291 / Tabs(t)) * math.sqrt(v['S']['c'][i]) +\
    #         0.0726483 * v['S']['c'][i] - 0.00574938 * v['S']['c'][i] * math.sqrt(v['S']['c'][i])
    # # Pressure correction
    # gasconst_bar_cm3_o_mol_k = 83.14510  # libthdyct
    # zds = v['S']['c'][i] - 34.8  # salinity-34.8
    # zt_degc = Tabs(t) - 273.15  # temperature in celsius
    # zrt = gasconst_bar_cm3_o_mol_k * Tabs(t)  # R*t_k, R in bar*cm3/(mol*K)
    # zdvi = -25.50 - 0.151 * zds + 0.1271 * zt_degc  # volume change for ionization
    # zdki = (-3.08 - 0.578 * zds + 0.0877 * zt_degc) * 1.0e-03  # compressibility change for ionization
    # K1_pp =(-zdvi + zdki * p_bar(t, i) / 2) * p_bar(t, i) / zrt
    # # Final K1 value
    # K1 = math.exp(K1_p0 + K1_pp)
    return K1

def K2_CO2(t, i, water_temp):
    """Estimate second dissociation constant of carbonic acid K2
! in mol/kg-SW on the SWS pH-scale (Miller 1995) from Darwin"""
    # terms used more than once for:
    # temperature
    tk = Tabs(t, water_temp)
    tk100 = tk / 100
    invtk = 1.0 / tk
    # salinity
    s = v['S']['c'][i]
    s2 = s * s
    # at constant pressure (zero)
    K2_p0 = 10**(-1 * (1394.7 * invtk + 4.777 - 0.0184 * s + 0.000118 * s2))
    # Pressure correction
    K2 = K2_p0 * math.exp((16.4 - 0.040 * water_temp) * p_bar(t, i, water_temp) / (83.143 * tk))
    # at constant pressure (zero)
    # K2_p0 = -0.84226 - 3741.1288 / Tabs(t) - 1.437139 * math.log(Tabs(t)) +\
    #         (-0.128417 - 24.41239 / Tabs(t)) * math.sqrt(v['S']['c'][i]) +\
    #         0.1195308 * v['S']['c'][i] - 0.00912840 * v['S']['c'][i] * math.sqrt(v['S']['c'][i])
    # # Pressure correction
    # gasconst_bar_cm3_o_mol_k = 83.14510  # libthdyct
    # zds = v['S']['c'][i] - 34.8  # salinity-34.8
    # zt_degc = Tabs(t) - 273.15  # temperature in celsius
    # zrt = gasconst_bar_cm3_o_mol_k * Tabs(t)  # R*t_k, R in bar*cm3/(mol*K)
    # zdvi = -15.82 + 0.321 * zds - 0.0219 * zt_degc  # volume change for ionization
    # zdki = (1.13 - 0.314 * zds + 0.1475 * zt_degc) * 1.0e-03  # compressibility change for ionization
    # K2_pp =(-zdvi + zdki * p_bar(t, i) / 2) * p_bar(t, i) / zrt
    # # Final K1 value
    # K2 = math.exp(K2_p0 + K2_pp)
    return K2

def KB(t, i, water_temp):
    """ Calculate KB coefficient of dissociation of boric acid in seawater. mol/kg of solution
    Millero, 1995"""
    # terms used more than once for:
    # temperature
    tk = Tabs(t, water_temp)
    invtk = 1.0 / tk
    dlogtk = math.log(tk)
    # salinity
    s = v['S']['c'][i]
    s2 = s * s
    sqrts = math.sqrt(s)
    s15 = s**1.5
    KB = math.exp((-8966.90 - 2890.53 * sqrts - 77.942 * s + 1.728 * s15 - 0.0996 * s2) *
                  invtk + (148.0248 + 137.1942 * sqrts + 1.62142 * s) +
                  (-24.4344 - 25.085 * sqrts - 0.2474 * s) * dlogtk + 0.053105 * sqrts * tk)
    # lnKB = (-8966.90 - 2890.51 * v['S']['c'][i]**0.5 - 77.942 *
    #         v['S']['c'][i] + 1.726 * v['S']['c'][i]**1.5 - 0.0993 *
    #         v['S']['c'][i]**2) / Tabs(t) + (148.0248 + 137.194 *
    #                                         v['S']['c'][i]**0.5 + 1.62247 *
    #                                         v['S']['c'][i]) + (-24.4344 - 25.085 *
    #                                                            v['S']['c'][i]**0.5 - 0.2474 *
    #                                                            v['S']['c'][i]) * math.log(Tabs(t)) +\
    #        0.053105 * v['S']['c'][i]**0.5 * Tabs(t)
    # KB = math.exp(lnKB)
    return KB

def wISS(t, i, SPM):
    """Inorganic-sediment settling velocity [m/s], wMAX*SPM/(SPM+kISS) (Clark et al. 2022). `t`, `i` are unused (kept for the original call signature)."""
    ws = wMAX * (SPM / (SPM + kISS))
    return ws


# ---------------------------------------------------------------------------
# Jitted scalar forms of the rate / carbonate functions, for biogeo's compiled
# loop. The Python versions above are kept unchanged (they read v['S']['c'][i]
# and DEPTH[i] internally, which cannot be jitted); these take the salinity and
# depth as plain arguments instead. Arithmetic is transcribed verbatim, so the
# two forms agree bit-for-bit -- verified by the output-hash check.
#
# Note Tabs(t, T) is just 273.15 + T; `t` was never used.
# ---------------------------------------------------------------------------

@njit(cache=True)
def _fhet(water_temp):
    """Heterotrophic-respiration temperature factor (Q10 = 2.75, ref 5 degC), evaluated at the local water temperature."""
    return 2.75 ** (((273.15 + water_temp) - 278.0) / 10.0)


@njit(cache=True)
def _fnit(water_temp):
    """Nitrification temperature factor (Q10 = 5, ref 5 degC), evaluated at the local water temperature."""
    return 5.0 ** (((273.15 + water_temp) - 278.0) / 10.0)


@njit(cache=True)
def _o2sat(water_temp, S):
    """O2 saturation concentration [mmol/m^3] from temperature and salinity (Garcia & Gordon-type fit)."""
    abstemp = 273.15 + water_temp
    lno2sat = (-1.3529996 * 100 + 157228.8 / abstemp - 66371490 / (abstemp**2)
               + 12436780000 / (abstemp**3) - 8621061 * 100000 / (abstemp**4))
    f = -0.020573 + 12.142 / abstemp - 2363.1 / (abstemp**2)
    return math.exp(lno2sat + f * S)


@njit(cache=True)
def _p_bar(S, water_temp, depth):
    """Mid-column pressure [bar] from salinity, temperature and depth (jitted; density at mid-depth via density.dens)."""
    density = dens(S, water_temp, depth / 2)
    p = depth / 2 * density * G
    return p * 0.1


@njit(cache=True)
def _k0(water_temp, S):
    """Henry's constant K0 for CO2 [mol kg^-1 atm^-1] (Weiss 1974), from temperature and salinity."""
    tk = 273.15 + water_temp
    tk100 = tk / 100
    tk1002 = tk100 * tk100
    return math.exp(93.4517 / tk100 - 60.2409 + 23.3585 * math.log(tk100) + S *
                    (0.023517 - 0.023656 * tk100 + 0.0047036 * tk1002))


@njit(cache=True)
def _k1(water_temp, S, pb):
    """First carbonic-acid dissociation constant K1 [mol/kg], pressure-corrected at `pb` bar (Millero 1995)."""
    tk = 273.15 + water_temp
    invtk = 1.0 / tk
    dlogtk = math.log(tk)
    s2 = S * S
    K1_p0 = 10**(-1 * (3670.7 * invtk - 62.008 + 9.7944 * dlogtk - 0.0118 * S + 0.000116 * s2))
    return K1_p0 * math.exp((24.2 - 0.085 * water_temp) * pb / (83.143 * tk))


@njit(cache=True)
def _k2(water_temp, S, pb):
    """Second carbonic-acid dissociation constant K2 [mol/kg], pressure-corrected at `pb` bar (Millero 1995)."""
    tk = 273.15 + water_temp
    invtk = 1.0 / tk
    s2 = S * S
    K2_p0 = 10**(-1 * (1394.7 * invtk + 4.777 - 0.0184 * S + 0.000118 * s2))
    return K2_p0 * math.exp((16.4 - 0.040 * water_temp) * pb / (83.143 * tk))


@njit(cache=True)
def _kb(water_temp, S):
    """Boric-acid dissociation constant KB [mol/kg] from temperature and salinity (Millero 1995)."""
    tk = 273.15 + water_temp
    invtk = 1.0 / tk
    dlogtk = math.log(tk)
    s2 = S * S
    sqrts = math.sqrt(S)
    s15 = S**1.5
    return math.exp((-8966.90 - 2890.53 * sqrts - 77.942 * S + 1.728 * s15 - 0.0996 * s2) *
                    invtk + (148.0248 + 137.1942 * sqrts + 1.62142 * S) +
                    (-24.4344 - 25.085 * sqrts - 0.2474 * S) * dlogtk + 0.053105 * sqrts * tk)


@njit(cache=True)
def _pbar_rho(S, water_temp, depth):
    """Mid-column pressure [bar] AND in-situ density [kg/m^3], from ONE density.dens
    call (the expensive part of the carbonate stack). Returns (pbar, rho). Replaces the
    separate _p_bar so the density needed for the mmol/m^3 <-> mol/kg conversion below
    is not computed twice."""
    rho = dens(S, water_temp, depth / 2.0)
    # Mid-column hydrostatic pressure rho*g*h [Pa] -> bar (the Millero 1995 K1/K2
    # pressure corrections take pressure in bar). Pa->bar is 1e-5; the original *0.1
    # was a 1e4x unit error that gave pb ~ 657 bar at ~1.3 m depth instead of ~0.066,
    # inflating K1/K2/KB ~1.9x and biasing the DIAGNOSED pH ~0.3-0.5 low. FCO2 is
    # nearly unaffected -- the K1/K2 inflation is absorbed by the compensating H+
    # shift that conserves alkalinity -- so this corrects the pH field, not the flux.
    return (depth / 2.0 * rho * G) * 1.0e-5, rho


@njit(cache=True)
def _h_solve_kg(s, H_init, dic_kg, alk_kg, KB, K1, K2):
    """H+ [mol/kg] from DIC and alkalinity, Follows et al. (2006), solved consistently
    in mol/kg. This is the UNIT-CORRECT carbonate solve (ported from upstream C-GEM v2):
    DIC/ALK must be passed in mol/kg -- biogeo converts the mmol/m^3 state via the local
    density -- so they share the mol/kg basis of K1/K2/KB, unlike the legacy `pH()` which
    mixed mmol/m^3 state with mol/kg constants and patched only the borate term
    (CARB_BT_SCALE). Total boron follows salinity (0.1336*S mg/kg; Lee et al. 2010),
    here in mol/kg with no scale factor.

    The three inherited-defect fixes from `pH()` are retained: the H+ guess is hoisted
    out of the loop so the scheme actually iterates; a degenerate (near-zero) carbonate
    state returns a neutral pH 7 rather than a huge/negative H+; and any H+ implying a
    pH outside [2, 12] is rejected as non-physical (a bounded ice-out transient).
    """
    NEUTRAL = 1.0e-7          # pH 7 fallback for an ill-posed solve
    H_MIN, H_MAX = 1.0e-12, 1.0e-2   # pH 12 / pH 2 -- physical H+ bounds
    # 1 mmol/m^3 ~ 1e-6 mol/kg: reject the tiny-positive DIC/ALK left by transport just
    # after ice-out, exactly as the mmol/m^3 `pH()` rejected values below 1 mmol/m^3.
    if dic_kg < 1.0e-6 or alk_kg < 1.0e-6:
        return NEUTRAL
    bt = 0.1336 * s * 1.0e-3 / mass_mol_B          # total boron [mol/kg]
    hg = H_init if H_init > 0.0 else 1.0e-8
    for _ in range(1, pH_ite):
        bohg = bt * KB / (hg + KB)                 # borate share of alkalinity
        cag = alk_kg - bohg                        # carbonate alkalinity [mol/kg]
        if cag <= 0.0:
            return NEUTRAL
        gamm = dic_kg / cag
        dummy = (1.0 - gamm) * (1.0 - gamm) * K1 * K1 - 4.0 * K1 * K2 * (1.0 - 2.0 * gamm)
        if dummy < 0.0:
            return NEUTRAL                         # no real root
        hg = 0.5 * ((gamm - 1.0) * K1 + math.sqrt(dummy))
        if not (H_MIN < hg < H_MAX):
            return NEUTRAL                         # non-physical pH -> neutral
    return hg


# Jitted scalar helpers for the piston-velocity loop. Verbatim transcriptions of the
# Python D_O2 / Sc above, which read v['S']['c'][i] internally; these take salinity as
# a plain argument. (Tabs(t, T) is just 273.15 + T; `t` was never used.)
@njit(cache=True)
def _d_o2(water_temp):
    """Molecular diffusivity of O2 in water [m^2/s] as a function of temperature (for the piston velocity)."""
    abstemp = 273.15 + water_temp
    return (6.35 * abstemp - 1664.0) * 1.0e-11


@njit(cache=True)
def _sc(water_temp, S):
    """Schmidt number for CO2/O2 gas transfer as a function of temperature (and salinity), for the piston velocity."""
    sc0 = 1800.6 - 120.1 * water_temp + 3.7818 * (water_temp**2) - 0.047608 * (water_temp**3)
    return sc0 * (1 + 3.14e-3 * S)


@njit(cache=True)
def _sc_ch4(water_temp):
    """Schmidt number of CH4 in freshwater (Wanninkhof 2014, Table 1)."""
    t = water_temp
    return 1909.4 - 120.78 * t + 4.1555 * t * t - 0.080578 * t**3 + 0.00065777 * t**4


@njit(cache=True)
def _sc_n2o(water_temp):
    """Schmidt number of N2O in freshwater (Wanninkhof 2014, Table 1)."""
    t = water_temp
    return 2141.2 - 152.56 * t + 5.8963 * t * t - 0.12411 * t**3 + 0.0010655 * t**4


@njit(cache=True)
def _ch4_eq(water_temp, S, p_ch4):
    """Dissolved CH4 [mmol m^-3] in equilibrium with a moist atmosphere of mole
    fraction `p_ch4` [atm], from Wiesenburg & Guinasso (1979). Their fit returns
    nmol L^-1 (= 1e-3 mmol m^-3)."""
    tk = 273.15 + water_temp
    A1, A2, A3, A4 = -415.2807, 596.8104, 379.2599, -62.0757
    B1, B2, B3 = -0.059160, 0.032174, -0.0048198
    lnC = (A1 + A2 * (100.0 / tk) + A3 * math.log(tk / 100.0) + A4 * (tk / 100.0)
           + S * (B1 + B2 * (tk / 100.0) + B3 * (tk / 100.0) ** 2))
    c_nmol_l = p_ch4 * math.exp(lnC)          # nmol L^-1
    return c_nmol_l * 1.0e-3                   # -> mmol m^-3


@njit(cache=True)
def _n2o_eq(water_temp, S, p_n2o):
    """Dissolved N2O [mmol m^-3] in equilibrium with a moist atmosphere of mole
    fraction `p_n2o` [atm], from Weiss & Price (1980). Their fit returns mol L^-1
    atm^-1 (= 1e6 mmol m^-3 atm^-1)."""
    tk = 273.15 + water_temp
    A1, A2, A3, A4 = -165.8806, 222.8743, 92.0792, -1.48425
    B1, B2, B3 = -0.056235, 0.031619, -0.0048472
    lnF = (A1 + A2 * (100.0 / tk) + A3 * math.log(tk / 100.0) + A4 * (tk / 100.0) ** 2
           + S * (B1 + B2 * (tk / 100.0) + B3 * (tk / 100.0) ** 2))
    c_mol_l = p_n2o * math.exp(lnF)           # mol L^-1
    return c_mol_l * 1.0e6                     # -> mmol m^-3


@njit(cache=True)
def _piston_velocity_loop(U, DEPTH, T, c_S, kflow, kwind, vp, Uw_sal, Uw_tid, M, distance):
    """Jitted piston (gas-transfer) velocity per cell: combine a wind component and a
    current-shear component into vp [m/s], used for the O2 and CO2 air-water exchange."""
    for i in range(M + 1):
        kflow[i] = math.sqrt(abs(U[i]) * _d_o2(T[i]) / DEPTH[i])
        if i <= distance:
            kwind[i] = (1.0 / 3.6e5) * 0.31 * (Uw_sal**2) * (_sc(T[i], c_S[i]) / 660)**-0.5
        else:
            kwind[i] = (1.0 / 3.6e5) * 0.31 * (Uw_tid**2) * (_sc(T[i], c_S[i]) / 660)**-0.5
        vp[i] = kflow[i] + kwind[i]


_disp_capped = [False]   # report the numerical cap once per run, not per timestep


@njit(cache=True)
def _river_dispersion_loop(B, B_thread, DEPTH, Chezy, dispersion, q, M, G, DISP_MAX):
    """Jitted longitudinal dispersion coefficient per cell, Seo & Cheong (1998) from local
    width/depth/velocity with u* from the Chezy friction, capped at DISP_MAX for scheme stability.

    TWO widths, because they play different roles. The velocity is a CONVEYANCE quantity
    -- total discharge through the total cross-section -- so it uses the total width B.
    The aspect ratio W/H in the Seo & Cheong regression is a WITHIN-channel shear
    property, so it uses the per-thread width B_thread. For a single-thread channel the
    two arrays are identical and this is the original expression exactly.

    (Note the two choices are consistent: if n threads share the depth, the per-thread
    velocity q_i/(W_i*H) equals the conveyance velocity q/(B*H), so there is no separate
    per-thread velocity to compute.)"""
    # Returns the FIRST K that exceeded DISP_MAX (0.0 if none), so the Python wrapper
    # can emit the once-per-run warning.
    capped_K = 0.0
    for i in range(M + 1):
        W = B[i] if B[i] > 1e-3 else 1e-3
        Wt = B_thread[i] if B_thread[i] > 1e-3 else 1e-3
        H = DEPTH[i] if DEPTH[i] > 1e-3 else 1e-3
        U = q / (W * H)
        if U < 1e-9:
            dispersion[i] = 0.0
            continue
        ustar = math.sqrt(G) * U / Chezy[i]
        K = 5.915 * (Wt / H) ** 0.620 * (U / ustar) ** 1.428 * H * ustar
        if K > DISP_MAX:
            if capped_K == 0.0:
                capped_K = K
            K = DISP_MAX
        dispersion[i] = K
    return capped_K


def river_dispersion(Qr):
    """
    Longitudinal dispersion coefficient [m^2/s] at every grid point, written in
    place into variables.dispersion.

    Seo & Cheong (1998), regressed on 59 river data sets:

        K = 5.915 * (W/H)^0.620 * (U/u*)^1.428 * H * u*

    Shear velocity comes from the model's own friction rather than a channel slope:
    Chezy gives U = C*sqrt(R*S), so u* = sqrt(g*R*S) = sqrt(g)*U/C. This keeps the
    dispersion consistent with the momentum closure and needs no new parameter --
    which matters here because SWORD's slopes at these delta reaches are unusable
    (zeros and values as absurd as 2.5 m/m).

    Note U/u* = Chezy/sqrt(g) identically, so the velocity ratio is set by friction
    alone; K then scales with U, and with the width-to-depth ratio.

    Replaces the Savenije/Van der Burgh estuary form, which collapsed to zero beyond
    the first few grid points under observation-based geometry. See
    config.DISPERSION_MODEL for the numbers and the reasoning.
    """
    if DISPERSION_MODEL != "seo":
        raise SystemExit(f"unknown DISPERSION_MODEL={DISPERSION_MODEL!r}")

    capped_K = _river_dispersion_loop(B, B_thread, DEPTH, Chezy, dispersion, abs(Qr), M, G, DISP_MAX)
    if capped_K > 0.0 and not _disp_capped[0]:
        print(f"  [dispersion] Seo-Cheong K={capped_K:.0f} exceeds the numerical "
              f"ceiling DISP_MAX={DISP_MAX:.0f} m2/s; capping. "
              f"Reduce DELTI or coarsen DELXI to relax it.")
        _disp_capped[0] = True