"""
Utility functions module (translated from fun.c)
t is always time vector in the following functions. Use water_temp or Tabs(t) for temperature
"""

import math
from config import PI, AMPL, pfun, Uw_sal, Uw_tid, water_temp, WARMUP, distance, M, pH_ite, G, mass_mol_B
from variables import v, U, DEPTH, kflow, kwind, vp, Hplus
from density import dens
from numba import njit
import random

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

def p_bar(t, i):
    # water density approx [kg/m**3] based on T, S, D/P in the middle of the water column
    density = dens(v['S']['c'][i], Tabs(t) - 273.15, DEPTH[i] / 2)  # (assuming depth=p_bar)
    p = DEPTH[i] / 2 * density * G  # N/m^2 (10e-4 dbar)
    p = p * 0.1  # bar
    return p

def K0_CO2(t, i):
    """ Calculate K0 Henry's constant for CO2 exchanges as a function of temperature and salinity mmol m^−3 atm^−1 or
    mol kg-1 atm-1 (Regnier et al, 2013)"""
    tk = Tabs(t)
    tk100 = tk / 100
    tk1002 = tk100 * tk100
    KO_TS = math.exp(93.4517 / tk100 - 60.2409 + 23.3585 * math.log(tk100) + v['S']['c'][i] *
                     (0.023517 - 0.023656 * tk100 + 0.0047036 * tk1002))
    return KO_TS

def K1_CO2(t, i):
    """Estimate first dissociation constant of carbonic acid K1
! in mol/kg-SW on the SWS pH-scale (Miller 1995) from solvesaphe (Munhoven, 2013)"""
    # terms used more than once for:
    # temperature
    tk = Tabs(t)
    tk100 = tk / 100
    invtk = 1.0 / tk
    dlogtk = math.log(tk)
    # salinity
    s = v['S']['c'][i]
    s2 = s * s
    # at constant pressure (zero)
    K1_p0 = 10**(-1 * (3670.7 * invtk - 62.008 + 9.7944 * dlogtk - 0.0118 * s + 0.000116 * s2))
    # Pressure correction
    K1 = K1_p0 * math.exp((24.2 - 0.085 * water_temp) * p_bar(t, i) / (83.143 * tk))
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

def K2_CO2(t, i):
    """Estimate second dissociation constant of carbonic acid K2
! in mol/kg-SW on the SWS pH-scale (Miller 1995) from Darwin"""
    # terms used more than once for:
    # temperature
    tk = Tabs(t)
    tk100 = tk / 100
    invtk = 1.0 / tk
    # salinity
    s = v['S']['c'][i]
    s2 = s * s
    # at constant pressure (zero)
    K2_p0 = 10**(-1 * (1394.7 * invtk + 4.777 - 0.0184 * s + 0.000118 * s2))
    # Pressure correction
    K2 = K2_p0 * math.exp((16.4 - 0.040 * t) * p_bar(t, i) / (83.143 * tk))
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

def KB(t, i):
    """ Calculate KB coefficient of dissociation of boric acid in seawater. mol/kg of solution
    Millero, 1995"""
    # terms used more than once for:
    # temperature
    tk = Tabs(t)
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

@njit
def pH(t, i, s, Hplus, dic, alk, KB, K1, K2):
    # Estimate pH (Follows et al., 2006)
    bt_g = 0.1336 * s  # Boron concentration mg kg-1 based on salinity ratio from Lee et al., 2010
    bt = bt_g * mass_mol_B / 1000  # Boron concentration mol kg-1
    for ite in range(1, pH_ite):
        # 1 from Follows et al., 2006, iterate through pH for finding CO2aq
        hg = Hplus  # mol kg-1
        # 2 estimate contributions to total alk from boron
        bohg = bt * KB / (hg + KB)  # mol/kg
        # 3 estimate carbonate alkalinity
        fg = - bohg  # mol/kg
        cag = alk + fg  # mol/kg or mmol m-3
        # 4 improved estimate of hydrogen ion conc
        gamm = dic / cag  # mol/kg or mmol m-3
        dummy = (1.0 - gamm) * (1.0 - gamm) * K1 * K1 - 4.0 * K1 * K2 *\
                (1.0 - 2.0 * gamm)
        hg = 0.5 * ((gamm - 1.0) * K1 + math.sqrt(dummy))  # mol/kg or mmol m-3
    return hg