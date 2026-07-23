"""
Surface heat budget for the transported temperature field (step (b) of
docs/ice_model_plan.md 4b).

Without this, v['T'] is a pure mixing tracer: it blends the sea and river boundary
values and nothing warms the channel interior, so mid-channel temperature is
unphysically cold in summer. This module supplies the missing atmospheric forcing.

    rho_w * cp_w * H * dT/dt  =  Q_sw + Q_lw + Q_sens + Q_lat        [W m^-2]

with, per grid cell,

    Q_sw   = (1 - albedo) * I_sw                                  absorbed shortwave
    Q_lw   = eps_w * sigma * (eps_a * Ta_K^4 - Tw_K^4)            net longwave
    Q_sens = rho_a * cp_a * C_sens * U * (Ta - Tw)                sensible
    Q_lat  = rho_a * L_v * C_lat * U * 0.622/P * (e_a - e_s(Tw))  latent

Clear-sky air emissivity follows Brutsaert (1975), eps_a = 1.24 (e_a/Ta_K)^(1/7),
and saturation vapour pressure the Magnus form.

WHAT TO BE SUSPICIOUS OF
------------------------
1. **Humidity is assumed, not observed.** PRDA2's DEWP field is empty for 2022, so
   e_a is built from config.REL_HUMIDITY (0.85). Latent heat is a leading term in an
   Arctic surface budget; this is the single largest source of uncertainty here.
2. **No cloud longwave.** Longwave uses clear-sky emissivity. Cloud raises incoming
   longwave, so net longwave loss is biased high. NOTE this is a WARMING correction
   if added, and the SHORTWAVE cloud effect is already present because the solar
   forcing is observed -- so adding clouds would NOT cure the model's warm bias; it
   would slightly worsen it. The warm bias is dominated by the inflow BOUNDARY
   temperature (the air-temp regression, +2.6 C warm on the Sagavanirktok), not by
   the surface budget. Verified by term decomposition; see CLAUDE.md.
3. **Shortwave is a daily mean** (solarradiation.csv), so there is no diurnal cycle.
   Tolerable at these latitudes in midsummer, poor in the shoulder seasons.
4. **Temperature is clamped at the freezing point.** Heat loss below T_FREEZE is booked
   to `ice_energy_deficit` (SET per step, not accumulated). With ICE_MODEL on, step (c)
   -- now implemented in ice_module.py -- consumes it that same step to grow ice via
   rho_ice * L_fusion; with ICE_MODEL off nothing consumes it and it is an unused
   per-step diagnostic while T is simply held at T_FREEZE.
5. **Ice-covered cells are skipped, not the whole budget.** With ICE_MODEL on the budget
   runs YEAR-ROUND but skips cells the ice model reports as covered (the slab insulates
   them; ice_module holds them at freezing). With ICE_MODEL off it falls back to the
   legacy previousdays gate.
"""

import numpy as np
from numba import njit

from config import (M, DELTI, rho_w, cp_water, albedo_water, emiss_water, sigma_SB,
                    rho_air, cp_air, L_vap, C_sens, C_lat, P_atm, REL_HUMIDITY,
                    T_FREEZE, H_MIN_HEAT, HEAT_BUDGET, ICE_MODEL, ICE_FORM_THRESH)
from variables import v, DEPTH, ice_thickness

# Heat removed by the freezing-point clamp THIS STEP [J m^-2], per cell. ice_module
# reads it each step and converts it into new ice via rho_ice * L_fusion, so it is
# SET (overwritten) each call, not accumulated. With ICE_MODEL off nothing consumes
# it and it is simply a per-step diagnostic of the freeze-clamp energy.
ice_energy_deficit = np.zeros(M + 1, dtype=np.float64)

# Diagnostic: net surface flux actually applied [W m^-2], for output/inspection.
Qnet = np.zeros(M + 1, dtype=np.float64)


@njit(cache=True)
def _esat(T):
    """Saturation vapour pressure over water [hPa], Magnus form. T in degC."""
    return 6.112 * np.exp(17.67 * T / (T + 243.5))


@njit(error_model='numpy')
def _heat_loop(T, DEPTH, Qnet, deficit, ice_thickness, M,
               T_air, U_wind, I_sw, rel_hum, dt, ice_on):
    """Jitted surface heat budget: apply Q_sw+Q_lw+Q_sens+Q_lat to each OPEN-WATER cell's
    temperature in place; skip ice-covered cells (insulated); and book the per-step
    freeze-clamp energy in `deficit` (J/m^2) for ice_module to turn into ice."""
    ea = rel_hum * _esat(T_air)               # air vapour pressure [hPa]
    Ta_K = T_air + 273.15
    # Brutsaert clear-sky emissivity; guard the fractional power against ea <= 0
    eps_a = 1.24 * (max(ea, 1e-6) / Ta_K) ** (1.0 / 7.0)
    if eps_a > 1.0:
        eps_a = 1.0

    for i in range(1, M + 1):
        deficit[i] = 0.0

        # Insulation: an ice-covered cell is shielded from the atmosphere. ice_module
        # conducts heat out through the slab and holds the water at T_FREEZE; the
        # surface budget must not touch it, or it would double-count the exchange.
        if ice_on and ice_thickness[i] > ICE_FORM_THRESH:
            Qnet[i] = 0.0
            continue

        Tw = T[i]
        H = DEPTH[i] if DEPTH[i] > H_MIN_HEAT else H_MIN_HEAT
        Tw_K = Tw + 273.15

        Q_sw = (1.0 - albedo_water) * I_sw
        Q_lw = emiss_water * sigma_SB * (eps_a * Ta_K**4 - Tw_K**4)
        Q_sens = rho_air * cp_air * C_sens * U_wind * (T_air - Tw)
        Q_lat = (rho_air * L_vap * C_lat * U_wind
                 * 0.622 / P_atm * (ea - _esat(Tw)))

        Q = Q_sw + Q_lw + Q_sens + Q_lat
        Qnet[i] = Q

        dT = Q * dt / (rho_w * cp_water * H)
        Tnew = Tw + dT

        if Tnew < T_FREEZE:
            # Energy that would have cooled the water below freezing. ice_module turns
            # this into ice (latent heat of fusion); booked per step, not accumulated.
            deficit[i] = (T_FREEZE - Tnew) * rho_w * cp_water * H
            Tnew = T_FREEZE
        T[i] = Tnew


def heat_budget(t, T_air, U_wind, I_sw, previousdays, rel_hum=REL_HUMIDITY):
    """
    Apply the surface heat budget to v['T'] in place.

    `rel_hum` is the relative humidity fraction. main.py now passes the OBSERVED
    daily value from Deadhorse Airport (relhum_2022_frac.csv); the config constant
    REL_HUMIDITY remains only as the fallback default.

    With ICE_MODEL on, the budget runs YEAR-ROUND: it must keep cooling open water to
    the freezing point through autumn so ice_module has freeze energy to grow a cover,
    and it skips ice-covered cells internally (the slab insulates them). With ICE_MODEL
    off it falls back to the crude previousdays gate -- under ice the atmosphere is
    assumed not to touch the water at all.
    """
    if not HEAT_BUDGET:
        return
    if not ICE_MODEL and previousdays <= 0:
        return
    _heat_loop(v['T']['c'], DEPTH, Qnet, ice_energy_deficit, ice_thickness, M,
               float(T_air), float(U_wind), float(I_sw), float(rel_hum), float(DELTI),
               ICE_MODEL)
