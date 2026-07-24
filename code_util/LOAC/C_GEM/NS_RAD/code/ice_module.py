"""
Prognostic river-ice model (step (c) of docs/ice_model_plan.md).

Replaces the crude `previousdays` on/off gate with a physical ice cover that grows,
insulates, and melts. One prognostic per grid cell: `ice_thickness` [m]. The areal
cover `ice_frac` is diagnosed from it (0 or 1 here -- a cell is either open or capped;
sub-grid leads are not resolved on a 200 m 1-D channel).

ENERGETICS
----------
Two growth/melt mechanisms and one mechanical one:

1. **Frazil / initial freeze-up.** heat_module cools open water with the full surface
   budget. When a cell would drop below T_FREEZE it is clamped and the removed heat is
   booked, per step, in heat_module.ice_energy_deficit [J m^-2]. Here that energy is
   turned into new ice via the latent heat of fusion:  dh = E / (rho_ice * L_fusion).

2. **Conductive (Stefan) growth under an existing cover.** Once a slab exists the water
   below is held at the freezing point and the atmosphere no longer touches it directly
   (heat_module skips ice-covered cells). Ice thickens by conducting the base-to-surface
   heat flux out through the slab:

       Q_cond = k_ice * (T_freeze - T_surf) / h ,   T_surf = min(T_air, T_freeze)

   dh = Q_cond * dt / (rho_ice * L_fusion). This is the mechanism that builds the ~1.5-2 m
   of winter ice typical of North Slope rivers, which the frazil term alone cannot.

3. **Surface melt.** When air is above freezing the slab loses mass from the top to
   sensible heat and absorbed shortwave (ice albedo), thinning the cover through spring.

4. **Hydraulic break-up.** North Slope rivers do not simply melt out -- the spring
   freshet surge mechanically shatters and flushes the ice. When discharge exceeds
   BREAKUP_Q_FACTOR x winter baseflow AND air is above freezing, the cover is cleared.

BOTTOM-FAST ICE
---------------
Shallow reaches freeze to the bed. The local depth limits how much NEW ice can form -- once
the slab reaches the bed there is no water left to freeze -- but ice already present is
never removed by a falling water level: grounded ice rests on the bed at unchanged
thickness. Where the bed is reached the channel is closed (ice_frac stays 1 and downstream
couplings treat the cell as sealed). This is the North-Slope-specific behaviour that
motivated a real ice model.

The distinction matters because DEPTH is the live, tide- and surge-varying depth: capping
h to it outright destroyed ice on every low-water excursion and thinned the cover through
midwinter (see (5) in _ice_loop). Note that a grounded slab may now be thicker than the
instantaneous water depth; nothing downstream divides by (DEPTH - h), and the ice draft
still does not reduce the hydrodynamic cross-section (docs/ice_model_plan.md Tier 3).

COUPLINGS (applied elsewhere, driven by ice_frac / ice_thickness set here)
- heat_module: skips ice-covered cells (insulation); holds under-ice water at T_FREEZE.
- biogeo_module: gas exchange (piston velocity) scaled by (1 - ice_frac) -> 0 under ice;
  under-ice PAR attenuated by the slab (k_ice_PAR).
"""

import numpy as np
from numba import njit

from config import (M, DELTI, rho_ice, L_fusion, k_ice, albedo_ice, H_MIN_ICE,
                    ICE_FORM_THRESH, BREAKUP_Q_FACTOR, T_FREEZE, ICE_MODEL,
                    rho_air, cp_air, C_sens)
from variables import ice_thickness, ice_frac, v, DEPTH
from heat_module import ice_energy_deficit


@njit(cache=True)
def _ice_loop(ice_thickness, ice_frac, T, DEPTH, deficit, M,
              T_air, U_wind, I_sw, Qr, q_ref, dt):
    """Jitted per-cell ice update: grow ice from the heat budget's freeze energy, thicken
    by Stefan conduction, melt at the surface in spring, clear it on the freshet surge,
    and cap at the local depth (bottom-fast). Holds under-ice water at T_FREEZE and
    diagnoses ice_frac. Mutates ice_thickness / ice_frac / T in place."""
    latent = rho_ice * L_fusion                     # J m^-3 to freeze/melt unit thickness
    T_surf = T_air if T_air < T_FREEZE else T_FREEZE
    # Hydraulic break-up: the freshet surge mechanically shatters and flushes the cover.
    # Triggered on discharge relative to the ANNUAL MEAN (q_ref) rather than winter
    # baseflow, because the ice-affected winter percentile is ~0 and would make the
    # threshold meaningless. North Slope freshets run many times the annual mean, so
    # this fires only during the freshet -- no temperature gate is needed (and one made
    # breakup revert to slow thermal melt when a post-freshet cold snap refroze the ice).
    surge = (q_ref > 0.0) and (abs(Qr) > BREAKUP_Q_FACTOR * q_ref)

    for i in range(1, M + 1):
        h = ice_thickness[i]
        h_prev = h                                  # thickness at the start of the step,
                                                    # for the grounding rule at (5) below

        if h <= ICE_FORM_THRESH:
            # (1) freeze-up: convert this step's booked freeze energy into new ice
            E = deficit[i]
            if E > 0.0:
                h += E / latent
        else:
            # (2) conductive growth: base at T_FREEZE, surface at min(T_air, T_FREEZE)
            hc = h if h > H_MIN_ICE else H_MIN_ICE
            Q_cond = k_ice * (T_FREEZE - T_surf) / hc          # >0 when air below freezing
            h += Q_cond * dt / latent
            # (3) surface melt when the air is warm
            if T_air > T_FREEZE:
                Q_melt = (rho_air * cp_air * C_sens * U_wind * (T_air - T_FREEZE)
                          + (1.0 - albedo_ice) * I_sw)
                if Q_melt > 0.0:
                    h -= Q_melt * dt / latent
            # under-ice water is held at the freezing point (atmosphere is insulated away)
            T[i] = T_FREEZE

        # (4) hydraulic break-up on the freshet surge
        if surge:
            h = 0.0

        if h < 0.0:
            h = 0.0
        # (5) bottom-fast: ice cannot GROW past the bed, but GROUNDING does not destroy it.
        #
        # DEPTH is the LIVE water depth, and on this coast it swings ~0.8 m over a day on
        # tide + storm surge -- comparable to the whole water column on the three shallow
        # rivers. Truncating h to DEPTH on every low-water excursion (the previous
        # behaviour) was therefore a ONE-WAY RATCHET: the clipped ice was destroyed, and
        # when the water rose again it could only return by slow conduction. The cover
        # ratcheted THINNER through midwinter instead of growing -- domain-mean thickness
        # on the Sagavanirktok fell 0.81 -> 0.55 m over days 10-110, reversing direction on
        # 35 of 99 days, and 98% of all midwinter thinning events sat exactly at this cap.
        #
        # Physically, bottom-fast ice GROUNDS: when the water drains from under it, it
        # rests on the bed at unchanged thickness. So the depth limits how much NEW ice can
        # form (there is no water left to freeze) but never removes ice already there.
        # Spring surface melt at (3) and break-up at (4) are applied before this and still
        # thin the cover normally.
        if DEPTH[i] > 0.0 and h > DEPTH[i]:
            h = DEPTH[i] if DEPTH[i] > h_prev else h_prev

        ice_thickness[i] = h
        ice_frac[i] = 1.0 if h > ICE_FORM_THRESH else 0.0


def ice_step(t, T_air, U_wind, I_sw, Qr, q_ref):
    """
    Advance the ice cover one timestep, in place on variables.ice_thickness / ice_frac.

    Call AFTER heat_budget (which books the per-step freeze energy in
    heat_module.ice_energy_deficit) and BEFORE biogeo (which reads ice_frac to gate
    gas exchange and under-ice light). `q_ref` is the annual-mean discharge [m^3 s^-1]
    for the break-up trigger, computed once in main() from the discharge series.
    """
    if not ICE_MODEL:
        return
    _ice_loop(ice_thickness, ice_frac, v['T']['c'], DEPTH, ice_energy_deficit, M,
              float(T_air), float(U_wind), float(I_sw), float(Qr), float(q_ref),
              float(DELTI))
