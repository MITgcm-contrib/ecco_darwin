"""
Biogeochemical reaction network module (translated from biogeo.c)

The per-grid-point loop is @njit-compiled. It was the largest remaining hotspot
(~25% of runtime) once the hydrodynamic kernels and the density stack were done.
The kernel is a verbatim transcription of the original Python loop -- same
expressions, same order, same branches -- so results are bit-identical, verified by
the output-hash check described in docs/performance.md.

Two things about the compiled form:

* `@njit` WITHOUT `cache=True`. The kernel reads its ~28 rate constants (Pbmax,
  kox, KTOC, ...) as module globals, and numba freezes globals as compile-time
  constants. Numba's on-disk cache keys on the function's bytecode, NOT on the
  values it captured -- so with caching enabled, editing config.py would silently
  keep using the previously compiled constants. Recompiling each process costs a
  second or two against a ~20-minute run, which is the right trade. If you ever add
  `cache=True` here, pass the constants as arguments instead.
* The scalar helpers in fun_module (`_fhet`, `_k1`, `_o2sat`, ...) are jitted
  variants of the Python functions of the same purpose; the originals are untouched
  and still used elsewhere.

One redundancy was removed in the process, arithmetically identical because the
inputs are unchanged within an iteration: `p_bar` is computed once per cell and
shared by K1 and K2 rather than each recomputing it, halving the remaining
`density.dens` calls. (`Fhet`/`Fnit` were also hoisted out of the loop while
temperature was a single scalar; they moved back inside once temperature became a
transported per-cell field.)

TEMPERATURE IS A TRANSPORTED FIELD. The kernel takes `c_T` (variables.v['T']['c'])
and evaluates every temperature-dependent rate at the LOCAL value `Ti = c_T[i]`.
It used to take a single `water_temp` scalar applied uniformly to all M grid points.
"""

import math

import numpy as np
from numba import njit

from variables import v, DEPTH, M, vp, GPP, NPP_NO3, NPP_NH4, NPP, phy_death, Si_consumption, aer_deg, denit, nit,\
    O2_ex, NEM, FCO2, Hplus, ice_frac, ice_thickness,\
    rdoc_ox, photo, ch4_ox, ch4_ex, n2o_prod, n2o_ex, sod
from config import (
    KD1, KD2, Pbmax, alpha, kexcr, kgrowth, kmaint, kmort, redsi, redn,
    redp, KdSi, KN, KPO4, Euler, kox, KTOC, KO2_ox, KO2_nit, KinO2, KNO3, knit, KNH4, kdenit, DELTI, TS, Chla2CMIN, KE,
    ICE_MODEL, k_ice_PAR,
    # Arctic biogeochemistry extension
    ARCTIC_BGC,
    krefr, K_RDOC, PHOTO_EFF, f_photo_lab, k_ch4ox, K_ch4, K_ch4O2, pCH4_atm,
    y_nit, y_denit, pN2O_atm, k_sod, K_sod, k_bdenit, k_methano,
)
from fun_module import (piston_velocity,
                        _fhet, _fnit, _o2sat, _pbar_rho, _k0, _k1, _k2, _kb, _h_solve_kg,
                        _sc, _sc_ch4, _sc_n2o, _ch4_eq, _n2o_eq)
from file_module import Rates


@njit(error_model='numpy')  # deliberately not cache=True -- see module docstring
# error_model='numpy' is REQUIRED, not cosmetic. Numba's default ('python') raises
# ZeroDivisionError on float division by zero; numpy returns +/-inf. The species
# arrays are numpy, so numpy semantics are what the uncompiled model actually does
# -- and this loop does divide by zero in the ice-gated branch, where the state has
# been zeroed and inf propagates harmlessly into terms that are then discarded.
# With the default model the kernel raises where the interpreter does not.
def _biogeo_loop(c_DIA, c_dSi, c_NO3, c_NH4, c_PO4, c_O2, c_TOC, c_S, c_SPM,
                 c_DIC, c_ALK, c_pH, c_T, c_RDOC, c_CH4, c_N2O, DEPTH, vp,
                 GPP, NPP_NO3, NPP_NH4, NPP, phy_death, Si_consumption,
                 aer_deg, denit, nit, O2_ex, NEM, FCO2, Hplus,
                 rdoc_ox, photo, ch4_ox, ch4_ex, n2o_prod, n2o_ex, sod,
                 ice_frac, ice_thickness, k_ice_PAR, ice_on, arctic_on,
                 M, t, pCO2, I0, previousdays):
    """Jitted per-cell biogeochemistry: light-limited primary production, aerobic
    degradation, denitrification, nitrification, O2 and CO2 gas exchange, the mol/kg
    carbonate solve (pH, CO2*, FCO2), and the species state update -- all evaluated at
    the LOCAL temperature c_T[i] and coupled to the ice cover (see the ICE COUPLING block
    below). Mutates the concentration and rate arrays in place; nothing is returned."""
    for i in range(1, M + 1):
        # ICE COUPLING. With the prognostic ice model on, the crude previousdays gate
        # is retired in favour of the real ice cover:
        #   * incident PAR is attenuated through the slab (Beer-Lambert, k_ice_PAR),
        #     so under thick ice the light integral collapses and GPP -> 0 on its own
        #     while thin spring ice still admits an under-ice bloom;
        #   * gas exchange (O2 and CO2) is scaled by the open-water fraction, so a
        #     sealed cell neither outgasses nor ventilates;
        #   * state is CONSERVED and keeps reacting (respiration builds under-ice DIC)
        #     rather than being zeroed, which also removes the near-zero DIC/ALK that
        #     used to blow the carbonate solve up at ice-out.
        # With ice_on False the block reproduces the legacy previousdays behaviour.
        if ice_on:
            I_use = I0 * math.exp(-k_ice_PAR * ice_thickness[i])
            openfac = 1.0 - ice_frac[i]
            react = True
        else:
            I_use = I0
            openfac = 1.0
            react = previousdays > 0
        # Temperature is now per-cell, so the Arrhenius factors vary along the
        # channel and can no longer be hoisted out of this loop (they were, while
        # water_temp was a single scalar applied to all M points).
        Ti = c_T[i]
        fhet = _fhet(Ti)
        fnit = _fnit(Ti)
        # Underwater light field and nutrient dependence. I_use is the incident PAR
        # reaching the water surface -- equal to I0 in open water, attenuated by the
        # ice slab under a cover (see the ICE COUPLING block above).
        KD = KD1 + KD2 * (1000.0 * c_SPM[i] + 100.0)
        Ebottom = max(1e-3, I_use * math.exp(-KD * DEPTH[i]))
        # chla to carbon ration g Chla(g C)-1
        Chla2C = max(Chla2CMIN, Chla2CMIN * (1 + 4 * math.exp((-1 / 2) * ((I_use * 1e6 / 86400) / KE)) *
                                ((c_NO3[i] + c_NH4[i]) / (c_NH4[i] + c_NO3[i] + KN))))
        # Photoacclimation parameter muE m^-2 s^-1
        Ek = (1 / Chla2C) * Pbmax / alpha
        psurf = (Pbmax * (1 - math.exp(-((I_use * 1e6 / 86400) / Ek))) *
                 ((c_NO3[i] + c_NH4[i]) / (c_NH4[i] + c_NO3[i] + KN)))  # s-1
        pbot = (Pbmax * (1 - math.exp(-((Ebottom * 1e6 / 86400) / Ek))) *
                ((c_NO3[i] + c_NH4[i]) / (c_NH4[i] + c_NO3[i] + KN)))  # s-1

        if (I_use > 0 and Ebottom > 0) and react:
            # Depth-integrated light-limited production. With Beer-Lambert light
            # I(z) = I_use * exp(-KD*z) and a Platt-style P-I curve P = Pbmax(1 - exp(-I/Ek)),
            # the analytic integral of P over the water column reduces to a difference of
            # exponential-integral (incomplete-gamma) terms evaluated at the surface (psurf)
            # and bottom (pbot) light levels. appGAMMAsurf/bot approximate that special
            # function in two regimes: p <= 1 uses a small-argument log-series, p > 1 a
            # large-argument continued fraction (standard E1 evaluation). `integral` (per
            # unit KD) is the resulting depth-integrated growth, fed into GPP below.
            # Gamma approximation for surface
            if psurf <= 1.0 and psurf > 0:
                appGAMMAsurf = -(
                    math.log(psurf) + Euler - psurf
                    + psurf**2 / 4.0 - psurf**3 / 18.0
                    + psurf**4 / 96.0 - psurf**5 / 600.0
                )
            else:
                appGAMMAsurf = math.exp(-psurf) * (
                    1.0 / (psurf + 1.0 - (1.0 / (psurf + 3.0 - (4.0 / (psurf + 5.0 - (9.0 / (psurf + 7.0 - (16.0 / (psurf + 9.0)))))))))
                )

            # Gamma approximation for bottom
            if pbot <= 1.0 and pbot > 0:
                appGAMMAbot = -(
                    math.log(pbot) + Euler + pbot
                    - pbot**2 / 4.0 + pbot**3 / 18.0
                    - pbot**4 / 96.0 + pbot**5 / 600.0
                )
            else:
                appGAMMAbot = math.exp(-pbot) * (
                    1.0 / (pbot + 1.0 - (1.0 / (pbot + 3.0 - (4.0 / (pbot + 5.0 - (9.0 / (pbot + 7.0 - (16.0 / (pbot + 9.0)))))))))
                )
            integral = (appGAMMAsurf - appGAMMAbot + math.log(I_use / Ebottom)) / KD
        else:
            integral = 0.0

        if c_dSi[i] <= 0.0:
            nlim = 0.0
        else:
            nlim = (c_dSi[i] / (c_dSi[i] + KdSi)
                    * (c_NO3[i] + c_NH4[i]) / (c_NH4[i] + c_NO3[i] + KN)
                    * (c_PO4[i] / (c_PO4[i] + KPO4)))
        fNH4 = 0.0 if c_NH4[i] <= 0.0 else c_NH4[i] / (10.0 + c_NH4[i])

        # Biogeochemical reaction rates [mmol m^-3 s^-1]
        GPP[i] = Pbmax * c_DIA[i] * nlim * integral
        NPP_NO3[i] = (1.0 - fNH4) * ((GPP[i] / DEPTH[i]) * (1.0 - kexcr) * (1.0 - kgrowth) - kmaint * c_DIA[i])
        NPP_NH4[i] = fNH4 * ((GPP[i] / DEPTH[i]) * (1.0 - kexcr) * (1.0 - kgrowth) - kmaint * c_DIA[i])
        NPP[i] = NPP_NO3[i] + NPP_NH4[i]
        phy_death[i] = kmort * c_DIA[i]
        Si_consumption[i] = max(0.0, redsi * NPP[i])
        aer_deg[i] = kox * fhet * c_TOC[i] / (c_TOC[i] + KTOC) * c_O2[i] / (c_O2[i] + KO2_ox)
        denit[i] = kdenit * fhet * c_TOC[i] / (c_TOC[i] + KTOC) * KinO2 / (c_O2[i] + KinO2) * c_NO3[i] / (c_NO3[i] + KNO3)
        nit[i] = knit * fnit * c_O2[i] / (c_O2[i] + KO2_nit) * c_NH4[i] / (c_NH4[i] + KNH4)
        # Air-water O2 exchange, shut off in proportion to ice cover.
        O2_ex[i] = openfac * (vp[i] / DEPTH[i]) * (_o2sat(Ti, c_S[i]) - c_O2[i])
        NEM[i] = NPP[i] - aer_deg[i] - denit[i]

        # Carbonate system, solved consistently in mol/kg (ported from C-GEM v2). The
        # dissociation constants K0/K1/K2/KB are on a mol/kg-SW basis, so the state must
        # be converted from mmol/m^3 to mol/kg via the local in-situ density before the
        # speciation -- the legacy path mixed the two and only patched the borate term.
        # density.dens is the expensive call; _pbar_rho returns pressure AND density from
        # one evaluation, and K1/K2 share that pressure.
        pb, rho_i = _pbar_rho(c_S[i], Ti, DEPTH[i])       # bar, kg/m^3
        k1 = _k1(Ti, c_S[i], pb)
        k2 = _k2(Ti, c_S[i], pb)
        dic_kg = c_DIC[i] * 1.0e-3 / rho_i                # mmol/m^3 -> mol/kg
        alk_kg = c_ALK[i] * 1.0e-3 / rho_i

        # H+ [mol/kg] from the mol/kg state (Follows et al., 2006). Guess = 10**(-pH)
        # from the previous step; the solver keeps the pH-fix guards (iterate, neutral
        # fallback, physical bound).
        Hplus[i] = _h_solve_kg(c_S[i], 10 ** (-c_pH[i]), dic_kg, alk_kg,
                               _kb(Ti, c_S[i]), k1, k2)
        # CO2* [mol/kg] from the speciation, then back to mmol/m^3 with the same density.
        co2s_kg = dic_kg / (1.0 + (k1 / Hplus[i]) + (k1 * k2 / (Hplus[i] * Hplus[i])))
        co2s = co2s_kg * rho_i * 1.0e3                    # mmol/m^3
        k0_v = _k0(Ti, c_S[i]) * rho_i * 1.0e3            # Henry, mmol m^-3 atm^-1
        # CO2 flux. SIGN CONVENTION: POSITIVE = outgassing (sea -> air), i.e. RCO2 > 0
        # when the water is CO2-supersaturated (co2s > K0*pCO2). Scaled by open-water
        # fraction (consistent with O2_ex): a sealed cell neither outgasses nor ventilates.
        vpCO2 = 0.915 * vp[i]  # O2 -> CO2 piston velocity (Regnier et al., 2002) [m/s]
        RCO2 = vpCO2 * (co2s - k0_v * pCO2)  # mmol C m^-2 s^-1, >0 outgassing
        FCO2[i] = openfac * RCO2 / DEPTH[i]  # mmol C m^-3 s^-1, >0 outgassing

        # ============================================================================
        # ARCTIC BIOGEOCHEMISTRY EXTENSION (docs/arctic_biogeochemistry.md). All rates
        # [mmol m^-3 s^-1]. In-water reactions (rdoc_ox, ch4_ox, benthic) run under ice
        # too -- only the GAS-EXCHANGE terms are openfac-scaled, exactly like O2_ex/FCO2,
        # so DIC/CH4 build up under the cover and vent at ice-out. Photolysis self-gates
        # on the (already ice-attenuated) light field.
        # ----------------------------------------------------------------------------
        # OPT-IN (config.ARCTIC_BGC). When off, every new rate is zero and the state
        # updates below reduce EXACTLY to the shipped network -> real-river fields are
        # bit-identical. The three new tracers still advect (as inert passive tracers).
        if arctic_on:
            # (1) refractory DOC: slow aerobic oxidation + CDOM photomineralisation -> DIC
            rdoc_ox[i] = krefr * fhet * c_RDOC[i] / (c_RDOC[i] + K_RDOC) * c_O2[i] / (c_O2[i] + KO2_ox)
            Iabs = I_use - Ebottom                      # PAR absorbed in the column [muE m^-2 s^-1]
            if Iabs < 0.0:
                Iabs = 0.0
            photo[i] = PHOTO_EFF * (c_RDOC[i] / (c_RDOC[i] + K_RDOC)) * Iabs / DEPTH[i]
            # (2a) methane: methanotrophy (CH4 + 2 O2 -> CO2) + air-water exchange
            ch4_ox[i] = k_ch4ox * fhet * c_CH4[i] / (c_CH4[i] + K_ch4) * c_O2[i] / (c_O2[i] + K_ch4O2)
            vp_ch4 = vp[i] * math.sqrt(_sc(Ti, c_S[i]) / _sc_ch4(Ti))   # Schmidt-rescaled piston vel
            ch4_ex[i] = openfac * (vp_ch4 / DEPTH[i]) * (c_CH4[i] - _ch4_eq(Ti, c_S[i], pCH4_atm))
            # (2b) nitrous oxide: yield from nitrification + denitrification, + air-water exchange
            n_denit = (94.4 / 106.0) * denit[i]         # N reduced by denitrification [mmol N m^-3 s^-1]
            n2o_prod[i] = 0.5 * (y_nit * nit[i] + y_denit * n_denit)   # 2 N atoms per N2O molecule
            vp_n2o = vp[i] * math.sqrt(_sc(Ti, c_S[i]) / _sc_n2o(Ti))
            n2o_ex[i] = openfac * (vp_n2o / DEPTH[i]) * (c_N2O[i] - _n2o_eq(Ti, c_S[i], pN2O_atm))
            # (3) benthic efflux (per unit volume = per-area flux / depth): sediment O2 demand
            #     -> DIC; benthic denitrification -> alkalinity + NO3 loss; methanogenesis -> CH4
            sod[i]  = k_sod * fhet * c_O2[i] / (c_O2[i] + K_sod) / DEPTH[i]
            b_denit = k_bdenit * fhet * c_NO3[i] / (c_NO3[i] + KNO3) / DEPTH[i]
            b_methano = k_methano * fhet * (1.0 - c_O2[i] / (c_O2[i] + K_sod)) / DEPTH[i]
        else:
            rdoc_ox[i] = 0.0; photo[i] = 0.0; ch4_ox[i] = 0.0; ch4_ex[i] = 0.0
            n2o_prod[i] = 0.0; n2o_ex[i] = 0.0; sod[i] = 0.0
            b_denit = 0.0; b_methano = 0.0
        # total heterotrophic + abiotic respiration for the metabolism diagnostic
        # (reduces to NPP - aer_deg - denit when the extension is off)
        NEM[i] = NPP[i] - aer_deg[i] - denit[i] - rdoc_ox[i] - ch4_ox[i] - sod[i]

        if react:
            # Update biogeochemical state variables [mmol m^-3]
            c_DIA[i] = c_DIA[i] + (NPP[i] - phy_death[i]) * DELTI
            c_dSi[i] = c_dSi[i] - Si_consumption[i] * DELTI
            c_NO3[i] = c_NO3[i] + (-94.4 / 106.0 * denit[i] + nit[i] - redn * NPP_NO3[i] - b_denit) * DELTI
            c_NH4[i] = c_NH4[i] + (redn * (aer_deg[i] - NPP_NH4[i]) - nit[i]) * DELTI
            c_PO4[i] = c_PO4[i] + redp * (aer_deg[i] + denit[i] - NPP[i]) * DELTI
            # O2 sinks now also: refractory-DOC oxidation, methanotrophy (2 O2 per CH4), and SOD
            c_O2[i] = c_O2[i] + (-aer_deg[i] + NPP_NH4[i] + (138.0 / 106.0) * NPP_NO3[i] - 2.0 * nit[i]
                                 + O2_ex[i] - rdoc_ox[i] - 2.0 * ch4_ox[i] - sod[i]) * DELTI
            # labile DOC gains the labile fraction of the photoproduct
            c_TOC[i] = c_TOC[i] + (-aer_deg[i] - denit[i] + phy_death[i] + f_photo_lab * photo[i]) * DELTI
            # refractory/chromophoric DOC: consumed by slow oxidation and photolysis
            c_RDOC[i] = c_RDOC[i] + (-rdoc_ox[i] - photo[i]) * DELTI
            # methane: methanotrophy + outgassing sinks, benthic methanogenesis source
            c_CH4[i] = c_CH4[i] + (-ch4_ox[i] - ch4_ex[i] + b_methano) * DELTI
            # nitrous oxide: production from N cycling, outgassing sink
            c_N2O[i] = c_N2O[i] + (n2o_prod[i] - n2o_ex[i]) * DELTI
            # -FCO2 because outgassing (positive) REMOVES DIC. DIC now also gains refractory
            # oxidation, the direct (non-labile) photoproduct, methane oxidation, and benthic
            # respiration (SOD) -- the Arctic land-to-ocean CO2 terms.
            c_DIC[i] = c_DIC[i] + (-FCO2[i] - NPP_NO3[i] - NPP_NH4[i] + aer_deg[i] + denit[i]
                                   + rdoc_ox[i] + (1.0 - f_photo_lab) * photo[i]
                                   + ch4_ox[i] + sod[i]) * DELTI
            # ALK: benthic denitrification raises alkalinity (~1 eq per mol NO3 reduced)
            c_ALK[i] = c_ALK[i] + (((15 / 106) * aer_deg[i]) + ((93.4 / 106) * denit[i]) - (2 * nit[i])
                                   - ((15 / 106) * NPP_NH4[i]) + ((17 / 106) * NPP_NO3[i])
                                   + b_denit) * DELTI

            # pH from the mol/kg H+. `_h_solve_kg` always returns a POSITIVE H+ in
            # [1e-12, 1e-2] (or the neutral 1e-7 for a degenerate carbonate state), so
            # the guard below is belt-and-braces: the negative/zero-H+ pathologies of the
            # old mmol/m^3 solve -- which produced math-domain errors at ice-out -- can no
            # longer occur (state is also conserved, not zeroed, under the ice model).
            if Hplus[i] > 0.0:
                c_pH[i] = -math.log10(Hplus[i])   # H+ in mol/kg -> pH
        else:
            GPP[i] = 0.0
            NPP_NO3[i] = 0.0
            NPP_NH4[i] = 0.0
            NPP[i] = 0.0
            phy_death[i] = 0.0
            Si_consumption[i] = 0.0
            aer_deg[i] = 0.0
            denit[i] = 0.0
            nit[i] = 0.0
            O2_ex[i] = 0.0
            NEM[i] = 0.0
            # Arctic-extension rates zeroed consistently with the legacy gate.
            rdoc_ox[i] = 0.0
            photo[i] = 0.0
            ch4_ox[i] = 0.0
            ch4_ex[i] = 0.0
            n2o_prod[i] = 0.0
            n2o_ex[i] = 0.0
            sod[i] = 0.0
            # TIER 0a: gate the CO2 flux consistently with O2_ex.
            # FCO2 is computed above this branch, so under ice it was evaluated from
            # DIC/ALK that this same branch had zeroed -- yielding first a spurious
            # wrong-signed uptake, then inf/NaN once Hplus collapsed. In the four-site
            # production run that put NaN into ~5% of FCO2.dat rows (918 of 18007 for
            # Colville), every one of them inside the ice-gated window. Oxygen
            # exchange was already gated here; carbon was not. That asymmetry was
            # unintended.
            FCO2[i] = 0.0

            # Update biogeochemical state variables [mmol m^-3]
            c_DIA[i] = 0.0
            c_dSi[i] = 0.0
            c_NO3[i] = 0.0
            c_NH4[i] = 0.0
            c_PO4[i] = 0.0
            c_O2[i] = 0.0
            c_TOC[i] = 0.0
            c_RDOC[i] = 0.0
            c_CH4[i] = 0.0
            c_N2O[i] = 0.0
            c_DIC[i] = 0.0
            c_ALK[i] = 0.0
            c_pH[i] = 0.0


def biogeo(t, Uw_sal, Uw_tid, T, pCO2, I0, previousdays):
    """Biogeochemical reaction network simulation."""
    piston_velocity(t, Uw_sal, Uw_tid, T)
    # convert solar radiation in W/m2 to muE m^-2 s^-1
    I0 = I0 * 4.57

    _biogeo_loop(v['DIA']['c'], v['dSi']['c'], v['NO3']['c'], v['NH4']['c'],
                 v['PO4']['c'], v['O2']['c'], v['TOC']['c'], v['S']['c'],
                 v['SPM']['c'], v['DIC']['c'], v['ALK']['c'], v['pH']['c'],
                 T, v['RDOC']['c'], v['CH4']['c'], v['N2O']['c'], DEPTH, vp,
                 GPP, NPP_NO3, NPP_NH4, NPP, phy_death, Si_consumption,
                 aer_deg, denit, nit, O2_ex, NEM, FCO2, Hplus,
                 rdoc_ox, photo, ch4_ox, ch4_ex, n2o_prod, n2o_ex, sod,
                 ice_frac, ice_thickness, k_ice_PAR, ICE_MODEL, ARCTIC_BGC,
                 M, t, pCO2, I0, previousdays)

    # Write biogeochemical process rates [mmol m^-3 s-1]
    if (float(t) / float(TS * DELTI)) % 1 == 0:
        Rates(NPP, "NPP.dat", t)
        Rates(aer_deg, "aer_deg.dat", t)
        Rates(denit, "denit.dat", t)
        Rates(nit, "nit.dat", t)
        Rates(O2_ex, "O2_ex.dat", t)
        Rates(NEM, "NEM.dat", t)
        Rates(phy_death, "phy_death.dat", t)
        Rates(FCO2, "FCO2.dat", t)
        # Arctic biogeochemistry extension rates
        Rates(rdoc_ox, "rdoc_ox.dat", t)
        Rates(photo, "photo.dat", t)
        Rates(ch4_ox, "ch4_ox.dat", t)
        Rates(ch4_ex, "ch4_ex.dat", t)
        Rates(n2o_prod, "n2o_prod.dat", t)
        Rates(n2o_ex, "n2o_ex.dat", t)
        Rates(sod, "sod.dat", t)
