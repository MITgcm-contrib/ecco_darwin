import numpy as np
import math
from variables import (
    v, DEPTH, M, vp, GPP, NPP_NO3, NPP_NH4, NPP, phy_death, Si_consumption, aer_deg,
    denit, nit, O2_ex, NEM, kflow, kwind, U
)
from config import (
    KD1, KD2, Pbmax, alpha, kexcr, kgrowth, kmaint, kmort, redsi, redn,
    redp, KdSi, KN, KPO4, Euler, kox, KTOC, KO2_ox, KO2_nit, KinO2, KNO3, knit, KNH4, kdenit, DELTI, TS,
    USE_CO2_FLUX
)
if USE_CO2_FLUX:
    from variables import FCO2

from fun_module import I0, Fhet, Fnit, O2sat, piston_velocity, co2_step_semi_implicit

def gamma_approx(p):
    """Accurate approximation of the incomplete gamma integral for primary production."""

    res = np.zeros_like(p)
    small = p <= 1.0
    large = ~small

    # Safe log
    with np.errstate(divide="ignore", invalid="ignore"):

        # --- For small p: Taylor expansion around 0
        ps = p[small]
        res[small] = -(
            np.log(ps + 1e-300) + Euler - ps
            + ps**2 / 4.0
            - ps**3 / 18.0
            + ps**4 / 96.0
            - ps**5 / 600.0
        )

        # --- For large p: Continued fraction approximation
        pl = p[large]
        denom = (
            pl + 1.0
            - (1.0 / (
                pl + 3.0
                - (4.0 / (
                    pl + 5.0
                    - (9.0 / (
                        pl + 7.0
                        - (16.0 / (pl + 9.0))
                    ))
                ))
            ))
        )
        res[large] = np.exp(-pl) / denom

    return res

def biogeo(t, io):
    """Biogeochemical reaction network simulation."""
    piston_velocity(t,v['S']['c'],U,DEPTH,kflow,kwind,vp)

    i_arr = np.arange(1, M + 1)

    # Light field
    I0_val = I0(t)
    Fhet_val = Fhet(t)
    Fnit_val = Fnit(t)

    KD = KD1 + KD2 * (1000.0 * v['SPM']['c'][1:M+1] + 100.0)
    Ebottom = I0_val * np.exp(-KD * DEPTH[1:M+1])
    psurf = I0_val * alpha / Pbmax
    pbot = Ebottom * alpha / Pbmax

    appGAMMAsurf = gamma_approx(np.full_like(pbot, psurf))
    appGAMMAbot = gamma_approx(pbot)
    with np.errstate(divide='ignore'):
        integral = (appGAMMAsurf - appGAMMAbot + np.log(I0_val / Ebottom)) / KD
    integral = np.where(Ebottom >= 1e-300, integral, 0.0)

    # Nutrient limitation
    dSi = v['dSi']['c'][1:M+1]
    NO3 = v['NO3']['c'][1:M+1]
    NH4 = v['NH4']['c'][1:M+1]
    PO4 = v['PO4']['c'][1:M+1]

    nlim = np.where(
        dSi <= 0,
        0.0,
        (dSi / (dSi + KdSi)) *
        ((NO3 + NH4) / (NO3 + NH4 + KN)) *
        (PO4 / (PO4 + KPO4))
    )

    fNH4 = np.where(NH4 <= 0.0, 0.0, NH4 / (10.0 + NH4))

    DIA = v['DIA']['c'][1:M+1]
    depth = DEPTH[1:M+1]

    # GPP & NPP
    GPP[1:M+1] = Pbmax * DIA * nlim * integral
    prod = (GPP[1:M+1] / depth) * (1.0 - kexcr) * (1.0 - kgrowth) - kmaint * DIA
    NPP_NO3[1:M+1] = (1.0 - fNH4) * prod
    NPP_NH4[1:M+1] = fNH4 * prod
    NPP[1:M+1] = NPP_NO3[1:M+1] + NPP_NH4[1:M+1]

    # Mortality & silicon use
    phy_death[1:M+1] = kmort * DIA
    Si_consumption[1:M+1] = np.maximum(0.0, redsi * NPP[1:M+1])

    TOC = v['TOC']['c'][1:M+1]
    O2 = v['O2']['c'][1:M+1]

    # Degradation
    aer_deg[1:M+1] = (
        kox * Fhet_val * TOC / (TOC + KTOC) *
        O2 / (O2 + KO2_ox)
    )
    denit[1:M+1] = (
        kdenit * Fhet_val * TOC / (TOC + KTOC) *
        KinO2 / (O2 + KinO2) * NO3 / (NO3 + KNO3)
    )
    nit[1:M+1] = (
        knit * Fnit_val * O2 / (O2 + KO2_nit) * NH4 / (NH4 + KNH4)
    )

    O2_eq = O2sat(t, v['S']['c'][1:M+1])
    O2_ex[1:M+1] = (vp[1:M+1] / depth) * (O2_eq - O2)

    NEM[1:M+1] = NPP[1:M+1] - aer_deg[1:M+1] - denit[1:M+1]

    # Update state variables
    v['DIA']['c'][1:M+1] += (NPP[1:M+1] - phy_death[1:M+1]) * DELTI
    v['dSi']['c'][1:M+1] -= Si_consumption[1:M+1] * DELTI
    v['NO3']['c'][1:M+1] += (-94.4 / 106.0 * denit[1:M+1] + nit[1:M+1] - redn * NPP_NO3[1:M+1]) * DELTI
    v['NH4']['c'][1:M+1] += (redn * (aer_deg[1:M+1] - NPP_NH4[1:M+1]) - nit[1:M+1]) * DELTI
    v['PO4']['c'][1:M+1] += redp * (aer_deg[1:M+1] + denit[1:M+1] - NPP[1:M+1]) * DELTI
    v['O2']['c'][1:M+1] += (-aer_deg[1:M+1] + NPP_NH4[1:M+1] + (138.0 / 106.0) * NPP_NO3[1:M+1] - 2.0 * nit[1:M+1] + O2_ex[1:M+1]) * DELTI
    v['TOC']['c'][1:M+1] += (-aer_deg[1:M+1] - denit[1:M+1] + phy_death[1:M+1]) * DELTI

    # ---------------- CO2 FLUX (optional) ----------------
    if USE_CO2_FLUX:

        co2_step_semi_implicit(
            t,
            v['S']['c'], DEPTH, vp,
            v['DIC']['c'], v['ALK']['c'], v['pH']['c'],
            NPP_NO3, NPP_NH4, aer_deg, denit, nit,
            FCO2
        )

        # optional output
        if t % (TS * DELTI) == 0:
            io.write_rates(FCO2, "FCO2.dat", t)
            io.write_tracer(v['pH']['c'], v['pH']['name'], t)

    # Optional output
    if t % (TS * DELTI) == 0:
        io.write_rates(NPP, "NPP.dat", t)
        io.write_rates(aer_deg, "aer_deg.dat", t)
        io.write_rates(denit, "denit.dat", t)
        io.write_rates(nit, "nit.dat", t)
        io.write_rates(O2_ex, "O2_ex.dat", t)
        io.write_rates(NEM, "NEM.dat", t)
        io.write_rates(phy_death, "phy_death.dat", t)
