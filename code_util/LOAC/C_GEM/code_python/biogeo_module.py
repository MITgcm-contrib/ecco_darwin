"""
Biogeochemical reaction network module (translated from biogeo.c)
"""

import math
from variables import v, DEPTH, M, vp, GPP, NPP_NO3, NPP_NH4, NPP, phy_death, Si_consumption, aer_deg, denit, nit, O2_ex, NEM
from config import (
    KD1, KD2, Pbmax, alpha, kexcr, kgrowth, kmaint, kmort, redsi, redn,
    redp, KdSi, KN, KPO4, Euler, kox, KTOC, KO2, KinO2, KNO3, knit, KNH4, kdenit, DELTI, TS
)
from fun_module import I0, Fhet, Fnit, O2sat, piston_velocity
from file_module import Rates

def biogeo(t):
    """Biogeochemical reaction network simulation."""
    piston_velocity(t)

    for i in range(1, M + 1):
        # Underwater light field and nutrient dependence
        KD = KD1 + KD2 * (1000.0 * v['SPM']['c'][i] + 100.0)
        Ebottom = I0(t) * math.exp(-KD * DEPTH[i])
        psurf = I0(t) * alpha / Pbmax
        pbot = Ebottom * alpha / Pbmax

        if I0(t) > 0 or Ebottom >= 1.e-300:
            # Gamma approximation for surface
            if psurf <= 1.0:
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
            if pbot <= 1.0:
                appGAMMAbot = -(
                    math.log(pbot) + Euler + pbot
                    - pbot**2 / 4.0 + pbot**3 / 18.0
                    - pbot**4 / 96.0 + pbot**5 / 600.0
                )
            else:
                appGAMMAbot = math.exp(-pbot) * (
                    1.0 / (pbot + 1.0 - (1.0 / (pbot + 3.0 - (4.0 / (pbot + 5.0 - (9.0 / (pbot + 7.0 - (16.0 / (pbot + 9.0)))))))))
                )
            integral = (appGAMMAsurf - appGAMMAbot + math.log(I0(t) / Ebottom)) / KD
        else:
            integral = 0.0

        nlim = (
            0.0 if v['dSi']['c'][i] <= 0.0
            else v['dSi']['c'][i] / (v['dSi']['c'][i] + KdSi)
            * (v['NO3']['c'][i] + v['NH4']['c'][i]) / (v['NH4']['c'][i] + v['NO3']['c'][i] + KN)
            * (v['PO4']['c'][i] / (v['PO4']['c'][i] + KPO4))
        )
        fNH4 = 0.0 if v['NH4']['c'][i] <= 0.0 else v['NH4']['c'][i] / (10.0 + v['NH4']['c'][i])

        # Biogeochemical reaction rates [mmol m^-3 s^-1]
        GPP[i] = Pbmax * v['DIA']['c'][i] * nlim * integral
        NPP_NO3[i] = (1.0 - fNH4) * ((GPP[i] / DEPTH[i]) * (1.0 - kexcr) * (1.0 - kgrowth) - kmaint * v['DIA']['c'][i])
        NPP_NH4[i] = fNH4 * ((GPP[i] / DEPTH[i]) * (1.0 - kexcr) * (1.0 - kgrowth) - kmaint * v['DIA']['c'][i])
        NPP[i] = NPP_NO3[i] + NPP_NH4[i]
        phy_death[i] = kmort * v['DIA']['c'][i]
        Si_consumption[i] = max(0.0, redsi * NPP[i])
        aer_deg[i] = kox * Fhet(t) * v['TOC']['c'][i] / (v['TOC']['c'][i] + KTOC) * v['O2']['c'][i] / (v['O2']['c'][i] + KO2)
        denit[i] = kdenit * Fhet(t) * v['TOC']['c'][i] / (v['TOC']['c'][i] + KTOC) * KinO2 / (v['O2']['c'][i] + KinO2) * v['NO3']['c'][i] / (v['NO3']['c'][i] + KNO3)
        nit[i] = knit * Fnit(t) * v['O2']['c'][i] / (v['O2']['c'][i] + KO2) * v['NH4']['c'][i] / (v['NH4']['c'][i] + KNH4)
        O2_ex[i] = (vp[i] / DEPTH[i]) * (O2sat(t, i) - v['O2']['c'][i])
        NEM[i] = NPP[i] - aer_deg[i] - denit[i]

        # Update biogeochemical state variables [mmol m^-3]
        v['DIA']['c'][i] = v['DIA']['c'][i] + (NPP[i] - phy_death[i]) * DELTI
        v['dSi']['c'][i] = v['dSi']['c'][i] - Si_consumption[i] * DELTI
        v['NO3']['c'][i] = v['NO3']['c'][i] + (-94.4 / 106.0 * denit[i] + nit[i] - redn * NPP_NO3[i]) * DELTI
        v['NH4']['c'][i] = v['NH4']['c'][i] + (redn * (aer_deg[i] - NPP_NH4[i]) - nit[i]) * DELTI
        v['PO4']['c'][i] = v['PO4']['c'][i] + redp * (aer_deg[i] + denit[i] - NPP[i]) * DELTI
        v['O2']['c'][i] = v['O2']['c'][i] + (-aer_deg[i] + NPP_NH4[i] + (138.0 / 106.0) * NPP_NO3[i] - 2.0 * nit[i] + O2_ex[i]) * DELTI
        v['TOC']['c'][i] = v['TOC']['c'][i] + (-aer_deg[i] - denit[i] + phy_death[i]) * DELTI

# Write biogeochemical process rates [mmol m^-3 s-1]
    if (float(t) / float(TS * DELTI)) % 1 == 0:
        Rates(NPP, "NPP.dat", t)
        Rates(aer_deg, "aer_deg.dat", t)
        Rates(denit, "denit.dat", t)
        Rates(nit, "nit.dat", t)
        Rates(O2_ex, "O2_ex.dat", t)
        Rates(NEM, "NEM.dat", t)
        Rates(phy_death, "phy_death.dat", t)
