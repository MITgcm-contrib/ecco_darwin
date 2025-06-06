"""
Main model routine (translated from main.c)
"""
import os
import math
from os import listdir
from typing import Union, Any

from file_module import exfread
from init_module import init
from hyd_module import hyd
from transport_module import transport
from biogeo_module import biogeo
from sed_module import sed
from config import M, MAXT, DELTI, WARMUP, TS, DEPTH_lb, B_lb, PI, G, EL, B_ub, DELXI
from variables import dispersion

def main():
    """Main model routine."""

    # Remove existing output files
    for fileName in listdir("."):
        # Check file extension
        if fileName.endswith('.dat'):
            # Remove File
            os.remove(fileName)

    init()  # Initialize the model

    # Main simulation loop
    for t in range(0, MAXT + 1, DELTI):
        # Print the number of simulated days
        print(f"t: {float(t) / (24.0 * 60.0 * 60.0):.2f} days")

        # Read forcing files
        # River discharge
        Qr, previousdays = exfread('/Users/rsavelli/Documents/CMS_LOAC/forcings_NorthSlope/colville_river_discharge_2022_m3sec.csv', t)
        Qr = - Qr

        # Convergence length [m]
        LC = 1.0 / -(1.0 / float(EL) * math.log(float(B_ub) / float(B_lb)))
        # Van der Burgh's coefficient [-]
        K = 4.38 * DEPTH_lb ** 0.36 * B_lb ** -0.21 * LC ** -0.14
        # Dispersion coefficient
        for i in range(M + 1):
            N = -PI * Qr / (DEPTH_lb * B_lb)
            D0 = 26 * math.sqrt(N * G) * DEPTH_lb ** 1.5
            beta = -(K * LC * Qr) / (D0 * B_lb * DEPTH_lb)
            d = D0 * (1.0 - beta * (math.exp((i * DELXI) / LC) - 1.0))
            dispersion[i] = max(d, 0.0)

        # Wind
        Uw_sal, previousdays = exfread('/Users/rsavelli/Documents/CMS_LOAC/forcings_NorthSlope/windspeed.csv', t)
        Uw_tid, previousdays = exfread('/Users/rsavelli/Documents/CMS_LOAC/forcings_NorthSlope/windspeed.csv', t)

        # Solar radiation
        I0, previousdays = exfread('/Users/rsavelli/Documents/CMS_LOAC/forcings_NorthSlope/solarradiation.csv', t)
        # Air Temperature
        # pCO2
        pCO2, previousdays = exfread('/Users/rsavelli/Documents/CMS_LOAC/forcings_NorthSlope/pCO2_barrow_2022.csv', t)  #microatm
        pCO2 = pCO2 * 1e-6  # atm
        # Water Temperature
        water_temp, previousdays = exfread('/Users/rsavelli/Documents/CMS_LOAC/forcings_NorthSlope/watertemp.csv', t)

        # Run hydrodynamics and transport processes
        hyd(t, Qr)
        transport(t, previousdays)

        # Run biogeochemical and sediment processes after warmup
        if t > WARMUP:
            biogeo(t, Uw_sal, Uw_tid, water_temp, pCO2, I0, previousdays)
            sed(t, previousdays)

if __name__ == "__main__":
    main()
