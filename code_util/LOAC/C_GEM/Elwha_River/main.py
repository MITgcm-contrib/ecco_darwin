"""
Main model routine (translated from main.c)
"""
import os
from os import listdir
from init_module import init
from hyd_module import hyd
from transport_module import transport
from biogeo_module import biogeo
from sed_module import sed
from config import MAXT, DELTI, WARMUP
from pre_model_stratup_check import check_output_writability

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

        # Run hydrodynamics and transport processes
        hyd(t)
        transport(t)

        # Run biogeochemical and sediment processes after warmup
        if t > WARMUP:
            biogeo(t)
            sed(t)

if __name__ == "__main__":

# list of all output files to write to
    output_files = [
    "depth.dat",     # water depth data
    "width.dat",     # estuarine width data
    "S.dat",         # salinity data
    "SPM.dat",       # suspended particulate matter concentrations
    "DIA.dat",       # diatoms concentrations
    "NH4.dat",       # ammonium concentrations
    "NO3.dat",       # nitrate concentrations
    "O2.dat",        # oxygen concentrations
    "PO4.dat",       # phosphate concentrations
    "DSi.dat",       # dissolved silica concentrations
    "TOC.dat",       # total organic carbon concentrations
    "NPP.dat",       # Net Primary Production rates
    "aer_deg.dat",   # aerobic degradation rates
    "denit.dat",     # denitrification rates
    "nit.dat",       # nitrification rates
    "O2_ex.dat",     # air-water O₂ exchange rates
    "NEM.dat",       # Net Ecosystem Metabolism rates

    # FCO2 version outputs:
    "DIC.dat",       # dissolved inorganic carbon concentration
    "ALK.dat",       # alkalinity
    "pH.dat",        # pH
    "FCO2.dat"       # air-water CO₂ exchange rates
    ]

    if not check_output_writability(output_files):
        exit("Aborting: one or more files are locked or not writable.")


    main()
