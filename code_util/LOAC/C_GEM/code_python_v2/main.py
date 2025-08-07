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
            #sed(t)

if __name__ == "__main__":
    main()
