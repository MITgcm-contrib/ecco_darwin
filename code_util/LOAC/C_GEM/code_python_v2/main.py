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
import numpy as np
from file_module import OutputManager
np.set_printoptions(precision=16)

def main():
    """Main model routine."""
    io = OutputManager(prefix="outputs/", overwrite_all=True)

    init()  # Initialize the model
    try:
        # Main simulation loop
        for t in range(0, MAXT + 1, DELTI):
            # Print the number of simulated days
            print(f"t: {float(t) / (24.0 * 60.0 * 60.0):.2f} days")

            # Run hydrodynamics and transport processes
            hyd(t, io)
            transport(t, io)

            # Run biogeochemical and sediment processes after warmup
            if t > WARMUP:
                biogeo(t, io)
                sed(t, io)
    finally:
        io.close()

if __name__ == "__main__":
    main()
