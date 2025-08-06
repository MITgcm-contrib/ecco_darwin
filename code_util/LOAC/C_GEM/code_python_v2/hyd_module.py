"""
Hydrodynamics module (translated from hyd.c)
"""

from config import M, DELTI, TS
from variables import D, Dold2
from uphyd_module import new_bc, new_uh, update
from tridag_module import tridag
from file_module import hydwrite

def hyd(t):
    """Main hydrodynamics routine."""
    # Store previous depths
    for i in range(1, M + 1):
        Dold2[i] = D[i]

    # Set new boundary conditions
    new_bc(t)
    rsum = 0.0
    i = 0

    # Directly solve the tridiagonal matrix system
    tridag()
    update()

    # Update U and H values
    new_uh()

    # Write hydrodynamic data to files at specified intervals
    if (float(t) / float(TS * DELTI)) % 1 == 0:
        hydwrite(t)

