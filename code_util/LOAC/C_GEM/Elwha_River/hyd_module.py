"""
Hydrodynamics module (translated from hyd.c)
"""

from config import M, M1, M2, TOL, TS, DELTI
from variables import D, Dold2, TH, E, TU, Y, DEPTH
from uphyd_module import new_bc, new_uh, update
from tridag_module import coeff_a, tridag, conv
from file_module import hydwrite

def hyd(t):
    """Main hydrodynamics routine."""
    #print("DEPTH id at start of hyd():", id(DEPTH), "first 5 values:", DEPTH[:5])
    #print("DEPTH values at start of hyd():", DEPTH)
    # Store previous depths
    for i in range(1, M + 1):
        Dold2[i] = D[i]

    # Set new boundary conditions
    new_bc(t)
    rsum = 0.0
    i = 0

    # Iterative process to converge solution
    while rsum != 2.:
        i += 1
        coeff_a(t)
        tridag()
        rsum = conv(3, M1, TOL, TH, E) + conv(2, M2, TOL, TU, Y)
        update()
        #print("DEPTH values at end of hyd() iteration:", DEPTH)

    # Update U and H values
    new_uh()

    # Write hydrodynamic data to files at specified intervals
    if (float(t) / float(TS * DELTI)) % 1 == 0:
        hydwrite(t)

