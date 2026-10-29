#!/usr/bin/env python3
"""Inputs for sediment_seiche: a closed 20 km x 5 m basin (x-z), 40 columns
of 500 m with land at both ends, a 0.5 m linear surface tilt at t=0 and a
uniform bed of mud and fine sand (kg/m^2). Writes big-endian float32 files
in the current directory."""
import numpy as np

nx, ny = 40, 1
H = 5.0           # basin depth (m)
eta0 = 0.5        # surface displacement at the end walls (m)
bedMud = 5.0      # initial bed mass of mud (kg/m^2)
bedSand = 5.0     # initial bed mass of sand (kg/m^2)


def write(name, a):
    a.astype('>f4').tofile(name)


bathy = -H*np.ones((ny, nx))
bathy[:, 0] = 0.0
bathy[:, -1] = 0.0
write('bathy.bin', bathy)

# linear tilt, high at the west wall, low at the east wall, zero mean
x = (np.arange(nx) - (nx - 1)/2.0)/((nx - 3)/2.0)
eta = -eta0*x*np.ones((ny, nx))
eta[bathy == 0] = 0.0
write('eta.bin', eta)

write('bed_mud.bin', bedMud*(bathy < 0))
write('bed_sand.bin', bedSand*(bathy < 0))
