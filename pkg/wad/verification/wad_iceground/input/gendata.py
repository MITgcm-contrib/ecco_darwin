#!/usr/bin/env python3
"""wad_iceground: grounded sea ice that melts and floats (pkg/wad +
pkg/seaice with mass loading, useRealFreshWaterFlux).

Closed basin, 60 x 10 cells of 50 m, 5 r* levels, rest level DATUM = 1 m
above MSL. Bed from -2 m MSL (west) to +0.3 m (east). 0.8 m of sea ice
everywhere (AREA 1), loading the surface: the ice bottom (the model
surface) sits a draft (0.71 m) below MSL, and where the water is shallower
than the draft the column is at the film (grounded ice). A warm constant
atmosphere melts the ice (15 C, 300 W/m2, 5 m/s); as it thins, water flows
back under it and the grounding line moves east. Writes bathy.bin,
eta.bin, heff.bin, area.bin (64-bit); heff_part.bin, area_part.bin,
bathy_deep.bin, eta_part.bin (deeper shelf, ice east of 1 km only, open
water at MSL, loaded ice in equilibrium under the ramped load of
wadLoadDepth = 0.3 m) for input.dyn.
"""
import numpy as np

NX, NY, DX = 60, 10, 50.0
DATUM = 1.0
H0 = 0.8
RHO_I, RHO_W = 910.0, 1025.0
x = (np.arange(NX) + 0.5) * DX
zb = np.repeat((-2.0 + 2.3 * x / (NX * DX))[None, :], NY, 0)     # m MSL
land = np.zeros((NY, NX), bool)
land[0, :] = land[-1, :] = True
land[:, 0] = land[:, -1] = True
draft = H0 * RHO_I / RHO_W
eta = np.maximum(-draft, zb + 0.05)                                # ice bottom or film
bathy = np.where(land, 0.0, zb - DATUM)
w = lambda n, a: a.astype('>f8').tofile(n)
w('bathy.bin', bathy)
w('eta.bin', np.where(land, 0.0, eta - DATUM))
w('heff.bin', np.where(land, 0.0, H0))
w('area.bin', np.where(land, 0.0, 1.0))
# input.dyn: ice only east of 1 km (open water to drift into)
part = (~land) & (x[None, :] > 1000.0)
w('heff_part.bin', np.where(part, H0, 0.0))
w('area_part.bin', np.where(part, 1.0, 0.0))
# input.dyn: a deeper shelf (-4 m to +0.3 m MSL) so that the floating ice
# has 0.7-2.4 m of water under it (Colville: 2-6 m), ice east of 1 km only,
# no side walls (periodic in y, a wide shelf: in a 400 m channel the wall
# friction holds even weak ice against the wind)
land2 = np.zeros((NY, NX), bool)
land2[:, 0] = land2[:, -1] = True
part2 = (~land2) & (x[None, :] > 1000.0)
zb2 = np.repeat((-4.0 + 4.3 * x / (NX * DX))[None, :], NY, 0)
# loaded ice in equilibrium with sea level 0 under the ramped load
# (data.wad wadLoadDepth = LOAD_D): a column with D of water under the ice
# passes the fraction f(D) = (D - CRIT)/(LOAD_D - CRIT) (0..1) of the load,
# so it floats (D >= LOAD_D) where zb + LOAD_D + draft <= 0, and otherwise
# settles where zb + D + f(D) draft = 0 (grounded ice partly on the bed),
# or at the film above sea level
CRIT, FILM, LOAD_D = 0.10, 0.05, 0.3
dfloat = -draft - zb2
dramp = (-zb2 + CRIT * draft / (LOAD_D - CRIT)) / (1.0 + draft / (LOAD_D - CRIT))
dice = np.where(dfloat >= LOAD_D, dfloat, np.clip(dramp, FILM, LOAD_D))
dice = np.where(zb2 + CRIT > 0, FILM, dice)
eta2 = np.where(part2, zb2 + dice, np.maximum(0.0, zb2 + 0.05))
w('bathy_deep.bin', np.where(land2, 0.0, zb2 - DATUM))
w('eta_part.bin', np.where(land2, 0.0, eta2 - DATUM))
w('heff_part.bin', np.where(part2, H0, 0.0))
w('area_part.bin', np.where(part2, 1.0, 0.0))
# no snow (pkg/seaice starts with 0.2 m x AREA without a snow file, which
# the equilibrium surface above does not include)
w('hsnow_part.bin', np.zeros((NY, NX)))
g = (~land) & (zb + draft > 0)
print(f'draft {draft:.2f} m; grounded at start: {g.sum()} of {(~land).sum()} cells')
