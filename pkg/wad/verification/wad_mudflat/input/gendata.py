#!/usr/bin/env python3
"""Inputs for wad_mudflat: an idealized macrotidal mudflat with a tidal-creek
network, a river and passive tracers.

Grid: 100 x 60 x 5, dx = 100 m (x: cross-shore, 10 km; y: alongshore,
6 km), 5 levels of 2.5 m. Column i=0 is the open boundary (Flather M2 tide,
2 m amplitude, code/obcs_calc.F); the last column (i=99) and the last row
(j=59) are land. Elevations are relative to mean sea level:
  x < 2 km     offshore basin, 10 m deep
  2 - 3 km     slope from -10 m to -3 m
  3 - 10 km    mudflat rising from -3 m to +2 m (slope ~1:1400)
  creeks       a main creek along y = 3 km from the flat edge to the
               eastern wall (2.8 m deep at the flat edge, shoaling to
               0.3 m at x = 9.5 km, 0.3 m from there to the wall; the
               river enters at its landward end) and 6 side branches
               (1.2 m deep at the junction, or the main-creek depth there
               if that is less, shoaling to 0.2 m), carved
               into the flat with Gaussian cross-sections; seaward of the
               flat edge the main creek keeps its 2.8 m depth below the
               flat plane extended seaward, until it meets the slope
               (x ~ 2.55 km), so the creek mouth has no sill
Model rest level DATUM = 2.5 m above mean sea level (every cell is an ocean
cell); initial state at high water (+2 m), at rest.
River: 20 m3/s of fresh water (addMass, salt_addMass = 0) into the top level
of the last ocean cell of the main creek, next to the eastern wall.
Passive tracers: 1 = river water (relaxed to 1 at the river cell by
pkg/rbcs); 2 = "flat water" (1 over the mudflat, x > 3 km, at t = 0).
Usage:  python3 gendata.py [OUTDIR]
With no OUTDIR the files go to input/ next to this script; with OUTDIR
they go there. The grid can be refined for runs outside
testreport with WAD_NX (default 100), WAD_NY (60), WAD_NR (5): the domain
stays 10 x 6 km and 12.5 m deep (dx = 10000/NX, dz = 12.5/NR).
"""
import os
import sys

import numpy as np

NX = int(os.environ.get('WAD_NX', 100))
NY = int(os.environ.get('WAD_NY', 60))
NR = int(os.environ.get('WAD_NR', 5))
DX, DZ = 10000.0 / NX, 12.5 / NR
DATUM = 2.5
RIVER_Q = 20.0                          # m3/s
RHO = 1000.0

here = os.path.dirname(os.path.abspath(__file__))
x = (np.arange(NX) - 0.5) * DX          # cell centres; i=0 is the OB cell
y = (np.arange(NY) + 0.5) * DX
X, Y = np.meshgrid(x, y)                 # (NY, NX)

# --- bathymetry (m, relative to mean sea level)
zb = np.where(X < 2000, -10.0,
     np.where(X < 3000, -10.0 + 7.0 * (X - 2000) / 1000,
              -3.0 + 5.0 * (X - 3000) / 7000))
flat = X >= 3000

# creeks: main creek + side branches (segments), Gaussian cross-sections
def seg_dist(px, py, ax, ay, bx, by):
    """distance to segment a-b and position s in [0,1] along it"""
    vx, vy = bx - ax, by - ay
    s = np.clip(((px - ax) * vx + (py - ay) * vy) / (vx * vx + vy * vy), 0, 1)
    return np.hypot(px - (ax + s * vx), py - (ay + s * vy)), s

carve = np.zeros_like(zb)
XE = x[NX - 2]                          # last ocean column (wall at NX-1)
d, s = seg_dist(X, Y, 3000, 3000, XE, 3000)
dep = 2.8 - 2.5 * np.clip(s * (XE - 3000) / 6500, 0, 1)
carve = np.maximum(carve, dep * np.exp(-0.5 * (d / 150.0) ** 2))
# main creek through the top of the slope: floor at the flat plane
# (extended seaward) minus 2.8 m, where that is below the slope
dm, _ = seg_dist(X, Y, 2000, 3000, 3000, 3000)
zmouth = -3.0 + 5.0 * (X - 3000) / 7000 - 2.8 * np.exp(-0.5 * (dm / 150.0) ** 2)
for x0, side in ((4200, 1), (5000, -1), (5800, 1), (6600, -1), (7400, 1), (8200, -1)):
    d, s = seg_dist(X, Y, x0, 3000, x0 + 1400, 3000 + side * 1900)
    # junction depth: 1.2 m, but not deeper than the main creek there
    d0 = min(1.2, 2.8 - 2.5 * min(1.0, (x0 - 3000) / 6500))
    carve = np.maximum(carve, (d0 - (d0 - 0.2) * s) * np.exp(-0.5 * (d / 100.0) ** 2))
zb = np.where(flat, zb - carve, zb)
zb = np.where((X >= 2000) & (X < 3000), np.minimum(zb, zmouth), zb)

ocean = np.ones((NY, NX), bool)
ocean[:, -1] = False
ocean[-1, :] = False
r_low = np.where(ocean, zb - DATUM, 0.0)

# --- initial free surface: high water (+2 m), film where the bed is higher
eta = np.where(ocean, np.maximum(2.0 - DATUM, r_low + 0.05), 0.0)

# --- river: top-level cell at the landward end of the main creek
ir = NX - 2
jr = int(np.argmin(np.abs(y - 3000)))
add = np.zeros((NR, NY, NX))
add[0, jr, ir] = RIVER_Q * RHO          # kg/s
rmask = np.zeros((NR, NY, NX))
rmask[:, jr, ir] = 1.0

# --- passive tracers
ptr1 = np.zeros((NR, NY, NX))                     # river water
ptr2 = np.where((flat & ocean)[None].repeat(NR, 0), 1.0, 0.0)   # flat water

OUT = sys.argv[1] if len(sys.argv) > 1 else None


def w(name, a, sub='input'):
    d = OUT or os.path.join(os.path.dirname(here), sub)
    a.astype('>f8').tofile(os.path.join(d, name))

w('bathy.bin', r_low)
w('eta.bin', eta)
w('addMass.bin', add)
w('rbcs_mask.bin', rmask)
w('rbcs_ptr1.bin', rmask)
w('ptr1.bin', ptr1)
w('ptr2.bin', ptr2)
print(f'bed range {zb[ocean].min():.2f} .. {zb[ocean].max():.2f} m;'
      f' river cell (i,j)=({ir},{jr}), bed {zb[jr, ir]:.2f} m;'
      f' intertidal (bed -2..+2 m): {100 * np.mean((zb > -2) & ocean):.0f} % of cells')
