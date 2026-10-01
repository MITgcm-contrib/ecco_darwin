#!/usr/bin/env python3
"""Inputs for wad_balzano (Balzano 1998-style tidal flat tests).

Geometry (x = 0 at the open boundary face, 100 m cells, 140 columns):
  column 1         open-boundary cell (x < 0)
  columns 2..139   the beach, 0 < x < 13.8 km
  column 140       land (wall)
Bed below mean sea level:  z_b(x) = -5 + 5 x / 13800   ("slope")
  step : z_b + 1 m for x > 6900 m (a 1 m bed step half way up the beach)
  pool : a 0.5 m sill on 9200 < x < 9600 m and a 1 m deep pool on
         9600 < x < 11400 m; at low tide the pool must hold water at
         the sill crest (+ at most wadCritDepth)
These follow the spirit of Balzano's tests 1, 3 and 5 (sloping beach,
beach with a discontinuity, beach with a pond); the exact published
geometries differ.

Model rest level is DATUM = 2.5 m above mean sea level (= tideDatum in
code/obcs_calc.F): R_low = z_b - DATUM, eta_model = eta - DATUM.
Initial state: high water (eta = +2 m), at rest.

Usage:  python3 gendata.py     (run from any directory; writes
        input/bathy.bin, input.step/bathy.bin, input.pool/bathy.bin,
        input/eta_init.bin)
"""
import os
import numpy as np

NX, DX = 140, 100.0
DATUM, AMP = 2.5, 2.0
L = 13800.0

x = (np.arange(NX) - 1 + 0.5) * DX        # column 1 centred at -50 m
beach = (np.arange(NX) >= 1) & (np.arange(NX) <= NX - 2)
ocean = np.arange(NX) <= NX - 2            # OB cell + beach


def bed(kind):
    zb = -5.0 + 5.0 * x / L
    if kind == 'step':
        zb = zb + np.where(x > 6900.0, 1.0, 0.0)
    elif kind == 'pool':
        zb = zb + np.where((x > 9200.0) & (x < 9600.0), 0.5, 0.0)
        zb = zb - np.where((x > 9600.0) & (x < 11400.0), 1.0, 0.0)
    return zb


def r_low(kind):
    return np.where(ocean, bed(kind) - DATUM, 0.0)


if __name__ == '__main__':
    top = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    for kind, sub in (('slope', 'input'), ('step', 'input.step'),
                      ('pool', 'input.pool')):
        r_low(kind).astype('>f8').tofile(os.path.join(top, sub, 'bathy.bin'))
    eta = np.where(ocean, AMP - DATUM, 0.0)
    eta.astype('>f8').tofile(os.path.join(top, 'input', 'eta_init.bin'))
    print('wrote bathy.bin (slope, step, pool) and eta_init.bin')
