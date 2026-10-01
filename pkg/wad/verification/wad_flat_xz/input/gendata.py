#!/usr/bin/env python3
"""Inputs for wad_flat_xz: a stratified ocean at rest over a 20 m shelf
(0 < x < 5 km) and a beach rising to +2 m at x = 10 km, so that the top
of the beach is dry (film) at rest. The model rest level is DATUM = 3 m
above sea level, so every cell is an ocean cell. T depends on z only
(20 C at z = 0 to 10 C at z = -20 m in model coordinates): any velocity
is spurious (r* pressure-gradient error + WAD). Phase-2 gate of
MITgcm_WAD_design.md: max|u| < 1 mm/s.
Run from the input directory:  python3 gendata.py
"""
import numpy as np
NX, DX, NR, DZ = 100, 100.0, 10, 2.5
DATUM = 3.0                      # model rest level above sea level
x = (np.arange(NX) + 0.5) * DX
zb = np.where(x < 5000, -20.0, -20.0 + 22.0 * (x - 5000) / 5000)   # 20 m shelf -> +2 m
ocean = np.arange(NX) < NX - 1
r_low = np.where(ocean, zb - DATUM, 0.0)
r_low.astype('>f8').tofile('bathy.bin')
eta = np.where(ocean, np.maximum(0.0 - DATUM, r_low + 0.05), 0.0)       # sea level 0, film above
eta.astype('>f8').tofile('eta.bin')
# linear stratification in z (not r*): T = 20 at z=0 .. 10 at z=-20, same T in every column
zc = -DATUM - (np.arange(NR) + 0.5) * DZ         # rest-level r centres; r* maps them per column
t = np.zeros((NR, 1, NX))
for i in range(NX):
    D = eta[i] - r_low[i] if ocean[i] else 0.0
    for k in range(NR):
        # r* level centre -> actual z of layer centre in this column
        z = eta[i] + (zc[k] + DATUM) * D / (-r_low[i] if ocean[i] else 1.0) if ocean[i] else 0.0
        t[k, 0, i] = 20.0 + 0.5 * max(z, -20.0)
t.astype('>f8').tofile('T.bin')
print('ok', zb.max(), (eta - r_low)[ocean].min())
