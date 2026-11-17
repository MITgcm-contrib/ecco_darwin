#!/usr/bin/env python3
"""Nares run_arch_ICT (2026-10-20, Dustin "stage IC thin"): thin the Jan 1 2015 sea-ice IC at the Robeson mouth.
In i 840-920 (0-based columns 840..919, the 'mouth' region of the arch analyses), where the floe thickness HEFF/AREA > 2.0 m,
set HEFF = AREA * 2.0 m (CryoSat-2 Jan 1-7 2015 median 1.98 m, see ic_vs_cs2smos.py). AREA, HSNOW, UICE, VICE unchanged.
Input/output: big-endian real*8 1120x360 (the pickup_seaice_si*.0000604872.bin format).
Usage: make_ic_thin_mouth.py <in HEFF> <in AREA> <out HEFF>
"""
import sys
import numpy as np

FH, FA, FO = sys.argv[1:4]
NX, NY, HMAX, I0, I1 = 1120, 360, 2.0, 840, 920
h = np.fromfile(FH, '>f8').reshape(NY, NX); a = np.fromfile(FA, '>f8').reshape(NY, NX)
o = h.copy()
m = np.zeros_like(h, bool); m[:, I0:I1] = (a[:, I0:I1] > 0) & (h[:, I0:I1] > HMAX * a[:, I0:I1])
o[m] = HMAX * a[m]
assert (o <= h).all() and np.isfinite(o).all()
print(f'{m.sum()} cells thinned; mouth volume {h[:, I0:I1].sum():.0f} -> {o[:, I0:I1].sum():.0f} m (sum of HEFF); '
      f'max floe thickness in mouth {np.nanmax(np.where(a[:, I0:I1] > 0, o[:, I0:I1] / np.maximum(a[:, I0:I1], 1e-9), 0)):.3f} m; '
      f'cells changed outside mouth: {(o != h)[:, np.r_[0:I0, I1:NX]].sum()}')
o.astype('>f8').tofile(FO)
