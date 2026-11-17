#!/usr/bin/env python3
"""Nares minimum depth (2026-10-15, Dustin D130: 20 m min depth for the tide tests; adapted from the WestAO
prep/deepen_min_depth.py used for the 3 m / 5 m / 10 m passes there).

Why: in the tide tests (run_arch_T, run_arch_TC2) spurious 2-5 m/s surface currents live in the 10 m min-depth cells (SE corner
of the Lincoln Sea, i 1104-1115 j 4-31; also i 494 j 326), start within hours of t=0 and, with tides, drive 10 m columns below
hFacInf (STOP in CALC_R_STAR). Every open-water wet cell shallower than HMIN is deepened to HMIN (mask unchanged).
Differences from the WestAO script:
  - wet = (bathy < 0) & (shelficeTopo > bathy) as in gen_ob_index.py / nares_ob_points.py; only cells with NO ice shelf
    (shelficeTopo == 0) are deepened (Petermann cavity untouched);
  - hFacMin = 0.3 (Nares data), 12 top-level thicknesses (same delR as WestAO), so a 20 m column gets 11 wet levels;
  - grid 1120 x 360 x 90.
Cells within 2 of a domain edge are left alone (the OBCS files carry no values for levels that were dry).
Newly wet levels in every 3-D pickup_*.data (90x360x1120 real*8) get the column's old bottom wet-level value.
Works IN PLACE on the files given (run it on copies). Backups: <bathy>.pre_min<H>m_<date>, <pickups>/min_depth_backup_<H>m.npz.
Usage: deepen_min_depth_nares.py <bathy.bin> <shelficeTopo.bin> <pickups dir> <run dir with hFacC.data> [HMIN=20] [iter]
"""
import datetime, glob, os, sys
import numpy as np

B, SH, P, G = sys.argv[1:5]
HMIN = float(sys.argv[5]) if len(sys.argv) > 5 else 20.0
IT = sys.argv[6] if len(sys.argv) > 6 else '0000604872'
NR, NY, NX = 90, 360, 1120
HFACMIN = 0.3
DELR_TOP = [1.00, 1.14, 1.30, 1.49, 1.70, 1.93, 2.20, 2.50, 2.84, 3.21, 3.63, 4.10]   # data delR, first 12 levels
bk_b = f'{B}.pre_min{HMIN:g}m_{datetime.date.today():%Y%m%d}'
bk_p = os.path.join(P, f'min_depth_backup_{HMIN:g}m.npz')
for f in [bk_b, bk_p] + glob.glob(f'{B}.pre_min{HMIN:g}m_*'):
    if os.path.exists(f):
        sys.exit(f'FATAL: {f} exists -- already applied')

d = np.fromfile(B, '>f8').reshape(NY, NX)
sh = np.fromfile(SH, '>f8').reshape(NY, NX)
hC = np.fromfile(f'{G}/hFacC.data', '>f4').reshape(NR, NY, NX)
nw_old = (hC > 0).sum(0)
wet = (d < 0) & (sh > d)
assert (wet == (nw_old > 0)).all(), f'wet mask differs from hFacC: {(wet != (nw_old > 0)).sum()} cells'
edge = np.zeros_like(wet); edge[:2] = edge[-2:] = True; edge[:, :2] = edge[:, -2:] = True
open_water = sh == 0
sel = wet & (-d < HMIN) & ~edge & open_water
print(f'wet cells shallower than {HMIN:g} m: {(wet & (-d < HMIN)).sum()}  (edge skipped: {(wet & (-d < HMIN) & edge).sum()}, '
      f'under ice shelf skipped: {(wet & (-d < HMIN) & ~edge & ~open_water).sum()}); deepened: {sel.sum()}')

rF = np.concatenate([[0.], np.cumsum(DELR_TOP)])
assert rF[-1] >= HMIN, 'extend DELR_TOP'
nw_new = 0
for k, dz in enumerate(DELR_TOP):
    h = min(max((HMIN - rF[k]) / dz, 0.), 1.)
    if h >= 0.5 * HFACMIN: nw_new = k + 1      # model/src/ini_masks_etc.F: dry if hFac < hFacMin/2
print(f'wet levels at {HMIN:g} m: {nw_new}; old wet levels in deepened cells: '
      + ', '.join(f'{k}:{(nw_old[sel] == k).sum()}' for k in range(1, nw_new + 1)))
assert (nw_old[sel] <= nw_new).all() and (nw_old[sel] >= 1).all()

jb, ib = np.nonzero(sel)
files = sorted(f for f in glob.glob(f'{P}/pickup_*.{IT}.data') if os.path.getsize(f) == NR * NY * NX * 8)
print(f'{len(files)} 3-D pickup files: {[os.path.basename(f) for f in files]}; {len(jb)} columns')
save = {'j': jb, 'i': ib, 'depth_old': d[sel], 'nw_old': nw_old[sel]}
for f in files:
    a = np.memmap(f, '>f8', 'r').reshape(NR, NY, NX)
    save[os.path.basename(f)] = np.array(a[:nw_new, jb, ib])
np.savez(bk_p, **save)
os.rename(B, bk_b)
dn = d.copy(); dn[sel] = -HMIN
dn.astype('>f8').tofile(B)
kb = nw_old[sel] - 1
for f in files:
    a = np.memmap(f, '>f8', 'r+').reshape(NR, NY, NX)
    col = a[:nw_new, jb, ib]
    bot = col[kb, np.arange(len(jb))]
    for k in range(1, nw_new):
        m = k > kb
        col[k, m] = bot[m]
    a[:nw_new, jb, ib] = col
    a.flush()
S = np.memmap(f'{P}/pickup_SALT.{IT}.data', '>f8', 'r').reshape(NR, NY, NX)[:nw_new, jb, ib]
print(f'new-depth columns: S range {S.min():.2f}-{S.max():.2f}; zero S values: {(S == 0).sum()}')
print(f'bathy: wet {wet.sum()} -> {((dn < 0) & (sh > dn)).sum()}, open-water min depth {-(dn[(dn < 0) & open_water & ~edge]).max():.2f} m')
print('MINDEPTH_DONE')
