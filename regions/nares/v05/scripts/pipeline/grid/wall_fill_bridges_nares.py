#!/usr/bin/env python3
"""Nares: wall off min-depth FILL bridges (2026-10-15, after run_arch_TC3; Dustin "Yes" to steps 1-2).

The TC3 hotspots (i1109-1117 j28-37: 2.5 m/s; i666-669 j186-189: 2.1 m/s; i493-494 j326: 1.6 m/s) are all thin channels of
cells that are exactly -10.0 m in Nares_bathymetry.bin, i.e. near-dry/land GEBCO cells that the deepen-only 10 m rule made wet,
which short-circuit two parts of the domain that are already connected elsewhere (sill jets). Rule:
  fill  = wet & (bathy == -10.0) & no ice shelf
  A. shortcut channels: a 4-connected fill cluster whose adjacent non-fill water splits into >= 2 groups that are not connected
     to each other inside a local box (cluster bbox + 5 cells), and whose removal keeps the wet mask ONE 4-connected component
     -> land. (First version required >= 2 GLOBAL components and so missed i666 j188, where both sides join elsewhere.)
  B. sole bridges: a fill cluster whose removal cuts a pocket off the main body -> cluster AND pocket to land, but only where
     TC3 had |u| > 1 m/s in the cluster (the evidence); other sole bridges are listed for review, not changed. A pocket with
     open-boundary cells is reported so data.obcs OB_I*/OB_J* can be zeroed there.
Edge band (2 cells) never touched (OBCS rows). Writes a NEW file (never edits the input):
  <out bathy> + <out>.changed.npz (j, i, depth_old) and prints every removed cluster.
Usage: wall_fill_bridges_nares.py <Nares_bathymetry.bin (10 m)> <Nares_shelficeTopo.bin> <out bathy.bin> [smax_hourly.npy]
"""
import os, sys
import numpy as np
from scipy import ndimage

B, SH, OUT = sys.argv[1:4]
SMAX = sys.argv[4] if len(sys.argv) > 4 else None
NX, NY = 1120, 360
if os.path.exists(OUT):
    sys.exit(f'FATAL: {OUT} exists')
b = np.fromfile(B, '>f8').reshape(NY, NX)
sh = np.fromfile(SH, '>f8').reshape(NY, NX)
wet = (b < 0) & (sh > b)
c4 = ndimage.generate_binary_structure(2, 1)
lab0, n0 = ndimage.label(wet, c4)
assert n0 == 1, f'input wet mask has {n0} 4-connected components'
edge = np.zeros_like(wet); edge[:2] = edge[-2:] = True; edge[:, :2] = edge[:, -2:] = True
fill = wet & (b == -10.0) & (sh == 0) & ~edge
rest, nr = ndimage.label(wet & ~fill, c4)
fl, nf = ndimage.label(fill, c4)
smax = np.load(SMAX) if SMAX else None
out = wet.copy()
removed, review, pockets = [], [], []
for p in range(1, nf + 1):
    m = fl == p
    ys, xs = np.nonzero(m)
    j0, j1, i0, i1 = max(ys.min() - 5, 0), min(ys.max() + 6, NY), max(xs.min() - 5, 0), min(xs.max() + 6, NX)
    loc = (wet & ~fill)[j0:j1, i0:i1]
    ring = (ndimage.binary_dilation(m, c4) & ~m & wet & ~fill)[j0:j1, i0:i1]
    if not ring.any():
        continue
    ll, _ = ndimage.label(loc, c4)
    ngrp = np.unique(ll[ring]).size
    trial = out & ~m
    lab, n = ndimage.label(trial, c4)
    sm = float(smax[m].max()) if smax is not None else np.nan
    if n == 1:
        if ngrp >= 2:
            out = trial
            removed.append(('A', m.sum(), xs.min(), xs.max(), ys.min(), ys.max(), ngrp, sm, 0, 0))
        continue
    sizes = ndimage.sum(np.ones_like(lab), lab, range(1, n + 1))
    main = sizes.argmax() + 1
    pocket = trial & (lab != main)
    nob = pocket[:, :2].sum() + pocket[:, -2:].sum() + pocket[:2, :].sum() + pocket[-2:, :].sum()
    rec = ('B', m.sum(), xs.min(), xs.max(), ys.min(), ys.max(), ngrp, sm, int(pocket.sum()), int(nob))
    if sm > 1.0:
        out = out & ~m & ~pocket
        removed.append(rec); pockets.append(pocket)
    else:
        review.append(rec)
print(f'wet {wet.sum()}, fill (exactly 10 m, open water, off edge) {fill.sum()} in {nf} clusters')
hdr = ' rule ncell   i range     j range   groups  TC3 max|u|  pocket cells  pocket OB cells'
pr = lambda r: print(f'  {r[0]}  {r[1]:5d}  {r[2]:4d}-{r[3]:4d}  {r[4]:3d}-{r[5]:3d}   {r[6]:3d}     {r[7]:5.2f}   {r[8]:8d}  {r[9]:8d}')
print(f'CHANGED ({len(removed)} clusters):'); print(hdr)
for r in sorted(removed, key=lambda r: -np.nan_to_num(r[7])): pr(r)
print(f'REVIEW: sole bridges left unchanged (TC3 |u| <= 1 m/s), {len(review)}:'); print(hdr)
for r in sorted(review, key=lambda r: -r[8]): pr(r)
lab, n = ndimage.label(out, c4)
assert n == 1
gone = wet & ~out
bn = b.copy(); bn[gone] = 0.0
jj, ii = np.nonzero(gone)
np.savez(OUT + '.changed.npz', j=jj, i=ii, depth_old=b[gone])
bn.astype('>f8').tofile(OUT)
if smax is not None:
    for name, (i0, i1, j0, j1) in {'SE corner': (1109, 1117, 28, 37), 'i666 j188': (666, 669, 186, 189),
                                   'i494 j326': (493, 494, 326, 326)}.items():
        box = np.zeros_like(wet); box[j0:j1 + 1, i0:i1 + 1] = True
        print(f'  hotspot {name}: fast (>1 m/s) cells now land {(gone & box & (smax > 1)).sum()} of {(wet & box & (smax > 1)).sum()}')
print(f'wrote {OUT}: wet {out.sum()} (one 4-connected component); min open-water depth {-(bn[out & (sh == 0)]).max():.1f} m')
print('WALL_DONE')
