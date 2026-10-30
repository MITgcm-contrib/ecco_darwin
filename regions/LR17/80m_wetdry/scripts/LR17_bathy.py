#!/usr/bin/env python3
"""Wetting-drying bathymetry for the LR17 80 m grid (Pertuis Charentais).

Usage:  python3 LR17_bathy.py [ZMAX] [--zmouth=Z] [--ztip=Z] [--out=NAME] [--nofig]
  --zmouth, --ztip: thalweg profile ends (m MSL, default -2.5 and -0.5; v4: -4 and -3);
  --out: output file name in data/LR17/work (default LR17_zb.bin); --nofig: no figure.
  LR17_ROOT (env): root that holds data/ (default: the parent of scripts/).
Grid: 936 x 875, 0.0008 deg (62 x 89 m), lon -1.7001 + ..., lat 45.7499 + ...
(Adeola's LR17 80 m Darwin P1). Source: data/LR17/bathymetry/bathy_LR17_936x875,
the 80 m mean of the SHOM HOMONIM 20 m DEM (WGS84, mean sea level 'NM',
negative down), which holds the intertidal flats and the land.
Floodable land (option A, 2026-09-29): tidal flats and salt marsh seaward of
the sea defenses only. OpenStreetMap natural=coastline (which follows the
dikes here), man_made=dyke and man_made=embankment (data/LR17/osm) are drawn
onto the grid as barriers; the model keeps every cell reached from the open
sea (the western and southern edges) through cells below ZMAX (default 5 m
MSL) without crossing a barrier, plus the main rivers upstream of the
defenses: the OSM centerlines of the Charente, Seudre, Sevre Niortaise and
Lay (data/LR17/osm/osm_mainrivers.json), 3 cells wide, as far as they stay
below ZMAX and connected. Only the river cells extend inland, so the fill
cannot spill out of a river into the marsh (letting channels below MSL, or
every OSM river area, cross the barriers leaked through polder drains and
marsh channels into 280-460 km2 of polder).
Polders behind the dikes stay land however low they are.
River thalwegs (2026-09-30): the DEM gives the carved rivers their banks or
water surface (1-2.6 m above MSL, with sills), so each carved river cell is
capped by a profile from Z_TIP (-0.5 m MSL) at the upstream end to Z_MOUTH
(-2.5 m) where it meets deeper water: the Lay and Sevre now reach the sea.
Writes data/LR17/work/LR17_zb.bin (float32 big-endian, m MSL, land = NaN
... written as +99) and figures/LR17/LR17_bathymetry.{pdf,png}.

Figure caption: LR17 80 m bathymetry. a, Model bed elevation (m MSL;
black line: 0 m MSL; gray: land, above ZMAX or behind the sea defenses).
b, Floodable cells (bed -4 to +5 m MSL, blue) and low polders kept dry
(-1 to +3 m MSL, orange), with the OSM coastline and dikes (black).
"""
import json
from collections import deque
import os
import sys

import numpy as np
from scipy import ndimage
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import wetdrystyle as ws  # noqa: E402

TOP = os.environ.get('LR17_ROOT', os.path.dirname(HERE))
D = os.path.join(TOP, 'data', 'LR17')
NX, NY, DG = 936, 875, 0.0008
LON0, LAT0 = -1.7001, 45.7499
OPT = dict(a[2:].split('=', 1) if '=' in a else (a[2:], '1') for a in sys.argv[1:] if a.startswith('--'))
POS = [a for a in sys.argv[1:] if not a.startswith('--')]
ZMAX = float(POS[0]) if POS else 5.0
LAND = 99.0
# m MSL: river thalweg at the sea end and the upstream end (v1-v3: -2.5, -0.5)
Z_MOUTH, Z_TIP = float(OPT.get('zmouth', -2.5)), float(OPT.get('ztip', -0.5))
Z_SEA = -2.5                    # m MSL: a river mouth is a river cell next to sea deeper than this
OUT = OPT.get('out', 'LR17_zb.bin')


def rasterize(elements):
    """barrier cells from OSM ways (sampled at a fifth of a cell)"""
    bar = np.zeros((NY, NX), bool)
    for e in elements:
        g = e.get('geometry', [])
        for p, q in zip(g[:-1], g[1:]):
            n = int(max(abs(q['lon'] - p['lon']), abs(q['lat'] - p['lat'])) / (DG / 5)) + 2
            lo = np.linspace(p['lon'], q['lon'], n)
            la = np.linspace(p['lat'], q['lat'], n)
            i = np.floor((lo - LON0) / DG).astype(int)
            j = np.floor((la - LAT0) / DG).astype(int)
            ok = (i >= 0) & (i < NX) & (j >= 0) & (j < NY)
            bar[j[ok], i[ok]] = True
    return bar


def fill_polys(elements):
    """cells inside OSM areas: each way, or all member ways of a relation,
    drawn as a boundary and its holes filled"""
    out = np.zeros((NY, NX), bool)
    for e in elements:
        if e['type'] == 'way':
            parts = [e]
        else:
            parts = [dict(geometry=m.get('geometry', [])) for m in e.get('members', [])
                     if m.get('role') in ('outer', '') and m.get('geometry')]
        edge = rasterize(parts)
        out |= ndimage.binary_fill_holes(edge)
    return out


def main(write=True, fig=True):
    b = np.fromfile(os.path.join(D, 'bathymetry', 'bathy_LR17_936x875'), '>f4') \
        .reshape(NY, NX).astype(float)
    osm = json.load(open(os.path.join(D, 'osm', 'osm_coast_dykes.json')))['elements']
    bar = rasterize(osm)
    rivers = [e for e in json.load(open(os.path.join(D, 'osm', 'osm_mainrivers.json')))['elements']
              if e['tags'].get('name') != 'Le Petit Lay']
    river = ndimage.binary_dilation(rasterize(rivers), np.ones((3, 3), bool)) & (b < ZMAX)
    lab, _ = ndimage.label((b < ZMAX) & ~bar)
    seeds = np.unique(np.r_[lab[:, 0][b[:, 0] < -2], lab[0, :][b[0, :] < -2]])
    wet0 = np.isin(lab, seeds[seeds > 0])
    lab2, _ = ndimage.label(wet0 | river)
    wet = np.isin(lab2, np.unique(lab2[wet0]))
    zb = np.where(wet, b, LAND)
    # river thalwegs: the DEM holds the water surface or banks, not the
    # channel bed, so the carved rivers sat 1-2.6 m above MSL with sills and
    # ran dry between source and sea. Each carved river cell's bed is capped
    # by a profile deepening linearly (in steps along the channel) from
    # Z_TIP at the upstream end to Z_MOUTH where it meets deeper water.
    rv = river & wet
    mouth = rv & ndimage.binary_dilation(wet & ~rv & (zb < Z_SEA),
                                         np.array([[0, 1, 0], [1, 1, 1], [0, 1, 0]], bool))
    dist = np.full(zb.shape, -1)
    q = deque()
    for j, i in zip(*np.where(mouth)):
        dist[j, i] = 0
        q.append((j, i))
    while q:
        j, i = q.popleft()
        for jj, ii in ((j + 1, i), (j - 1, i), (j, i + 1), (j, i - 1)):
            if 0 <= jj < NY and 0 <= ii < NX and rv[jj, ii] and dist[jj, ii] < 0:
                dist[jj, ii] = dist[j, i] + 1
                q.append((jj, ii))
    lab_r, nr = ndimage.label(rv)
    carved = 0
    for k in range(1, nr + 1):
        c = (lab_r == k) & (dist >= 0)
        if not c.any():
            continue
        dmax = max(dist[c].max(), 1)
        prof = Z_MOUTH + (Z_TIP - Z_MOUTH) * dist / dmax
        low = c & (zb > prof)
        zb = np.where(low, prof, zb)
        carved += low.sum()
    print(f'river thalwegs: {carved} carved river cells lowered to the '
          f'{Z_MOUTH:+.1f} .. {Z_TIP:+.1f} m MSL profile')
    os.makedirs(os.path.join(D, 'work'), exist_ok=True)
    if write:
        zb.astype('>f4').tofile(os.path.join(D, 'work', OUT))
        print('wrote', os.path.join(D, 'work', OUT))
    cell = DG * 111.2e3 * np.cos(np.radians(46.1)) * DG * 111.2e3 / 1e6
    inter = wet & (b > -3) & (b < 3)
    polder = (b > -1) & (b < 3) & ~wet
    print(f'model cells {wet.sum()} ({wet.mean():.2f}); bed {b[wet].min():.1f} .. {b[wet].max():.2f} m; '
          f'-3..+3 m in the model {inter.sum() * cell:.0f} km2; low land (-1..+3 m) kept dry '
          f'(polders) {polder.sum() * cell:.0f} km2; barrier cells {bar.sum()}; '
          f'river-area cells {river.sum()} ({(river & wet).sum()} in the model)')
    if not fig or 'nofig' in OPT:
        return dict(b=b, zb=zb, wet=wet, rv=rv, dist=dist, lab_r=lab_r)

    ws.use()
    from matplotlib.colors import TwoSlopeNorm
    from matplotlib.lines import Line2D
    from matplotlib.patches import Patch
    from matplotlib.ticker import FuncFormatter, MultipleLocator
    lon = LON0 + (np.arange(NX) + 0.5) * DG
    lat = LAT0 + (np.arange(NY) + 0.5) * DG
    fig, ax = plt.subplots(1, 2, figsize=ws.size(183, 88), constrained_layout=True)
    for a in ax:
        a.set_facecolor('#CFCAB8')
        a.set_aspect(1 / np.cos(np.radians(46.1)))
        for sp in a.spines.values():
            sp.set_visible(True)
        a.xaxis.set_major_locator(MultipleLocator(0.2))
        a.yaxis.set_major_locator(MultipleLocator(0.2))
        a.xaxis.set_major_formatter(FuncFormatter(lambda x, _: f'{abs(x):.1f}°{"W" if x < 0 else "E"}'))
        a.yaxis.set_major_formatter(FuncFormatter(lambda y, _: f'{y:.1f}°N'))
    m = ax[0].pcolormesh(lon, lat, np.ma.masked_where(~wet, b), cmap='wetdry_topo',
                         norm=TwoSlopeNorm(0, -40, 5), shading='auto', rasterized=True)
    ax[0].contour(lon, lat, np.where(wet, b, np.nan), levels=[0], colors='k', linewidths=0.3)
    cb = ws.cbar(fig, m, ax[0], 'Bed elevation (m MSL)', shrink=0.8, aspect=25)
    cb.set_ticks([-40, -30, -20, -10, 0, 5])
    ax[0].legend(handles=[Patch(fc='#CFCAB8', ec='none',
                                label=f'Land (above {ZMAX:g} m MSL\nor behind sea defenses)'),
                          Line2D([], [], color='k', lw=0.5, label='0 m MSL')],
                 loc='upper right', borderaxespad=0.3)
    show = np.where(wet & (b > -4), 1, np.where(polder, 2, 0)).astype(float)
    cm = matplotlib.colors.ListedColormap(['#CFCAB8', '#56B4E9', '#E69F00'])
    ax[1].pcolormesh(lon, lat, show, cmap=cm, vmin=-0.5, vmax=2.5, shading='auto',
                     rasterized=True)
    for e in osm:
        g = e.get('geometry', [])
        ax[1].plot([p['lon'] for p in g], [p['lat'] for p in g], color='k', lw=0.4)
    ax[1].legend(handles=[Patch(fc='#56B4E9', label='Floodable, bed −4 to +5 m MSL'),
                          Patch(fc='#E69F00', label='Polder kept dry, −1 to +3 m MSL'),
                          Line2D([], [], color='k', lw=0.5, label='OSM coastline and dikes')],
                 loc='center left', borderaxespad=0.3)
    ax[1].set_xlim(LON0, LON0 + NX * DG)
    ax[1].set_ylim(LAT0, LAT0 + NY * DG)
    ws.letter(ax[0], 'a')
    ws.letter(ax[1], 'b')
    ws.save(fig, 'LR17_bathymetry')
    return dict(b=b, zb=zb, wet=wet, rv=rv, dist=dist, lab_r=lab_r)

if __name__ == '__main__':
    main()
