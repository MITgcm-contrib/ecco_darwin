#!/usr/bin/env python3
"""Check and plot a wad_balzano run.

Usage:  python3 balzano.py RUNDIR [--plot]
Reads Eta.*.data (and data for deltaT), the bathymetry (bathy.bin) and
the %WAD_MON lines of RUNDIR/output.txt. Reports:
  - minimum water depth (must stay >= wadMinDepth = 0.05 m),
  - volume budget residual (from %WAD_MON),
  - shoreline position at each low/high water of the last cycle vs the
    "bathtub" estimate (first bed point above the boundary elevation),
  - for the pool geometry: lowest pool level of the last cycle vs the
    sill crest (should be crest + at most wadCritDepth).
"""
import glob
import os
import re
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import gendata as G   # noqa: E402

DMIN, DCRIT = 0.05, 0.10
PERIOD = 43200.0


def run(rundir, plot=False):
    dt = None
    with open(os.path.join(rundir, 'data')) as f:
        for line in f:
            m = re.match(r'\s*deltaT\s*=\s*([0-9.eE+-]+)', line)
            if m:
                dt = float(m.group(1))
    r_low = np.fromfile(os.path.join(rundir, 'bathy.bin'), '>f8')
    zb = r_low + G.DATUM
    kind = min(('slope', 'step', 'pool'),
               key=lambda k: np.abs(G.r_low(k) - r_low).max())
    snaps = []
    for fn in sorted(glob.glob(os.path.join(rundir, 'Eta.*.data'))):
        it = int(fn.split('.')[-2])
        snaps.append((it * dt, np.fromfile(fn, '>f8')[:G.NX] + G.DATUM))
    dmin = min((e - zb)[G.beach].min() for _, e in snaps)
    print(f'{rundir}: geometry={kind}  {len(snaps)} snapshots')
    print(f'  min depth over beach, all times: {dmin:.4f} m '
          f'({"OK" if dmin >= DMIN - 1e-9 else "BELOW FILM"})')
    with open(os.path.join(rundir, 'output.txt')) as f:
        bud = [l.split() for l in f if 'budgetErr/vol0=' in l]
    if bud:
        print(f'  volume budget: cum. inflow {float(bud[-1][4]):.4e} m3, '
              f'residual {float(bud[-1][5]):.3e} m3 '
              f'(rel. {float(bud[-1][6]):.2e})')
    # shoreline at low / high water of the last cycle
    last = [s for s in snaps if s[0] >= 2 * PERIOD]
    for phase, name in ((0.0, 'high'), (0.5, 'low')):
        t_target = 2 * PERIOD + phase * PERIOD
        t, e = min(last, key=lambda s: abs(s[0] - t_target))
        eta_b = G.AMP * np.cos(2 * np.pi * t / PERIOD)
        dep = e - zb
        wet = np.nonzero((dep > DCRIT + 1e-6) & G.beach)[0]
        # connected wet region from the open boundary
        k = 1
        while k + 1 < G.NX - 1 and dep[k + 1] > DCRIT + 1e-6:
            k += 1
        above = np.nonzero(G.beach & (zb > eta_b))[0]
        x_bath = G.x[above[0]] - G.DX / 2 if above.size else G.x[-2] + G.DX / 2
        print(f'  {name} water t={t/3600:5.1f} h: boundary eta {eta_b:+.2f} m,'
              f' model shoreline {G.x[k] + G.DX/2:7.0f} m,'
              f' bathtub {x_bath:7.0f} m')
    if kind == 'pool':
        pool = G.beach & (G.x > 9600) & (G.x < 11400)
        crest = zb[G.beach & (G.x > 9200) & (G.x < 9600)].max()
        tl, lvl = min(((t, e[pool].mean()) for t, e in last),
                      key=lambda p: p[1])
        print(f'  pool: lowest level {lvl:+.3f} m at t={tl/3600:.1f} h,'
              f' sill crest {crest:+.3f} m, excess {lvl - crest:+.3f} m'
              f' (expect 0 .. {DCRIT})')
    if plot:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(figsize=(9, 4.5))
        xs = G.x[G.beach] / 1000
        ax.fill_between(xs, zb[G.beach], -5.5, color='0.8', lw=0)
        cyc = [s for s in snaps if s[0] >= 2 * PERIOD]
        for t, e in cyc[::3]:
            dep = (e - zb)[G.beach]
            ax.plot(xs, np.where(dep > DCRIT, e[G.beach], np.nan), lw=1,
                    label=f'{(t - 2*PERIOD)/3600:4.1f} h')
        ax.set_xlabel('x (km)')
        ax.set_ylabel('elevation (m)')
        ax.set_ylim(-5.5, 2.5)
        ax.set_title(f'wad_balzano ({kind}): water level over the 3rd cycle')
        ax.legend(fontsize=7, ncol=3, loc='lower left')
        fig.tight_layout()
        fn = os.path.join(rundir, f'balzano_{kind}.png')
        fig.savefig(fn, dpi=130)
        print('  wrote', fn)


if __name__ == '__main__':
    if len(sys.argv) < 2:
        print(__doc__)
        sys.exit(1)
    run(sys.argv[1], plot='--plot' in sys.argv)
