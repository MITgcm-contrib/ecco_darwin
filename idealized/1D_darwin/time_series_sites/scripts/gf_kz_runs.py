#!/usr/bin/env python3
"""
Green's-function perturbation runs for a 20th control, KZ_BG: background vertical diffusivity
x(1 + delta) at 75-250 m (the band of the runs/exp kz3 experiments), one run per site.

  python3 scripts/gf_kz_runs.py [SITE ...]   -> runs/gf_kz/KZ_BG/<SITE>   (default: the 5 original sites)

Each run is a copy of runs/baseline/<SITE> (Mac build, corrected rbcs timing, control parameters),
so the sensitivity is equiv(KZ_BG run) - equiv(baseline run) with the same compiler (gf_kz_cache.py).
"""
import os, sys
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from make_exp import copy_run
from gf_solve import ZC
from gf_controls import KZ_CONTROL

NAME, DELTA, (ZTOP, ZBOT) = KZ_CONTROL[0], KZ_CONTROL[1], KZ_CONTROL[2]
SITES = ['HOT', 'BATS', 'HydroS', 'PAPA', 'PAP']
BASE = os.environ.get('GF_BASE', 'runs/baseline')     # 2026-10-18: overridable like gf_kz_cache.py
KZ = os.environ.get('GF_KZ', 'runs/gf_kz')

if __name__ == '__main__':
    for s in sys.argv[1:] or SITES:
        d = '%s/%s/%s' % (KZ, NAME, s)
        os.makedirs(os.path.dirname(d), exist_ok=True)
        copy_run(os.path.realpath('%s/%s' % (BASE, s)), d)
        p = os.path.join(d, 'diffkr_1x1x50')
        k = np.fromfile(p, '>f4').copy()
        os.remove(p)                                   # may be a symlink into the baseline inputs
        m = (ZC >= ZTOP) & (ZC <= ZBOT)
        k[m] *= 1. + DELTA
        k.astype('>f4').tofile(p)
        for f in os.listdir(d):                        # drop pid/log files copied from the baseline
            if f.endswith('.pid') or f == 'experiment.txt':
                os.remove(os.path.join(d, f))
        open(os.path.join(d, 'experiment.txt'), 'w').write(
            '%s: %s/%s with background diffusivity x%.2f at %g-%g m (%d levels)\n'
            % (NAME, BASE, s, 1. + DELTA, ZTOP, ZBOT, m.sum()))
        print(d, 'levels', m.sum(), 'kz %.2e -> %.2e' % (k[m][0] / (1 + DELTA), k[m][0]))
