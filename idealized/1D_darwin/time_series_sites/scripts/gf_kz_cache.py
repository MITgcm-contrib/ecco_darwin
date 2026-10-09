#!/usr/bin/env python3
"""
Append the KZ_BG sensitivity column to the Green's-function cache.

  python3 scripts/gf_kz_cache.py [cache_in] [cache_out] [bins]   (default gf_results/cache -> gf_results/cache_kz,
                                                                    bins gf_results/bins)

For each site: g = equiv(runs/gf_kz/KZ_BG/<SITE>) - equiv(runs/baseline/<SITE>), on the cached bins;
writes <cache_out>/<SITE>.npz with the original m0 (Pleiades control) and G = [G_19 | g], and copies
<SITE>.pkl. Both runs are Mac gfortran builds, so compiler differences cancel in g.
Then: python3 scripts/gf_solve.py <ctrl> <gf> obs <out> --cache gf_results/cache_kz --kz [--only SITE]
"""
import os, shutil, sys
import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from gf_solve import Model, equivalents
from gf_controls import KZ_CONTROL
from gf_bins import load_bins

SITES = ['HOT', 'BATS', 'HydroS', 'PAPA', 'PAP']

if __name__ == '__main__':
    cin = sys.argv[1] if len(sys.argv) > 1 else 'gf_results/cache'
    cout = sys.argv[2] if len(sys.argv) > 2 else 'gf_results/cache_kz'
    bdir = sys.argv[3] if len(sys.argv) > 3 else 'gf_results/bins'
    os.makedirs(cout, exist_ok=True)
    for s in SITES:
        z = np.load(os.path.join(cin, s + '.npz'))
        b = load_bins(bdir, s)              # full bins (times/depths); same order as <SITE>.pkl
        assert len(b) == len(pd.read_pickle(os.path.join(cin, s + '.pkl'))) == len(z['m0'])
        mb = equivalents(Model('runs/baseline/%s' % s), b, 432000.)
        mk = equivalents(Model('runs/gf_kz/%s/%s' % (KZ_CONTROL[0], s)), b, 432000.)
        g = mk - mb
        np.savez(os.path.join(cout, s + '.npz'), m0=z['m0'], G=np.column_stack([z['G'], g]))
        shutil.copy(os.path.join(cin, s + '.pkl'), os.path.join(cout, s + '.pkl'))
        v = b['variable'].values
        print('%-6s %d bins; mean |g| by variable: %s' % (s, len(b), ', '.join(
            '%s %.3g' % (k, np.nanmean(np.abs(g[v == k]))) for k in sorted(set(v)))), flush=True)
