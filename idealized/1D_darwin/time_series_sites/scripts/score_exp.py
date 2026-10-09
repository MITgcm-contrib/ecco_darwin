#!/usr/bin/env python3
"""
Score follow-up experiments (runs/exp) with the same cost as gf_solve.py / gf_confirm.py.

  python3 scripts/score_exp.py <solve_out_dir> <bins_dir> <cache_dir> <run_dir>:<SITE> [...]

Prints the mean normalized misfit per (site, variable) of each forward run, (d - m)^2 / r, using
r and the ok mask from <solve_out_dir>/gf_system.npz. Writes exp_scores.csv in <solve_out_dir>
(does not touch confirm.csv).
"""
import os, sys
import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from gf_solve import Model, equivalents
from gf_bins import load_bins, SITES


def main():
    sol, bdir, cache = sys.argv[1:4]
    z = np.load(os.path.join(sol, 'gf_system.npz'))
    order = sorted(SITES)
    start, sl = 0, {}
    for s in order:
        n = len(pd.read_pickle(os.path.join(cache, s + '.pkl')))
        sl[s] = slice(start, start + n); start += n
    rows = []
    for spec in sys.argv[4:]:
        run, s = spec.split(':')
        b = load_bins(bdir, s)
        m = equivalents(Model(run), b, 432000.)
        i = sl[s]
        d, r, ok = z['d'][i], z['r'][i], z['ok'][i]
        t = pd.DataFrame({'variable': b['variable'].values, 'J': (d - m) ** 2 / r})[ok].groupby('variable').mean()
        t['run'] = '/'.join(os.path.normpath(run).split('/')[-2:]); t['site'] = s
        rows.append(t.reset_index())
        print('%-24s %-6s mean J %.3f' % (t['run'].iloc[0], s, t['J'].mean()), flush=True)
    out = pd.concat(rows)
    out.to_csv(os.path.join(sol, 'exp_scores.csv'), index=False)
    pd.set_option('display.width', 250)
    print(out.pivot_table(index='variable', columns='run', values='J').round(2).to_string())


if __name__ == '__main__':
    main()
