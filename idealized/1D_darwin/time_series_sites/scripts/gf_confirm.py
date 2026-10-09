#!/usr/bin/env python3
"""
Score confirmation runs (optimized parameter sets) with the same cost as gf_solve.py and compare
with the control and the Green's-function linear prediction.

  python3 gf_confirm.py <solve_out_dir> <bins_dir> <cache_dir> <run_dir>:<SITE> [...]

<solve_out_dir> holds gf_system.npz from the gf_solve.py fit whose parameters the runs used.
Prints mean normalized misfit per (site, variable): control (J0), linear prediction (J1),
actual forward run (Jrun).
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
    eta = pd.read_csv(os.path.join(sol, 'gf_parameters.csv'))['eta'].values
    # site order in the solve = sorted run-directory names; slice the stacked system per site
    order = sorted(SITES)
    sizes = {s: len(pd.read_pickle(os.path.join(cache, s + '.pkl'))) for s in order}
    start, sl = 0, {}
    for s in order:
        sl[s] = slice(start, start + sizes[s]); start += sizes[s]
    rows = []
    for spec in sys.argv[4:]:
        run, s = spec.split(':')
        b = load_bins(bdir, s)
        m = equivalents(Model(run), b, 432000.)
        i = sl[s]
        d, m0, G, r, ok = z['d'][i], z['m0'][i], z['G'][i], z['r'][i], z['ok'][i]
        # the confirmation run may fix some controls at 0 (e.g. KPOM): use its own eta
        man = open(os.path.join(run, 'opt_manifest.txt')).read()
        e = eta.copy()
        if 'fixed at control:' in man and 'none' not in man.split('fixed at control:')[1].splitlines()[0]:
            from gf_controls import CONTROLS
            fixed = man.split('fixed at control:')[1].splitlines()[0]
            for k, c in enumerate(CONTROLS):
                if "'%s'" % c[0] in fixed:
                    e[k] = 0.
        m1 = m0 + G @ e
        df = pd.DataFrame({'variable': b['variable'].values, 'J0': (d - m0) ** 2 / r,
                           'J1': (d - m1) ** 2 / r, 'Jrun': (d - m) ** 2 / r})[ok]
        t = df.groupby('variable').mean()
        t['run'] = os.path.basename(os.path.normpath(run)); t['site'] = s
        rows.append(t.reset_index())
        print('%-16s %-6s mean J: control %.3f | linear %.3f | actual %.3f' % (
            t['run'].iloc[0], s, t['J0'].mean(), t['J1'].mean(), t['Jrun'].mean()), flush=True)
    out = pd.concat(rows)
    out.to_csv(os.path.join(sol, 'confirm.csv'), index=False)
    pd.set_option('display.width', 200)
    print(out.pivot_table(index='variable', columns=['run'], values=['J0', 'J1', 'Jrun']).round(2).to_string())


if __name__ == '__main__':
    main()
