#!/usr/bin/env python3
"""
Which controls matter in each single-site Green's-function fit.

  python3 gf_importance.py <gf_results_dir> <cache_dir>   ->  <gf_results_dir>/importance.csv

For site S (fit only_S), with J = mean over the site's variables of the mean normalized misfit
(same cost as gf_solve.py, linear GF model):
  J0      control,  Jopt  optimized,
  dJ_out(i) = J(eta with control i reverted to 0) - Jopt   (> 0: the fit needs control i)
  dJ_in(i)  = J(only control i applied) - J0               (< 0: control i alone helps)
"""
import os, sys
import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from gf_controls import CONTROLS

SITES = ['BATS', 'HOT', 'HydroS', 'PAP', 'PAPA']        # stacking order in gf_solve (sorted)


def main():
    root, cache = sys.argv[1:3]
    names = [c[0] for c in CONTROLS]
    meta = pd.concat([pd.read_pickle(os.path.join(cache, s + '.pkl')) for s in SITES], ignore_index=True)
    rows = []
    for s in ['HOT', 'BATS', 'HydroS', 'PAPA', 'PAP']:
        z = np.load(os.path.join(root, 'only_%s' % s, 'out', 'gf_system.npz'))
        par = pd.read_csv(os.path.join(root, 'only_%s' % s, 'out', 'gf_parameters.csv')).set_index('control')
        eta = par.loc[names, 'eta'].values
        d, m0, G, r, ok = z['d'], z['m0'], z['G'], z['r'], z['ok'] & (meta['site'].values == s)
        var = meta['variable'].values

        def J(e):
            m = m0 + G @ e
            j = pd.Series(((d - m) ** 2 / r)[ok]).groupby(var[ok]).mean()
            return j.mean()
        J0, Jopt = J(np.zeros_like(eta)), J(eta)
        for i, n in enumerate(names):
            e_out = eta.copy(); e_out[i] = 0.
            e_in = np.zeros_like(eta); e_in[i] = eta[i]
            rows.append(dict(site=s, control=n, factor=par.loc[n, 'factor'],
                             sd=abs(par.loc[n, 'delta']) * par.loc[n, 'eta_std'],
                             dJ_out=J(e_out) - Jopt, dJ_in=J(e_in) - J0, J0=J0, Jopt=Jopt))
    t = pd.DataFrame(rows)
    t.to_csv(os.path.join(root, 'importance.csv'), index=False)
    for s, g in t.groupby('site', sort=False):
        g = g.sort_values('dJ_out', ascending=False)
        print('\n%s: J %.2f -> %.2f (linear)' % (s, g.J0.iloc[0], g.Jopt.iloc[0]))
        print(g.head(5)[['control', 'factor', 'sd', 'dJ_out', 'dJ_in']].round(3).to_string(index=False))


if __name__ == '__main__':
    main()
