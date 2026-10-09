#!/usr/bin/env python3
"""Summarize global / hold-out / per-site GF solves:  python3 gf_summary.py gf_results"""
import os, sys, glob
import pandas as pd

root = sys.argv[1]
runs = ['global'] + sorted(d for d in os.listdir(root) if d.startswith(('holdout_', 'only_')))
fac, cost = {}, {}
for r in runs:
    p = os.path.join(root, r, 'out')
    if not os.path.exists(os.path.join(p, 'gf_parameters.csv')):
        continue
    t = pd.read_csv(os.path.join(p, 'gf_parameters.csv'))
    fac[r] = t.set_index('control').apply(lambda x: '%.2f±%.2f' % (x['factor'], abs(x['delta']) * x['eta_std']), axis=1)
    c = pd.read_csv(os.path.join(p, 'cost_by_site_variable.csv'))
    cost[r] = c.groupby('site')[['J0', 'J1']].mean().apply(lambda x: '%.2f→%.2f' % (x['J0'], x['J1']), axis=1)
pd.set_option('display.width', 250); pd.set_option('display.max_columns', 30)
print('Optimized multiplicative factor on each control (±1 sd), by fit:\n')
print(pd.DataFrame(fac).to_string())
print('\nMean normalized misfit per site, control → optimized (linear prediction), by fit:\n')
print(pd.DataFrame(cost).to_string())
for r in runs:
    lg = os.path.join(root, r + '.log')
    if os.path.exists(lg):
        for l in open(lg):
            if l.startswith(('mean normalized cost', 'out-of-sample')):
                print('%-15s %s' % (r, l.strip()))
