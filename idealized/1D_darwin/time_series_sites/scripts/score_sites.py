#!/usr/bin/env python3
"""
Score forward runs at any site (incl. the new ones without a GF system) with the same cost as gf_solve.py:
bins from bin_obs (monthly x layer x variable), r = max(sd, 0.1|mean|)^2 per (site, variable, layer),
>= 20 bins per (site, variable) and >= 3 bins per layer (gf_solve --min_layer_bins 3), mean (d - m)^2 / r
per variable, then the mean over variables.

  python3 scripts/score_sites.py <out_csv> <run_dir>:<SITE> [...]

MAXDEPTH (below) limits a site's bins, e.g. CARIACO (anoxic below ~250 m) to the top 200 m.
"""
import os, sys
import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import obs_utils
from gf_solve import Model, bin_obs, equivalents, LAYERS

# Extra unit conversions for units that only occur at the new sites (2026-10-09). Applied here, NOT in
# obs_utils, so the existing GF cost is unchanged (POC ug -> mmol/m3 is in obs_utils since 2026-10-09).
_base = obs_utils._to_model_units
def _to_model_units(var, units, v):
    out = _base(var, units, v)
    u = units.fillna('').str.lower()
    if var == 'POC_flux':
        out[u.str.startswith('gc/m^2')] *= 1000.                                    # g -> mg C m-2 d-1
    elif var == 'PP':
        out[u.str.contains('/hr')] *= 12.                                           # CARIACO: hourly x 12 h
    return out
obs_utils._to_model_units = _to_model_units

OBS = 'obs'
MAXDEPTH = {'CARIACO': 200.}
_bins = {}


_load = obs_utils.load_obs
def load_obs(*a, **k):
    """Cost-consistency filters for the new sites: pCO2/fCO2 only in the top 10 m (CARIACO has calculated
    fCO2 at depth, OWSM a 500 m sensor; the model pCO2 is surface only); drop uncalibrated CTD fluorescence Chl."""
    o = _load(*a, **k)
    o = o[(o['variable'] != 'pCO2') | (o['depth_m'] <= 10.)]
    return o[~((o['variable'] == 'Chl') & o['units'].fillna('').str.lower().str.contains('ctd fluorescence'))]
import gf_solve
gf_solve.load_obs = load_obs


def site_bins(s):
    if s not in _bins:
        b = bin_obs(s, OBS)
        if s in MAXDEPTH:
            b = b[b['depth'] <= MAXDEPTH[s]].reset_index(drop=True)
        g = b.groupby(['variable', 'layer'])['value']
        sd = g.transform('std').fillna(0.).values
        r = np.maximum(sd, 0.1 * np.abs(g.transform('mean').values)) ** 2
        ok = (r > 0) & (b.groupby('variable')['value'].transform('size').values >= 20) & \
             (g.transform('size').values >= 3)
        _bins[s] = (b, r, ok)
    return _bins[s]


def main():
    out, rows = sys.argv[1], []
    for spec in sys.argv[2:]:
        run, s = spec.split(':')
        b, r, ok = site_bins(s)
        m = equivalents(Model(run), b, 432000.)
        k = ok & np.isfinite(m) & np.isfinite(b['value'].values)
        t = pd.DataFrame({'variable': b['variable'].values, 'J': (b['value'].values - m) ** 2 / r})[k]
        t = t.groupby('variable')['J'].agg(['mean', 'size']).rename(columns={'mean': 'J', 'size': 'nbins'})
        t['run'] = '/'.join(os.path.normpath(run).split('/')[-2:]); t['site'] = s
        rows.append(t.reset_index())
        print('%-28s %-8s mean J %.3f (%d vars)' % (t['run'].iloc[0], s, t['J'].mean(), len(t)), flush=True)
    a = pd.concat(rows)
    a.to_csv(out, index=False)
    pd.set_option('display.width', 250)
    print(a.pivot_table(index='variable', columns='run', values='J').round(2).to_string())


if __name__ == '__main__':
    main()
