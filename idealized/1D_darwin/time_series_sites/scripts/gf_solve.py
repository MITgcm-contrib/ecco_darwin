#!/usr/bin/env python3
"""
Green's-function parameter estimation for the 1-D ECCO-Darwin columns
(after Menemenlis et al. 2005, MWR; Carroll et al. 2020, JAMES, doi:10.1029/2019MS001888).

  python3 gf_solve.py <ctrl_runs> <gf_runs> <obs_dir> <out_dir> [--holdout SITE | --only SITE] [--prior 1.0]

  --only SITE fits the parameters to that site's observations alone (per-site optimization,
  same Green's-function runs); the other sites are then scored as out-of-sample.

Data vector: observations binned by (site, variable, depth layer, year-month);
model equivalents are the same bins of the model's daily (2-D) or 5-day (3-D)
means, sampled on the observation days.  For each control i (perturbation
factor 1+delta_i, gf_controls.py) the Green's function is
    G[:, i] = (m_i - m_0)        (unit eta_i == the tested perturbation)
and the estimate is
    eta = (P^-1 + G' R^-1 G)^-1 G' R^-1 (d - m_0),   P = prior^2 I,
with R diagonal: per (site, variable) variance of the binned observations
(floor: 10% of the mean) / n_obs-in-bin + model-representation error.
Optimized parameter values: p = p0 * (1 + delta_i * eta_i).

Outputs: eta, posterior std, cost before/after per (site, variable),
and the parameter table, in <out_dir>.
"""
import os, sys, glob, argparse, datetime as dt
import numpy as np
import pandas as pd
import netCDF4

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from gf_controls import CONTROLS, RUNDIR, KZ_CONTROL
from obs_utils import load_obs

T0 = dt.datetime(1992, 1, 1)
DRF = np.array([
    10.00, 10.00, 10.00, 10.00, 10.00, 10.00, 10.00, 10.01,
    10.03, 10.11, 10.32, 10.80, 11.76, 13.42, 16.04, 19.82, 24.85,
    31.10, 38.42, 46.50, 55.00, 63.50, 71.58, 78.90, 85.15, 90.18,
    93.96, 96.58, 98.25, 99.25, 100.01, 101.33, 104.56, 111.33, 122.83,
    139.09, 158.94, 180.83, 203.55, 226.50, 249.50, 272.50, 295.50, 318.50,
    341.50, 364.50, 387.50, 410.50, 433.50, 456.50])
ZC = np.cumsum(DRF) - DRF / 2.
LAYERS = [0, 50, 100, 200, 500, 1000, 2000, 6000]
RHO = 1.025          # kg/L, umol/kg -> mmol/m3 for bottle data reported per kg

# obs variable -> (model diagnostic, list?, scale applied to MODEL to obs units)
# Model ptracers are mmol/m3 (Chl mg/m3); PP diag is mmol C m-3 s-1; pCO2 atm.
VARMAP = {
    'NO3':  ('TRAC02', 1.), 'PO4': ('TRAC05', 1.), 'SiO2': ('TRAC07', 1.),
    'DIC':  ('TRAC01', 1.), 'ALK': ('TRAC18', 1.), 'O2':   ('TRAC19', 1.),
    'Chl':  ('CHL', 1.),                                     # sum TRAC27..31
    'PP':   ('PP', 12.011 * 86400.),                         # -> mg C m-3 d-1
    # obs POC is all particles: model = detrital POC + plankton C (c01-c07). DOC/PON not used:
    # Darwin DOC is the labile pool only (no ~40 uM refractory DOC), PON adds little beyond POC.
    'POC':  ('POCTOT', 1.),
    'pCO2': ('pCO2', 1e6),                                   # -> uatm
    'POC_flux': ('POCFLUX', 12.011 * 86400.),                # wC_sink*POC -> mg C m-2 d-1
}


class Model:
    """Model output for one run dir: daily 2-D and daily/5-day 3-D NetCDF."""
    def __init__(self, run):
        self.d = {}
        self.wC = 1.1574074074074075E-004                    # darwin default 10 m/d
        for line in open(os.path.join(run, 'darwin_params.txt')):
            if line.strip().upper().startswith('WC_SINK'):
                self.wC = float(line.split('=')[1].strip().rstrip(',').replace('D', 'E'))
        for lst in ['diags3D', 'diags2D']:
            f = sorted(glob.glob(os.path.join(run, 'mnc_out', lst + '.*.nc')))
            if not f:
                raise FileNotFoundError(run)
            ds = netCDF4.Dataset(f[0])
            t = np.array(ds['T'][:])
            freq = np.median(np.diff(t)) if len(t) > 1 else 86400.
            self.d[lst] = (ds, t, freq)

    def series(self, name, avg_to=None):
        """Return (t_center_days, array[ntime, nlev]) for diag `name`."""
        for lst, (ds, t, freq) in self.d.items():
            if name == 'POCTOT' and 'TRAC12' in ds.variables:
                a = sum(np.array(ds['TRAC%02d' % n][:]) for n in [12] + list(range(20, 27)))
            elif name == 'POCFLUX' and 'TRAC12' in ds.variables:
                a = np.array(ds['TRAC12'][:]) * self.wC         # mmol C m-2 s-1
            elif name == 'CHL' and 'TRAC27' in ds.variables:
                a = sum(np.array(ds['TRAC%02d' % n][:]) for n in range(27, 32))
            elif name in ds.variables:
                a = np.array(ds[name][:])
            else:
                continue
            a = a.reshape(len(t), -1)
            tc = (t - freq / 2.) / 86400.
            if avg_to and avg_to > freq:          # block-average daily -> avg_to
                n = int(round(avg_to / freq))
                m = (len(t) // n) * n
                a = a[:m].reshape(-1, n, a.shape[1]).mean(1)
                tc = tc[:m].reshape(-1, n).mean(1)
            return tc, a
        raise KeyError(name)


def bin_obs(site, obs_root):
    o = load_obs(site, obs_root, variables=list(VARMAP))
    # POC flux: only shallow traps (100-500 m), where wC_sink x POC is a fair equivalent
    o = o[(o['variable'] != 'POC_flux') | o['depth_m'].between(100., 500.)]
    o['tday'] = (o['t'] - pd.Timestamp(T0)).dt.total_seconds() / 86400.
    o['layer'] = np.digitize(o['depth_m'], LAYERS) - 1
    o['ym'] = o['t'].dt.year * 100 + o['t'].dt.month
    g = o.groupby(['variable', 'layer', 'ym'])
    b = g.agg(value=('value_m', 'mean'), n=('value_m', 'size'),
              depth=('depth_m', 'mean')).reset_index()
    b['site'] = site
    b['days'] = g['tday'].apply(list).values
    b['depths'] = g['depth_m'].apply(list).values
    return b


def equivalents(model, bins, avg3d):
    cache, out = {}, np.zeros(len(bins))
    for j, r in enumerate(bins.itertuples()):
        name, scale = VARMAP[r.variable]
        if name not in cache:
            cache[name] = model.series(name, avg_to=avg3d)
        tc, a = cache[name]
        vals = []
        for day, z in zip(r.days, r.depths):
            it = int(np.argmin(np.abs(tc - day)))
            prof = a[it]
            vals.append(prof[0] if prof.size == 1 else np.interp(z, ZC, prof))
        out[j] = np.mean(vals) * scale
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('ctrl'); ap.add_argument('gf'); ap.add_argument('obs'); ap.add_argument('out')
    ap.add_argument('--holdout', default=None)
    ap.add_argument('--cache', default='gf_results/cache', help='model-equivalent cache (shared by all solves)')
    ap.add_argument('--recompute', action='store_true', help='ignore cached model equivalents')
    ap.add_argument('--only', default=None, help='fit to this site only, or a comma list of sites (regional fit)')
    ap.add_argument('--prior', type=float, default=1.0)
    ap.add_argument('--min_bins', type=int, default=20, help='min bins per (site, variable) pair')
    ap.add_argument('--min_layer_bins', type=int, default=1,
                    help='min bins per (site, variable, layer) cell: its error variance comes from those bins')
    ap.add_argument('--avg3d', type=float, default=432000.)
    ap.add_argument('--kz', action='store_true', help='add the KZ_BG control (needs --cache from gf_kz_cache.py)')
    A = ap.parse_args()
    os.makedirs(A.out, exist_ok=True)
    sites = sorted(d for d in os.listdir(A.ctrl) if os.path.isfile(os.path.join(A.ctrl, d, 'data')))
    ctrls = list(CONTROLS) + ([(KZ_CONTROL[0], KZ_CONTROL[1], [], KZ_CONTROL[3])] if A.kz else [])
    names = [c[0] for c in ctrls]

    D, M0, Gc, rows = [], [], [], []
    cache_dir = A.cache
    os.makedirs(cache_dir, exist_ok=True)
    for s in sites:
        cf = os.path.join(cache_dir, s)
        if os.path.exists(cf + '.npz') and not A.recompute:
            z = np.load(cf + '.npz'); b = pd.read_pickle(cf + '.pkl')
            m0, Gs = z['m0'], z['G']
            if Gs.shape[1] != len(names):
                sys.exit('%s: cache has %d controls, solve wants %d (use --kz w/ gf_results/cache_kz)' % (cf, Gs.shape[1], len(names)))
        else:
            b = bin_obs(s, A.obs)
            m0 = equivalents(Model(os.path.join(A.ctrl, s)), b, A.avg3d)
            Gs = np.array([equivalents(Model(os.path.join(A.gf, RUNDIR.get(n, n), s)), b, A.avg3d) - m0
                           for n in names]).T
            np.savez(cf + '.npz', m0=m0, G=Gs)
            b[['site', 'variable', 'layer', 'ym', 'n', 'value']].to_pickle(cf + '.pkl')
        D.append(b['value'].values); M0.append(m0); Gc.append(Gs)
        rows.append(b[['site', 'variable', 'layer', 'ym', 'n', 'value']])
        print('%s: %d bins' % (s, len(b)), flush=True)
    meta = pd.concat(rows, ignore_index=True)
    d, m0, G = np.concatenate(D), np.concatenate(M0), np.vstack(Gc)

    # observation-error variance per (site, variable, layer)
    sd = meta.groupby(['site', 'variable', 'layer'])['value'].transform('std').fillna(0.).values
    floor = 0.1 * np.abs(meta.groupby(['site', 'variable', 'layer'])['value'].transform('mean').values)
    r = np.maximum(sd, floor) ** 2
    ok = np.isfinite(d) & np.isfinite(m0) & np.all(np.isfinite(G), 1) & (r > 0)
    # require >= MIN_BINS bins per (site, variable) so a few samples cannot steer the fit
    nb = meta.groupby(['site', 'variable'])['value'].transform('size').values
    ok &= nb >= A.min_bins
    nl = meta.groupby(['site', 'variable', 'layer'])['value'].transform('size').values
    ok &= nl >= A.min_layer_bins
    train = ok & (meta['site'].values != A.holdout)
    if A.only:
        train = ok & np.isin(meta['site'].values, A.only.split(','))   # one site or a comma list (regional fit)
    # equal weight per (site, variable): each pair contributes its mean normalized misfit
    ngrp = meta.groupby(['site', 'variable'])['value'].transform('size').values
    W = 1. / r / ngrp
    GtW = G[train].T * W[train]
    Pinv = np.eye(len(names)) / A.prior ** 2
    H = Pinv + GtW @ G[train]
    eta = np.linalg.solve(H, GtW @ (d - m0)[train])
    post = np.sqrt(np.diag(np.linalg.inv(H)))
    m1 = m0 + G @ eta                      # linear prediction of optimized run

    meta['J0'] = np.where(ok, (d - m0) ** 2 / r, np.nan)      # normalized misfit per bin (report)
    meta['J1'] = np.where(ok, (d - m1) ** 2 / r, np.nan)
    cost = meta.groupby(['site', 'variable'])[['J0', 'J1']].mean()
    cost.to_csv(os.path.join(A.out, 'cost_by_site_variable.csv'))
    tab = pd.DataFrame({'control': names, 'delta': [c[1] for c in ctrls],
                        'eta': eta, 'eta_std': post,
                        'factor': [1. + c[1] * e for c, e in zip(ctrls, eta)],
                        'description': [c[3] for c in ctrls]})
    tab.to_csv(os.path.join(A.out, 'gf_parameters.csv'), index=False)
    np.savez(os.path.join(A.out, 'gf_system.npz'), d=d, m0=m0, G=G, r=r, ok=ok, train=train)
    print(tab.to_string(index=False))
    def pairmean(mask, col):
        return meta[mask].groupby(['site', 'variable'])[col].mean().mean()
    print('\nmean normalized cost over site-variable pairs (train%s): %.3f -> %.3f (linear prediction)' % (
        ', holdout ' + A.holdout if A.holdout else (', fit to ' + A.only + ' only' if A.only else ''),
        pairmean(train, 'J0'), pairmean(train, 'J1')))
    if A.holdout or A.only:
        h = ok & ~train
        print('out-of-sample %s: %.3f -> %.3f' % (A.holdout or 'other sites', pairmean(h, 'J0'), pairmean(h, 'J1')))
    print(cost.round(2).to_string())


if __name__ == '__main__':
    main()
