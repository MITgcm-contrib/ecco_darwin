#!/usr/bin/env python3
"""
1-D column vs ECCO-Darwin v5r1 vs observations at one site, Nature style (nature-figures skill).

  python3 compare_obs.py <run_dir> <site_dir> <obs_root> <out_path_without_ext> [--label 1-D]

a-f: 1992-2025 mean profiles (0-1000 m; Chl and PP 0-250 m) of NO3, DIC, ALK, O2, Chl, PP,
     with observations binned by model level (median; 25-75 % band).
g-j: monthly surface (0-20 m) DIC, pCO2 (0-10 m obs), total Chl, and POC flux at ~150 m
     (model = wC_sink x POC; obs 100-200 m).
1-D model vermillion, v5r1 gray dashed, observations blue.  Writes <out>.pdf + <out>.png.
"""
import os, sys, argparse, datetime as dt
import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))   # nature_style.py ships alongside
from nature_style import *          # noqa: F401,F403
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from obs_utils import load_obs
from gf_solve import Model, ZC, DRF, T0

ZF = np.concatenate([[0.], np.cumsum(DRF)])
LABEL = {'HOT': 'HOT (ALOHA)', 'BATS': 'BATS', 'HydroS': 'Hydrostation S', 'PAPA': 'OWS Papa', 'PAP': 'PAP-SO'}
PROF = [('NO3', 'TRAC02', 1., 'NO3', 'NO$_3$ (mmol m$^{-3}$)'), ('DIC', 'TRAC01', 1., 'DIC', 'DIC (mmol m$^{-3}$)'),
        ('ALK', 'TRAC18', 1., 'ALK', 'ALK (mmol m$^{-3}$)'), ('O2', 'TRAC19', 1., 'O2', 'O$_2$ (mmol m$^{-3}$)'),
        ('Chl', 'CHL', 1., None, 'Chl (mg m$^{-3}$)'),
        ('PP', 'PP', 12.011 * 86400., 'primProd', 'PP (mg C m$^{-3}$ d$^{-1}$)')]


def v5r1(site_dir, name, freq='monthly'):
    d = os.path.join(site_dir, 'v5r1_' + freq)
    f = os.path.join(d, name + '.bin')
    if not os.path.exists(f):
        return None, None
    steps = np.loadtxt(os.path.join(d, 'steps.txt'), dtype=int)
    a = np.fromfile(f, '>f4').reshape(len(steps), -1)
    per = 86400. if freq == 'daily' else 30.4 * 86400.
    t = np.array([T0 + dt.timedelta(seconds=s * 1200. - per / 2.) for s in steps])
    return t, a


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('run'); ap.add_argument('site_dir'); ap.add_argument('obs'); ap.add_argument('out')
    ap.add_argument('--label', default='1-D')
    A = ap.parse_args()
    site = os.path.basename(os.path.normpath(A.run))
    M = Model(A.run)
    obs = load_obs(site, A.obs)
    obs['k'] = np.clip(np.searchsorted(ZF, obs['depth_m'].values, side='right') - 1, 0, 49)
    kmax = int(np.searchsorted(ZC, 1000.))

    fig = plt.figure(figsize=size(DOUBLE, 130), constrained_layout=True)
    gs = fig.add_gridspec(3, 6, height_ratios=[1.35, 1, 1])
    letters = iter('abcdefghij')
    for j, (ov, mv, sc, vv, lab) in enumerate(PROF):
        ax = fig.add_subplot(gs[0, j])
        o = obs[obs['variable'] == ov]
        if len(o):
            q = o[o['k'] < kmax].groupby('k')['value_m'].quantile([.25, .5, .75]).unstack()
            ax.fill_betweenx(ZC[q.index], q[0.25], q[0.75], color=OBS, alpha=0.18, lw=0)
            ax.plot(q[0.5], ZC[q.index], 'o', color=OBS, ms=1.8, label='Observations')
        if vv:
            tv, b = v5r1(A.site_dir, vv)
            if b is not None:
                sc5 = 12.011 * 86400. if ov == 'PP' else 1.
                ax.plot(b[:, :kmax].mean(0) * sc5, ZC[:kmax], '--', color=GRAY, lw=1.0, label='v5r1')
        try:
            tc, a = M.series(mv)
            ax.plot(a[:, :kmax].mean(0) * sc, ZC[:kmax], color=MOD, lw=1.0, label=A.label)
        except KeyError:
            pass
        ax.set_ylim(1000 if ov not in ('Chl', 'PP') else 250, 0)
        ax.set_xlabel(lab)
        if j == 0:
            ax.set_ylabel('Depth (m)')
        panel(ax, next(letters))

    def ts(ax, ov, mv, sc, vv, lab, kz=(0, 20), v5sc=1., v5freq='monthly', obsdepth=(0, 20)):
        o = obs[(obs['variable'] == ov) & obs['depth_m'].between(*obsdepth)]
        if len(o):
            om = o.set_index('t')['value_m'].resample('MS').mean().dropna()
            ax.plot(om.index, om.values, 'o', color=OBS, ms=1.5, label='Observations')
        if vv:
            tv, b = v5r1(A.site_dir, vv, v5freq)
            if b is not None:
                kk = [0] if (b.shape[1] == 1 or v5freq == 'daily') else (ZC >= kz[0]) & (ZC <= kz[1])
                yv = pd.Series(b[:, kk].mean(1) * v5sc, index=tv).resample('MS').mean()
                ax.plot(yv.index, yv.values, '--', color=GRAY, lw=0.7, label='v5r1')
        tc, a = M.series(mv)
        k = (ZC >= kz[0]) & (ZC <= kz[1]) if a.shape[1] > 1 else [0]
        y = pd.Series(a[:, k].mean(1) * sc, index=[T0 + dt.timedelta(days=float(x)) for x in tc])
        y = y.resample('MS').mean()
        ax.plot(y.index, y.values, color=MOD, lw=0.8, label=A.label)
        ax.set_ylabel(lab); ax.set_xlim(T0, dt.datetime(2026, 1, 1))
        panel(ax, next(letters))

    ax = fig.add_subplot(gs[1, 0:3])
    ts(ax, 'DIC', 'TRAC01', 1., 'DIC', 'Surface DIC (mmol m$^{-3}$)')
    hs, ls = ax.get_legend_handles_labels()
    fig.legend(hs, ls, loc='outside lower center', ncol=3, handlelength=1.8)
    ts(fig.add_subplot(gs[1, 3:6]), 'pCO2', 'pCO2', 1e6, 'pCO2', 'Surface pCO$_2$ (µatm)', v5sc=1e6,
       v5freq='daily', obsdepth=(0, 10))
    ts(fig.add_subplot(gs[2, 0:3]), 'Chl', 'CHL', 1., None, 'Surface Chl (mg m$^{-3}$)')
    ts(fig.add_subplot(gs[2, 3:6]), 'POC_flux', 'POCFLUX', 12.011 * 86400., None,
       'POC flux, ~150 m (mg C m$^{-2}$ d$^{-1}$)', kz=(140, 160), obsdepth=(100, 200))
    save(fig, A.out)
    print('wrote', A.out + '.png')


if __name__ == '__main__':
    main()
