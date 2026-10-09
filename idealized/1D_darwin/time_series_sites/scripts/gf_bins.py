#!/usr/bin/env python3
"""
Observation bins for the Green's-function solve, in a version-independent format so the
model-equivalent step can run on Pleiades without the raw obs files.

  python3 gf_bins.py export <obs_root> <out_dir>          (Mac)   -> <out_dir>/<SITE>_bins.npz + _meta.csv
  python3 gf_bins.py equiv  <bins_dir> <ctrl> <gf> <cache> [SITE ...]   (PBS)  -> <cache>/<SITE>_mG.npz
  python3 gf_bins.py cache  <bins_dir> <cache>            (Mac)   -> gf_solve cache (<SITE>.npz + .pkl)
"""
import os, sys
import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
SITES = ['HOT', 'BATS', 'HydroS', 'PAPA', 'PAP']
COLS = ['site', 'variable', 'layer', 'ym', 'n', 'value']


def export(obs_root, out):
    from gf_solve import bin_obs
    os.makedirs(out, exist_ok=True)
    for s in SITES:
        b = bin_obs(s, obs_root)
        lens = np.array([len(x) for x in b['days']])
        np.savez(os.path.join(out, '%s_bins.npz' % s), lens=lens,
                 days=np.concatenate(b['days'].values).astype('f8'),
                 depths=np.concatenate(b['depths'].values).astype('f8'))
        b[COLS].to_csv(os.path.join(out, '%s_meta.csv' % s), index=False)
        print(s, len(b), 'bins', lens.sum(), 'obs')


def load_bins(bdir, s):
    meta = pd.read_csv(os.path.join(bdir, '%s_meta.csv' % s))
    z = np.load(os.path.join(bdir, '%s_bins.npz' % s))
    idx = np.concatenate([[0], np.cumsum(z['lens'])])
    meta['days'] = [list(z['days'][idx[i]:idx[i + 1]]) for i in range(len(meta))]
    meta['depths'] = [list(z['depths'][idx[i]:idx[i + 1]]) for i in range(len(meta))]
    return meta


def equiv_site(args):
    bdir, ctrl, gf, cache, s = args
    from gf_solve import Model, equivalents
    from gf_controls import CONTROLS, RUNDIR
    b = load_bins(bdir, s)
    m0 = equivalents(Model(os.path.join(ctrl, s)), b, 432000.)
    G = np.array([equivalents(Model(os.path.join(gf, RUNDIR.get(c[0], c[0]), s)), b, 432000.) - m0
                  for c in CONTROLS]).T
    np.savez(os.path.join(cache, '%s_mG.npz' % s), m0=m0, G=G)
    return s, G.shape


if __name__ == '__main__':
    mode = sys.argv[1]
    if mode == 'export':
        export(sys.argv[2], sys.argv[3])
    elif mode == 'equiv':
        from multiprocessing import Pool
        bdir, ctrl, gf, cache = sys.argv[2:6]
        sites = sys.argv[6:] or SITES
        os.makedirs(cache, exist_ok=True)
        with Pool(len(sites)) as p:
            for s, sh in p.imap_unordered(equiv_site, [(bdir, ctrl, gf, cache, s) for s in sites]):
                print(s, 'G', sh, flush=True)
    elif mode == 'cache':
        bdir, cache = sys.argv[2:4]
        for s in SITES:
            z = np.load(os.path.join(cache, '%s_mG.npz' % s))
            np.savez(os.path.join(cache, s + '.npz'), m0=z['m0'], G=z['G'])
            pd.read_csv(os.path.join(bdir, '%s_meta.csv' % s)).to_pickle(os.path.join(cache, s + '.pkl'))
            print(s, 'cached')
