#!/usr/bin/env python3
"""
Export Green's-function obs bins for the new sites (2026-10-09) in gf_bins.py format, using the
new-site cost filters of score_sites.py (CARIACO top 200 m; pCO2 top 10 m; no uncalibrated CTD
fluorescence Chl; g/hr unit fixes), so the GF cost matches the forward-run scores.

  python3 scripts/gf_bins_new.py <out_dir> [SITE ...]    (default: the six new sites)
"""
import os, sys
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import score_sites
from gf_bins import COLS

NEW = ['SOTS', 'K2', 'CARIACO', 'DYFAMED', 'MUNIDA', 'OWSM']

if __name__ == '__main__':
    out = sys.argv[1]
    os.makedirs(out, exist_ok=True)
    for s in sys.argv[2:] or NEW:
        b, r, ok = score_sites.site_bins(s)
        lens = np.array([len(x) for x in b['days']])
        np.savez(os.path.join(out, '%s_bins.npz' % s), lens=lens,
                 days=np.concatenate(b['days'].values).astype('f8'),
                 depths=np.concatenate(b['depths'].values).astype('f8'))
        b[COLS].to_csv(os.path.join(out, '%s_meta.csv' % s), index=False)
        print(s, len(b), 'bins', lens.sum(), 'obs', flush=True)
