#!/usr/bin/env python3
"""
Build rbcs relaxation files from extracted v5r1 monthly THETA / SALTanom (+35)
and, for deep BGC relaxation, v5r1 monthly ptracers (relax_ptrNN_1x1x50_mon).
Same logic as the last step of extract_v5r1_sites.py:
  python3 make_relax.py <site_dir> [<site_dir> ...]
Records: Jan1992 (pad), Jan1992 ... Dec2025, Dec2025 (pad) = 410; below-bottom
levels filled with the deepest wet value.
"""
import os, sys
import numpy as np

NR = 50
# v5r1 monthly diagnostic name -> ptracer index (data.ptracers order)
BGC = {'DIC': 1, 'NO3': 2, 'PO4': 5, 'FeT': 6, 'SiO2': 7, 'ALK': 18, 'O2': 19}
for d in sys.argv[1:]:
    nwet = int([l for l in open(os.path.join(d, 'site_info.txt')) if l.startswith('nwet')][0].split('=')[1])
    m = os.path.join(d, 'v5r1_monthly')
    T = np.fromfile(os.path.join(m, 'THETA.bin'), '>f4').reshape(-1, NR).astype('f8')
    S = np.fromfile(os.path.join(m, 'SALTanom.bin'), '>f4').reshape(-1, NR) + 35.
    for nm, a in [('T', T), ('S', S)]:
        a = a.copy()
        a[:, nwet:] = a[:, nwet - 1:nwet]
        a = np.concatenate([a[:1], a, a[-1:]], axis=0)
        a.astype('>f4').tofile(os.path.join(d, 'relax_%s_1x1x50_mon' % nm))
        print(d, nm, a.shape, 'surface range %.2f..%.2f' % (a[:, 0].min(), a[:, 0].max()))
    for nm, it in BGC.items():
        a = np.fromfile(os.path.join(m, nm + '.bin'), '>f4').reshape(-1, NR).astype('f8')
        a[:, nwet:] = a[:, nwet - 1:nwet]
        a = np.concatenate([a[:1], a, a[-1:]], axis=0)
        a.astype('>f4').tofile(os.path.join(d, 'relax_ptr%02d_1x1x50_mon' % it))
        print(d, nm, 'ptr%02d' % it, a.shape, 'k30 mean %.4g' % a[:, 29].mean())
