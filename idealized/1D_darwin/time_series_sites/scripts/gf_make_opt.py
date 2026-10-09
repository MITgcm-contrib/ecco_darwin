#!/usr/bin/env python3
"""
Build forward runs with an optimized Green's-function parameter set (all controls at once).

  python3 gf_make_opt.py <gf_parameters.csv> <ctrl_runs> <out_dir> <tag> SITE [SITE ...] [--fix KPOM ...]

Each control's factor (1 + delta * eta, from gf_solve.py) is applied to its namelist entries,
as in gf_make_runs.py. Controls named after --fix stay at their control value (eta = 0).
Creates <out_dir>/<tag>_<SITE>/ linked to the control run's inputs, with daily diagnostics.
"""
import os, sys, shutil
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from gf_controls import CONTROLS, KZ_CONTROL
import numpy as np
from gf_solve import ZC
from gf_make_runs import NAMELISTS, scale_vector, scale_scalar, darwin_defaults


def main():
    args = sys.argv[1:]
    fix = []
    if '--fix' in args:
        i = args.index('--fix'); fix = args[i + 1:]; args = args[:i]
    pcsv, ctrl, out, tag = args[:4]
    sites = args[4:]
    tab = pd.read_csv(pcsv).set_index('control')
    here = os.path.dirname(os.path.abspath(__file__))
    for s in sites:
        src = os.path.abspath(os.path.join(ctrl, s))
        dst = os.path.join(out, '%s_%s' % (tag, s))
        if os.path.exists(dst):
            shutil.rmtree(dst)
        shutil.copytree(src, dst, symlinks=True,
                        ignore=shutil.ignore_patterns('mnc_out', 'output.txt', 'STD*', '*.pid', 'pickup*',
                                                      'darwin_*.txt', '*.data', '*.meta'))
        os.makedirs(os.path.join(dst, 'mnc_out'), exist_ok=True)
        defaults = darwin_defaults(os.path.join(here, 'darwin_params_default.txt'))
        texts = {f: open(os.path.join(dst, f)).read() for f in NAMELISTS}
        log = ['optimized set %s from %s; fixed at control: %s' % (tag, pcsv, fix or 'none')]
        for name, delta, entries, desc in CONTROLS:
            fac = 1. if name in fix else float(tab.loc[name, 'factor'])
            if abs(fac - 1.) < 1e-6:
                continue
            log.append('control %s factor %.4f' % (name, fac))
            for f, var, idx in entries:
                if idx:
                    texts[f] = scale_vector(texts[f], var, idx, fac, log)
                else:
                    texts[f] = scale_scalar(texts[f], var, fac, defaults[var.upper()], log)
        for f, t in texts.items():
            open(os.path.join(dst, f), 'w').write(t)
        kn = KZ_CONTROL[0]                     # optional 20th control: scale diffkr_1x1x50 at 75-250 m
        if kn in tab.index and kn not in fix and abs(float(tab.loc[kn, 'factor']) - 1.) > 1e-6:
            fac = float(tab.loc[kn, 'factor'])
            p = os.path.join(dst, 'diffkr_1x1x50')
            k = np.fromfile(p, '>f4').copy()
            os.remove(p)                       # a symlink into the control inputs: never write through it
            m = (ZC >= KZ_CONTROL[2][0]) & (ZC <= KZ_CONTROL[2][1])
            k[m] *= fac
            k.astype('>f4').tofile(p)
            log.append('control %s factor %.4f' % (kn, fac))
            log.append('  diffkr_1x1x50 x%.4f at %g-%g m' % (fac, KZ_CONTROL[2][0], KZ_CONTROL[2][1]))
        open(os.path.join(dst, 'opt_manifest.txt'), 'w').write('\n'.join(log) + '\n')
        print(dst, '|', '; '.join(l for l in log[1:] if l.startswith('control')))


if __name__ == '__main__':
    main()
