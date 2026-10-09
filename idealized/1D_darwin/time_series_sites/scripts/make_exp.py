#!/usr/bin/env python3
"""
Targeted follow-up experiments on top of the site-optimized parameter sets (runs/exp/).

  python3 scripts/make_exp.py

Each experiment copies a source run directory (inputs symlinked) and applies one or more
modifications; experiment.txt records what changed.
  ksatfe3   : KSATFET x3 for diatoms + other large eukaryotes (types 1-2)  -> stronger Fe limitation
  upwell    : BGC relaxation mask 0 above 50 m -> 1 at 100 m, tau 90 d    -> stand-in for Ekman upwelling
  picogrow2 : PCMAX x2 for Synechococcus + Prochlorococcus (types 3-5)
  kz3       : background diffusivity x3 at 75-250 m
"""
import os, re, shutil, sys
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from gf_make_runs import scale_vector
from gf_solve import ZC

IGN = shutil.ignore_patterns('mnc_out', 'output.txt', 'STD*', '*.pid', 'pickup*', 'darwin_*.txt', '*.data', '*.meta')


def copy_run(src, dst):
    if os.path.exists(dst):
        shutil.rmtree(dst)
    shutil.copytree(src, dst, symlinks=True, ignore=IGN)
    os.makedirs(os.path.join(dst, 'mnc_out'), exist_ok=True)
    exe = os.path.join(dst, 'mitgcmuv')
    if os.path.lexists(exe):
        os.remove(exe)
    os.symlink(os.path.abspath('build/mitgcmuv'), exe)


def traits(d, var, idx, fac, log):
    p = os.path.join(d, 'data.traits')
    t = open(p).read()                       # read first: open(p, 'w') would truncate it
    open(p, 'w').write(scale_vector(t, var, idx, fac, log))


def upwell(d, log):
    m = os.path.join(d, 'relax_mask_bgc_1x1x50')
    os.remove(m)
    np.clip((ZC - 50.) / 50., 0, 1).astype('>f4').tofile(m)
    p = os.path.join(d, 'data.rbcs')
    t = re.sub(r'tauRelaxPTR\((\d+)\)=31536000\.', r'tauRelaxPTR(\1)=7776000.', open(p).read())
    open(p, 'w').write(t)
    log.append('BGC relaxation mask 0 above 50 m, ramp to 1 at 100 m; tau 90 d (was 150-250 m, 1 yr)')


def kz3(d, log):
    p = os.path.join(d, 'diffkr_1x1x50')
    k = np.fromfile(p, '>f4').copy()
    os.remove(p)
    m = (ZC >= 75) & (ZC <= 250)
    k[m] *= 3
    k.astype('>f4').tofile(p)
    log.append('background diffusivity x3 at 75-250 m')


EXPS = [
    ('PAPA_ksatfe3', 'runs/opt/onlyPAPA_PAPA', ['ksatfe3']),
    ('PAPA_upwell', 'runs/opt/onlyPAPA_PAPA', ['upwell']),
    ('PAPA_ksatfe3_upwell', 'runs/opt/onlyPAPA_PAPA', ['ksatfe3', 'upwell']),
    ('HOT_picogrow2', 'runs/opt/onlyHOT_HOT', ['picogrow2']),
    ('BATS_picogrow2', 'runs/opt/onlyBATS_BATS', ['picogrow2']),
    ('HOT_kz3', 'runs/opt/onlyHOT_HOT', ['kz3']),
    ('BATS_kz3', 'runs/opt/onlyBATS_BATS', ['kz3']),
    ('HOT_kz3_picogrow2', 'runs/opt/onlyHOT_HOT', ['kz3', 'picogrow2']),
]

if __name__ == '__main__':
    os.makedirs('runs/exp', exist_ok=True)
    for name, src, mods in EXPS:
        d = 'runs/exp/' + name
        copy_run(src, d)
        log = ['%s: based on %s' % (name, src)]
        for m in mods:
            if m == 'ksatfe3':
                traits(d, 'KSATFET', [1, 2], 3., log)
            elif m == 'picogrow2':
                traits(d, 'PCMAX', [3, 4, 5], 2., log)
            elif m == 'upwell':
                upwell(d, log)
            elif m == 'kz3':
                kz3(d, log)
        open(os.path.join(d, 'experiment.txt'), 'w').write('\n'.join(log) + '\n')
        print(' | '.join(log))
