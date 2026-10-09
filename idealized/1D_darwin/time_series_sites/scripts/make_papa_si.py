#!/usr/bin/env python3
"""
OWS Papa surface-silicate experiments (runs/exp), on top of the Papa-only 20-control set (runs/opt_kz/onlyPAPA_PAPA).

  python3 scripts/make_papa_si.py

  PAPA_Si_upwell        SiO2-only upwelling proxy: SiO2 (ptracer 7) relaxed toward v5r1 w/ mask 0 above 50 m -> 1 at
                        100 m and below, tau 90 d; all other BGC tracers keep the standard deep mask (150-250 m) + 1 yr.
                        Needs per-tracer masks: build_mask21 (RBCS_SIZE.h maskLEN = 21; bit-identical to build/ when
                        all masks are equal, 10-day test runs/test_mask21).
  PAPA_RSiC05           diatom Si:C (R_SIC, type 1) x0.5 -> half the silicate uptake per unit diatom growth
  PAPA_Si_upwell_RSiC05 both
"""
import os, sys
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from make_exp import copy_run, traits
from gf_solve import ZC

SRC = 'runs/opt_kz/onlyPAPA_PAPA'
SI = 7                                                   # ptracer index of SiO2 (data.ptracers)


def si_upwell(d, log):
    np.clip((ZC - 50.) / 50., 0, 1).astype('>f4').tofile(os.path.join(d, 'relax_mask_si_1x1x50'))
    p = os.path.join(d, 'data.rbcs')
    t = open(p).read()
    masks = ''.join(" relaxMaskFile(%d)='%s',\n" % (2 + i, 'relax_mask_si_1x1x50' if i == SI else 'relax_mask_bgc_1x1x50')
                    for i in range(2, 20))          # ptracers 2-19 -> masks 4-21 (ptracer 1 keeps mask 3)
    t = t.replace(" relaxMaskFile(3)='relax_mask_bgc_1x1x50',\n", " relaxMaskFile(3)='relax_mask_bgc_1x1x50',\n" + masks)
    t = t.replace('tauRelaxPTR(%d)=31536000.' % SI, 'tauRelaxPTR(%d)=7776000.' % SI)
    os.remove(p); open(p, 'w').write(t)                   # remove first: never write through a symlink
    exe = os.path.join(d, 'mitgcmuv')
    os.remove(exe); os.symlink(os.path.abspath('build_mask21/mitgcmuv'), exe)
    log.append('SiO2-only relaxation: mask 0 above 50 m -> 1 at 100 m, tau 90 d (build_mask21)')


EXPS = [('PAPA_Si_upwell', ['si_upwell']), ('PAPA_RSiC05', ['rsic05']),
        ('PAPA_Si_upwell_RSiC05', ['si_upwell', 'rsic05'])]

if __name__ == '__main__':
    for name, mods in EXPS:
        d = 'runs/exp/' + name
        copy_run(SRC, d)
        log = ['%s: based on %s' % (name, SRC)]
        for m in mods:
            if m == 'si_upwell':
                si_upwell(d, log)
            elif m == 'rsic05':
                traits(d, 'R_SIC', [1], 0.5, log)
        open(os.path.join(d, 'experiment.txt'), 'w').write('\n'.join(log) + '\n')
        print(' | '.join(log))
