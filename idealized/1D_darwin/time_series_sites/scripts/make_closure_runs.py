#!/usr/bin/env python3
"""
Grazer + diazotroph closure at any site (2026-10-14; master watch D105; gf_results/zoo/README.md).

  python3 scripts/make_closure_runs.py <dst_parent> <src_run_dir>:<name> [...]
  python3 scripts/make_closure_runs.py --inplace <run_dir> [...]      (edits existing control dirs, e.g. the GF ctrl)

Copies each source run (outputs/pickups left out) to <dst_parent>/<name> (or edits <run_dir> in place) and applies,
in one go, the three edits tested at BATS + HOT:
  q01 (make_zoo_runs.py)    MORT(6:7) = 0, MORT2(6:7) = MORT_lin(6:7) / 0.1 mmol C m-3 (from the source's own MORT)
  d2  (make_diazo_runs.py)  type 2 -> diazotroph: DIAZO 1, USENH4/NO2/NO3 0 (INTEGER traits); c02 IC >= 1e-3 mmol C m-3,
                            Chl02 IC >= c02 x (type-5 Chl:C)
  s2  (make_diazo_slow.py)  PCMAX(2) = 4.70e-6 s-1 (0.41 d-1, v06 diazotroph)
"""
import os, sys, shutil
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from make_zoo_runs import read_arr, set_arr

J, ZOO, ZREF, SEED, PCMAX_DIAZO = 1, [5, 6], 0.1, 1e-3, 4.70e-6


def set_int(txt, name, vals):
    m, _ = read_arr(txt, name, len(vals))
    return txt[:m.start()] + ' %s = %s,\n' % (name, ', '.join('%d' % x for x in vals)) + txt[m.end():]


def apply(d, src):
    """Apply q01 + d2 + s2 to run dir d (namelists + c02/Chl02 ICs). Returns the log line."""
    p = os.path.join(d, 'data.traits'); txt = open(p).read(); log = []
    _, dz = read_arr(txt, 'DIAZO', 7)
    assert not any(dz), d + ': already has a diazotroph (closure applied twice?)'
    _, mort = read_arr(txt, 'MORT', 7); _, mort2 = read_arr(txt, 'MORT2', 7)
    assert mort[ZOO[0]] > 0 and not any(mort2[j] for j in ZOO), d + ': zoo MORT/MORT2 not in the v05 form'
    for j in ZOO:
        mort2[j] = mort[j] / ZREF; mort[j] = 0.
    txt = set_arr(set_arr(txt, 'MORT', mort), 'MORT2', mort2)
    log.append('q01: MORT(6:7)=0, MORT2(6:7)=%s' % ', '.join('%.3e' % mort2[j] for j in ZOO))
    for nm, val in [('DIAZO', 1), ('USENH4', 0), ('USENO2', 0), ('USENO3', 0)]:
        _, a = read_arr(txt, nm, 7); a = [int(x) for x in a]; a[J] = val; txt = set_int(txt, nm, a)
    _, a = read_arr(txt, 'PCMAX', 7); log.append('s2: PCMAX(2) %.3e -> %.3e' % (a[J], PCMAX_DIAZO))
    a[J] = PCMAX_DIAZO; txt = set_arr(txt, 'PCMAX', a)
    if os.path.islink(p):
        os.remove(p)
    open(p, 'w').write(txt)
    fc, fchl, fref, fcref = [os.path.join(d, 'ptracer%02d_ini' % n) for n in (20 + J, 27 + J, 31, 24)]
    c = np.fromfile(fc, '>f4'); chl5 = np.fromfile(fref, '>f4'); c5 = np.fromfile(fcref, '>f4')
    assert c.size == 50
    c0 = c[0]; c = np.maximum(c, SEED)
    ratio = np.where(c5 > 0, chl5 / np.maximum(c5, 1e-30), np.nan)
    ratio = np.where(np.isfinite(ratio), ratio, np.nanmean(ratio))
    chl = np.maximum(np.fromfile(fchl, '>f4'), c * ratio)
    for f in (fc, fchl):
        os.remove(f)                                   # break any link to the source before writing
    c.astype('>f4').tofile(fc); chl.astype('>f4').tofile(fchl)
    log.append('d2: type 2 diazotroph, c02 IC surface %.3g -> %.3g' % (c0, c[0]))
    with open(os.path.join(d, 'experiment.txt'), 'a') as f:
        f.write('closure q01+d2+s2 (scripts/make_closure_runs.py, from %s): %s\n' % (src, '; '.join(log)))
    return '; '.join(log)


def main():
    if sys.argv[1] == '--inplace':
        for d in sys.argv[2:]:
            print('%-40s %s' % (d, apply(d, d + ' (in place)')))
        return
    dst_p = sys.argv[1]
    for spec in sys.argv[2:]:
        src, name = spec.rsplit(':', 1)
        d = os.path.join(dst_p, name)
        assert not os.path.exists(d), d + ' exists'
        shutil.copytree(src, d, symlinks=True, ignore=shutil.ignore_patterns(
            'mnc_out', 'output.txt', 'STD*', '*.pid', 'pickup*', '*.log', 'darwin_*.txt', 'scratch*'))
        os.makedirs(os.path.join(d, 'mnc_out'))
        print('%-28s %s' % (d, apply(d, src)))


if __name__ == '__main__':
    main()
