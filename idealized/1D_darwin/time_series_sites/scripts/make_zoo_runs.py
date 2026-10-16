#!/usr/bin/env python3
"""
Zooplankton-closure test runs (2026-10-13; zooplankton c06/c07 die out in the 1-D columns, gf_results/zoo/README.md).

  python3 scripts/make_zoo_runs.py <src_parent> <dst_parent> <SITE> [SITE ...]

Copies <src_parent>/gyre_<SITE> (outputs left out; namelists, inputs, executable as in the source, like
make_nonneg_runs.py) to <dst_parent>/<variant>_<SITE> and edits only the zooplankton (types 6-7) entries of
data.traits:
  q1     MORT = 0, MORT2 = MORT_lin / 1.0 mmol C m-3   quadratic loss = today's linear loss at Z = 1 mmol C m-3
  q01    MORT = 0, MORT2 = MORT_lin / 0.1 mmol C m-3   ... at Z = 0.1 (obs mesozooplankton 0-150 m ~0.05-0.1)
  k3     KGRAZESAT 10.2 -> 3 mmol C m-3 (linear loss kept): break-even prey 1.22 -> 0.36
  q01k3  q01 + k3
TEMPMORT2 = 0 and EXPORTFRACMORT2 (0.7 for types 6-7) stay as in data.traits.
"""
import os, re, sys, shutil

ZOO = [5, 6]                                     # 0-based type indices of c06, c07
VARIANTS = {'q1': dict(zref=1.0), 'q01': dict(zref=0.1), 'k3': dict(ksat=3.0), 'q01k3': dict(zref=0.1, ksat=3.0)}


def read_arr(txt, name, n):
    m = re.search(r'^(\s*%s\s*=)(.*?)(?=^\s*[A-Z_0-9]+\s*=|^\s*&|^\s*/|\Z)' % name, txt, re.M | re.S)
    vals = []
    for tok in m.group(2).replace('\n', ' ').split(','):
        tok = tok.strip()
        if not tok:
            continue
        if '*' in tok:
            k, v = tok.split('*'); vals += [float(v)] * int(k)
        else:
            vals.append(float(tok))
    assert len(vals) == n, (name, len(vals))
    return m, vals


def set_arr(txt, name, vals):
    m, _ = read_arr(txt, name, len(vals))
    return txt[:m.start()] + ' %s = %s,\n' % (name, ', '.join('%.15E' % v for v in vals)) + txt[m.end():]


def main():
    src_p, dst_p = sys.argv[1:3]
    for s in sys.argv[3:]:
        src = os.path.join(src_p, 'gyre_' + s)
        for v, opt in VARIANTS.items():
            d = os.path.join(dst_p, '%s_%s' % (v, s))
            assert not os.path.exists(d), d + ' exists'
            shutil.copytree(src, d, symlinks=True, ignore=shutil.ignore_patterns(
                'mnc_out', 'output.txt', 'STD*', '*.pid', 'pickup*', '*.log', 'run.log'))
            os.makedirs(os.path.join(d, 'mnc_out'))
            p = os.path.join(d, 'data.traits')
            txt = open(p).read()
            _, mort = read_arr(txt, 'MORT', 7)
            _, mort2 = read_arr(txt, 'MORT2', 7)
            _, ksat = read_arr(txt, 'KGRAZESAT', 7)
            log = []
            if 'zref' in opt:
                for j in ZOO:
                    mort2[j] = mort[j] / opt['zref']; mort[j] = 0.
                txt = set_arr(set_arr(txt, 'MORT', mort), 'MORT2', mort2)
                log.append('MORT(6:7) = 0, MORT2(6:7) = %s (s-1 per mmol C m-3; Zref %.1f)'
                           % (', '.join('%.4E' % mort2[j] for j in ZOO), opt['zref']))
            if 'ksat' in opt:
                for j in ZOO:
                    ksat[j] = opt['ksat']
                txt = set_arr(txt, 'KGRAZESAT', ksat)
                log.append('KGRAZESAT(6:7) = %.1f' % opt['ksat'])
            open(p, 'w').write(txt)
            with open(os.path.join(d, 'experiment.txt'), 'a') as f:
                f.write('zoo closure test %s (scripts/make_zoo_runs.py, from %s): %s\n' % (v, src, '; '.join(log)))
            print('%-24s %s' % (d, '; '.join(log)))


if __name__ == '__main__':
    main()
