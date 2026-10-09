#!/usr/bin/env python3
"""
Create Green's-function perturbation run directories from the control runs.

  python3 gf_make_runs.py <ctrl_runs_dir> <gf_runs_dir> [control ...]

For each control in gf_controls.CONTROLS (or those named) and each site in
<ctrl_runs_dir>, makes <gf_runs_dir>/<control>/<site>/ with:
  - the control's namelists, with data.traits / data.darwin entries scaled by (1+delta)
  - symlinks to the control's input files and executable
  - data.diagnostics: 2-D fields daily, 3-D fields every 5 days (smaller output)
A gf_manifest.txt records every edit (old -> new values).
"""
import os, re, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from gf_controls import CONTROLS
import make_runs

NAMELISTS = ['data', 'data.cal', 'data.darwin', 'data.exf', 'data.gchem', 'data.ggl90', 'data.mnc',
             'data.pkg', 'data.ptracers', 'data.rbcs', 'data.seaice', 'data.traits', 'eedata']
SKIP = set(NAMELISTS) | {'data.diagnostics', 'output.txt', 'mnc_out', 'mitgcmuv'}


def expand(vals):
    out = []
    for tok in [t.strip() for t in vals.split(',') if t.strip()]:
        if '*' in tok and not tok.startswith("'"):
            n, v = tok.split('*', 1)
            out += [v.strip()] * int(n)
        else:
            out.append(tok)
    return out


def fnum(x):
    return float(x.replace('D', 'E').replace('d', 'e'))


def scale_vector(text, var, idx, fac, log):
    """Scale 1-based entries idx of namelist array var (handles n*value)."""
    pat = re.compile(r'(^[ \t]*%s[ \t]*=)(.*?)(?=^[ \t]*[A-Za-z_][A-Za-z_0-9]*[ \t]*=|^[ \t]*/|^[ \t]*&)'
                     % re.escape(var), re.M | re.S | re.I)
    m = pat.search(text)
    if not m:
        raise KeyError('%s not found' % var)
    vals = expand(m.group(2).replace('\n', ' '))
    old = [vals[i - 1] for i in idx]
    for i in idx:
        vals[i - 1] = '%.15E' % (fnum(vals[i - 1]) * fac)
    log.append('%s(%s): %s -> %s' % (var, ','.join(map(str, idx)), old, [vals[i - 1] for i in idx]))
    # wrap: Fortran namelist records must stay short (< ~200 chars)
    new = ' ' + ',\n   '.join(', '.join(vals[i:i + 5]) for i in range(0, len(vals), 5)) + ',\n'
    return text[:m.start(2)] + new + text[m.end(2):]


def scale_scalar(text, var, fac, default, log, group='&DARWIN_PARAMS'):
    pat = re.compile(r'^([ \t]*%s[ \t]*=[ \t]*)([^,\n]+)(,?)' % re.escape(var), re.M | re.I)
    m = pat.search(text)
    if m:
        old = fnum(m.group(2).strip())
        new = old * fac
        text = text[:m.start(2)] + '%.15E' % new + text[m.end(2):]
    else:   # absent: insert after the group header with the darwin default
        old = default
        new = old * fac
        g = text.index(group) + len(group)
        text = text[:g] + '\n %s = %.15E,' % (var, new) + text[g:]
    log.append('%s: %.6E -> %.6E%s' % (var, old, new, '' if m else ' (added; default)'))
    return text


def darwin_defaults(path):
    d = {}
    for line in open(path):
        m = re.match(r'^\s*([A-Za-z_0-9]+)\s*=\s*([-+0-9.EeDd]+)\s*,', line)
        if m:
            d[m.group(1).upper()] = fnum(m.group(2))
    return d


def diag_text():
    s = make_runs.diagnostics(86400.)
    return s.replace('  frequency(1) = 86400.0,', '  frequency(1) = 432000.0,', 1)


def main():
    ctrl, gfdir = sys.argv[1:3]
    want = sys.argv[3:]
    sites = sorted(d for d in os.listdir(ctrl) if os.path.isfile(os.path.join(ctrl, d, 'data')))
    for name, delta, entries, desc in CONTROLS:
        if want and name not in want:
            continue
        fac = 1. + delta
        for site in sites:
            src = os.path.abspath(os.path.join(ctrl, site))
            dst = os.path.join(gfdir, name, site)
            os.makedirs(os.path.join(dst, 'mnc_out'), exist_ok=True)
            dp = os.path.join(src, 'darwin_params.txt')
            if not os.path.exists(dp):     # control not run yet: same darwin build defaults
                dp = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'darwin_params_default.txt')
            defaults = darwin_defaults(dp)
            texts = {f: open(os.path.join(src, f)).read() for f in NAMELISTS}
            log = ['control %s: %s, factor %.3f' % (name, desc, fac)]
            for f, var, idx in entries:
                if idx:
                    texts[f] = scale_vector(texts[f], var, idx, fac, log)
                else:
                    texts[f] = scale_scalar(texts[f], var, fac, defaults[var.upper()], log)
            for f, t in texts.items():
                open(os.path.join(dst, f), 'w').write(t)
            open(os.path.join(dst, 'data.diagnostics'), 'w').write(diag_text())
            open(os.path.join(dst, 'gf_manifest.txt'), 'w').write('\n'.join(log) + '\n')
            for f in os.listdir(src):
                if f in SKIP or f.endswith(('.pid', '.data', '.meta', '.txt', '.nc')) or f.startswith(('STD', 'pickup')):
                    continue
                p = os.path.join(dst, f)
                if os.path.lexists(p):
                    os.remove(p)
                os.symlink(os.path.realpath(os.path.join(src, f)), p)
            p = os.path.join(dst, 'mitgcmuv')
            if os.path.lexists(p):
                os.remove(p)
            os.symlink(os.path.realpath(os.path.join(src, 'mitgcmuv')), p)
        print('%-12s x%.2f  %s' % (name, fac, '; '.join(log[1:])))


if __name__ == '__main__':
    main()
