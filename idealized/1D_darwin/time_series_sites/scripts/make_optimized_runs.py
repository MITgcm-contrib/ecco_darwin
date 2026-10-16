#!/usr/bin/env python3
"""
Build ready-to-run 1-D columns with the current best (optimized) parameters, in one command (2026-10-15).

  python3 make_optimized_runs.py <sites_dir> <input_template_dir> <runs_dir> <mitgcmuv> [SITE ...]

Writes <runs_dir>/<SITE>/ for each site (default: all sites in <sites_dir>); run with ./mitgcmuv > output.txt.
The parameters are applied as the chain that produced them (each step is an existing script, so a set-up built here
is identical to the runs it was scored on):
  1. make_runs.py              v05 set-up from the site inputs (negative IC/relaxation values clipped to zero)
  2. gf_make_opt.py, step 1    regional set from the clipped-IC Green's-function fit (19 controls + mixing):
                               results/optimized/{gyre,highlat}_step1.csv
  3. make_closure_runs.py      grazer + nitrogen-fixer closure: quadratic zooplankton mortality, plankton type 2 ->
                               slow non-quota diazotroph (0.41 d-1)
  4. gf_make_opt.py, step 2    second Green's-function iteration around 2-3 (20 controls incl. MORT2_zoo, PCMAX_Syn,
                               PCMAX_diazo; GF_CONTROLS=closure): results/optimized/{gyre,highlat}_step2.csv
Regions: gyre = HOT, BATS, HydroS; high latitude = PAPA, PAP, SOTS, K2, MUNIDA, OWSM. CARIACO + DYFAMED have no
regional set (structural misfit): they get steps 1 + 3 only (v05 parameters + closure).
Each run dir gets opt_manifest.txt (steps 2 + 4: every namelist change) and experiment.txt (step 3).
"""
import os, sys, shutil, subprocess

HERE = os.path.dirname(os.path.abspath(__file__))
PARAMS = os.environ.get('OPT_PARAMS', os.path.join(HERE, '..', 'results', 'optimized'))
REGION = {'HOT': 'gyre', 'BATS': 'gyre', 'HydroS': 'gyre',
          'PAPA': 'highlat', 'PAP': 'highlat', 'SOTS': 'highlat', 'K2': 'highlat', 'MUNIDA': 'highlat',
          'OWSM': 'highlat'}


def run(cmd, closure=False):
    env = dict(os.environ)
    env.pop('GF_CONTROLS', None)
    if closure:
        env['GF_CONTROLS'] = 'closure'
    subprocess.run([sys.executable] + cmd, check=True, env=env)


def main():
    if len(sys.argv) < 5:
        sys.exit(__doc__)
    sites_dir, tmpl, runs, exe = sys.argv[1:5]
    sites = sys.argv[5:] or sorted(d for d in os.listdir(sites_dir) if os.path.isdir(os.path.join(sites_dir, d)))
    work = os.path.join(runs, '.build')
    if os.path.exists(work):
        shutil.rmtree(work)
    py = lambda f: os.path.join(HERE, f)
    run([py('make_runs.py'), sites_dir, tmpl, os.path.join(work, 'base'), exe] + sites)
    for s in sites:
        reg = REGION.get(s)
        if reg:
            run([py('gf_make_opt.py'), os.path.join(PARAMS, reg + '_step1.csv'), os.path.join(work, 'base'),
                 os.path.join(work, 'step1'), reg, s])
            d1 = os.path.join(work, 'step1', '%s_%s' % (reg, s))
        else:
            d1 = os.path.join(work, 'base', s)
        m1 = open(os.path.join(d1, 'opt_manifest.txt')).read() if reg else ''
        run([py('make_closure_runs.py'), '--inplace', d1])
        dst = os.path.join(runs, s)
        if os.path.exists(dst):
            shutil.rmtree(dst)
        if reg:
            run([py('gf_make_opt.py'), os.path.join(PARAMS, reg + '_step2.csv'), os.path.dirname(d1),
                 os.path.join(work, 'step2'), 'opt', os.path.basename(d1)], closure=True)
            shutil.move(os.path.join(work, 'step2', 'opt_' + os.path.basename(d1)), dst)
        else:
            shutil.move(d1, dst)
        clean = lambda t: t.replace(os.path.abspath(PARAMS), 'results/optimized').replace(PARAMS, 'results/optimized') \
                           .replace(os.path.abspath(d1), s).replace(d1, s)
        ex = os.path.join(dst, 'experiment.txt'); xt = clean(open(ex).read()); open(ex, 'w').write(xt)
        mf = os.path.join(dst, 'opt_manifest.txt')
        m2 = open(mf).read() if reg else ''
        open(mf, 'w').write(clean(
            '%s: optimized 1-D set-up (scripts/make_optimized_runs.py)\n' % s +
            '# step 1: v05 set-up from the site inputs (make_runs.py; negative ICs + relaxation values clipped to 0)\n' +
            ('# step 2: %s set, first Green\'s-function fit\n%s' % (reg, m1) if reg else
             '# step 2: none (no regional set; v05 parameters)\n') +
            '# step 3: closure\n%s' % xt +
            ('# step 4: %s set, second Green\'s-function iteration\n%s' % (reg, m2) if reg else '')))
        print('ready:', dst, '(%s)' % (reg + ' set + closure' if reg else 'v05 parameters + closure'))
    shutil.rmtree(work)


if __name__ == '__main__':
    main()
