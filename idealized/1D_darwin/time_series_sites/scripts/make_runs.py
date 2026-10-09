#!/usr/bin/env python3
"""
Build per-site 1-D run directories from extracted site inputs + namelist templates.

  python3 make_runs.py <sites_dir> <input_template_dir> <runs_dir> <mitgcmuv> [site ...]

For each site: copies namelists (filling @SITE@, @F0@, @XG@, @YG@ in `data`),
writes data.diagnostics and relax_mask_1x1x50, and links the extracted inputs.
"""
import os, sys, math, shutil
import numpy as np

DIAG3D = ['THETA', 'SALT', 'UVEL', 'VVEL', 'GGL90Kr', 'GGL90TKE', 'PP', 'PAR'] + \
         ['TRAC%02d' % n for n in range(1, 32)]
DIAG2D = ['MXLDEPTH', 'ETAN', 'oceQnet', 'oceQsw', 'TFLUX', 'SFLUX', 'EXFwspee',
          'SIarea', 'SIheff', 'fluxCO2', 'pCO2', 'pH', 'apCO2', 'gO2surf', 'sfcSolFe']
INPUTS = ['bathy_1x1', 'diffkr_1x1x50', 'T_ini_1x1x50', 'S_ini_1x1x50', 'iron_dust_1x1_mon',
          'relax_T_1x1x50_mon', 'relax_S_1x1x50_mon', 'EXF', 'apCO2'] + \
         ['ptracer%02d_ini' % n for n in range(1, 32)] + \
         ['relax_ptr%02d_1x1x50_mon' % n for n in (1, 2, 5, 6, 7, 18, 19)]
DRF = [10.00, 10.00, 10.00, 10.00, 10.00, 10.00, 10.00, 10.01, 10.03, 10.11, 10.32, 10.80, 11.76,
       13.42, 16.04, 19.82, 24.85, 31.10, 38.42, 46.50, 55.00, 63.50, 71.58, 78.90, 85.15, 90.18,
       93.96, 96.58, 98.25, 99.25, 100.01, 101.33, 104.56, 111.33, 122.83, 139.09, 158.94, 180.83,
       203.55, 226.50, 249.50, 272.50, 295.50, 318.50, 341.50, 364.50, 387.50, 410.50, 433.50, 456.50]
BGC_Z0, BGC_Z1 = 150., 250.      # BGC relaxation mask: 0 above Z0, linear ramp, 1 below Z1
OMEGA = 2. * math.pi / 86164.
DX = 0.3333333


def site_info(path):
    d = {}
    for line in open(path):
        k, v = line.split(' = ', 1)
        d[k.strip()] = v.strip()
    return d


def diagnostics(freq=86400.):
    """Two lists (3-D, 2-D) written to NetCDF via pkg/mnc: diags3D.*.nc, diags2D.*.nc."""
    out = ['# Diagnostics: freq %.0f s, NetCDF (mnc) => one file per list, time appended' % freq,
           ' &DIAGNOSTICS_LIST', ' diag_mnc = .TRUE.,', ' dumpAtLast = .TRUE.,']
    for n, (fn, names) in enumerate([('diags3D', DIAG3D), ('diags2D', DIAG2D)], 1):
        out += ["  frequency(%d) = %.1f," % (n, freq), "   fileName(%d) = '%s'," % (n, fn)]
        out += ["   fields(%d,%d) = '%-8s'," % (m, n, v) for m, v in enumerate(names, 1)]
    out += [' /', ' &DIAG_STATIS_PARMS', ' /']
    return '\n'.join(out) + '\n'


def main():
    sites_dir, tmpl, runs, exe = sys.argv[1:5]
    names = sys.argv[5:] or sorted(os.listdir(sites_dir))
    for name in names:
        sd = os.path.abspath(os.path.join(sites_dir, name))
        info = site_info(os.path.join(sd, 'site_info.txt'))
        xc, yc = float(info['xc']), float(info['yc'])
        rd = os.path.join(runs, name)
        os.makedirs(os.path.join(rd, 'mnc_out'), exist_ok=True)
        for f in os.listdir(tmpl):
            s = open(os.path.join(tmpl, f)).read()
            s = (s.replace('@SITE@', '%s (%.3fN, %.3fE, %s m)' % (name, yc, xc, info['depth']))
                  .replace('@F0@', '%.6e' % (2. * OMEGA * math.sin(math.radians(yc))))
                  .replace('@XG@', '%.6f' % (xc - DX / 2.))
                  .replace('@YG@', '%.6f' % (yc - DX / 2.)))
            open(os.path.join(rd, f), 'w').write(s)
        open(os.path.join(rd, 'data.diagnostics'), 'w').write(diagnostics())
        np.ones(50, '>f4').tofile(os.path.join(rd, 'relax_mask_1x1x50'))
        zc = np.cumsum(DRF) - np.array(DRF) / 2.
        np.clip((zc - BGC_Z0) / (BGC_Z1 - BGC_Z0), 0., 1.).astype('>f4').tofile(
            os.path.join(rd, 'relax_mask_bgc_1x1x50'))
        for f in INPUTS:
            dst = os.path.join(rd, f)
            if os.path.lexists(dst):
                os.remove(dst)
            os.symlink(os.path.join(sd, f), dst)
        dst = os.path.join(rd, 'mitgcmuv')
        if os.path.lexists(dst):
            os.remove(dst)
        os.symlink(os.path.abspath(exe), dst)
        print('%-7s -> %s  (xc %.3f, yc %.3f, f0 %.3e)' % (
            name, rd, xc, yc, 2. * OMEGA * math.sin(math.radians(yc))))


if __name__ == '__main__':
    main()
