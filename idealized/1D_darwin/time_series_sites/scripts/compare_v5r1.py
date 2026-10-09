#!/usr/bin/env python3
"""
Compare a 1-D column run with the v5r1 3-D solution at the same llc270 point.

  python3 compare_v5r1.py <run_dir> <site_dir> <out.png>

Reads daily-mean NetCDF diagnostics from <run_dir>/mnc_out and the v5r1 monthly/daily
series extracted by extract_v5r1_sites.py from <site_dir>.
"""
import os, sys, glob, datetime as dt
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

NR = 50
DT = 1200.
T0 = dt.datetime(1992, 1, 1)
DRF = np.array([
    10.00, 10.00, 10.00, 10.00, 10.00, 10.00, 10.00, 10.01,
    10.03, 10.11, 10.32, 10.80, 11.76, 13.42, 16.04, 19.82, 24.85,
    31.10, 38.42, 46.50, 55.00, 63.50, 71.58, 78.90, 85.15, 90.18,
    93.96, 96.58, 98.25, 99.25, 100.01, 101.33, 104.56, 111.33, 122.83,
    139.09, 158.94, 180.83, 203.55, 226.50, 249.50, 272.50, 295.50, 318.50,
    341.50, 364.50, 387.50, 410.50, 433.50, 456.50])
ZC = np.cumsum(DRF) - DRF / 2.


def steps2time(steps, period):
    # record time stamp = end of averaging window; plot at its centre
    return [T0 + dt.timedelta(seconds=s * DT - period / 2.) for s in steps]


def read_run(run, name):
    """Daily-mean field from the run's NetCDF diagnostics (mnc_out/diags3D|diags2D)."""
    import netCDF4
    for lst in ['diags3D', 'diags2D']:
        for f in sorted(glob.glob(os.path.join(run, 'mnc_out', lst + '.*.nc'))):
            d = netCDF4.Dataset(f)
            if name in d.variables:
                a = np.array(d[name][:]).reshape(len(d['T']), -1)
                t = np.array([T0 + dt.timedelta(seconds=float(x) - 43200.) for x in d['T'][:]])
                return t, a
    raise KeyError(name)


def read_v5r1(site, freq, name, nlev):
    d = os.path.join(site, 'v5r1_' + freq)
    f = os.path.join(d, name + '.bin')
    if not os.path.exists(f):
        return None, None
    steps = np.loadtxt(os.path.join(d, 'steps.txt'), dtype=int)
    a = np.fromfile(f, '>f4')
    a = a.reshape(len(steps), -1)          # nlev inferred (v5r1 daily pCO2/pH are 3-D)
    per = 86400. if freq == 'daily' else 30.4 * 86400.
    return np.array(steps2time(steps, per)), a


def main():
    run, site, out = sys.argv[1:4]
    name = os.path.basename(os.path.normpath(run))
    nwet = int([l for l in open(os.path.join(site, 'site_info.txt')) if l.startswith('nwet')][0].split('=')[1])
    zmax = 500.
    kz = np.searchsorted(ZC, zmax)

    fig, ax = plt.subplots(4, 2, figsize=(14, 13), constrained_layout=True)
    # --- T Hovmoller: 1-D vs v5r1
    t1, T1 = read_run(run, 'THETA')
    tm, Tm = read_v5r1(site, 'monthly', 'THETA', NR)
    lev = np.linspace(np.nanmin(Tm[:, :kz]), np.nanmax(Tm[:, :kz]), 21)
    for a, t, T, ttl in [(ax[0, 0], t1, T1, '1-D'), (ax[0, 1], tm, Tm, 'v5r1 (monthly)')]:
        cs = a.contourf(t, ZC[:kz], T[:, :kz].T, lev, cmap='RdYlBu_r', extend='both')
        a.set_ylim(zmax, 0); a.set_ylabel('depth (m)'); a.set_title('%s  THETA (degC), %s' % (name, ttl))
    fig.colorbar(cs, ax=ax[0, :], shrink=0.8)

    # --- NO3 Hovmoller
    t1, N1 = read_run(run, 'TRAC02')
    tm, Nm = read_v5r1(site, 'monthly', 'NO3', NR)
    lev = np.linspace(0, np.nanmax(Nm[:, :kz]), 21)
    for a, t, N, ttl in [(ax[1, 0], t1, N1, '1-D'), (ax[1, 1], tm, Nm, 'v5r1 (monthly)')]:
        cs = a.contourf(t, ZC[:kz], N[:, :kz].T, lev, cmap='viridis', extend='max')
        a.set_ylim(zmax, 0); a.set_ylabel('depth (m)'); a.set_title('NO3 (mmol N m$^{-3}$), %s' % ttl)
    fig.colorbar(cs, ax=ax[1, :], shrink=0.8)

    # --- surface time series
    def line(a, t1, y1, tv, yv, ttl):
        a.plot(t1, y1, lw=0.6, label='1-D (daily)')
        if tv is not None:
            a.plot(tv, yv, lw=0.9, label='v5r1')
        a.set_title(ttl); a.legend(fontsize=8)

    t1, M1 = read_run(run, 'MXLDEPTH'); td, Md = read_v5r1(site, 'daily', 'mldDepth', 1)
    line(ax[2, 0], t1, M1[:, 0], td, Md[:, 0] if Md is not None else None, 'Mixed-layer depth (m)')
    ax[2, 0].invert_yaxis()
    chl1 = sum(read_run(run, 'TRAC%02d' % n)[1][:, 0] for n in range(27, 32))
    chld = None
    for n in range(1, 6):
        tv, c = read_v5r1(site, 'daily', 'surfChl%d' % n, 1)
        chld = c[:, 0] if chld is None else chld + c[:, 0]
    line(ax[2, 1], t1, chl1, tv, chld, 'Surface total Chl (mg m$^{-3}$)')
    t1, P1 = read_run(run, 'pCO2'); td, Pd = read_v5r1(site, 'daily', 'pCO2', 1)
    line(ax[3, 0], t1, P1[:, 0] * 1e6, td, Pd[:, 0] * 1e6 if Pd is not None else None, 'Surface pCO2 (uatm)')
    t1, D1 = read_run(run, 'TRAC01'); tm, Dm = read_v5r1(site, 'monthly', 'DIC', NR)
    line(ax[3, 1], t1, D1[:, 0], tm, Dm[:, 0], 'Surface DIC (mmol C m$^{-3}$)')
    fig.savefig(out, dpi=130)
    print('wrote', out)


if __name__ == '__main__':
    main()
