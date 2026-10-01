#!/usr/bin/env python3
"""Thacker (1981) planar surface in a frictionless parabolic basin (1-D).

Basin depth below Thacker's mean level:  H(x) = h0 (1 - x^2/a^2)
Exact solution (u uniform in space):
    u(t)   = U sin(w t)
    eta(x,t) = -(U w / g) x cos(w t) - (U^2 / 2g) cos^2(w t)
    w = sqrt(2 g h0) / a
valid where H + eta > 0. The water body is a rigid parabolic lens of
half-width a whose centre moves as  x_c(t) = -(U/w) cos(w t).

MITgcm has no land above its rest level, so the model rest level is put
DATUM metres above Thacker's mean level:  R_low = -(H + DATUM),
eta_model = eta - DATUM, and cells beyond |x| > XLAND are land.

Usage:
    python3 analytic.py gen                 # write bathy + initial eta
    python3 analytic.py cmp RUNDIR [--plot] [--exact]
        compare Eta.*.data in RUNDIR with the time-discrete solution
        (default) or the continuous one (--exact)
"""
import glob
import os
import re
import sys

import numpy as np

# --- parameters (must match data: delR, dXspacing, deltaT, gravity)
G = 9.81
A = 3000.0          # basin half-width at mean level [m]
H0 = 10.0           # central depth [m]
XAMP = float(os.environ.get('WAD_XAMP', 1000.0))  # lens-centre excursion [m]
DATUM = 10.0        # model rest level above Thacker's mean level [m]
XLAND = 4200.0      # cells with |x| > XLAND are land
# grid; WAD_NX overrides (the domain stays 10 km: DX = 10000/NX)
NX = int(os.environ.get('WAD_NX', 250))
DX = 10000.0 / NX
DMIN = 0.05         # wadMinDepth

W = np.sqrt(2.0 * G * H0) / A
T = 2.0 * np.pi / W
U = XAMP * W

x = (np.arange(NX) + 0.5) * DX - NX * DX / 2.0
H = H0 * (1.0 - x**2 / A**2)
ocean = np.abs(x) <= XLAND
r_low = np.where(ocean, -(H + DATUM), 0.0)


def eta_exact(t):
    """Thacker eta relative to the mean level, and wet mask."""
    c = np.cos(W * t)
    eta = -(U * W / G) * x * c - (U**2 / (2.0 * G)) * c**2
    return eta, (H + eta) > 0.0


def eta_discrete(t, dt):
    """Same planar solution, but time-stepped exactly as MITgcm does with
    implicSurfPress = implicDiv2DFlow = 1 (backward Euler for the surface
    pressure gradient and the divergence, depth at time n):
        u' = u - g dt s'          s' = s + (2 h0/a^2) dt u'
        c' = c - dt u' s
    so the difference to the model is spatial / wetting-drying error only.
    """
    eta, wet, _ = _discrete(t, dt)
    return eta, wet


def amp_discrete(t, dt):
    """Amplitude of the time-discrete oscillation relative to the exact
    one: sqrt(s^2 + (u w/g)^2) / (U w/g)  (= 1 for the exact solution)."""
    return _discrete(t, dt)[2]


def _discrete(t, dt):
    k = 2.0 * H0 / A**2
    u, sl, c = 0.0, -U * W / G, -U**2 / (2.0 * G)
    for _ in range(int(round(t / dt))):
        un = (u - G * dt * sl) / (1.0 + G * k * dt * dt)
        c = c - dt * un * sl
        sl = sl + k * dt * un
        u = un
    eta = sl * x + c
    amp = np.hypot(sl, u * W / G) / (U * W / G)
    return eta, (H + eta) > 0.0, amp


def shoreline_of(eta_wet):
    """Exact shoreline (H + eta = 0 roots) of a planar surface."""
    eta, _ = eta_wet
    sl = (eta[1] - eta[0]) / DX
    c = eta[0] - sl * x[0]
    # h0 (1 - x^2/a^2) + sl x + c = 0
    qa, qb, qc = -H0 / A**2, sl, H0 + c
    disc = np.sqrt(qb * qb - 4 * qa * qc)
    r = sorted([(-qb + disc) / (2 * qa), (-qb - disc) / (2 * qa)])
    return r[0], r[1]


def shoreline(t):
    xc = -(U / W) * np.cos(W * t)
    return xc - A, xc + A


def model_eta(eta):
    """Thacker eta -> model eta, with the film where dry."""
    e, wet = eta
    return np.where(ocean, np.maximum(e - DATUM, r_low + DMIN), 0.0)


def gen():
    here = os.path.dirname(os.path.abspath(__file__))
    r_low.astype('>f8').tofile(os.path.join(here, 'bathy_thacker.bin'))
    model_eta(eta_exact(0.0)).astype('>f8').tofile(
        os.path.join(here, 'eta_thacker.bin'))
    print(f'w = {W:.9e} 1/s  T = {T:.6f} s  U = {U:.4f} m/s')
    print(f'max |eta| over the lens at t=0: {np.abs(eta_exact(0)[0][eta_exact(0)[1]]).max():.4f} m')
    print(f'wrote bathy_thacker.bin, eta_thacker.bin ({NX} x 1)')


def read_run(rundir):
    """Return list of (time, eta_model) for all Eta.*.data files."""
    dt = None
    with open(os.path.join(rundir, 'data')) as f:
        for line in f:
            m = re.match(r'\s*deltaT\s*=\s*([0-9.eE+-]+)', line)
            if m:
                dt = float(m.group(1))
    out = []
    for fn in sorted(glob.glob(os.path.join(rundir, 'Eta.*.data'))):
        it = int(fn.split('.')[-2])
        out.append((it * dt, np.fromfile(fn, '>f8').reshape(-1)[:NX]))
    return out, dt


def wet_edges(depth, thresh):
    """x edges of the wet region (depth > thresh) connected to the
    deepest point (isolated puddles on the beach are ignored)."""
    wet = depth > thresh
    k = int(np.argmax(depth))
    if not wet[k]:
        return np.nan, np.nan
    lo = k
    while lo > 0 and wet[lo - 1]:
        lo -= 1
    hi = k
    while hi < NX - 1 and wet[hi + 1]:
        hi += 1
    return x[lo] - DX / 2, x[hi] + DX / 2


def cmp(rundir, plot=False, ref='discrete'):
    runs, dt = read_run(rundir)
    amp = U * W / G * A          # surface elevation amplitude at the shoreline
    print(f'eta amplitude scale A = (U w/g) a = {amp:.4f} m ; T = {T:.4f} s'
          f' ; reference = {ref}')
    print(f'{"t/T":>7} {"L2(eta)/A":>10} {"L2/Aref":>9} {"Linf/A":>9} '
          f'{"shoreL err":>10} {"shoreR err":>10} {"dV/V0":>11}')
    vol0 = None
    rows = []
    for t, em in runs:
        ex = eta_discrete(t, dt) if ref == 'discrete' else eta_exact(t)
        e_ex, wet_ex = ex
        depth = np.where(ocean, em - r_low, 0.0)
        vol = depth[ocean].sum() * DX
        if vol0 is None:
            vol0 = vol
        # compare where both analytic and model are wet (away from film)
        both = wet_ex & (depth > 0.1) & ocean
        err = (em + DATUM - e_ex)[both]
        l2 = np.sqrt(np.mean(err**2)) / amp if err.size else np.nan
        linf = np.abs(err).max() / amp if err.size else np.nan
        sl, sr = shoreline_of(ex)
        ml, mr = wet_edges(depth, 0.1)
        aref = amp_discrete(t, dt) if ref == 'discrete' else 1.0
        rows.append((t / T, l2, l2 / aref, linf, ml - sl, mr - sr,
                     (vol - vol0) / vol0))
        print(f'{t/T:7.3f} {l2:10.4e} {l2/aref:9.3e} {linf:9.3e} '
              f'{ml-sl:10.1f} {mr-sr:10.1f} {(vol-vol0)/vol0:11.3e}')
    if plot:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(figsize=(9, 4.5))
        bed = -H
        ax.fill_between(x, bed, bed.min() - 1, color='0.8', lw=0)
        sel = [r for r in runs if abs(r[0] / T - round(r[0] / T * 4) / 4)
               < 1e-3][-5:]
        for t, em in sel:
            e_ex, wet_ex = (eta_discrete(t, dt) if ref == 'discrete'
                            else eta_exact(t))
            l, = ax.plot(x, np.where(wet_ex, e_ex, np.nan), '-', lw=1,
                         label=f't/T={t/T:.2f}')
            dep = em - r_low
            ax.plot(x, np.where((dep > 0.1) & ocean, em + DATUM, np.nan),
                    '.', ms=3, color=l.get_color())
        ax.set_xlim(-4500, 4500)
        ax.set_ylim(-11, 4)
        ax.set_xlabel('x (m)')
        ax.set_ylabel('elevation (m)')
        ax.set_title(f'Thacker 1-D: {ref} reference (lines) vs pkg/wad (dots)')
        ax.legend(fontsize=8, loc='lower center', ncol=3)
        fig.tight_layout()
        fn = os.path.join(rundir, 'thacker_cmp.png')
        fig.savefig(fn, dpi=130)
        print('wrote', fn)
    return rows


if __name__ == '__main__':
    if len(sys.argv) >= 2 and sys.argv[1] == 'gen':
        gen()
    elif len(sys.argv) >= 3 and sys.argv[1] == 'cmp':
        cmp(sys.argv[2], plot='--plot' in sys.argv,
            ref='exact' if '--exact' in sys.argv else 'discrete')
    else:
        print(__doc__)
        sys.exit(1)
