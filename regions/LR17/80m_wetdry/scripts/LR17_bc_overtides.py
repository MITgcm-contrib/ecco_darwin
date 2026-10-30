#!/usr/bin/env python3
"""LR17 open-boundary test: IBI overtides replaced by FES M4, MS4 and M6.

Usage:  python3 LR17_bc_overtides.py MAIN_DIR OBCS_IN OBCS_OUT [START]
                                    [--cons=M4,MS4,M6] [--bands=3.6-4.2,5.6-6.2]
                                    [--fig=LR17_bc_overtides]
MAIN_DIR: the wetting_drying checkout that holds data/LR17 (IBI downloads,
FES atlases, LR17_zb.bin) and scripts/fes_tide.py. OBCS_IN: the v2 boundary
files OB[WS]{eta,u,v}.bin of a run from START (default 2022-03-01, hourly
from START 00 UTC, as gendata.py writes them). Writes the same six files to
OBCS_OUT (copy OB[WS]{t,s}.bin from OBCS_IN unchanged) and
figures/LR17/LR17_bc_overtides.{pdf,png}.

At every open-boundary point the IBI series that gendata.py interpolates
(sea level zos and depth-mean current ubar, vbar, hourly means placed at the
half hour as in v2) is band-passed over the 50-day download (20 Feb - 10 Apr
2022, FFT, linear trend removed; until 2026-10-20 the 22-day download was
used, which leaked up to 0.27 m rms of the band into the run window) in the quarter-diurnal (3.6-4.2 cpd: MN4,
M4, MS4) and sixth-diurnal (5.6-6.2 cpd: M6, 2MS6) bands. That band signal is
subtracted from the v2 files and the FES prediction of M4, MS4 and M6 is
added (sea level: FES2022b, current: FES2014; nodal corrections, as
fes_tide.py). --cons and --bands change the FES constituents and the IBI
bands (e.g. --cons=M4,MS4 --bands=3.6-4.2 keeps the IBI sixth-diurnal band). The current change is depth-uniform, as gendata.py adds the
depth-mean current to every level. Everything else stays as in v2.
The script first rebuilds the v2 boundary sea level from IBI and stops if it
differs from OBCS_IN by more than 1 mm.

Figure caption: open-boundary overtides for the LR17 boundary test. a, Sea
level in the quarter- and sixth-diurnal bands at the middle of the west
boundary, 1-5 March 2022: IBI (removed) and FES2022b M4 + MS4 + M6 (added).
b, Root-mean-square of that band signal along the two open boundaries.
"""
import os
import sys

import numpy as np
import xarray as xr
from scipy import ndimage
from scipy.interpolate import RegularGridInterpolator
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import wetdrystyle as ws  # noqa: E402

NX, NY, NR, DG = 936, 875, 14, 0.0008
LON0, LAT0 = -1.7001, 45.7499
lon = LON0 + (np.arange(NX) + 0.5) * DG
lat = LAT0 + (np.arange(NY) + 0.5) * DG
DATUM, ZOB = 5.5, -2.0
TAG = '2022-02-24_2022-03-18'                    # the download gendata.py used (v2 check)
TAG_BAND = '2022-02-20_2022-04-10'               # 50-day download for the band-pass
BANDS = [(3.6, 4.2), (5.6, 6.2)]                 # cycles per day (--bands)
OVERTIDES = ['M4', 'MS4', 'M6']                  # FES constituents added (--cons)


def fill_nearest(a):
    bad = ~np.isfinite(a)
    if bad.all() or not bad.any():
        return a
    idx = ndimage.distance_transform_edt(bad, return_distances=False, return_indices=True)
    return a[tuple(idx)]


def along(field3d, glat, glon, pts):
    """(ntime, npt): gendata.py's bilinear interpolation after filling land"""
    return np.array([RegularGridInterpolator((glat, glon), fill_nearest(f), bounds_error=False,
                                             fill_value=None)(pts) for f in field3d])


def half_hour(a):
    """gendata.at_hours (stamp_times False): value at full hour T = mean of the
    hourly means stamped T-1h and T; returns records for times[1:]"""
    return 0.5 * (a[:-1] + a[1:])


def band(a):
    """quarter- and sixth-diurnal band of hourly series a (ntime, npt)"""
    n = a.shape[0]
    t = np.arange(n)
    good = np.all(np.isfinite(a), axis=0)
    a = np.where(good, a, 0.0)
    p = np.polyfit(t, a, 1)
    r = a - (p[0][None] * t[:, None] + p[1][None])
    F = np.fft.rfft(r, axis=0)
    f = np.fft.rfftfreq(n, d=1 / 24.0)
    keep = np.zeros(f.size, bool)
    for lo, hi in BANDS:
        keep |= (f >= lo) & (f <= hi)
    return np.fft.irfft(np.where(keep[:, None], F, 0.0), n=n, axis=0)


V3_CONS = ['M2', 'S2', 'N2', 'K2', '2N2', 'MU2', 'NU2', 'L2', 'T2', 'M4', 'MS4']   # LR17 v3
V3_BANDS = [(1.80, 2.06), (3.6, 4.2)]


def swap_terms(D, pts, hts, cons=None, bands=None, tag_band=None, min_edge_days=8.0):
    """IBI band signal to remove and FES tide to add at boundary points
    pts = (lat, lon) pairs, hourly full-hour times hts. Returns dict of
    (ntime, npt) arrays: b_eta, b_u, b_v (IBI zos - its mean, ubar, vbar,
    half-hour convention, band-passed over the 50-day download) and f_eta,
    f_u, f_v (FES2022b sea level, FES2014 depth-mean current). Used by main()
    and by gendata.py --fes-boundary. tag_band names the IBI download used for the
    band-pass (default TAG_BAND); the run window must lie at least min_edge_days
    from both ends of it, or the FFT leaks into the window."""
    import fes_tide
    fes_tide.D = D
    cons = cons or OVERTIDES
    old = BANDS[:]
    if bands:
        BANDS[:] = bands
    try:
        tag = tag_band or TAG_BAND
        sshL = xr.open_dataset(os.path.join(D, 'ibi', f'ibi_ssh_{tag}.nc'))
        curL = xr.open_dataset(os.path.join(D, 'ibi', f'ibi_cur_{tag}.nc'))
        glatL, glonL = sshL.latitude.values, sshL.longitude.values
        assert np.array_equal(sshL.time.values, curL.time.values)
        tL = sshL.time.values[1:]
        kL = np.searchsorted(tL, hts)
        assert np.array_equal(tL[kL], hts), f'run window outside the IBI download {tag}'
        edge = min((hts[0] - tL[0]), (tL[-1] - hts[-1])) / np.timedelta64(1, 'D')
        if edge < min_edge_days:
            sys.exit(f'band-pass download {tag} ends {edge:.1f} d from the run window '
                     f'(need >= {min_edge_days:g} d): download a longer IBI ssh+cur series')
        zl = sshL.zos.values - np.nanmean(sshL.zos.values, axis=0)
        out = {'b_eta': band(half_hour(along(zl, glatL, glonL, pts)))[kL],
               'b_u': band(half_hour(along(curL.ubar.values, glatL, glonL, pts)))[kL],
               'b_v': band(half_hour(along(curL.vbar.values, glatL, glonL, pts)))[kL]}
    finally:
        BANDS[:] = old
    plat, plon = pts[:, 0], pts[:, 1]
    for kind, src in (('eta', 'FES2022b'), ('u', 'FES2014'), ('v', 'FES2014')):
        A = fes_tide.atlas(src, kind, cons)
        if sorted(A) != sorted(cons):
            sys.exit(f'{src} {kind}: found only {sorted(A)}')
        amp, pha, names = fes_tide.constants(A, plat, plon)
        out['f_' + kind] = fes_tide.predict_fast(hts, amp, pha, names, float(np.mean(plat))).T
    return out


def main(top, obin, obout, start='2022-03-01', fig_name='LR17_bc_overtides'):
    D = os.path.join(top, 'data', 'LR17')
    sys.path.insert(0, os.path.join(top, 'scripts'))
    import fes_tide
    zb = np.fromfile(os.path.join(D, 'work', 'LR17_zb.bin'), '>f4').reshape(NY, NX)
    model = zb < 90
    model[:, 0] &= zb[:, 0] < ZOB
    model[0, :] &= zb[0, :] < ZOB
    ssh = xr.open_dataset(os.path.join(D, 'ibi', f'ibi_ssh_{TAG}.nc'))
    cur = xr.open_dataset(os.path.join(D, 'ibi', f'ibi_cur_{TAG}.nc'))
    glat, glon = ssh.latitude.values, ssh.longitude.values
    zos = ssh.zos.values - np.nanmean(ssh.zos.values, axis=0)
    stime = ssh.time.values
    assert np.array_equal(stime, cur.time.values)
    t_full = stime[1:]                                   # hours with a half-hour value
    t0 = np.datetime64(start + 'T00:00')
    nrec = os.path.getsize(os.path.join(obin, 'OBWeta.bin')) // (4 * NY)
    hts = t0 + np.arange(nrec) * np.timedelta64(1, 'h')
    k = np.searchsorted(t_full, hts)
    assert np.array_equal(t_full[k], hts)
    print(f'band-pass: 50-day IBI download {TAG_BAND}; window {hts[0]} .. {hts[-1]}')
    sides = {'W': (np.full(NY, lon[0] - DG / 2), lat, model[:, 0], NY),
             'S': (lon, np.full(NX, lat[0] - DG / 2), model[0, :], NX)}
    os.makedirs(obout, exist_ok=True)
    show = {}
    for s, (plon, plat, open_, n) in sides.items():
        pts = np.c_[plat, plon]
        rd = lambda name, shape: np.fromfile(os.path.join(obin, f'OB{s}{name}.bin'), '>f4') \
            .reshape(shape).astype(float)
        eta_v2 = rd('eta', (nrec, n))
        u_v2, v_v2 = rd('u', (nrec, NR, n)), rd('v', (nrec, NR, n))
        eta_i = half_hour(along(zos, glat, glon, pts)) - DATUM
        err = np.abs(eta_i[k][:, open_] - eta_v2[:, open_]).max()
        print(f'{s}: v2 boundary sea level rebuilt from IBI, max difference {err * 1e3:.3f} mm')
        if err > 1e-3:
            sys.exit('rebuilt IBI boundary does not match OBCS_IN')
        # band-pass over the 50-day download: the run window sits >= 9 days from
        # both ends (the 22-day series leaked up to 0.27 m rms into the window)
        T = swap_terms(D, pts, hts)
        b_eta, b_u, b_v = T['b_eta'], T['b_u'], T['b_v']
        f = {'eta': T['f_eta'], 'u': T['f_u'], 'v': T['f_v']}
        m2 = open_[None]
        m3 = open_[None, None]
        eta = np.where(m2, eta_v2 - b_eta + f['eta'], eta_v2)
        u = np.where(m3, u_v2 + (f['u'] - b_u)[:, None, :], u_v2)
        v = np.where(m3, v_v2 + (f['v'] - b_v)[:, None, :], v_v2)
        for name, a in (('eta', eta), ('u', u), ('v', v)):
            a.astype('>f4').tofile(os.path.join(obout, f'OB{s}{name}.bin'))
        o = open_
        print(f'{s}: band rms on open points, IBI removed / FES added: sea level '
              f'{np.sqrt((b_eta[:, o] ** 2).mean()):.3f} / {np.sqrt((f["eta"][:, o] ** 2).mean()):.3f} m, '
              f'u {np.sqrt((b_u[:, o] ** 2).mean()):.3f} / {np.sqrt((f["u"][:, o] ** 2).mean()):.3f}, '
              f'v {np.sqrt((b_v[:, o] ** 2).mean()):.3f} / {np.sqrt((f["v"][:, o] ** 2).mean()):.3f} m/s')
        show[s] = (plat if s == 'W' else plon, o, b_eta, f['eta'])
    print('wrote', os.path.join(obout, 'OB[WS]{eta,u,v}.bin'))

    ws.use()
    fig, axs = plt.subplots(1, 2, figsize=ws.size(183, 70), constrained_layout=True)
    pos, o, b_eta, f_eta = show['W']
    jm = np.where(o)[0][o.sum() // 2]
    th = (hts - t0) / np.timedelta64(1, 'h') / 24.0
    a = axs[0]
    a.plot(th, b_eta[:, jm], color=ws.OKABE[0], lw=0.8, label='IBI band (removed)')
    a.plot(th, f_eta[:, jm], color=ws.OKABE[1], lw=0.8,
           label=f'FES2022b {" + ".join(OVERTIDES)} (added)')
    a.set_xlabel(f'Days from {start} 00 UTC')
    a.set_ylabel(f'Sea level, {lat[jm]:.2f}°N west boundary (m)')
    a.legend(loc='upper center', bbox_to_anchor=(0.5, -0.25), ncol=2, frameon=False)
    b = axs[1]
    for s, ls in (('W', '-'), ('S', '--')):
        pos, o, b_eta, f_eta = show[s]
        rb = np.where(o, np.sqrt((b_eta ** 2).mean(0)), np.nan)
        rf = np.where(o, np.sqrt((f_eta ** 2).mean(0)), np.nan)
        x = np.arange(pos.size)
        b.plot(x, rb, color=ws.OKABE[0], ls=ls, lw=0.8, label=f'IBI, {"west" if s == "W" else "south"}')
        b.plot(x, rf, color=ws.OKABE[1], ls=ls, lw=0.8, label=f'FES, {"west" if s == "W" else "south"}')
    b.set_xlabel('Boundary cell (west: south to north; south: west to east)')
    b.set_ylabel('Band rms sea level (m)')
    b.legend(loc='upper center', bbox_to_anchor=(0.5, -0.25), ncol=2, frameon=False)
    ws.letter(a, 'a')
    ws.letter(b, 'b')
    ws.save(fig, fig_name)


if __name__ == '__main__':
    opts = dict(a[2:].split('=', 1) for a in sys.argv[1:] if a.startswith('--'))
    args = [a for a in sys.argv[1:] if not a.startswith('--')]
    if 'cons' in opts:
        OVERTIDES[:] = opts['cons'].split(',')
    if 'bands' in opts:
        BANDS[:] = [tuple(float(x) for x in b.split('-')) for b in opts['bands'].split(',')]
    print(f'FES constituents {OVERTIDES}; IBI bands removed {BANDS} cpd')
    main(*args[:4], fig_name=opts.get('fig', 'LR17_bc_overtides'))
