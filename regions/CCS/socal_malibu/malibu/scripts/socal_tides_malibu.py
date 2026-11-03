#!/usr/bin/env python3
"""Tidal open-boundary forcing for the Malibu nest from the socal 2-km tidal run (2026-10-13).

The socal nest run_v3_prod (FES2014 boundary currents calibrated against FES2022b, 2-day ramp; see
downscaling/socal/DECISIONS.md) wrote hourly-mean 3-D UVEL/VVEL (diags/hourly/state_3d: THETA, SALT,
UVEL, VVEL; 168 x 101 x 60, deltaT 120 s, startDate 2025-01-01). Here, at each Malibu open-boundary
face, we take the socal depth-integrated normal transport per unit width q = sum_k U hFacW drF
(m2/s), fit tidal harmonics to it and convert to the depth-uniform normal velocity that pkg/obcs
adds after OBCS_BALANCE_FLOW (obcs_add_tides.F: OB?u += am cos(2 pi (myTime - ph) / T)):
  am = q_amplitude / H_malibu (H = Malibu wet depth at the face), ph in s.
Using transport keeps the socal tidal volume flux through each Malibu boundary (the tidal sea level
in the nest comes from that flux, since the balance step removes only the daily-mean net inflow).

Fit: ordinary least squares on hours after the 2-day socal ramp, t = myTime of both runs
(startDate 2025-01-01 00:00, the hourly mean centered 1800 s before each stamp); unknowns: mean, trend,
M2, N2, S2 (+K2), K1 (+P1), O1, Q1. A 179-day record cannot separate K2 from S2 or P1 from K1
(182.6-d beat), so K2 and P1 ride on S2 and K1 with the amplitude ratio and phase offset of the socal
boundary forcing itself (forcings/OBCS_tides_cal1, median over all OB points). Amplitudes are
divided by the hourly-mean damping sin(w h/2)/(w h/2) (M2 0.989). The fixed amplitude/phase holds
the nodal factors of mid-2025 (fitted over Jan-Jun 2025), so it is valid for this window only.

Output (<out>/): OB{W,S,E}{am,ph}_tides.bin, big-endian float32, 10 records (tidalComponents) of
the boundary length (72 W/E, 128 S), constituents M2,S2,N2,K2,K1,O1,P1,Q1 (+2 zero records), as the
socal data.obcs; data.obcs_tides.snippet; tides_malibu.json (fit skill, amplitudes, transports).

Usage: python3 socal_tides_malibu.py <socal run dir> <socal grid npz> <socal OBCS_tides dir> <out dir>
"""
import glob
import json
import os
import sys

import numpy as np
from scipy.interpolate import RegularGridInterpolator

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from make_ic_obcs import nest_grid, hfac_c, boundary_points   # noqa: E402

SNX, SNY, SNR, SDT = 168, 101, 60, 120.0
CONS = ['M2', 'S2', 'N2', 'K2', 'K1', 'O1', 'P1', 'Q1']
PER = dict(M2=44714.164313, S2=43200.000017, N2=45570.053572, K2=43082.045289,
           K1=86164.090578, O1=92949.630100, P1=86637.204835, Q1=96726.083762)   # socal data.obcs
FIT = ['M2', 'N2', 'S2', 'K1', 'O1', 'Q1']     # K2 rides on S2, P1 on K1
RIDE = {'S2': 'K2', 'K1': 'P1'}
NTC = 10
RAMP = 2 * 86400.0
SIDES = {'west': 'W', 'south': 'S', 'east': 'E'}


def socal_ratios(tdir):
    """K2/S2 and P1/K1 amplitude ratio + phase offset (rad) of the socal boundary forcing."""
    out = {}
    for main, rid in RIDE.items():
        r, d = [], []
        for s, n in (('W', SNY), ('E', SNY), ('S', SNX), ('N', SNX)):
            am = np.fromfile(os.path.join(tdir, f'OB{s}am_tides.bin'), '>f8').reshape(NTC, n)
            ph = np.fromfile(os.path.join(tdir, f'OB{s}ph_tides.bin'), '>f8').reshape(NTC, n)
            a0, a1 = am[CONS.index(main)], am[CONS.index(rid)]
            ok = a0 > 1e-3
            r.append(a1[ok] / a0[ok])
            th0 = 2 * np.pi * ph[CONS.index(main)][ok] / PER[main]
            th1 = 2 * np.pi * ph[CONS.index(rid)][ok] / PER[rid]
            d.append(np.angle(np.exp(1j * (th1 - th0))))
        r, d = np.concatenate(r), np.concatenate(d)
        out[rid] = dict(ratio=float(np.median(r)), dphase=float(np.angle(np.mean(np.exp(1j * d)))),
                        ratio_iqr=[float(np.percentile(r, 25)), float(np.percentile(r, 75))])
    return out


def design(t, rides):
    cols, names = [np.ones_like(t), (t - t.mean()) / 86400.0], ['mean', 'trend']
    for c in FIT:
        w = 2 * np.pi / PER[c]
        cc, ss = np.cos(w * t), np.sin(w * t)
        if c in RIDE:
            rd = RIDE[c]
            wr, r, d = 2 * np.pi / PER[rd], rides[rd]['ratio'], rides[rd]['dphase']
            cc = cc + r * np.cos(wr * t - d)
            ss = ss + r * np.sin(wr * t - d)
        cols += [cc, ss]
        names += [c + '_c', c + '_s']
    return np.stack(cols, 1), names


def main():
    run, gnpz, tdir, out = sys.argv[1:5]
    os.makedirs(out, exist_ok=True)
    g = nest_grid()
    hf = hfac_c(g)
    dz = np.asarray(g['dz'])
    # Malibu wet depth at the normal-velocity faces (hFacW/S = min of the two cells, as ini_masks_etc)
    Hm = {'west': (np.minimum(hf[:, :, 0], hf[:, :, 1]) * dz[:, None]).sum(0),
          'east': (np.minimum(hf[:, :, -2], hf[:, :, -1]) * dz[:, None]).sum(0),
          'south': (np.minimum(hf[:, 0, :], hf[:, 1, :]) * dz[:, None]).sum(0)}
    pts = {s: boundary_points(g, s, 'u' if s != 'south' else 'v') for s in SIDES}

    sg = np.load(gnpz)
    xc, yc, hfs, drf = sg['XC'][0], sg['YC'][:, 0], sg['HFacC'], sg['drF']
    dx, dy = np.diff(xc).mean(), np.diff(yc).mean()
    xg, ygr = xc - dx / 2, yc - dy / 2                  # socal is a regular 1/48 deg lon/lat grid
    i0 = np.searchsorted(xc, g['xg'][0] - 3 * dx); i1 = np.searchsorted(xc, g['xg'][-1] + 3 * dx)
    j0 = np.searchsorted(yc, g['yg'][0] - 3 * dy); j1 = np.searchsorted(yc, g['yg'][-1] + 3 * dy)
    hW = np.minimum(hfs[:, j0:j1, i0:i1], hfs[:, j0:j1, i0 - 1:i1 - 1]) * drf[:, None, None]
    hS = np.minimum(hfs[:, j0:j1, i0:i1], hfs[:, j0 - 1:j1 - 1, i0:i1]) * drf[:, None, None]
    Hs_u, Hs_v = hW.sum(0), hS.sum(0)

    fs = sorted(glob.glob(os.path.join(run, 'diags', 'hourly', 'state_3d.*.data')))
    t = np.array([int(f.split('.')[-2]) * SDT - 1800.0 for f in fs])
    keep = t >= RAMP
    fs, t = [f for f, k in zip(fs, keep) if k], t[keep]
    q = {s: np.empty((t.size, pts[s][0].size)) for s in SIDES}
    for n, f in enumerate(fs):
        a = np.memmap(f, '>f4', 'r', shape=(4, SNR, SNY, SNX))
        qu = (np.asarray(a[2, :, j0:j1, i0:i1]) * hW).sum(0)
        qv = (np.asarray(a[3, :, j0:j1, i0:i1]) * hS).sum(0)
        for s in SIDES:
            x, y = pts[s]
            if s == 'south':
                f2 = RegularGridInterpolator((ygr[j0:j1], xc[i0:i1]), qv)
            else:
                f2 = RegularGridInterpolator((yc[j0:j1], xg[i0:i1]), qu)
            q[s][n] = f2(np.c_[y, x])
        if n % 500 == 0:
            print(f'hour {n}/{t.size}', flush=True)

    rides = socal_ratios(tdir)
    X, names = design(t, rides)
    res = {'hours': int(t.size), 't0_h': float(t[0] / 3600), 't1_h': float(t[-1] / 3600), 'rides': rides,
           'periods_s': PER, 'sides': {}}
    snippet = [' useOBCStides = .TRUE.,',
               '# M2,S2,N2,K2,K1,O1,P1,Q1 from socal run_v3_prod hourly transports (socal_tides_malibu.py);',
               '# phase zero = startDate 2025-01-01, nodal factors of Jan-Jun 2025 (valid for this window)']
    snippet += [f' tidalPeriod({k + 1}) = {PER[c]:.6f},' for k, c in enumerate(CONS)]
    qin = {}
    for s, L in SIDES.items():
        Q = q[s]
        coef, *_ = np.linalg.lstsq(X, Q, rcond=None)
        resid = Q - X @ coef
        detr = Q - X[:, :2] @ coef[:2]
        skill = 1 - resid.var(0) / np.maximum(detr.var(0), 1e-12)
        am = np.zeros((NTC, Q.shape[1]))
        ph = np.zeros((NTC, Q.shape[1]))
        H = Hm[s]
        for c in FIT:
            a, b = coef[names.index(c + '_c')], coef[names.index(c + '_s')]
            w = 2 * np.pi / PER[c]
            damp = np.sin(w * 1800.0) / (w * 1800.0)
            A, phi = np.hypot(a, b) / damp, np.arctan2(b, a)      # q = A cos(w t - phi)
            ride = [(RIDE[c], A * rides[RIDE[c]]['ratio'], phi + rides[RIDE[c]]['dphase'])] if c in RIDE else []
            for cc, AA, pp in [(c, A, phi)] + ride:
                k = CONS.index(cc)
                wk = 2 * np.pi / PER[cc]
                am[k] = np.where(H > 1.0, AA / np.maximum(H, 1.0), 0.0)
                ph[k] = np.mod(pp / wk, PER[cc])
        am.astype('>f4').tofile(os.path.join(out, f'OB{L}am_tides.bin'))
        ph.astype('>f4').tofile(os.path.join(out, f'OB{L}ph_tides.bin'))
        snippet += [f" OB{L}amFile = 'tides/OB{L}am_tides.bin',", f" OB{L}phFile = 'tides/OB{L}ph_tides.bin',"]
        ds = (np.diff(g['yg']) * 111.32e3 if s != 'south'
              else np.diff(g['xg']) * 111.32e3 * np.cos(np.deg2rad(g['yg'][1])))
        wet = H > 1.0
        trc = {c: (am[CONS.index(c)] * H * ds * np.exp(1j * 2 * np.pi * ph[CONS.index(c)] / PER[c]))[wet].sum()
               for c in CONS}                                   # complex transport, east/north positive
        tr = {c: float(abs(v)) for c, v in trc.items()}
        for c in CONS:
            qin[c] = qin.get(c, 0) + (-trc[c] if s == 'east' else trc[c])   # inflow positive
        res['sides'][s] = dict(n=int(Q.shape[1]), wet=int(wet.sum()),
                               skill_median=float(np.median(skill[wet])), skill_p10=float(np.percentile(skill[wet], 10)),
                               am_median_cm_s={c: round(100 * float(np.median(am[CONS.index(c)][wet])), 2) for c in CONS},
                               am_max_cm_s={c: round(100 * float(am[CONS.index(c)][wet].max()), 2) for c in CONS},
                               net_transport_amp_m3_s=tr)
        print(s, json.dumps(res['sides'][s]))
    # implied domain-mean M2 elevation if the three boundaries' M2 transports converge (sanity check)
    A_dom = float((np.diff(g['xg'])[None, :] * 111.32e3 * np.cos(np.deg2rad(g['yc']))[:, None] *
                   (np.diff(g['yg'])[:, None] * 111.32e3) * (hf[0] > 0)).sum())
    res['domain_wet_area_m2'] = A_dom
    res['implied_mean_eta_cm'] = {c: round(100 * float(abs(qin[c])) / (2 * np.pi / PER[c]) / A_dom, 2) for c in CONS}
    print('implied domain-mean elevation amplitude (cm):', res['implied_mean_eta_cm'])
    open(os.path.join(out, 'data.obcs_tides.snippet'), 'w').write('\n'.join(snippet) + '\n')
    json.dump(res, open(os.path.join(out, 'tides_malibu.json'), 'w'), indent=1)
    print('wrote', out)


if __name__ == '__main__':
    main()
