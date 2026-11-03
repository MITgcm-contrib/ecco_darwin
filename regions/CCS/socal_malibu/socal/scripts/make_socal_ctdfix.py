#!/usr/bin/env python3
"""CalCOFI-based temperature correction of the SoCal 2-km nest IC and OBCS (2026-10-19).

The parent (LLC270 V5r1) is too warm in the upper ~60 m vs CalCOFI CTD casts
(socal_ctd_bias.py -> socal_ctd_bias.json):
  IC  (parent pickup 2025-01-01 vs 36 casts of 2501RL, Jan 4-15): shelf +1.3-2.0 C, offshore
      +0.9-1.4 C at 0-40 m, ~0 below 70 m.
  OBCS (hourly parent boundary values at the cast hour vs casts within 60 km of a boundary):
      2501RL (n=15) +0.5-0.9 C and 2504SH (n=14) +0.6-1.0 C at 0-35 m, ~0 below 60 m.
Correction = minus the bias, as a 5-level running median of the per-level median, weighted
1 down to 55 m and tapered linearly to 0 at 75 m (below, the medians are within +-0.15 C and
change sign between cruises).
  IC:   shelf profile where the bottom is < 300 m, offshore profile > 700 m, linear in between;
        applied to the RAW parent pickup (replaces socal_ic_murfix.py).
  OBCS: one profile for all four boundaries, the 2501RL one until 2025-01-10 (cruise center),
        linear in time to the 2504SH one at 2025-04-04, held after. Dry (zero) values untouched.
Writes forcings/pickups/pickup_THETA_ctdfix.0000867888.data (+ .meta copy) and
forcings/OBCS_ctdfix/THETA_{west,east,south,north}.bin, and a summary json.

usage: make_socal_ctdfix.py <socal_ctd_bias.json> <summary.json>
"""
import json, os, shutil, sys
import numpy as np

S = '/nobackup/dcarrol2/downscaling/socal'
R = os.path.join(S, 'regional_run', 'run_v3_prod')
NX, NY, NR = 168, 101, 60


def prof(med, rc):
    p = np.array([np.nan if v is None else v for v in med], float)
    sm = np.array([np.nanmedian(p[max(0, k - 2):k + 3]) for k in range(len(p))])
    w = np.clip((75.0 - rc) / 20.0, 0.0, 1.0)
    return np.nan_to_num(sm) * w


def main():
    bj, out = sys.argv[1:3]
    b = json.load(open(bj))
    rc = -np.fromfile(os.path.join(R, 'RC.data'), '>f4').astype(float)
    dep = np.fromfile(os.path.join(R, 'Depth.data'), '>f4').reshape(NY, NX).astype(float)
    ic_sh, ic_off = prof(b['ic_median_2501_shelf'], rc), prof(b['ic_median_2501_offshore'], rc)
    ob_jan, ob_apr = prof(b['obcs_median_2501_all'], rc), prof(b['obcs_median_2504_all'], rc)
    # IC
    f = np.clip((dep - 300.0) / 400.0, 0.0, 1.0)                      # 0 shelf -> 1 offshore
    dT = (1 - f)[None] * ic_sh[:, None, None] + f[None] * ic_off[:, None, None]
    src = os.path.join(S, 'forcings', 'pickups', 'pickup_THETA.0000867888')
    th = np.fromfile(src + '.data', '>f8').reshape(NR, NY, NX)
    new = np.where(th != 0, th - dT, th)
    dst = os.path.join(S, 'forcings', 'pickups', 'pickup_THETA_ctdfix.0000867888')
    assert not os.path.exists(dst + '.data'), dst
    new.astype('>f8').tofile(dst + '.data')
    shutil.copy(src + '.meta', dst + '.meta')
    summ = dict(levels_m=rc[:35].round(1).tolist(), ic_shelf=ic_sh[:35].round(3).tolist(),
                ic_offshore=ic_off[:35].round(3).tolist(), ob_jan=ob_jan[:35].round(3).tolist(),
                ob_apr=ob_apr[:35].round(3).tolist(),
                ic_surface_mean_change=float((new[0] - th[0])[th[0] != 0].mean()))
    # OBCS
    t0 = np.datetime64('2025-01-01T00:00')
    tj, ta = np.datetime64('2025-01-10T00:00'), np.datetime64('2025-04-04T00:00')
    od = os.path.join(S, 'forcings', 'OBCS_ctdfix')
    os.makedirs(od, exist_ok=True)
    for side, n in (('west', NY), ('east', NY), ('south', NX), ('north', NX)):
        a = np.fromfile(os.path.join(S, 'forcings', 'OBCS', 'THETA_%s.bin' % side), '>f4').reshape(-1, NR, n)
        nrec = a.shape[0]
        t = t0 + np.arange(nrec) * np.timedelta64(1, 'h') + np.timedelta64(30, 'm')
        al = np.clip((t - tj) / (ta - tj), 0.0, 1.0)
        dd = (1 - al)[:, None] * ob_jan[None, :] + al[:, None] * ob_apr[None, :]   # (nrec, NR)
        b2 = np.where(a != 0, a - dd[:, :, None], a).astype('>f4')
        o = os.path.join(od, 'THETA_%s.bin' % side)
        assert not os.path.exists(o), o
        b2.tofile(o)
        summ['ob_%s' % side] = dict(nrec=nrec, mean_change_top=float((b2[:, 0] - a[:, 0])[a[:, 0] != 0].mean()),
                                    min_new=float(b2[b2 != 0].min()), max_new=float(b2.max()))
        print(side, summ['ob_%s' % side])
    json.dump(summ, open(out, 'w'), indent=1)
    print(json.dumps({k: v for k, v in summ.items() if not k.startswith('ob_')}))


if __name__ == '__main__':
    main()
