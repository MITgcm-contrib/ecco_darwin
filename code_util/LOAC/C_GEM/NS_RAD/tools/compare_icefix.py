"""
Before/after for the bottom-fast ice grounding fix, all four rivers.

  before  runs/multichannel_preicefix/<site>   ice truncated to the LIVE water depth every
                                               timestep -- a one-way ratchet, because the
                                               clipped ice was destroyed and could only
                                               return by slow conduction
  after   runs/definitive/<site>               depth limits how much NEW ice can form, but
                                               grounded ice keeps its thickness

Both sides carry the same (multi-channel) geometry, so this isolates the ice change.

THE SYMPTOM this fixes: DEPTH swings ~0.8 m per day here on tide + storm surge --
comparable to the entire water column on the three shallow rivers -- so the old cap chopped
the cover twice a day. Domain-mean thickness fell through midwinter instead of growing
(Sagavanirktok 0.81 -> 0.53 m over days 10-110), reversing direction on 35 of 99 days, and
98% of all midwinter thinning events sat exactly at the cap.

Usage: python tools/compare_icefix.py
"""
import sys
from pathlib import Path

import numpy as np
from netCDF4 import Dataset

ROOT = Path(__file__).resolve().parent.parent
SITES = ["colville", "kuparuk", "sagavanirktok", "canning"]
BEFORE = ROOT / "runs/multichannel_preicefix"
AFTER = ROOT / "runs/definitive"

Y2 = 365 * 86400.0
WINTER = (10, 110)      # deep winter of year 2: should be monotonic conductive growth


def ice_stats(path):
    d = Dataset(path)
    t = d.variables["time"][:]
    y2 = t >= Y2
    tt = (t[y2] - Y2) / 86400.0
    # column 0 is the seaward boundary and is fill/NaN in every run -- drop it
    th = np.asarray(d.variables["ice_thickness"][y2, 1:])
    H = np.asarray(d.variables["depth"][y2, 1:])
    fr = np.asarray(d.variables["ice_frac"][y2, 1:])

    w = (tt >= WINTER[0]) & (tt <= WINTER[1])
    x, Hw, tw = th[w], H[w], tt[w]
    dx = np.diff(x, axis=0)
    dec = dx < -1e-9
    capped = x[1:] >= Hw[1:] - 1e-9

    day = np.floor(tw).astype(int)
    dm = np.array([x[day == k].mean() for k in sorted(set(day))])

    return {
        "thin_pct": 100 * dec.mean(),
        "thin_at_cap": 100 * (dec & capped).sum() / max(dec.sum(), 1),
        "reversals": int(np.sum(np.diff(dm) < 0)),
        "n_days": len(dm) - 1,
        "h_start": dm[0], "h_end": dm[-1],
        "h_peak": float(th.max()),
        "ice_destroyed": float(-dx[dec].sum() / x.shape[1]),
        # season timing must NOT move: ice_frac should still flip twice per cell per year
        "frac_flips": int(np.abs(np.diff(fr, axis=0)).sum()),
        "n_cells": x.shape[1],
        "open_frac": float((fr < 0.5).mean()),
    }


def main():
    rows = []
    for s in SITES:
        b_p, a_p = BEFORE / s / "output.nc", AFTER / s / "output.nc"
        if not (b_p.exists() and a_p.exists()):
            print(f"skip {s}: missing {'before' if not b_p.exists() else 'after'} run")
            continue
        rows.append((s, ice_stats(b_p), ice_stats(a_p)))
    if not rows:
        raise SystemExit("no runs to compare")

    hdr = (f"{'river':15s}{'thinning steps %':>19s}{'day-to-day reversals':>23s}"
           f"{'winter growth (m)':>26s}{'peak ice':>17s}")
    sub = (f"{'':15s}{'before   after':>19s}{'before      after':>23s}"
           f"{'before        after':>26s}{'before  after':>17s}")
    print("\nBottom-fast ice grounding fix — deep winter (days 10-110), year 2\n" + "=" * len(hdr))
    print(hdr); print(sub); print("-" * len(hdr))
    for s, b, a in rows:
        print(f"{s:15s}"
              f"{b['thin_pct']:9.2f}{a['thin_pct']:10.2f}"
              f"{b['reversals']:14d}/{b['n_days']:<3d}{a['reversals']:5d}/{a['n_days']:<3d}"
              f"{b['h_start']:11.2f}->{b['h_end']:.2f}{a['h_start']:9.2f}->{a['h_end']:.2f}"
              f"{b['h_peak']:11.2f}{a['h_peak']:7.2f}")

    print("\nIce destroyed by the cap over days 10-110 (m per cell):")
    for s, b, a in rows:
        print(f"   {s:15s} {b['ice_destroyed']:6.2f}  ->  {a['ice_destroyed']:.2f}")

    # The test is that the fix did not MOVE the season -- i.e. after == before. The nominal
    # 2 flips/cell (one break-up, one freeze-up) is a separate, informational check: Canning
    # legitimately has one cell that freezes, briefly reopens in a late-October warm spell,
    # and refreezes (days 277.2 / 279.5 / 281.8), giving it 274. That is physical, and
    # identical in both runs -- so it must not be reported as a regression.
    print("\nSeason timing must NOT move (test: after == before):")
    for s, b, a in rows:
        nominal = 2 * b["n_cells"]
        ok = "OK" if a["frac_flips"] == b["frac_flips"] else "CHANGED"
        note = "" if a["frac_flips"] == nominal else f", {a['frac_flips']-nominal:+d} vs nominal {nominal} (freeze-thaw-freeze cell)"
        print(f"   {s:15s} ice_frac flips {b['frac_flips']} -> {a['frac_flips']} [{ok}{note}]"
              f"   open-water fraction {100*b['open_frac']:.2f}% -> {100*a['open_frac']:.2f}%")


if __name__ == "__main__":
    main()
