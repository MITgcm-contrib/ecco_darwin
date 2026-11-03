#!/usr/bin/env python3
"""Diatom-rich boundary phytoplankton community for Feb-Apr 2025 (2026-10-25, v4d twin of v4b).

calcofi_bgc.py splits the boundary total Chl over the 5 phytoplankton types with the parent's fixed
domain-mean class fractions (v4b surface Chl: 0.103 diatoms, 0.018 other large eukaryotes,
0.299 Synechococcus, 0 and 0.580 Prochlorococcus), so ~7% of the boundary phytoplankton C is diatoms all
year. CalHABMAP Santa Monica Pier counts show the Mar 10-17 and Apr 7-14 blooms at 89-99% diatoms
(cells, not biomass).

At every boundary point, level and day, with total Chl T = sum Chl1-5 and fractions f_i = Chl_i / T:
  target: diatoms (type 1) FDIA, types 2-5 share 1 - FDIA in their own proportions at that point;
  f_new = (1 - w) f_i + w target_i, w(t) = 1 over Feb 1 - Apr 30, ramped linearly from 0 over the first
  and last RAMP days (as make_pace_obcs_phyto.py), 0 outside;
  Chl_i_new = T f_new_i (total Chl unchanged);
  c_i_new = Chl_i_new / theta_i, theta_i = Chl_i / c_i at that point (each type keeps its own Chl:C);
  where a type has no Chl or C at a point, theta_i is the median of that type, level and boundary.
Grazers (c6, c7), nutrients and all other tracers are unchanged.

usage: make_obcs_community.py <v4b obcs dir (c1-c5, Chl1-Chl5)> <out dir> <json out> [FDIA]
"""
import datetime as dt
import json
import os
import sys
import numpy as np

NR, NY, NX, NREC = 50, 72, 128, 182
T0 = dt.date(2025, 1, 1)
D1, D2, RAMP = dt.date(2025, 2, 1), dt.date(2025, 4, 30), 5


def weight():
    w = np.zeros(NREC)
    k1, k2 = (D1 - T0).days, (D2 - T0).days
    w[k1:k2 + 1] = 1.0
    for i in range(RAMP):
        w[k1 + i] = w[k2 - i] = (i + 1) / (RAMP + 1)
    return w


def main():
    src, out, jout = sys.argv[1:4]
    fdia = float(sys.argv[4]) if len(sys.argv) > 4 else 0.7
    w = weight()
    os.makedirs(out, exist_ok=True)
    res = dict(fdia=fdia, window=[str(D1), str(D2)], ramp_days=RAMP, sides={})
    rd = lambda v, s, n: np.fromfile(os.path.join(src, 'malibu_128x72x50_%s_%s' % (v, s)), '>f4').reshape(NREC, NR, n).astype(float)
    for s, n in (('West', NY), ('East', NY), ('South', NX)):
        chl = np.stack([rd('Chl%d' % i, s, n) for i in range(1, 6)])        # (5, t, k, n)
        c = np.stack([rd('c%d' % i, s, n) for i in range(1, 6)])
        T = chl.sum(0)
        f = chl / np.maximum(T, 1e-30)
        th = chl / np.where(c > 0, c, np.nan)
        thm = np.nanmedian(np.where((chl > 0) & (c > 0), th, np.nan), axis=(1, 3), keepdims=True)  # (5, 1, k, 1)
        thm = np.where(np.isfinite(thm), thm, np.nanmedian(thm))
        th = np.where(np.isfinite(th) & (th > 0), th, thm)
        rest = f[1:].sum(0)
        tgt = np.empty_like(f)
        tgt[0] = fdia
        tgt[1:] = (1 - fdia) * np.where(rest > 0, f[1:] / np.maximum(rest, 1e-30), 0.25)
        ww = w[None, :, None, None]
        fn = (1 - ww) * f + ww * tgt
        chl_n = T[None] * fn
        c_n = np.where(chl_n > 0, chl_n / th, 0.0)
        for i in range(5):
            chl_n[i].astype('>f4').tofile(os.path.join(out, 'malibu_128x72x50_Chl%d_%s' % (i + 1, s)))
            c_n[i].astype('>f4').tofile(os.path.join(out, 'malibu_128x72x50_c%d_%s' % (i + 1, s)))
        k = (dt.date(2025, 3, 15) - T0).days
        sel = T[k, 0] > 0
        res['sides'][s] = dict(
            chl_total_maxabsdiff=float(np.abs(chl_n.sum(0) - T).max()),
            mar15_surface_chl_frac_old=[round(float(np.median(f[i, k, 0][sel])), 3) for i in range(5)],
            mar15_surface_chl_frac_new=[round(float(np.median(fn[i, k, 0][sel])), 3) for i in range(5)],
            mar15_surface_phytoC_old=round(float(np.median(c[:, k, 0].sum(0)[sel])), 3),
            mar15_surface_phytoC_new=round(float(np.median(c_n[:, k, 0].sum(0)[sel])), 3),
            jan15_unchanged=bool(np.allclose(c_n[:, 14], c[:, 14]) and np.allclose(chl_n[:, 14], chl[:, 14])))
        print(s, res['sides'][s])
    json.dump(res, open(jout, 'w'), indent=1)


if __name__ == '__main__':
    main()
