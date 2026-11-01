"""Scale down the fresh river plume (Gironde/Maumusson water in IBI) on the LR17 south open boundary.

  python3 LR17_plume_bc.py OBSs_IN OBSs_OUT BETA [--nref=100]
OBSs (south-boundary salinity, big-endian float32, records x 14 levels x 936 cells, hourly) as gendata.py writes it.
For each record and level, S_off = mean IBI salinity over the westernmost NREF open cells of the south boundary
(the offshore water); where IBI is fresher than S_off, S_new = S_off + BETA * (S - S_off), elsewhere S is kept.
BETA = 1 returns the input; BETA = 0 removes the fresh anomaly. Dry/closed cells (S = 0) are left at 0.
Why: in the v4.1 43-day run the freshening of the Pertuis in mid-March 2022 entered through the near-coast part of
the south boundary (LR17_fwbudget.py), where IBI surface salinity fell to 21-25 (3-12 below IBI offshore), while the
CMEMS multi-observation SSS shows a near-coast anomaly of only about 0.2-3 (0.125 deg grid).
"""
import sys
import numpy as np
NR, NX = 14, 936


def main(fin, fout, beta, nref=100):
    beta = float(beta)
    s = np.fromfile(fin, '>f4').reshape(-1, NR, NX).astype(np.float64)
    wet = s > 0
    ref = np.where(wet[:, :, :nref], s[:, :, :nref], np.nan)
    soff = np.nanmean(ref, axis=2, keepdims=True)                     # (nrec, NR, 1)
    soff = np.where(np.isfinite(soff), soff, 0.0)
    fresh = wet & (s < soff)
    out = np.where(fresh, soff + beta * (s - soff), s)
    out.astype('>f4').tofile(fout)
    d = (out - s)[wet]
    print('records %d; cells changed %.1f%%; salinity raised by mean %.3f, max %.2f; surface minimum %.2f -> %.2f' % (
        s.shape[0], 100 * fresh.sum() / wet.sum(), d.mean(), d.max(), s[:, 0][wet[:, 0]].min(), out[:, 0][wet[:, 0]].min()))


if __name__ == '__main__':
    a = [x for x in sys.argv[1:] if not x.startswith('--')]
    n = next((int(x.split('=')[1]) for x in sys.argv[1:] if x.startswith('--nref=')), 100)
    main(a[0], a[1], a[2], n)
