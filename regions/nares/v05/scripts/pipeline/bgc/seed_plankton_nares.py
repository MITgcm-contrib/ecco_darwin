#!/usr/bin/env python3
"""NARES COPY (2026-10-29) of west_AO prep/seed_plankton_westao.py: Nares grid via clip_bgc_negatives_nares,
and OB mode seeds ALL FOUR edges (Dustin, 2026-10-21 Nares note: IC and all OBs, esp. the Lincoln Sea east edge).
Seed plankton on the WestAO NORTH open boundary (2026-10-21, Dustin "1" via the Session watch; pan_AO method,
seed_plankton_panao.py). The V5r1 parent has c01-c07 ~0 under the central pack and the clip set them to exactly 0; darwin3
a_Xmin only stops mortality below Xmin, so water entering at 0 can never bloom (west_AO DECISIONS 2026-10-21).
Rule, upper 100 m only (levels whose top is < 100 m), north edge only, every record:
  c_k  -> Xmin_k where c_k < 1e-5 (option B, Dustin 2026-10-21; v2 raised every c_k < Xmin)   Xmin from darwin_traits.txt: c01-c05 0.0125, c06-c07 0.001 mmol C m-3
  Chl_k (k = 1..5) += dc_k * r_k    r_k = median Chl_k/c_k (mg Chl / mmol C) over all four edges of the same segment,
                                    upper 100 m, cells with c_k >= Xmin_k, every 48th record (printed; CHL2CMIN_k if no cell)
  Added C, N, P, Fe, Si are taken from DIC; NO3 -> NH4 -> NO2; PO4; FeT; SiO2 with the clip_bgc_negatives.py quotas and rules
  (pools floored at 0; any shortfall printed as unconserved_<element>).
Never edits its input: <dst> gets new *_north.bin copies of the changed tracers and symlinks to every other file of <src>.
Usage: seed_plankton_westao.py <src OB dir, e.g. OBCSc_seg_TAG> <dst dir, e.g. OBCSs_seg_TAG>
       seed_plankton_westao.py ic <src pickup dir> <dst dir> [iter=0000604872]   (2026-10-24, Dustin "Seed A": same rule on the
           initial condition, every wet cell of levels 1-KTOP; Chl:C medians from the IC itself; float64 pickups 90x680x860)
       seed_plankton_westao.py selftest
"""
import os, sys
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import clip_bgc_negatives_nares as C

DELR = np.array([1.00, 1.14, 1.30, 1.49, 1.70, 1.93, 2.20, 2.50, 2.84, 3.21, 3.63, 4.10, 4.61, 5.18, 5.79, 6.47, 7.20, 7.98, 8.83,
                 9.73, 10.69, 11.70, 12.76, 13.87, 15.03, 16.22, 17.45, 18.70, 19.97, 21.27, 22.56, 23.87])   # data delR (first 32)
KTOP = int((np.r_[0, np.cumsum(DELR)][:-1] < 100.).sum())           # levels with top < 100 m
XMIN = [0.0125] * 5 + [0.001] * 2
CSEED = 1e-5                                                  # 2026-10-21 Dustin "seed all B": seed only where c < 1e-5
CHL2CMIN = [1.738589830975993e-2, 1.658495418381498e-2, 2.619354007654867e-2, 1.456366814713903e-2, None]
PH = [f'c{k + 1:02d}' for k in range(7)]
CHL = [f'Chl{k + 1:02d}' for k in range(5)]
NAMES = PH + CHL + C.POOLS
stats = {}


def seed(f, ratio):
    """f: name -> array (nrec, KTOP, n), modified in place."""
    dc = []
    wet = f['DIC'] > 0                                         # v2 2026-10-21: dry OB cells (all tracers 0) are not seeded
    for k, nm in enumerate(PH):
        a = f[nm]
        add = np.where(wet & (a < CSEED), XMIN[k] - a, 0.)      # v3 (option B): only near-dead cells, c < 1e-5
        a += add; dc.append(add)
        s = stats.setdefault(nm, [0, 0.]); s[0] += int((add > 0).sum()); s[1] += float(add.sum())
    for k, nm in enumerate(CHL):
        f[nm] += dc[k] * ratio[k]
        s = stats.setdefault(nm, [0, 0.]); s[0] += int((dc[k] > 0).sum()); s[1] += float((dc[k] * ratio[k]).sum())
    comp = {'C': sum(dc), 'N': C.R_NC * sum(dc), 'P': C.R_PC * sum(dc), 'Fe': C.R_FEC * sum(dc),
            'Si': sum(r * c for r, c in zip(C.R_SIC, dc))}
    for el, pools in C.TAKE.items():
        need = comp[el].copy()
        for p in pools:
            take = np.minimum(need, np.maximum(f[p], 0.))
            f[p] -= take; need -= take
            s = stats.setdefault(p, [0, 0.]); s[0] += int((take > 0).sum()); s[1] += float(take.sum())
        s = stats.setdefault('unconserved_' + el, [0, 0.]); s[0] += int((need > 0).sum()); s[1] += float(need.sum())


def chl_ratio(src):
    rat = []
    for k in range(5):
        vals = []
        for edge, n in [('north', C.NX), ('south', C.NX), ('east', C.NY), ('west', C.NY)]:
            pc = os.path.realpath(os.path.join(src, f'{C.OB[PH[k]]}_{edge}.bin'))
            pq = os.path.realpath(os.path.join(src, f'{C.OB[CHL[k]]}_{edge}.bin'))
            nrec = os.path.getsize(pc) // (4 * C.NR * n)
            c = np.memmap(pc, '>f4', 'r', shape=(nrec, C.NR, n))[::48, :KTOP]
            q = np.memmap(pq, '>f4', 'r', shape=(nrec, C.NR, n))[::48, :KTOP]
            m = c >= XMIN[k]
            vals.append((np.asarray(q)[m] / np.asarray(c)[m]))
        v = np.concatenate(vals)
        r = float(np.median(v)) if v.size else (CHL2CMIN[k] or 0.02)
        print(f'  Chl{k + 1:02d}/c{k + 1:02d} median {r:.4f} mg/mmol C from {v.size} cells', flush=True)
        rat.append(r)
    return rat


def run(src, dst):
    if os.path.exists(dst) and os.listdir(dst):
        sys.exit(f'FATAL: {dst} exists and is not empty')
    os.makedirs(dst, exist_ok=True)
    print(f'upper 100 m = levels 1..{KTOP}', flush=True)
    ratio = chl_ratio(src)
    changed = set()
    for edge, n in [('north', C.NX), ('south', C.NX), ('east', C.NY), ('west', C.NY)]:
        names = {nm: f'{C.OB[nm]}_{edge}.bin' for nm in NAMES}
        paths = {nm: os.path.realpath(os.path.join(src, fn)) for nm, fn in names.items()}
        nrec = os.path.getsize(paths['DIC']) // (4 * C.NR * n)
        assert all(os.path.getsize(p) == nrec * 4 * C.NR * n for p in paths.values()), edge
        inp = {nm: np.memmap(p, '>f4', 'r', shape=(nrec, C.NR, n)) for nm, p in paths.items()}
        fh = {nm: open(os.path.join(dst, fn), 'wb') for nm, fn in names.items()}
        for r0 in range(0, nrec, 240):
            r1 = min(nrec, r0 + 240)
            f = {nm: np.array(inp[nm][r0:r1], 'f8') for nm in names}
            top = {nm: f[nm][:, :KTOP] for nm in names}            # views: seeding edits f in place
            seed(top, ratio)
            for nm in names:
                f[nm].astype('>f4').tofile(fh[nm])
        for h in fh.values():
            h.close()
        changed |= set(names.values())
        print(f'{edge}: {nrec} records', flush=True)
    C.link_rest(src, dst, changed)
    for nm in PH + CHL + C.POOLS + ['unconserved_' + e for e in C.TAKE]:
        if nm in stats:
            print(f'{nm:16s} {stats[nm][0]:12d} {stats[nm][1]:14.4e}')
    print('SEED_DONE')


def run_ic(src, dst, it='0000604872'):
    if os.path.exists(dst) and os.listdir(dst):
        sys.exit(f'FATAL: {dst} exists and is not empty')
    os.makedirs(dst, exist_ok=True)
    names = {nm: f'pickup_{C.IC[nm]}.{it}.data' for nm in NAMES}
    shp = (C.NR, C.NY, C.NX)
    inp = {nm: np.memmap(os.path.join(src, fn), '>f8', 'r', shape=shp) for nm, fn in names.items()}
    ratio = []
    for k in range(5):
        c = np.asarray(inp[PH[k]][:KTOP]); q = np.asarray(inp[CHL[k]][:KTOP]); m = c >= XMIN[k]
        r = float(np.median(q[m] / c[m])) if m.any() else (CHL2CMIN[k] or 0.02)
        print(f'  IC Chl{k + 1:02d}/c{k + 1:02d} median {r:.4f} mg/mmol C from {int(m.sum())} cells', flush=True); ratio.append(r)
    f = {nm: np.array(inp[nm][:KTOP], 'f8') for nm in names}
    seed(f, ratio)
    for nm, fn in names.items():
        a = np.array(inp[nm], 'f8'); a[:KTOP] = f[nm]
        assert np.isfinite(a).all(), nm
        a.astype('>f8').tofile(os.path.join(dst, fn))
        m = os.path.join(src, fn.replace('.data', '.meta'))
        if os.path.exists(m):
            os.symlink(os.path.realpath(m), os.path.join(dst, fn.replace('.data', '.meta')))
    C.link_rest(src, dst, set(names.values()) | {v.replace('.data', '.meta') for v in names.values()})
    print('IC levels 1..%d seeded' % KTOP)
    for nm in PH + CHL + C.POOLS + ['unconserved_' + e for e in C.TAKE]:
        if nm in stats:
            print(f'{nm:16s} {stats[nm][0]:12d} {stats[nm][1]:14.4e}')
    print('SEED_DONE')


def selftest():
    rng = np.random.default_rng(0)
    f = {nm: rng.uniform(0, 0.02, (5, KTOP, 7)) * (rng.uniform(size=(5, KTOP, 7)) > 0.5) for nm in PH + CHL}   # half exactly 0
    f.update({p: rng.uniform(0, 5, (5, KTOP, 7)) for p in C.POOLS})
    f['NO3'][0] = 0.                                           # force the NO3 -> NH4 -> NO2 cascade
    for nm in f:
        f[nm][:, :, -1] = 0.                                   # one dry column: must stay 0
    q = {'C': (1, 0, 0, 0, 0), 'N': (C.R_NC, 0, 0, 0, 0), 'P': (C.R_PC,), 'Fe': (C.R_FEC,)}
    def tot(f):
        ph = sum(f[nm] for nm in PH)
        return {'C': (ph + f['DIC']).sum(), 'N': (C.R_NC * ph + f['NO3'] + f['NH4'] + f['NO2']).sum(),
                'P': (C.R_PC * ph + f['PO4']).sum(), 'Fe': (C.R_FEC * ph + f['FeT']).sum(),
                'Si': (sum(r * f[nm] for r, nm in zip(C.R_SIC, PH)) + f['SiO2']).sum()}
    orig = {nm: f[nm].copy() for nm in PH}
    t0 = tot(f); seed(f, [0.1] * 5); t1 = tot(f)
    for el in t0:
        unc = stats.get('unconserved_' + el, [0, 0.])[1]
        assert abs(t1[el] - t0[el] - unc) < 1e-9 * max(1, abs(t0[el])), (el, t1[el] - t0[el], unc)
    for k, nm in enumerate(PH):
        lo = orig[nm][:, :, :-1] < CSEED
        assert (f[nm][:, :, :-1][lo] >= XMIN[k] - 1e-15).all() and (f[nm][:, :, :-1][~lo] == orig[nm][:, :, :-1][~lo]).all(), nm
    assert all((f[nm][:, :, -1] == 0).all() for nm in f), 'dry column seeded'
    assert all((f[p] >= 0).all() for p in C.POOLS)
    print('selftest OK: budget change == reported unconserved for C/N/P/Fe/Si; wet c_k < 1e-5 -> Xmin, others unchanged; dry cells untouched; pools >= 0')


if __name__ == '__main__':
    if sys.argv[1] == 'selftest':
        selftest()
    elif sys.argv[1] == 'ic':
        run_ic(sys.argv[2], sys.argv[3], *(sys.argv[4:5]))
    else:
        run(sys.argv[1], sys.argv[2])
