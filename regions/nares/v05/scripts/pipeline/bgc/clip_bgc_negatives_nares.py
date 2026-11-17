#!/usr/bin/env python3
"""NARES COPY (2026-10-29) of west_AO prep/clip_bgc_negatives.py: only NR/NY/NX changed.
Clip inherited negative BGC tracers at 0 and compensate the added mass (WestAO production restart, 2026-10-18, Dustin
"restart with fixes"; same rule as the pan_AO Jan 1 restart plan).

Clipped at 0: DOC, rDOC, DON, plankton c01-c07, Chl01-Chl05 (negatives inherited from the V5r1 parent; DECISIONS 2026-10-18).
Added mass is removed from the inorganic pools in the same cell/record, so C, N, P, Fe, Si stay conserved:
  dC  = dDOC + drDOC + sum_k dc_k                 -> DIC  -= dC
  dN  = dDON + sum_k R_NC_k dc_k                  -> NO3, then NH4, then NO2 (each floored at 0; 2026-10-18 v2:
                                                     v1 took all N from NO3 and left NO3 < 0 in H1 OB files)
  dP  = sum_k R_PC_k dc_k                         -> PO4  -= dP
  dFe = sum_k R_FeC_k dc_k                        -> FeT  -= dFe
  dSi = sum_k R_SiC_k dc_k                        -> SiO2 -= dSi
  Every pool is floored at 0; any amount that cannot be taken is printed as unconserved_<element>.
  Chl carries no mass. Quotas from run_prod_2015H1/darwin_traits.txt: R_NC 0.1333 (all 7), R_PC 1/120, R_FeC 8.333e-6,
  R_SiC 0.2667 for groups 3-4 (diatoms), 0 otherwise.
Files (model tracer -> IC pickup / OBCS file stem; parent index, rDOC inserted as tracer 9, data.ptracers / data.obcs):
  DIC pTr01/TRAC01, NO3 pTr02/TRAC02, PO4 pTr05/TRAC05, FeT pTr06/TRAC06, SiO2 pTr07/TRAC07, DOC pTr08/TRAC08, rDOC rDOC/rDOC,
  DON pTr09/TRAC09, c01-c07 pTr20-26/TRAC20-26, Chl01-05 pTr27-31/TRAC27-31.
Modes (never edits its inputs; the output dir gets symlinks to every unchanged file and new copies of the changed ones):
  ic   <src pickup dir> <dst dir> [iter=0000604872]   float64 (90, 680, 860) per file
  obcs <src OB dir> <dst dir>                         float32 (nrec, 90, n) per <stem>_<edge>.bin
Prints per-tracer clip counts/amounts and the compensation totals, and WARNS if compensation leaves a pool below -1e-6
(reported, not silently floored).
"""
import os, sys
import numpy as np

NR, NY, NX = 90, 360, 1120   # Nares copy 2026-10-29 (west_AO: 90, 680, 860)
R_NC, R_PC, R_FEC = 0.133333333333333, 8.333333333333333e-3, 8.333333333333334e-6
R_SIC = [0., 0., 0.266666666666667, 0.266666666666667, 0., 0., 0.]
IC = {'DIC': 'pTr01', 'NO3': 'pTr02', 'NO2': 'pTr03', 'NH4': 'pTr04', 'PO4': 'pTr05', 'FeT': 'pTr06', 'SiO2': 'pTr07', 'DOC': 'pTr08', 'rDOC': 'rDOC',
      'DON': 'pTr09', **{f'c{k + 1:02d}': f'pTr{20 + k:02d}' for k in range(7)},
      **{f'Chl{k + 1:02d}': f'pTr{27 + k:02d}' for k in range(5)}}
OB = {'DIC': 'TRAC01', 'NO3': 'TRAC02', 'NO2': 'TRAC03', 'NH4': 'TRAC04', 'PO4': 'TRAC05', 'FeT': 'TRAC06', 'SiO2': 'TRAC07', 'DOC': 'TRAC08', 'rDOC': 'rDOC',
      'DON': 'TRAC09', **{f'c{k + 1:02d}': f'TRAC{20 + k:02d}' for k in range(7)},
      **{f'Chl{k + 1:02d}': f'TRAC{27 + k:02d}' for k in range(5)}}
CLIP = ['DOC', 'rDOC', 'DON'] + [f'c{k + 1:02d}' for k in range(7)] + [f'Chl{k + 1:02d}' for k in range(5)]
TAKE = {'C': ['DIC'], 'N': ['NO3', 'NH4', 'NO2'], 'P': ['PO4'], 'Fe': ['FeT'], 'Si': ['SiO2']}
POOLS = ['DIC', 'NO3', 'NH4', 'NO2', 'PO4', 'FeT', 'SiO2']
stats = {}


def process(f):
    """f: dict name -> array (modified in place). Returns nothing; accumulates stats."""
    d = {}
    for nm in CLIP:
        a = f[nm]
        neg = a < 0
        d[nm] = np.where(neg, -a, 0.)
        s = stats.setdefault(nm, [0, 0., 0.])
        s[0] += int(neg.sum()); s[1] += float(d[nm].sum()); s[2] = min(s[2], float(a.min()))
        a[neg] = 0.
    dc = [d[f'c{k + 1:02d}'] for k in range(7)]
    comp = {'C': d['DOC'] + d['rDOC'] + sum(dc), 'N': d['DON'] + R_NC * sum(dc), 'P': R_PC * sum(dc),
            'Fe': R_FEC * sum(dc), 'Si': sum(r * c for r, c in zip(R_SIC, dc))}
    # take each element from its pools in order, never below 0; whatever cannot be taken is reported as unconserved
    for el, pools in TAKE.items():
        need = comp[el].copy()
        for p in pools:
            take = np.minimum(need, np.maximum(f[p], 0.))
            f[p] -= take; need -= take
            s = stats.setdefault(p, [0, 0., 0.])
            s[0] += int((take > 0).sum()); s[1] += float(take.sum()); s[2] = min(s[2], float(f[p].min()))
        s = stats.setdefault('unconserved_' + el, [0, 0., 0.])
        s[0] += int((need > 0).sum()); s[1] += float(need.sum())


def link_rest(src, dst, changed):
    os.makedirs(dst, exist_ok=True)
    for fn in sorted(os.listdir(src)):
        if fn in changed or os.path.exists(os.path.join(dst, fn)):
            continue
        os.symlink(os.path.realpath(os.path.join(src, fn)), os.path.join(dst, fn))


def mode_ic(src, dst, it):
    names = {nm: f'pickup_{st}.{it}.data' for nm, st in IC.items()}
    os.makedirs(dst, exist_ok=True)
    inp = {nm: np.memmap(os.path.join(src, fn), '>f8', 'r', shape=(NR, NY, NX)) for nm, fn in names.items()}
    out = {nm: np.memmap(os.path.join(dst, fn), '>f8', 'w+', shape=(NR, NY, NX)) for nm, fn in names.items()}
    for k in range(NR):
        f = {nm: np.array(inp[nm][k], 'f8') for nm in names}
        process(f)
        for nm in names:
            out[nm][k] = f[nm]
    for nm in names:
        out[nm].flush()
    for nm, fn in names.items():                       # the .meta files are unchanged: link them
        m = fn.replace('.data', '.meta')
        if os.path.exists(os.path.join(src, m)):
            os.symlink(os.path.realpath(os.path.join(src, m)), os.path.join(dst, m))
    link_rest(src, dst, set(names.values()))


def mode_obcs(src, dst):
    os.makedirs(dst, exist_ok=True)
    changed = set()
    for edge, n in [('north', NX), ('south', NX), ('east', NY), ('west', NY)]:
        names = {nm: f'{st}_{edge}.bin' for nm, st in OB.items()}
        paths = {nm: os.path.realpath(os.path.join(src, fn)) for nm, fn in names.items()}
        nrec = os.path.getsize(paths['DOC']) // (4 * NR * n)
        assert all(os.path.getsize(p) == nrec * 4 * NR * n for p in paths.values()), edge
        inp = {nm: np.memmap(p, '>f4', 'r', shape=(nrec, NR, n)) for nm, p in paths.items()}
        fh = {nm: open(os.path.join(dst, fn), 'wb') for nm, fn in names.items()}
        for r0 in range(0, nrec, 240):
            r1 = min(nrec, r0 + 240)
            f = {nm: np.array(inp[nm][r0:r1], 'f8') for nm in names}
            process(f)
            for nm in names:
                f[nm].astype('>f4').tofile(fh[nm])
        for h in fh.values():
            h.close()
        changed |= set(names.values())
        print(f'  {edge}: {nrec} records', flush=True)
    link_rest(src, dst, changed)


if __name__ == '__main__':
    mode = sys.argv[1]
    if os.path.exists(sys.argv[3]) and os.listdir(sys.argv[3]):
        sys.exit(f'FATAL: {sys.argv[3]} exists and is not empty')
    if mode == 'ic':
        mode_ic(sys.argv[2], sys.argv[3], sys.argv[4] if len(sys.argv) > 4 else '0000604872')
    elif mode == 'obcs':
        mode_obcs(sys.argv[2], sys.argv[3])
    else:
        sys.exit('mode: ic | obcs')
    print(f'{"tracer":6s} {"n changed":>12s} {"sum added/removed":>18s} {"min (before clip / pool after)":>30s}')
    for nm in CLIP + POOLS + ['unconserved_' + e for e in TAKE]:
        if nm in stats:
            s = stats[nm]; print(f'{nm:6s} {s[0]:12d} {s[1]:18.4e} {s[2]:30.4e}')
    bad = [p for p in POOLS if p in stats and stats[p][2] < -1e-6]
    print('WARNING: pool(s) below 0 after compensation: ' + ', '.join(bad) if bad else 'pools >= 0 (within 1e-6) after compensation')
    print('CLIP_DONE')
