#!/usr/bin/env python3
"""Health check for a Nares arch run (run on a compute node, not pfe).

1. Along-strait volume transport through x-sections at several i, from the newest pickup
   (U * hFacW * dyG * drF; rStar scaling of hFacW ignored, |eta| << depth). Positive = +i
   (toward the Lincoln Sea, east boundary), so the expected southward throughflow is negative;
   target ~ -0.85 Sv (DECISIONS 2026-10-06).
2. Sea ice from the newest 3-hourly SIdyn record: area-weighted mean SIarea and SIheff, and
   ice speed percentiles (stationary ice upstream of an arch shows up as speed ~ 0).
3. Seaice solver cost from STDOUT.0000: JFNK Newton/Krylov iteration counts and the
   timing table (SEAICE_* and total wall clock).
Usage: python3 check_arch_run.py <run dir>
"""
import glob, os, re, sys
import numpy as np

R = sys.argv[1]
NX, NY = 1120, 360


def meta(f):
    s = open(f).read()
    prec = 'f8' if 'float64' in s else 'f4'
    m = re.search(r"fldList = \{(.*?)\}", s, re.S)
    flds = re.findall(r"'\s*(\w+)\s*'", m.group(1)) if m else []
    return prec, flds


def grid(name):
    a = np.fromfile(f'{R}/{name}.data', '>' + meta(f'{R}/{name}.meta')[0])
    return a.reshape(-1, NY, NX) if a.size > NY * NX else a.reshape(NY, NX)


drf = np.fromfile(f'{R}/DRF.data', '>' + meta(f'{R}/DRF.meta')[0])
NR = drf.size
hfw = grid('hFacW').reshape(NR, NY, NX)
dyg = grid('DYG')

# ---- 1. transport
pk = glob.glob(f'{R}/pickup.0*.data') + glob.glob(f'{R}/pickup.ckpt*.data')
st = sorted(glob.glob(f'{R}/diags/monthly/state.*.data'))
area_u = hfw * dyg[None] * drf[:, None, None]
u = None
if pk:
    p = max(pk, key=os.path.getmtime)
    prec, flds = meta(p[:-5] + '.meta')
    assert flds[0] == 'Uvel', flds[:3]
    u = np.fromfile(p, '>' + prec, count=NR * NY * NX).reshape(NR, NY, NX)
    print(f'U from pickup {os.path.basename(p)}')
elif st:
    p = st[-1]
    prec, flds = meta(p[:-5] + '.meta')
    if 'UVEL' in flds:
        k = flds.index('UVEL')
        u = np.fromfile(p, '>' + prec, count=NR * NY * NX, offset=k * NR * NY * NX * int(prec[1])).reshape(NR, NY, NX)
        print(f'U = monthly-mean UVEL from {os.path.basename(p)}')
if u is None:
    print('no pickup or monthly state yet: transport skipped')
else:
    print(f'  |U| max {np.abs(u).max():.2f} m/s')
    assert np.isfinite(u).all()
    for i in [60, 200, 400, 560, 700, 900, 1050]:
        q = (u[:, :, i] * area_u[:, :, i]).sum() / 1e6
        print(f'  section i={i:4d}: transport {q:+.3f} Sv  (wet area {area_u[:, :, i].sum()/1e6:.2f} km2)')

# ---- 2. sea ice
sid = sorted(glob.glob(f'{R}/diags/3hourly/SIdyn.*.data'))
if sid:
    f = sid[-1]
    prec, flds = meta(f[:-5] + '.meta')
    a = np.fromfile(f, '>' + prec).reshape(len(flds), NY, NX)
    d = dict(zip(flds, a))
    ra = grid('RAC')
    wet = hfw[0] > 0
    A, H = d['SIarea'], d['SIheff']
    print(f'{os.path.basename(f)}: SIarea mean {(A*ra)[wet].sum()/ra[wet].sum():.3f}, '
          f'SIheff mean {(H*ra)[wet].sum()/ra[wet].sum():.2f} m, max {H.max():.2f} m')
    sp = np.hypot(d['SIuice'], d['SIvice'])[(A > 0.15) & wet]
    print('  ice speed (A>0.15) percentiles 5/50/95/max: ' +
          ' '.join(f'{v:.3f}' for v in np.percentile(sp, [5, 50, 95, 100])) + ' m/s; '
          f'fraction < 1 cm/s: {(sp < 0.01).mean():.2f}')
else:
    print('no SIdyn output yet')

# ---- 3. solver cost
s = open(f'{R}/STDOUT.0000').read()
nl = [int(x) for x in re.findall(r'JFNK.*?Nb\. of Newton iterations\s*=?\s*(\d+)', s)]
kr = [int(x) for x in re.findall(r'Nb\. of (?:linear|Krylov|FGMRES) iterations\s*=?\s*(\d+)', s)]
if nl:
    print(f'JFNK Newton its: n={len(nl)} mean {np.mean(nl):.1f} max {max(nl)}')
if kr:
    print(f'JFNK Krylov its: mean {np.mean(kr):.1f} max {max(kr)}')
fail = len(re.findall(r'JFNK.*(?:not converged|did not converge)', s, re.I))
print(f'JFNK non-converged messages: {fail}')
for sec in ['SEAICE_MODEL', 'SEAICE_DYNSOLVER', 'DYNAMICS', 'SOLVE_FOR_PRESSURE', 'ALL']:
    m = re.findall(r'Seconds in section "' + sec + r'.*?Wall clock time:\s*([\d.E+-]+)', s, re.S)
    if m:
        print(f'  wall clock {sec:20s} {float(m[0]):9.1f} s')
print('normal end:', 'Execution ended Normally' in s)
