"""Pass checks for Nares run_bgc_test2d (BGC_PLAN.md section 7). Light: reads STDOUT.0000 and 2 daily diag records."""
import os, re, sys, glob, numpy as np
R = sys.argv[1]; NX, NY, NR = 1120, 360, 90
so = open(f'{R}/STDOUT.0000', errors='replace').read()
print('NORMAL END:', 'Execution ended Normally' in so)
for k in ('DARWIN_ALLOW_RADIv2_EMUL', 'DARWIN_RADI_LIMITER'):
    print(k, 'in STDOUT:', k in so)
def timer(name):
    i = so.find(f'Seconds in section "{name}')
    if i < 0: return None
    m = re.search(r'Wall clock time:\s+([0-9.Ee+-]+)', so[i:i + 600]); return float(m.group(1)) if m else None
nt = 2880
for n in ('ALL ', 'MAIN_DO_LOOP', 'THERMODYNAMICS', 'GCHEM_FORCING_SEP  [F', 'DARWIN_LIGHT', 'DARWIN_PLANKTON', 'DYNAMICS  ', 'SEAICE_MODEL    [', 'BLOCKING_EXCH'):
    t = timer(n); print(f'timer {n:24s}', 'n/a' if t is None else f'{t:9.1f} s  = {t / nt:.3f} s/step')
t = timer('MAIN_DO_LOOP')
if t: print(f'=> {t / nt:.3f} s/step; {86400 / 60 * t / nt / 60:.1f} wall min per model day; {60 / (1440 * t / nt / 60):.2f} model days per wall hour (physics Pe16: 0.334 s/step, 7.5 d/h)')
# monitor: last darwin/ptracer min/max lines
for key in ('trcstat_ptracer01_min', 'trcstat_ptracer02_min', 'trcstat_ptracer07_min', 'trcstat_ptracer20_min', 'trcstat_ptracer22_max', 'dynstat_theta_max', 'seaice_heff_max'):
    v = re.findall(rf'%MON {key}\s+=\s+([-0-9.E+]+)', so)
    if v: print(f'monitor {key:24s} first {v[0]:>14s} last {v[-1]:>14s} (n={len(v)})')
nan = re.findall(r'NaN', so); print('NaN in STDOUT.0000:', len(nan))
err = [l for f in glob.glob(f'{R}/STDERR.0*') for l in open(f, errors='replace') if re.search('ERROR|ABNORMAL|STOP', l)]
print('STDERR error lines:', len(err)); [print('  ', l.rstrip()) for l in err[:5]]
# daily surface BGC, last record
def last(stem, nfld, k):
    fs = sorted(glob.glob(f'{R}/diags/daily/{stem}.*.data'))
    if not fs: return None, None
    a = np.fromfile(fs[-1], '>f4').reshape(nfld, -1, NY, NX); return a[k], fs[-1]
names = ['Chl', 'PP', 'plankC', 'NO3', 'PARF', 'pH', 'pCO2']
A, f = last('bgcsurf', 7, slice(None))
if A is not None:
    print('bgcsurf', os.path.basename(f))
    wet = A[5, 0] > 0
    for i, n in enumerate(names):
        v = A[i, 0][wet]; print(f'  {n:7s} min {v.min():.3g} median {np.median(v):.3g} max {v.max():.3g} frac>0 {np.mean(v > 0):.2f}')
# shelf light mask: PARF and RTtrans must be 0 under the Petermann shelf (all levels)
B, f = last('bgc3d', 3, slice(None))
topo = np.fromfile(f'{R}/Nares_shelficeTopo.bin', '>f8').reshape(NY, NX) if os.path.exists(f'{R}/Nares_shelficeTopo.bin') else None
if B is not None and topo is not None:
    cav = topo < 0
    print('shelf cavity columns:', int(cav.sum()))
    for i, n in enumerate(['PARF', 'RTtrans', 'PP']):
        v = B[i][:, cav]; print(f'  {n:7s} under shelf: max {np.abs(v).max():.3g}  (open water max {np.abs(B[i][:, ~cav]).max():.3g})')
