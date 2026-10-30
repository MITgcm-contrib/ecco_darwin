#!/usr/bin/env python3
"""
Closed 0-H m N and P budgets for a 1-D column run with the 'budget' diagnostics list
(ADVr, DFrE, DFrI, Forc, gDAR for every N and P tracer, plus Nfix/Denit).

  python3 scripts/bats_n_budget.py runs/budget/BATS [H=207]

Budget for the box from the surface to the interface nearest H (mmol m-2 yr-1):
  d(inventory)/dt = ADVr(H) + DFrE(H) + DFrI(H)       transport across the bottom face
                    + sum Forc                        RBCS relaxation
                    + sum gDAR                        Darwin (Nfix - denit - sinking out)
MITgcm stores r-fluxes at the top face of each level, positive upward (r = z), so an
upward flux at the box's bottom face is a gain for the box. Plankton N, P = C x R_NC, R_PC from darwin_traits.txt.
"""
import sys, re, numpy as np, xarray as xr

run = sys.argv[1]
H = float(sys.argv[2]) if len(sys.argv) > 2 else 207.
b = xr.open_dataset(f'{run}/mnc_out/budget.0000000000.t001.nc')
g = xr.open_dataset(f'{run}/mnc_out/grid.t001.nc')
d3 = xr.open_dataset(f'{run}/mnc_out/diags3D.0000000000.t001.nc')
dz = np.squeeze(g['drF'].values); zf = -np.squeeze(g['RF'].values)
rA = float(np.squeeze(g['rA'].values))
kb = int(np.argmin(abs(zf - H)))            # interface index; box = levels 0..kb-1
H = zf[kb]
tr = open(f'{run}/darwin_traits.txt').read()
rnc = float(re.search(r'R_NC=\s*\d+\*([0-9.Ee+-]+)', tr).group(1))
rpc = float(re.search(r'R_PC=\s*\d+\*([0-9.Ee+-]+)', tr).group(1))
PL = ['20', '21', '22', '23', '24', '25', '26']
SP = {'N': [('02', 1), ('03', 1), ('04', 1), ('09', 1), ('13', 1)] + [(p, rnc) for p in PL],
      'P': [('05', 1), ('10', 1), ('14', 1)] + [(p, rpc) for p in PL]}
YR = 365.25 * 86400
t = b['T'].values; yr = t.astype('datetime64[Y]').astype(int) + 1970
t3 = d3['T'].values; yr3 = t3.astype('datetime64[Y]').astype(int) + 1970

def v(name):
    return np.squeeze(b[name].values)

for el, sp in SP.items():
    terms = {}
    adv = dfe = dfi = frc = gd = 0.
    for s, f in sp:
        # flux across interface kb: MITgcm stores the flux at the top face of level k,
        # so the bottom face of the box is the top face of level kb
        adv = adv + f * v(f'ADVrTr{s}')[:, kb] / rA
        dfe = dfe + f * v(f'DFrETr{s}')[:, kb] / rA
        dfi = dfi + f * v(f'DFrITr{s}')[:, kb] / rA
        frc = frc + f * (v(f'ForcTr{s}')[:, :kb] * dz[:kb]).sum(1)
        gd = gd + f * (v(f'gDAR{s}')[:, :kb] * dz[:kb]).sum(1)
    terms['vert. advection in'] = adv * YR
    terms['diffusion in (expl.)'] = dfe * YR
    terms['diffusion in (impl.)'] = dfi * YR
    terms['RBCS relaxation'] = frc * YR
    terms['Darwin net (gDAR)'] = gd * YR
    if el == 'N':
        terms['  of which Nfix'] = (v('Nfix')[:, :kb] * dz[:kb]).sum(1) * YR
        terms['  of which denit'] = -(v('DenitN')[:, :kb] * dz[:kb]).sum(1) * YR
    inv = 0.
    for s, f in sp:
        inv = inv + f * (np.squeeze(d3[f'TRAC{s}'].values)[:, :kb] * dz[:kb]).sum(1)
    print(f'\n{el} budget, 0-{H:.0f} m (mmol {el} m-2 yr-1)')
    print(f'{"term":24s} {"1998-2010":>10s} {"2011-2023":>10s}')
    tot = {}
    for a, c in [(1998, 2010), (2011, 2023)]:
        m = (yr >= a) & (yr <= c)
        tot[(a, c)] = {k: x[m].mean() for k, x in terms.items()}
        # inventory change = (last daily mean of year c - last daily mean of year a-1) / n
        i0 = np.where(yr3 == a - 1)[0][-1]; i1 = np.where(yr3 == c)[0][-1]
        tot[(a, c)]['inventory change'] = (inv[i1] - inv[i0]) / (c - a + 1)
    keys = [k for k in terms] + ['inventory change']
    for k in keys:
        print(f'{k:24s} {tot[(1998,2010)][k]:10.2f} {tot[(2011,2023)][k]:10.2f}')
    for p in tot:
        s = sum(val for k, val in tot[p].items() if not k.startswith('  ') and k != 'inventory change')
        print(f'sum of sources {p}: {s:8.2f}  residual (inventory - sources): {tot[p]["inventory change"] - s:8.2f}')
