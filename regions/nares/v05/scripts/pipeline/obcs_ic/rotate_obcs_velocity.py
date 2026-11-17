#!/usr/bin/env python3
"""Rotate Nares OBCS velocities from geographic east/north into the Nares grid's own axes.

gen_obcs.py (and our gen_obcs_fast.py) rotate parent LLC (u, v) only to geographic east/north
(parent AngleCS/SN; gen_obcs.py:253-267). That is correct for lat-lon children (WestAO, Kerguelen)
but not for Nares, whose grid is rotated 54-83 deg: the east/north components were written as if
they were grid u/v. MITgcm: uE = CS*u - SN*v, vN = SN*u + CS*v, so the inverse is
    u_grid =  CS*uE + SN*vN
    v_grid = -SN*uE + CS*vN
with the CHILD AngleCS/SN at the boundary row/column (cell-centre values; the u and v points of a
boundary are within half a 1-km cell of them).

For each open boundary the existing UVEL_<b>.bin holds uE at that boundary's u points and VVEL_<b>.bin
holds vN at its v points (same (nt, Nr, n) layout), so both are available and are rotated pairwise.
Same for SIuice_<b>.bin / SIvice_<b>.bin (nt, n). Reads forcings/OBCS (never modified), writes rotated
files to forcings/OBCS_rot and symlinks everything else. Prints the grid-normal transport at E/W before
and after, so the throughflow offset can be re-derived from the corrected parent transport.
Usage: rotate_obcs_velocity.py <nares dir>
"""
import glob, os, sys
import numpy as np
import xarray as xr

N = sys.argv[1]
src, dst = f'{N}/forcings/OBCS', f'{N}/forcings/OBCS_rot'
os.makedirs(dst, exist_ok=True)
g = xr.open_dataset(f'{N}/Nares_ncgrid.nc')
CS, SN = g['AngleCS'].values, g['AngleSN'].values          # (Ny, Nx)
NY, NX = CS.shape
drF = g['drF'].values
NR = drF.size
hW, hS, dyG, dxG = g['HFacW'].values, g['HFacS'].values, g['dyG'].values, g['dxG'].values
assert np.isfinite(CS).all() and np.isfinite(SN).all()
assert np.allclose(CS**2 + SN**2, 1, atol=1e-6)

# boundary -> (CS, SN) along it (cell centres of the boundary row/column)
ang = {'east': (CS[:, -1], SN[:, -1]), 'west': (CS[:, 0], SN[:, 0]),
       'north': (CS[-1, :], SN[-1, :]), 'south': (CS[0, :], SN[0, :])}


def rd(name, n, nz):
    a = np.fromfile(f'{src}/{name}', '>f4')
    return a.reshape(-1, nz, n) if nz > 1 else a.reshape(-1, n)


def normal_transport(u, b):
    """Grid-normal volume transport through E/W faces (Sv), positive = +i."""
    i = NX - 1 if b == 'east' else 1        # MITgcm OBCS: east u at face Ieast, west at Iwest+1 (0-based)
    i = min(i, hW.shape[-1] - 1)
    area = hW[:, :, i] * dyG[:, i][None, :] * drF[:, None]
    return (u * area[None]).sum(axis=(1, 2)) / 1e6


for b, (cs, sn) in ang.items():
    n = cs.size
    for un, vn, nz in (('UVEL', 'VVEL', NR), ('SIuice', 'SIvice', 1)):
        fu, fv = f'{un}_{b}.bin', f'{vn}_{b}.bin'
        if not (os.path.exists(f'{src}/{fu}') and os.path.exists(f'{src}/{fv}')):
            print(f'{b}: {un}/{vn} missing, skipped'); continue
        uE, vN = rd(fu, n, nz), rd(fv, n, nz)
        assert uE.shape == vN.shape, (fu, uE.shape, vN.shape)
        c, s = (cs[None, None, :], sn[None, None, :]) if nz > 1 else (cs[None, :], sn[None, :])
        ug = c * uE + s * vN
        vg = -s * uE + c * vN
        # keep exact zeros (land / unset) as zeros
        ug = np.where((uE == 0) & (vN == 0), 0, ug)
        vg = np.where((uE == 0) & (vN == 0), 0, vg)
        for f, a in ((fu, ug), (fv, vg)):
            out = f'{dst}/{f}'
            if os.path.islink(out): os.remove(out)
            a.astype('>f4').tofile(out)
        if un == 'UVEL' and b in ('east', 'west'):
            t0, t1 = normal_transport(uE, b), normal_transport(ug, b)
            print(f'{b}: grid-normal transport (positive = +i, toward Lincoln Sea) '
                  f'as written {t0.mean():+.3f} Sv -> rotated {t1.mean():+.3f} Sv '
                  f'(range {t1.min():+.3f} to {t1.max():+.3f}) over {t1.size} records; '
                  f'angle {np.degrees(np.arctan2(sn, cs)).min():.1f}-{np.degrees(np.arctan2(sn, cs)).max():.1f} deg')
        else:
            print(f'{b}: rotated {un}/{vn} ({uE.shape[0]} records)')
for f in glob.glob(f'{src}/*.bin'):
    t = f'{dst}/{os.path.basename(f)}'
    if not os.path.exists(t):
        os.symlink(f, t)
print('rotation done ->', dst)
