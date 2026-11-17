#!/usr/bin/env python3
"""Nares throughflow correction (user decision 2026-10-06): add a constant barotropic -0.46 Sv
(southward = toward i=1) so the 2015-2025 mean through-strait transport is ~0.85 Sv
(Muenchow 2016): uniform normal velocity over the wet East (Lincoln Sea) OB faces and the same
volume over the wet West (Baffin) OB faces, so the net is zero and useOBCSbalance leaves it alone.
Reads forcings/OBCS/UVEL_{east,west}.bin (never written) and writes forcings/OBCS_tf/ with every
other file symlinked. MITgcm OBCS: east normal u at face i = Ieast (0-based NX-1), west at
i = Iwest+1 (0-based 1). File layout (nt, Nr, Ny), big-endian float32.
Usage: nares_throughflow_offset.py <nares dir> [Q_Sv=0.46] [src=OBCS] [dst=OBCS_tf]
2026-10-08: rerun on the rotated files (src OBCS_rot, dst OBCS_rot_tf); Q from obcs_boundary_budget.py."""
import sys, os, glob, numpy as np, xarray as xr
N = sys.argv[1]; Q = float(sys.argv[2]) * 1e6 if len(sys.argv) > 2 else 0.46e6
src = os.path.join(N, 'forcings', sys.argv[3] if len(sys.argv) > 3 else 'OBCS')
dst = os.path.join(N, 'forcings', sys.argv[4] if len(sys.argv) > 4 else 'OBCS_tf')
os.makedirs(dst, exist_ok=True)
g = xr.open_dataset(os.path.join(N, 'Nares_ncgrid.nc'))
hW, dyG, drF = g['HFacW'].values, g['dyG'].values, g['drF'].values      # (Nr,Ny,Nx+1?) check below
NR = drF.size
for b, i in (('east', hW.shape[-1] - 2 if hW.shape[-1] == 1121 else hW.shape[-1] - 1), ('west', 1)):
    area_k_j = hW[:, :, i] * dyG[:, i][None, :] * drF[:, None]           # (Nr, Ny) m^2
    A = area_k_j.sum()
    du = -Q / A                                                            # m/s, negative = southward
    u = np.fromfile(os.path.join(src, f'UVEL_{b}.bin'), '>f4').reshape(-1, NR, hW.shape[1])
    T0 = (u * area_k_j[None]).sum(axis=(1, 2))
    u2 = np.where(area_k_j[None] > 0, u + du, u).astype('>f4')
    T1 = (u2 * area_k_j[None]).sum(axis=(1, 2))
    u2.tofile(os.path.join(dst, f'UVEL_{b}.bin'))
    print(f'{b}: face i={i}, open area {A:.4e} m2, du {du*100:.3f} cm/s, '
          f'mean normal transport {T0.mean()/1e6:+.3f} -> {T1.mean()/1e6:+.3f} Sv over {u.shape[0]} records')
for f in glob.glob(os.path.join(src, '*.bin')):
    t = os.path.join(dst, os.path.basename(f))
    if not os.path.exists(t):
        os.symlink(os.path.realpath(f), t)
print('offset done ->', dst)
