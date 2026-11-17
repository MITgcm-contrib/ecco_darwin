#!/usr/bin/env python3
"""SHELFICEloadAnomalyFile for the static Petermann cavity, ported from MITgcm
verification/isomip/input/gendata.m (jmd95z branch): hydrostatic geopotential anomaly of the
reference profile (tRef/sRef from nares_tref_sref.py) relative to rhoConst, evaluated at the ice
base with the partial-cell correction drloc = 1 - hFacC(k). Units Pa (MITgcm multiplies by
recip_gravity for the ice mass). rhoConst = rhonil = 1027.5, gravity = 9.81 (input/data).
Usage: nares_shelfice_load.py <nares dir> <tref_sref.npz> <out .bin>"""
import sys, os, numpy as np, xarray as xr
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import jmd95   # copied from darwin3 utils/python/MITgcmutils (the downscaling env's MITgcmutils lacks it)
N, ts, out = sys.argv[1:4]
RHO0, G = 1027.5, 9.81
g = xr.open_dataset(os.path.join(N, 'Nares_ncgrid.nc'))
hC, drF = g['HFacC'].values, g['drF'].values
rF = -np.r_[0, np.cumsum(drF)]                  # face depths (m, <= 0)
rC = 0.5 * (rF[:-1] + rF[1:])
z = np.load(ts); t, s = z['tRef'], z['sRef']
p = np.abs(rC) * G * RHO0 * 1e-4               # dbar
drho = jmd95.dens(s, t, p) - RHO0
dzm = np.abs(np.r_[rF[0] - rC[0], 0.5 * np.diff(rC)])
dzp = np.abs(np.r_[0.5 * np.diff(rC), rC[-1] - rF[-1]])
phiF = np.zeros(drF.size + 1)
for k in range(drF.size):
    phiC = phiF[k] + dzm[k] * G * drho[k] / RHO0
    phiF[k + 1] = phiC + dzp[k] * G * drho[k] / RHO0
topo = np.fromfile(os.path.join(N, 'grid/Nares_shelficeTopo.bin'), '>f8').reshape(hC.shape[1:])
msk = (hC.sum(0) > 0)
load = np.zeros_like(topo)
for j, i in zip(*np.where((topo < 0) & msk)):
    ks = np.where(np.abs(rF[:-1]) < abs(topo[j, i]))[0]
    if ks.size == 0:
        continue
    k = ks.max(); kp1 = min(k + 1, drF.size)
    drloc = 1 - hC[k, j, i]
    load[j, i] = (phiF[k] + drloc * (phiF[kp1] - phiF[k])) * RHO0
load.astype('>f8').tofile(out)
c = topo < 0
print(f'shelf cells {c.sum()}, load anomaly min {load[c].min():.1f} max {load[c].max():.1f} Pa -> {out}')
