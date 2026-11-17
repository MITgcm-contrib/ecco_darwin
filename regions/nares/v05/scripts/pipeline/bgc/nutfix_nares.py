"""NARES COPY (2026-10-29) of west_AO nutfix_westao.py (python/west_AO/nutfix): pan_AO nutrient fix, factor v3.
Dustin 2026-10-29: "Check, then fix". Check (BGC_PLAN.md section 10): parent 2020 annual vs WOA18 annual + GLODAP in the
Nares boxes -> SiO2 ~2x low, PO4 ~1.5x low, NO3 25-40% low vs pre-bloom GLODAP; so apply the same rule as pan_AO/west_AO.
Rules (same as pan_AO nutfix_apply_panao_v3.py / west_AO):
  IC : NO3/PO4/SiO2 pickups (pTr02/05/07) x min(FOB, 3) = the WOA18-annual part of FIC. Nares lies in pan_AO's 'Baffin'
       sector, which has NO GLODAP shelf target (targets only Chukchi/Beaufort/Barents/Greenland), so pan_AO's own FIC
       there is the annual part only: the same value as here.
  OB : TRAC02/05/07_<edge>.bin x min(FOB, WIN / Jan-Apr-2015 max), >= 1 (time-invariant per cell and level).
Only values > 0 are scaled. Factor: bilinear in lat/lon on 1-degree cells, linear in depth.
Nares grid: curvilinear (rotated), so XC/YC are read from grid/Nares_ncgrid.nc; RC = cumsum(drF) - drF/2 (no RC in the file).
Edges (MITgcm naming, i = 1..1120 runs S->N along the strait): north = row j=360 (1120 pts), south = row j=1,
west = column i=1 (360 pts, North Water), east = column i=1120 (Lincoln Sea).
Writes NEW directories only (dst must not exist / be empty); symlinks for every unchanged entry.
Usage: nutfix_nares.py FACTOR_NPZ NCGRID ic  <src pickup dir> <dst dir>
       nutfix_nares.py FACTOR_NPZ NCGRID ob  <src OB dir> <dst dir>
       nutfix_nares.py FACTOR_NPZ NCGRID limit <src OB dir>          (print effective OB factors only)"""
import numpy as np, os, sys, shutil
import xarray as xr
from scipy.interpolate import RegularGridInterpolator
NX, NY, NR = 1120, 360, 90; TR = {'NO3': 2, 'PO4': 5, 'SiO2': 7}; IT = '0000604872'
EDGE = {'north': NX, 'south': NX, 'west': NY, 'east': NY}; NJA = 120 * 24     # Jan 1 - Apr 30 2015, hourly
FAC, NCG = sys.argv[1], sys.argv[2]; F = np.load(FAC); flat, flon, frc = F['lat'], F['lon'], F['rc']
G = xr.open_dataset(NCG)
XC = np.asarray(G['XC'].values, 'f8'); YC = np.asarray(G['YC'].values, 'f8'); drF = np.asarray(G['drF'].values, 'f8').ravel(); rc = np.cumsum(drF) - 0.5 * drF   # ncgrid has no RC
assert XC.shape == (NY, NX) and YC.shape == (NY, NX) and rc.size == NR and rc[0] > 0, (XC.shape, rc.size)
assert 70 < np.nanmedian(YC) < 85, 'Nares lat check'
def lonlat(edge):
    if edge == 'north': return XC[-1], YC[-1]
    if edge == 'south': return XC[0], YC[0]
    if edge == 'west': return XC[:, 0], YC[:, 0]
    return XC[:, -1], YC[:, -1]
def factor_at(A, lat, lon, fill=1.0):
    """A (50, nlat, nlon) on the 1-degree factor grid -> (NR, npts) (copied from nutfix_westao.py)"""
    A = A.astype('f8'); A = np.concatenate([A[:, :, -1:], A, A[:, :, :1]], 2); lo = np.r_[flon[0] - 1, flon, flon[-1] + 1]
    la = np.r_[flat[0] - 1, flat]; A = np.concatenate([np.full_like(A[:, :1], fill), A], 1)
    lat = np.clip(lat, la[0], la[-1]); lon = ((((lon + 180) % 360) - 180) - lo[0]) % 360 + lo[0]
    H = np.empty((frc.size, lat.size))
    for k in range(frc.size):
        H[k] = RegularGridInterpolator((la, lo), A[k], bounds_error=False, fill_value=fill)(np.c_[lat, lon])
    out = np.empty((NR, lat.size))
    for k in range(NR):
        z = min(max(rc[k], frc[0]), frc[-1]); j = min(max(np.searchsorted(frc, z), 1), frc.size - 1); w = (z - frc[j-1]) / (frc[j] - frc[j-1])
        out[k] = (1 - w) * H[j-1] + w * H[j]
    return out
def ob_factor(tr, edge, src):
    lon, lat = lonlat(edge); n = EDGE[edge]
    fob = factor_at(F['FOB_' + tr], lat, lon)
    U = factor_at(np.nan_to_num(F['WIN_' + tr], nan=1e9), lat, lon, fill=1e9); U[U > 1e8] = np.nan
    f = os.path.realpath('%s/TRAC%02d_%s.bin' % (src, TR[tr], edge))
    wmax = np.fromfile(f, '>f4', count=NJA * NR * n).reshape(NJA, NR, n).max(0)
    with np.errstate(invalid='ignore', divide='ignore'):
        lim = np.where(np.isfinite(U) & (wmax > 0), U / wmax, np.inf)
    return np.maximum(np.minimum(fob, lim), 1.0).astype('f4'), fob, wmax
def link_rest(src, dst, skip):
    for e in sorted(os.listdir(src)):
        if e in skip or os.path.lexists(dst + '/' + e): continue
        os.symlink(os.path.realpath(src + '/' + e), dst + '/' + e)
def newdir(d):
    if os.path.exists(d) and os.listdir(d): sys.exit('FATAL: %s exists and is not empty' % d)
    os.makedirs(d, exist_ok=True)
if __name__ == '__main__':
    mode = sys.argv[3]
    if mode == 'limit':
        for tr in TR:
            for edge in EDGE:
                fe, fob, wmax = ob_factor(tr, edge, sys.argv[4]); s = wmax[0] > 0
                print('%-4s %-5s FOB 0-100 m median %.2f | effective median %.2f (k1 %.2f, k17 %.2f, k30 %.2f) | Jan-Apr max k1 %.2f -> %.2f' % (
                      tr, edge, np.median(fob[:22][:, s]), np.median(fe[:22][:, s]), np.median(fe[0, s]), np.median(fe[16, s]), np.median(fe[29, s]),
                      np.median(wmax[0, s]), np.median((wmax[0] * fe[0])[s])), flush=True)
    elif mode == 'ic':
        src, dst = sys.argv[4], sys.argv[5]; newdir(dst)
        lat, lon = YC.ravel(), XC.ravel()
        for tr, t in TR.items():
            nm = 'pickup_pTr%02d.%s' % (t, IT); s = os.path.realpath(src + '/' + nm + '.data'); d = dst + '/' + nm + '.data'
            assert os.path.getsize(s) == NR * NX * NY * 8, s
            a = np.fromfile(s, '>f8').reshape(NR, NX * NY)
            fac = np.minimum(factor_at(F['FOB_' + tr], lat, lon), 3.0)
            b = np.where(a > 0, a * fac, a); assert np.isfinite(b).all(); b.astype('>f8').tofile(d + '.part'); os.replace(d + '.part', d)
            shutil.copyfile(os.path.realpath(src + '/' + nm + '.meta'), dst + '/' + nm + '.meta')
            for z0, z1 in ((0, 50), (50, 200), (200, 1000)):
                kk = (rc >= z0) & (rc < z1); w = a[kk] > 0
                print('IC %-4s %4d-%4d m: factor median %.2f p90 %.2f | mean %.2f -> %.2f' % (tr, z0, z1, np.median(fac[kk][w]), np.percentile(fac[kk][w], 90), a[kk][w].mean(), b[kk][w].mean()), flush=True)
        link_rest(src, dst, {'pickup_pTr%02d.%s.%s' % (t, IT, x) for t in TR.values() for x in ('data', 'meta')})
        print('IC done ->', dst, len(os.listdir(dst)), 'entries (source', len(os.listdir(src)), ')', flush=True)
    elif mode == 'ob':
        src, dst = sys.argv[4], sys.argv[5]; newdir(dst)
        FE = {(tr, e): ob_factor(tr, e, src)[0] for tr in TR for e in EDGE}
        names = set()
        for tr, t in TR.items():
            for e, n in EDGE.items():
                nm = 'TRAC%02d_%s.bin' % (t, e); names.add(nm); s = os.path.realpath(src + '/' + nm); d = dst + '/' + nm
                sz = os.path.getsize(s); nrec = sz // (NR * n * 4); assert nrec * NR * n * 4 == sz
                fe = FE[(tr, e)]; s_in = s_out = 0.0
                with open(s, 'rb') as fi, open(d + '.part', 'wb') as fo:
                    for r in range(0, nrec, 240):
                        a = np.fromfile(fi, '>f4', count=min(240, nrec - r) * NR * n).reshape(-1, NR, n)
                        b = np.where(a > 0, a * fe[None], a).astype('>f4'); b.tofile(fo)
                        s_in += float(a.clip(0).sum()); s_out += float(b.clip(0).sum())
                os.replace(d + '.part', d)
                print('%s: %d records, sum(out)/sum(in) %.3f' % (nm, nrec, s_out / max(s_in, 1e-30)), flush=True)
        link_rest(src, dst, names)
        print('OB done ->', dst, len(os.listdir(dst)), 'entries (source', len(os.listdir(src)), ')', flush=True)
