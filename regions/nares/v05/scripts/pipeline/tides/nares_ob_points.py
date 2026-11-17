#!/usr/bin/env python3
"""Open-boundary points file for gen_obcs_tides.py (downscaling/python/tides) on the Nares 1-km grid
(2026-10-14, tides_plan.md step 1).

Nares is a rotated lat-lon grid (i runs S -> N along the strait), so unlike WestAO the child-grid
AngleCS/AngleSN vary along each edge (east 43-69 deg, west 71-83 deg; DECISIONS 2026-10-08).
gen_obcs_tides.py rotates the FES (E, N) harmonic constants into grid (x, y) with these angles.

Edges, as in data.obcs / gen_ob_index.py:
    west  = i=1    (Ny=360 entries, one per row j)    North Water / N Baffin
    east  = i=1120 (Ny entries)                       Lincoln Sea
    south = j=1    (Nx=1120 entries, one per column i) Greenland side (Baffin segment)
    north = j=360  (Nx entries)                       Ellesmere side (Baffin + Lincoln segments)
b_has = wet edge cell, wet = (bathy < 0) & (shelficeTopo > bathy), exactly gen_ob_index.py.
lat/lon/angles are taken at the OB cell centre (C point); the u/v faces the tide is added to are
within 0.5 km of it.

Usage: nares_ob_points.py <grid dir> <Nares_ncgrid.nc> <out.npz> [--check-obcs data.obcs]
       [--check-angles <run dir with AngleCS.data/AngleSN.data (float32, model output)>]
"""
import argparse
import re

import numpy as np
import xarray as xr

NX, NY = 1120, 360
EXPECT = {'west': 337, 'east': 342, 'south': 185, 'north': 279}


def read_ob_index(path):
    """OB_Iwest/OB_Ieast (Ny) and OB_Jsouth/OB_Jnorth (Nx) from a data.obcs file"""
    txt = '\n'.join(l.split('#')[0].split('!')[0] for l in open(path).read().splitlines())
    out = {}
    for key, n in (('OB_Iwest', NY), ('OB_Ieast', NY), ('OB_Jsouth', NX), ('OB_Jnorth', NX)):
        m = re.search(rf'^\s*{key}\s*=(.*?)(?=^\s*[A-Za-z&/])', txt, flags=re.S | re.M)
        vals = [int(v) for v in re.findall(r'-?\d+', m.group(1))]
        assert len(vals) == n, (key, len(vals))
        out[key] = np.array(vals)
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('grid_dir')
    ap.add_argument('ncgrid')
    ap.add_argument('out')
    ap.add_argument('--check-obcs')
    ap.add_argument('--check-angles')
    a = ap.parse_args()

    b = np.fromfile(f'{a.grid_dir}/Nares_bathymetry.bin', '>f8').reshape(NY, NX)
    s = np.fromfile(f'{a.grid_dir}/Nares_shelficeTopo.bin', '>f8').reshape(NY, NX)
    wet = (b < 0) & (s > b)
    g = xr.open_dataset(a.ncgrid)
    xc, yc = g.XC.values, g.YC.values
    cs, sn = g.AngleCS.values, g.AngleSN.values
    assert xc.shape == (NY, NX)
    r = np.hypot(cs, sn)
    assert np.abs(r[wet] - 1).max() < 1e-6, np.abs(r[wet] - 1).max()

    if a.check_angles:   # model-written angles (float32) vs the ncgrid ones
        # Found 2026-10-14: the model's own AngleCS/SN (run_arch6m_D) agree with the ncgrid to < 0.09 deg
        # in the interior (small differences on tile-edge rows/columns), but are WRONG on the whole
        # east (i=1120) and north (j=360) edge rows: up to 167 deg, e.g. -90 deg instead of 54 deg at
        # (j=180, i=1120) -- the model has no neighbor beyond the facet edge there. The ncgrid angles
        # equal the geometric x-axis direction from XC/YC (checked), so the tide rotation uses them.
        mcs = np.fromfile(f'{a.check_angles}/AngleCS.data', '>f4').reshape(NY, NX)
        msn = np.fromfile(f'{a.check_angles}/AngleSN.data', '>f4').reshape(NY, NX)
        dang = np.degrees(np.abs(np.angle((mcs + 1j * msn) / (cs + 1j * sn))))
        inner = wet.copy()
        inner[0, :] = inner[-1, :] = inner[:, 0] = inner[:, -1] = False
        print(f'model vs ncgrid grid angle, wet interior: max {dang[inner].max():.3f} deg')
        assert dang[inner].max() < 0.2
        for bnd, ix in (('west', (slice(None), 0)), ('east', (slice(None), -1)),
                        ('south', (0, slice(None))), ('north', (-1, slice(None)))):
            w = wet[ix]
            print(f'  {bnd:5s} edge (wet): model angle off by up to {dang[ix][w].max():.2f} deg '
                  f'({(dang[ix][w] > 1).sum()} of {w.sum()} cells > 1 deg)')

    sl = {'west': (slice(None), 0), 'east': (slice(None), -1),
          'south': (0, slice(None)), 'north': (-1, slice(None))}
    out = {}
    for bnd, ix in sl.items():
        out[f'{bnd}_lat'] = yc[ix].astype('f8')
        out[f'{bnd}_lon'] = (xc[ix] % 360).astype('f8')
        out[f'{bnd}_cs'] = cs[ix].astype('f8')
        out[f'{bnd}_sn'] = sn[ix].astype('f8')
        out[f'{bnd}_has'] = wet[ix].copy()
        h = out[f'{bnd}_has']
        ang = np.degrees(np.arctan2(sn[ix][h], cs[ix][h]))
        print(f'{bnd:5s}: {h.sum():4d} of {h.size} open | lat {yc[ix][h].min():.2f}..{yc[ix][h].max():.2f} '
              f'| grid-x angle {ang.min():.1f}..{ang.max():.1f} deg')
        assert h.sum() == EXPECT[bnd], (bnd, h.sum())

    if a.check_obcs:     # the OB index arrays the model actually uses
        ob = read_ob_index(a.check_obcs)
        for bnd, key, val in (('west', 'OB_Iwest', 1), ('east', 'OB_Ieast', NX),
                              ('south', 'OB_Jsouth', 1), ('north', 'OB_Jnorth', NY)):
            used = ob[key] == val
            assert np.all(ob[key][~used] == 0), key
            assert np.array_equal(used, out[f'{bnd}_has']), (bnd, int(used.sum()))
        print('b_has == data.obcs OB index arrays on all four edges')

    np.savez(a.out, **out)
    print('->', a.out)


if __name__ == '__main__':
    main()
