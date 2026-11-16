#!/usr/bin/env python3
"""discharge2grid river temperature, real data (2026-11-16). Run in a PBS job (run_validate_v7.pbs).
V7a legacy: package (runoff map + temperature map 2015) vs glofas_to_mitgcm.py output (v1_prod, no extras)
V7b climatology rebuilt from the raw wattmp files (same file set as make_glofas_river_temp_forcing.py),
    sampled at the runoff-map source cells, vs the production temperature map at those cells
V7c raw path (coast rule, duplicates removed + WestAO hand list) with climatology temperature:
    heat flux and flux-weighted temperature vs the production forcing (map + extras + donors)
"""
import csv
import glob
import os
import sys
import time

import numpy as np

PKG = "$DATA/LOAC/discharge2grid_pkg"
sys.path.insert(0, PKG)
from discharge2grid import ModelGrid, Mapper, write_mitgcm_runoff, write_mitgcm_temperature  # noqa
from discharge2grid.readers import glofas, temperature as rt  # noqa

OUT = os.environ.get("OUT", "$DATA/LOAC/d2g_validation")
RAW = "$DATA/GloFas/raw/glofas_original_{year}.nc"
MAP = "$DATA/LOAC/glofas/river_runoff/GloFas_{year}"
TMAP = "$DATA/LOAC/glofas/river_temp/GloFas_temp_{year}"
WATTMP = "$DATA/LOAC/river_temp/wattmp*.bin"
GRID = "$DATA/LOAC/westAO/WestAO_grid.nc"
PROD_R = "$DATA/LOAC/westAO/runoff/GloFas_WestAO_2015"
PROD_T = "$DATA/LOAC/westAO/runoff/GloFas_temp_WestAO_2015"
HAND = os.path.join(PKG, "examples/glofas_extra_rivers_arctic_alaska.csv")
Y = 2015
T0 = time.time()


def say(*a):
    print(f"[{time.time() - T0:6.0f} s]", *a, flush=True)


grid = ModelGrid.from_netcdf(GRID, mask="mask", verbose=False)
mapper = Mapper(grid, sigma_km=18, max_snap_km=95, exclude_outside=True)
shape = grid.shape

# ---------------------------------------------------------------- V7a
src = glofas.runoff_map_sources(MAP.format(year=Y), mapper, year=Y, temp_path=TMAP.format(year=Y))
d = mapper.map_sources(src)
d.report(max_list=3)
pk = os.path.join(OUT, f"v7_pkg_GloFas_temp_WestAO_{Y}")
write_mitgcm_temperature(pk, d)
A = np.memmap(os.path.join(OUT, "v1_prod", f"GloFas_temp_WestAO_{Y}"), dtype=">f4", mode="r").reshape(-1, *shape)
B = np.memmap(pk, dtype=">f4", mode="r").reshape(-1, *shape)
nd = sum(int((np.asarray(A[k]) != np.asarray(B[k])).sum()) for k in range(A.shape[0]))
say(f"V7a RESULT: temperature {A.shape[0]} vs {B.shape[0]} records, values differing {nd} of {A.size:,} -> "
    f"{'BITWISE IDENTICAL' if nd == 0 and A.shape == B.shape else 'DIFFERENT'}")

# ---------------------------------------------------------------- V7b
files = sorted(glob.glob(WATTMP))
say(f"V7b: climatology from {len(files)} files ({os.path.basename(files[0])} .. {os.path.basename(files[-1])})")
clim = rt.TemperatureClimatology.from_binary(files, lat_range=(35.0, 89.9), **rt.WATTMP)
clim.save(os.path.join(OUT, "v7_wattmp_clim_35N.npz"))
ts, dkm = clim.sample(src.lon, src.lat, src.time)
wet = src.q > 0                                       # production assigned temperature where runoff
diff = np.abs(ts.astype(np.float32) - src.temp.astype(np.float32))
eq = (diff == 0) & wet
say(f"V7b RESULT: {int(wet.sum()):,} source-days with runoff; package climatology == production map "
    f"on {100 * eq.sum() / wet.sum():.2f} %; |diff| median {np.median(diff[wet]):.4f}, p99 "
    f"{np.percentile(diff[wet], 99):.4f}, max {diff[wet].max():.4f} C; lookup distance median "
    f"{np.median(dkm):.1f} km, max {dkm.max():.1f} km")
bad = wet & (diff > 0.01)
if bad.any():
    jj = np.unique(np.nonzero(bad)[1])
    say(f"   {jj.size} sources differ by > 0.01 C on some day (production averages the values of "
        f"neighbouring runoff cells over each cell's 3x3); largest:")
    worst = jj[np.argsort(-np.nanmax(np.where(bad, diff, 0), axis=0)[jj])][:6]
    for j in worst:
        k = int(np.argmax(np.where(bad[:, j], diff[:, j], 0)))
        say(f"   {src.lon[j]:9.3f} {src.lat[j]:7.3f} day {k}: package {ts[k, j]:.3f} production "
            f"{src.temp[k, j]:.3f} (Q {src.q[k, j]:.1f} m3/s)")

# ---------------------------------------------------------------- V7c
scan = glofas.scan_raw(RAW.format(year=Y), lat_range=(40.0, 90.0), verbose=False)
out = glofas.coastal_outlets(scan, lat_range=(40.0, 90.0))
keep = glofas.flag_coastal_duplicates(out, scan) < 0
sn = mapper.snap(out["lon"], out["lat"])
near = np.flatnonzero(keep & ((sn.status == 0) | (sn.snap_km <= 3 * mapper.max_snap_km)))
hand = [h for h in csv.DictReader(open(HAND)) if "WestAO" in h["domains"]]
s = glofas.sources_from_raw(RAW, Y, out["lon"][near], out["lat"][near])
ex = glofas.sources_from_raw(RAW, Y, [float(h["lon"]) for h in hand], [float(h["lat"]) for h in hand],
                             name=[h["name"] for h in hand], kind=["extra"] * len(hand),
                             deliver_lon=[float(h["deliver_lon"]) for h in hand],
                             deliver_lat=[float(h["deliver_lat"]) for h in hand])
s = rt.add_temperature(s.concat(ex), clim)
s, donors = rt.fill_missing(s, max_km=150)
say(f"V7c: raw path {s.n:,} sources; {s.attrs.get('temperature')}; {s.attrs.get('temperature_fill')}")
dr = mapper.map_sources(s)
dr.report(max_list=3)
PR = np.memmap(PROD_R, dtype=">f4", mode="r").reshape(-1, *shape)
PT = np.memmap(PROD_T, dtype=">f4", mode="r").reshape(-1, *shape)
hp, hr, qp, qr = [], [], [], []
for k in range(s.nt):
    r = np.asarray(PR[k], float)
    t = np.asarray(PT[k], float)
    qp.append(float((r * grid.area).sum()))
    hp.append(float((r * t * grid.area).sum()))
    t2, r2 = dr.temperature_points(k, return_runoff=True)
    g = np.isfinite(t2)
    qr.append(float((r2 * dr.point_area).sum()))
    hr.append(float((r2[g] * t2[g] * dr.point_area[g]).sum()))
hp, hr, qp, qr = map(np.array, (hp, hr, qp, qr))
tp, tr = hp / qp, hr / qr                                # domain flux-weighted temperature per day
say(f"V7c RESULT: flux-weighted river temperature, annual mean: production {np.sum(hp) / np.sum(qp):.3f} C, "
    f"raw path {np.sum(hr) / np.sum(qr):.3f} C; daily mean |diff| {np.mean(np.abs(tp - tr)):.3f} C, "
    f"max {np.max(np.abs(tp - tr)):.3f} C, corr {np.corrcoef(tp, tr)[0, 1]:.4f}")
say(f"   heat flux (sum Q T) raw/production: {np.sum(hr) / np.sum(hp):.4f}; water raw/production "
    f"{np.sum(qr) / np.sum(qp):.4f}")
for m0, m1, lab in ((0, 90, "Jan-Mar"), (90, 181, "Apr-Jun"), (181, 273, "Jul-Sep"), (273, 365, "Oct-Dec")):
    say(f"   {lab}: production {np.sum(hp[m0:m1]) / np.sum(qp[m0:m1]):6.3f} C, raw path "
        f"{np.sum(hr[m0:m1]) / np.sum(qr[m0:m1]):6.3f} C")
say("done")
