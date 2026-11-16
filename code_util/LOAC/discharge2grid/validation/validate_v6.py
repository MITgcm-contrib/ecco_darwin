#!/usr/bin/env python3
"""discharge2grid validation part 2 (2026-11-16), after V1-V5 (validate_pfe.py):
V6a cell areas: GloFAS pixarea vs the spherical formula of glofas_to_mitgcm.py (explains map/raw 0.993?)
V6b where the coast rule (duplicates removed) differs from the runoff map, 2015, 40-90N, by 10x5 deg box
V6c WestAO 2015, raw path with duplicates removed + hand list vs production: 1-deg boxes
V6d legacy runoff map converted with pixarea instead of the formula: WestAO total vs raw
Run in a PBS job (run_validate_v6.pbs)."""
import csv
import os
import sys
import time

import numpy as np

PKG = "$DATA/LOAC/discharge2grid_pkg"
sys.path.insert(0, PKG)
from discharge2grid import ModelGrid, Mapper  # noqa
from discharge2grid.readers import glofas  # noqa

OUT = os.environ.get("OUT", "$DATA/LOAC/d2g_validation")
RAW = "$DATA/GloFas/raw/glofas_original_{year}.nc"
MAP = "$DATA/LOAC/glofas/river_runoff/GloFas_{year}"
PIX = "$DATA/GloFas/GLOFAS_pixarea_Global_03min.nc"
GRID = "$DATA/LOAC/westAO/WestAO_grid.nc"
PROD = "$DATA/LOAC/westAO/runoff/GloFas_WestAO_2015"
HAND = os.path.join(PKG, "examples/glofas_extra_rivers_arctic_alaska.csv")
Y = 2015
T0 = time.time()


def say(*a):
    print(f"[{time.time() - T0:6.0f} s]", *a, flush=True)


# ---------------------------------------------------------------- V6a
lats = np.array([42.5, 50.025, 55.025, 60.025, 65.025, 70.025, 75.025, 80.025, 85.025])
pix = glofas.cell_area_from_file(PIX, lats, np.zeros(lats.size))
form = glofas.glofas_cell_area(lats)
for la, p, f in zip(lats, pix, form):
    say(f"V6a lat {la:6.2f}: pixarea {p / 1e6:8.3f} km2, formula {f / 1e6:8.3f} km2, formula/pixarea "
        f"{f / p:.5f}")

# ---------------------------------------------------------------- V6b
scan = glofas.scan_raw(RAW.format(year=Y), lat_range=(40.0, 90.0), verbose=False)
out = glofas.coastal_outlets(scan, lat_range=(40.0, 90.0))
dup = glofas.flag_coastal_duplicates(out, scan)
keep = dup < 0
mm = np.memmap(MAP.format(year=Y), dtype=">f4", mode="r")
nd = mm.size // (3000 * 7200)
mm = mm.reshape(nd, 3000, 7200)
mrows = (3000 - 1 - np.rint((scan.lat - glofas.GLO_LAT0) / glofas.GLO_INC)).astype(int)
m0, m1 = int(mrows.min()), int(mrows.max()) + 1
acc = np.zeros((scan.lat.size, 7200))
for k in range(nd):
    acc += np.nan_to_num(np.asarray(mm[k, m0:m1, :], dtype=np.float64))[mrows - m0]
pixrow = glofas.cell_area_from_file(PIX, scan.lat, np.zeros(scan.lat.size))
qmap = acc / nd * pixrow[:, None]                       # raw-equivalent m3/s (pixarea)
inb = ((scan.lat >= 40) & (scan.lat <= 90))[:, None]
qmap = np.where(inb, qmap, 0.0)
tot_map = qmap.sum()
tot_raw = out["mean_q"][keep].sum()
say(f"V6b 2015 40-90N: map x pixarea {tot_map:,.0f} m3/s; coast rule, duplicates removed {tot_raw:,.0f} "
    f"({100 * (tot_raw / tot_map - 1):+.1f} %); with duplicates {out['mean_q'].sum():,.0f}")
LON2 = np.broadcast_to(scan.lon[None, :], qmap.shape)
LAT2 = np.broadcast_to(scan.lat[:, None], qmap.shape)
bm = {}
for lo_, la_, q in zip(LON2[qmap > 0], LAT2[qmap > 0], qmap[qmap > 0]):
    key = (int(np.floor(lo_ / 10) * 10), int(np.floor(la_ / 5) * 5))
    bm[key] = bm.get(key, 0.0) + q
br = {}
for lo_, la_, q in zip(out["lon"][keep], out["lat"][keep], out["mean_q"][keep]):
    key = (int(np.floor(lo_ / 10) * 10), int(np.floor(la_ / 5) * 5))
    br[key] = br.get(key, 0.0) + q
keys = sorted(set(bm) | set(br), key=lambda k: -abs(br.get(k, 0) - bm.get(k, 0)))
rows = []
for k in keys[:20]:
    r, m = br.get(k, 0.0), bm.get(k, 0.0)
    rows.append(dict(lon0=k[0], lat0=k[1], coast_rule=f"{r:.0f}", runoff_map=f"{m:.0f}", diff=f"{r - m:.0f}"))
    say(f"   box lon {k[0]:5d}..{k[0] + 10:<5d} lat {k[1]:3d}..{k[1] + 5:<3d}: coast rule {r:10,.0f}  map "
        f"{m:10,.0f}  diff {r - m:+9,.0f} m3/s")
with open(os.path.join(OUT, "v6b_boxes.csv"), "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=list(rows[0]))
    w.writeheader()
    w.writerows(rows)
# small cells: how much of the difference is cells < 10 m3/s?
small = keep & (out["mean_q"] < 10)
qms = qmap[(qmap > 0) & (qmap < 10)].sum()
say(f"   cells < 10 m3/s: coast rule {out['mean_q'][small].sum():,.0f} ({int(small.sum()):,} cells) vs map "
    f"{qms:,.0f} ({int(((qmap > 0) & (qmap < 10)).sum()):,} cells)")
big = keep & (out["mean_q"] >= 10)
say(f"   cells >= 10 m3/s: coast rule {out['mean_q'][big].sum():,.0f} ({int(big.sum()):,}) vs map "
    f"{qmap[qmap >= 10].sum():,.0f} ({int((qmap >= 10).sum()):,})")

# ---------------------------------------------------------------- V6c
grid = ModelGrid.from_netcdf(GRID, mask="mask", verbose=False)
mapper = Mapper(grid, sigma_km=18, max_snap_km=95, exclude_outside=True)
sn = mapper.snap(out["lon"], out["lat"])
near = np.flatnonzero(keep & ((sn.status == 0) | (sn.snap_km <= 3 * mapper.max_snap_km)))
hand = [h for h in csv.DictReader(open(HAND)) if "WestAO" in h["domains"]]
s = glofas.sources_from_raw(RAW, Y, out["lon"][near], out["lat"][near])
ex = glofas.sources_from_raw(RAW, Y, [float(h["lon"]) for h in hand], [float(h["lat"]) for h in hand],
                             name=[h["name"] for h in hand], kind=["extra"] * len(hand),
                             deliver_lon=[float(h["deliver_lon"]) for h in hand],
                             deliver_lat=[float(h["deliver_lat"]) for h in hand])
d = mapper.map_sources(s.concat(ex))
mean_new = np.zeros(grid.size)
for k in range(d.sources.nt):
    mean_new[d.points] += d.runoff_points(k)
mean_new = (mean_new / d.sources.nt).reshape(grid.shape)
P = np.memmap(PROD, dtype=">f4", mode="r").reshape(-1, *grid.shape)
mean_prod = np.zeros(grid.shape)
for k in range(P.shape[0]):
    mean_prod += np.asarray(P[k], float)
mean_prod /= P.shape[0]
tn, tp = float((mean_new * grid.area).sum()), float((mean_prod * grid.area).sum())
say(f"V6c WestAO 2015: raw path (dups removed + hand list) {tn:,.0f} m3/s vs production {tp:,.0f} "
    f"({100 * (tn / tp - 1):+.1f} %)")
dq = (mean_new - mean_prod) * grid.area
lon180 = np.where(grid.lon > 180, grid.lon - 360, grid.lon)
keys = np.stack([np.floor(lon180).ravel(), np.floor(grid.lat).ravel()], axis=1)
ub, inv = np.unique(keys, axis=0, return_inverse=True)
sd = np.bincount(inv.ravel(), weights=dq.ravel())
sp = np.bincount(inv.ravel(), weights=(mean_prod * grid.area).ravel())
for t in np.argsort(-np.abs(sd))[:15]:
    say(f"   1-deg box lon {int(ub[t, 0]):5d} lat {int(ub[t, 1]):3d}: production {sp[t]:9,.1f}  "
        f"raw path - production {sd[t]:+9,.1f} m3/s")
pos, neg = sd[sd > 0].sum(), sd[sd < 0].sum()
say(f"   sum of positive box differences {pos:,.0f}, negative {neg:,.0f} m3/s")

# ---------------------------------------------------------------- V6d
sm = glofas.runoff_map_sources(MAP.format(year=Y), mapper, year=Y, area_file=PIX, verbose=False)
sf = glofas.runoff_map_sources(MAP.format(year=Y), mapper, year=Y, verbose=False)
km = mapper.snap(sm.lon, sm.lat).kept
say(f"V6d WestAO 2015 runoff map, kept sources: x formula area {sf.q[:, km].sum(axis=1).mean():,.1f} m3/s, "
    f"x pixarea {sm.q[:, km].sum(axis=1).mean():,.1f} m3/s "
    f"(ratio {sf.q[:, km].sum() / sm.q[:, km].sum():.5f})")
say("done")
