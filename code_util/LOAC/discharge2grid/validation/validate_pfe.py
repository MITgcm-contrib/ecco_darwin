#!/usr/bin/env python3
"""Real-data validation of the shareable package discharge2grid (2026-11-16). Run in a PBS job
(run_validate_d2g.pbs), never on a front end. Writes $OUT/*.csv and prints a report.

V1 legacy reproduction: package (GloFAS runoff map 2015) vs glofas_to_mitgcm.py (md5 16d3b87c),
   WestAO grid, sigma 18, snap 95, exclude outside, no extra outlets -> must be bitwise identical.
V2 outlet rule: coastal cells of raw avg_dis 2015 (40-90N) vs the cells the runoff map 2015 holds.
V3 dead ends: automatic candidates (40-90N) vs the 27 hand-triaged extra rivers.
V4 raw path, WestAO 2015: raw coastal outlets + hand list -> runoff; totals vs the production
   forcing file (runoff map + extras + fills) and vs V1; dense vs sparse size.
V5 dropouts: is the outlet cell coastal in the raw file in a good year and a dropout year?
"""
import csv
import os
import subprocess
import sys
import time

import numpy as np

PKG = "$DATA/LOAC/discharge2grid_pkg"
sys.path.insert(0, PKG)
from discharge2grid import ModelGrid, Mapper, write_mitgcm_runoff, write_sparse_netcdf  # noqa
from discharge2grid.readers import glofas  # noqa
from discharge2grid.geometry import gc_km  # noqa

OUT = os.environ.get("OUT", "$DATA/LOAC/d2g_validation")
RAW = "$DATA/GloFas/raw/glofas_original_{year}.nc"
MAP = "$DATA/LOAC/glofas/river_runoff/GloFas_{year}"
GRID = "$DATA/LOAC/westAO/WestAO_grid.nc"
PROD_SCRIPT = "$DATA/LOAC/pan_AO/glofas_to_mitgcm.py"
PROD_FORCING = "$DATA/LOAC/westAO/runoff/GloFas_WestAO_2015"
HAND = os.path.join(PKG, "examples/glofas_extra_rivers_arctic_alaska.csv")
YEAR = 2015
T0 = time.time()


def say(*a):
    print(f"[{time.time() - T0:7.0f} s]", *a, flush=True)


def write_csv(path, rows, fields):
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields, extrasaction="ignore")
        w.writeheader()
        for r in rows:
            w.writerow(r)


def records(path, shape):
    n = os.path.getsize(path) // (4 * int(np.prod(shape)))
    return np.memmap(path, dtype=">f4", mode="r", shape=(n,) + tuple(shape))


os.makedirs(OUT, exist_ok=True)
grid = ModelGrid.from_netcdf(GRID, mask="mask")
say(grid.summary())
mapper = Mapper(grid, sigma_km=18, max_snap_km=95, exclude_outside=True)

import traceback

STATE = {}


def part_v1():
    say("V1: production script, WestAO 2015, water only, no extras")
    v1p = os.path.join(OUT, "v1_prod")
    os.makedirs(v1p, exist_ok=True)
    with open(os.path.join(OUT, "v1_prod.log"), "w") as fh:
        subprocess.run([sys.executable, "-u", PROD_SCRIPT, "--grid-nc", GRID, "--nc-mask", "mask",
                        "--tag", "WestAO", "--out-dir", v1p, "--years", str(YEAR), "--sigma-km", "18",
                        "--max-snap-km", "95", "--exclude-outside-domain"], stdout=fh, stderr=subprocess.STDOUT,
                       check=True)
    say("V1: package, legacy runoff map")
    src_map = glofas.runoff_map_sources(MAP.format(year=YEAR), mapper, year=YEAR)
    d_map = mapper.map_sources(src_map)
    info1 = d_map.report()
    v1k = os.path.join(OUT, f"v1_pkg_GloFas_WestAO_{YEAR}")
    write_mitgcm_runoff(v1k, d_map)
    A = records(os.path.join(v1p, f"GloFas_WestAO_{YEAR}"), grid.shape)
    B = records(v1k, grid.shape)
    ndiff = sum(int((np.asarray(A[k]) != np.asarray(B[k])).sum()) for k in range(A.shape[0]))
    say(f"V1 RESULT: {A.shape[0]} vs {B.shape[0]} records; values differing: {ndiff} of {A.size:,} "
        f"-> {'BITWISE IDENTICAL' if ndiff == 0 and A.shape == B.shape else 'DIFFERENT'}")

    STATE['B'] = B


def part_v2_v4():
    B = STATE.get('B')
    say("V2: scan raw 2015, 40-90N")
    scan = glofas.scan_raw(RAW.format(year=YEAR), lat_range=(40.0, 90.0))
    say(f"   raw lat {scan.lat[0]:.3f} .. {scan.lat[-1]:.3f} (rows {scan.rows[0]}..{scan.rows[-1]}), "
        f"lon {scan.lon[0]:.3f} .. {scan.lon[-1]:.3f}, {scan.ndays} days")
    out = glofas.coastal_outlets(scan, lat_range=(40.0, 90.0))
    dup = glofas.flag_coastal_duplicates(out, scan)
    # runoff map 2015 mean on the same rows (map is north-to-south = raw row order, checked below)
    mm = np.memmap(MAP.format(year=YEAR), dtype=">f4", mode="r")
    nd = mm.size // (3000 * 7200)
    mm = mm.reshape(nd, 3000, 7200)
    r0, r1 = int(scan.rows[0]), int(scan.rows[-1]) + 1
    # map is north-to-south on the GloFAS grid; index it by LATITUDE/LONGITUDE, not by assumption
    mrows = (3000 - 1 - np.rint((scan.lat - glofas.GLO_LAT0) / glofas.GLO_INC)).astype(int)
    mcols = np.rint((scan.lon - glofas.GLO_LON0) / glofas.GLO_INC).astype(int)
    say(f"   map rows for the raw rows: {mrows[0]}..{mrows[-1]} (raw {r0}..{r1 - 1}); map cols = raw cols: "
        f"{np.array_equal(mcols, np.arange(7200))}")
    m0, m1 = int(mrows.min()), int(mrows.max()) + 1
    acc = np.zeros((r1 - r0, 7200))
    for k in range(nd):
        blk = np.asarray(mm[k, m0:m1, :], dtype=np.float64)
        acc += np.nan_to_num(blk[mrows - m0][:, mcols])
    area = glofas.glofas_cell_area(scan.lat)[:, None]
    qmap = acc / nd * area                                     # m3/s on raw rows/cols
    inband = ((scan.lat >= 40) & (scan.lat <= 90))[:, None] & np.ones((1, 7200), bool)
    mapcells = (qmap > 0) & inband
    raw_mask = np.zeros_like(mapcells)
    raw_mask[out["row"] - r0, out["col"]] = True
    both, only_map, only_raw = mapcells & raw_mask, mapcells & ~raw_mask, raw_mask & ~mapcells
    qr = scan.qmean
    say(f"V2 RESULT: map cells {int(mapcells.sum()):,} ({qmap[mapcells].sum():,.0f} m3/s); raw coastal "
        f"{int(raw_mask.sum()):,} ({qr[raw_mask].sum():,.0f} m3/s)")
    say(f"   both {int(both.sum()):,}: map {qmap[both].sum():,.0f} vs raw {qr[both].sum():,.0f} m3/s; "
        f"ratio map/raw at shared cells: median {np.median(qmap[both] / np.maximum(qr[both], 1e-12)):.6f}, "
        f"p1 {np.percentile(qmap[both] / np.maximum(qr[both], 1e-12), 1):.4f}, "
        f"p99 {np.percentile(qmap[both] / np.maximum(qr[both], 1e-12), 99):.4f}")
    say(f"   only in map {int(only_map.sum()):,} ({qmap[only_map].sum():,.0f} m3/s); only raw-coastal "
        f"{int(only_raw.sum()):,} ({qr[only_raw].sum():,.0f} m3/s)")
    for lab, sel, val in (("only-map", only_map, qmap), ("only-raw", only_raw, qr)):
        jj, ii = np.nonzero(sel)
        o = np.argsort(-val[jj, ii])[:12]
        rows = []
        for t in o:
            j, i = jj[t], ii[t]
            nb = scan.ocean[max(j - 1, 0):j + 2, max(i - 1, 0):i + 2]
            rows.append(dict(lon=f"{scan.lon[i]:.3f}", lat=f"{scan.lat[j]:.3f}", q_map=f"{qmap[j, i]:.1f}",
                             q_raw=f"{qr[j, i]:.1f}", raw_is_sea=int(scan.ocean[j, i]), sea_nbrs=int(nb.sum())))
            say(f"   {lab}: {rows[-1]}")
        write_csv(os.path.join(OUT, f"v2_{lab}.csv"), rows, list(rows[0]) if rows else ["lon"])
    isdup = dup >= 0
    dup_in_map = int(mapcells[out["row"][isdup] - r0, out["col"][isdup]].sum())
    say(f"   along-coast duplicates flagged: {int(isdup.sum()):,} cells, {out['mean_q'][isdup].sum():,.0f} m3/s; "
        f"{dup_in_map} of them are also in the runoff map (= double counted there)")
    for lim in (1, 10, 100, 1000):
        s = out["mean_q"] >= lim
        say(f"   raw coastal cells >= {lim:>4} m3/s: {int(s.sum()):>7,} carrying "
            f"{100 * out['mean_q'][s].sum() / out['mean_q'].sum():5.1f} %")
    write_csv(os.path.join(OUT, "v2_outlets_2015_40N.csv"),
              [dict(lon=f"{out['lon'][k]:.3f}", lat=f"{out['lat'][k]:.3f}", row=int(out['row'][k]),
                    col=int(out['col'][k]), mean_q=f"{out['mean_q'][k]:.4g}", n_sea=int(out['n_sea'][k]),
                    dup_of=int(dup[k]), in_map=int(mapcells[out['row'][k] - r0, out['col'][k]]))
               for k in range(out["lon"].size)],
              ["lon", "lat", "row", "col", "mean_q", "n_sea", "dup_of", "in_map"])

    # ============================================================== V3
    say("V3: dead-end candidates vs the hand list")
    dead = glofas.dead_end_candidates(scan, min_q=100, max_sea_km=250)
    hand = list(csv.DictReader(open(HAND)))
    dl = np.array([d["lon"] for d in dead]); dt = np.array([d["lat"] for d in dead])
    found = 0
    for h in hand:
        dd = gc_km(float(h["lon"]), float(h["lat"]), dl, dt) if dl.size else np.array([np.inf])
        k = int(np.argmin(dd))
        hit = dd[k] <= 15
        found += hit
        say(f"   {h['name']:<15} {'FOUND' if hit else 'missed'}  nearest candidate {dd[k]:7.1f} km "
            f"({dead[k]['mean_q']:,.0f} m3/s)" if dl.size else f"   {h['name']} missed")
    say(f"V3 RESULT: {len(dead)} candidates >= 100 m3/s within 250 km of the sea; {found} of {len(hand)} "
        f"hand-triaged rivers among them")
    write_csv(os.path.join(OUT, "v3_dead_end_candidates_40N.csv"), dead,
              ["lon", "lat", "mean_q", "channel_cells", "sea_km", "sea_lon", "sea_lat"])

    # ============================================================== V4
    say("V4: raw path on WestAO 2015")
    sn = mapper.snap(out["lon"], out["lat"])
    near = np.flatnonzero((sn.status == 0) | (sn.snap_km <= 3 * mapper.max_snap_km))
    hw = [h for h in hand if "WestAO" in h["domains"]]
    res = {}
    for lab, drop_dup in (("coastal+hand", False), ("coastal-dups+hand", True)):
        sel = near[~isdup[near]] if drop_dup else near
        s = glofas.sources_from_raw(RAW, YEAR, out["lon"][sel], out["lat"][sel])
        ex = glofas.sources_from_raw(RAW, YEAR, [float(h["lon"]) for h in hw], [float(h["lat"]) for h in hw],
                                     name=[h["name"] for h in hw], kind=["extra"] * len(hw),
                                     deliver_lon=[float(h["deliver_lon"]) for h in hw],
                                     deliver_lat=[float(h["deliver_lat"]) for h in hw])
        s = s.concat(ex)
        d = mapper.map_sources(s)
        info = d.report(max_list=5)
        tot = np.array([float((d.runoff_points(k) * d.point_area).sum()) for k in range(s.nt)])
        res[lab] = (d, tot)
        say(f"   {lab}: {s.n:,} sources, {d.n_point:,} receiving cells, mean {tot.mean():,.1f} m3/s")
    area_w = grid.area
    P = records(PROD_FORCING, grid.shape)
    tot_prod = np.array([float((np.asarray(P[k], float) * area_w).sum()) for k in range(P.shape[0])])
    tot_v1 = (np.array([float((np.asarray(B[k], float) * area_w).sum()) for k in range(B.shape[0])])
              if B is not None else np.array([np.nan]))
    for lab, (d, tot) in res.items():
        say(f"V4 RESULT {lab}: {tot.mean():,.1f} m3/s vs production forcing {tot_prod.mean():,.1f} "
            f"({100 * (tot.mean() / tot_prod.mean() - 1):+.2f} %) vs V1 runoff map only {tot_v1.mean():,.1f}; "
            f"daily corr with production {np.corrcoef(tot, tot_prod)[0, 1]:.5f}")
    d, tot = res["coastal+hand"]
    p = os.path.join(OUT, f"v4_WestAO_runoff_{YEAR}.nc")
    write_sparse_netcdf(p, d, title="WestAO 2015 raw-path validation")
    dense = 4 * grid.size * d.sources.nt
    say(f"   sparse file {os.path.getsize(p) / 1e6:,.1f} MB ({d.n_point:,} points) vs dense MITgcm file "
        f"{dense / 1e6:,.1f} MB ({grid.size:,} cells)")
    # where do raw path and production differ? annual-mean runoff difference aggregated in 1-deg boxes
    mean_new = np.zeros(grid.size)
    for k in range(d.sources.nt):
        mean_new[d.points] += d.runoff_points(k)
    mean_new = (mean_new / d.sources.nt).reshape(grid.shape)
    mean_prod = np.zeros(grid.shape)
    for k in range(P.shape[0]):
        mean_prod += np.asarray(P[k], float)
    mean_prod /= P.shape[0]
    dq = (mean_new - mean_prod) * area_w                      # m3/s per cell
    keys = np.stack([np.floor(np.where(grid.lon > 180, grid.lon - 360, grid.lon)).ravel(),
                     np.floor(grid.lat).ravel()], axis=1)
    ub, inv = np.unique(keys, axis=0, return_inverse=True)
    s_ = np.bincount(inv.ravel(), weights=dq.ravel())
    o = np.argsort(-np.abs(s_))[:15]
    rows = []
    for t in o:
        lo_, la_ = int(ub[t, 0]), int(ub[t, 1])
        rows.append(dict(box_lon=lo_, box_lat=la_, diff_m3s=f"{s_[t]:.1f}"))
        say(f"   1-deg box lon {lo_:5d} lat {la_:3d}: raw path - production = {s_[t]:+10.1f} m3/s")
    write_csv(os.path.join(OUT, "v4_box_differences.csv"), rows, ["box_lon", "box_lat", "diff_m3s"])



def part_v5():
    say("V5: dropout cells, good year vs dropout year")
    cases = [("Olfusa", -21.225, 63.875, (1990, 1991)), ("Pur", 77.725, 67.575, (1998, 1999)),
             ("Chukotka -173.1", -173.125, 66.975, (1991, 1992)), ("GoA -141.3", -141.275, 59.975, (2007, 2008))]
    for name, lo_, la_, yrs in cases:
        for y in yrs:
            sc = glofas.scan_raw(RAW.format(year=y), lat_range=(la_ - 0.2, la_ + 0.2), verbose=False)
            j = int(np.argmin(np.abs(sc.lat - la_)))
            i = int(np.argmin(np.abs(sc.lon - lo_)))
            nb = sc.ocean[j - 1:j + 2, i - 1:i + 2]
            mm = np.memmap(MAP.format(year=y), dtype=">f4", mode="r")
            ndy = mm.size // (3000 * 7200)
            mm = mm.reshape(ndy, 3000, 7200)
            row = int(3000 - 1 - round((sc.lat[j] - glofas.GLO_LAT0) / glofas.GLO_INC))
            col = int(round((sc.lon[i] - glofas.GLO_LON0) / glofas.GLO_INC))
            mv = np.nanmean(np.asarray(mm[:, row - 1:row + 2, col - 1:col + 2], float), axis=0) * \
                glofas.glofas_cell_area(sc.lat[j])
            say(f"   {name:<16} {y}: raw Q {sc.qmean[j, i]:8.1f} m3/s, cell is sea {bool(sc.ocean[j, i])}, "
                f"sea neighbours {int(nb.sum()) - int(sc.ocean[j, i])}; raw 3x3 Q "
                f"{np.array2string(np.round(sc.qmean[j - 1:j + 2, i - 1:i + 2], 0), separator=',').replace(chr(10), '')}; "
                f"map 3x3 {np.array2string(np.round(mv, 0), separator=',').replace(chr(10), '')}")
    say("done")


for f in (part_v1, part_v2_v4, part_v5):
    try:
        f()
    except Exception:
        say(f"{f.__name__} FAILED:")
        traceback.print_exc()
        sys.stdout.flush()
