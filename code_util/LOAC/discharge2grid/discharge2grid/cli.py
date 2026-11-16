"""Command line:  python -m discharge2grid {outlets,temp-clim,run,dropouts} ...

  outlets   find GloFAS outlets (drainage map or coast rule) and dead-end candidates
  temp-clim build a daily river-temperature climatology from gridded files (once)
  run       sources -> grid -> MITgcm binaries and/or sparse netCDF, year by year
            (runoff, plus river temperature when a temperature source is given)
  dropouts  scan per-year source tables for rivers that vanish in some years

Run any subcommand with -h for its options; docs/INSTRUCTIONS.md walks through a full case.
"""
import argparse
import csv
import glob
import os
import re
import sys
import time

import numpy as np


# ---------------------------------------------------------------------- helpers
def _years(spec):
    out = []
    for part in str(spec).replace(" ", "").split(","):
        if not part:
            continue
        if "-" in part:
            a, b = part.split("-")
            out += list(range(int(a), int(b) + 1))
        else:
            out.append(int(part))
    return out


def _read_csv(path):
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh))


def _write_csv(path, rows, fields):
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields, extrasaction="ignore")
        w.writeheader()
        for r in rows:
            w.writerow(r)


def _grid(a):
    from .grid import ModelGrid
    if a.grid_mds:
        return ModelGrid.from_mds(a.grid_mds)
    return ModelGrid.from_netcdf(a.grid, lon=a.grid_lon, lat=a.grid_lat, area=a.grid_area,
                                 mask=a.grid_mask, level=a.mask_level, depth=a.grid_depth)


# ---------------------------------------------------------------------- outlets
def cmd_outlets(a):
    from .readers import glofas
    years = _years(a.years) if a.years else []
    t0 = time.time()
    if a.ldd:
        out = glofas.outlets_from_ldd(a.ldd, lat_range=a.lat_range, var=a.ldd_var,
                                      max_sea_km=a.dead_max_sea_km)
        mq = np.full(out["lon"].size, np.nan)
        if a.raw and years:                            # mean discharge, for information
            q, _ = glofas.read_outlet_series(a.raw.format(year=years[0]), out["lon"], out["lat"],
                                             var=a.var)
            mq = q.mean(axis=0)
        rows = [dict(name=f"glofas_{out['lon'][k]:.3f}_{out['lat'][k]:.3f}",
                     lon=f"{out['lon'][k]:.3f}", lat=f"{out['lat'][k]:.3f}", row=int(out["row"][k]),
                     col=int(out["col"][k]), mean_q=f"{mq[k]:.4g}", kind=out["kind"][k],
                     sea_km=f"{out['sea_km'][k]:.1f}", use=int(out["use"][k]))
                for k in range(out["lon"].size)]
        _write_csv(a.out, rows, ["name", "lon", "lat", "row", "col", "mean_q", "kind", "sea_km", "use"])
        kinds = {k: int((out["kind"] == k).sum()) for k in ("mouth", "near_sea", "inland")}
        print(f"LDD outlets: {kinds} (inland pits get use=0); wrote {a.out}")
        if np.isfinite(mq).any():
            for k in ("mouth", "near_sea", "inland"):
                s = out["kind"] == k
                print(f"   {k:<9} {np.nansum(mq[s]):12,.0f} m3/s ({years[0]})")
        print(f"({time.time() - t0:.0f} s)")
        return
    if not (a.raw and years):
        raise SystemExit("outlets: give --raw and --years (or --ldd)")
    files = [a.raw.format(year=y) for y in years]
    scan = glofas.scan_raw(files, var=a.var, lat_range=a.lat_range)
    out = glofas.coastal_outlets(scan, min_q=a.min_q, lat_range=a.lat_range)
    dup = glofas.flag_coastal_duplicates(out, scan, rel_tol=a.dup_tol, min_q=a.dup_min_q)
    drop = not a.keep_duplicates
    rows = []
    for k in range(out["lon"].size):
        rows.append(dict(name=f"glofas_{out['lon'][k]:.3f}_{out['lat'][k]:.3f}",
                         lon=f"{out['lon'][k]:.3f}", lat=f"{out['lat'][k]:.3f}",
                         row=int(out["row"][k]), col=int(out["col"][k]),
                         mean_q=f"{out['mean_q'][k]:.4g}", n_sea=int(out["n_sea"][k]),
                         kind="coastal", dup_of=int(dup[k]),
                         use=0 if (drop and dup[k] >= 0) else 1))
    _write_csv(a.out, rows, ["name", "lon", "lat", "row", "col", "mean_q", "n_sea", "kind",
                             "dup_of", "use"])
    q = out["mean_q"]
    print(f"{len(rows):,} coastal outlets (mean Q > {a.min_q:g}) carrying {q.sum():,.0f} m3/s "
          f"over {scan.ndays} days; {int((dup >= 0).sum())} along-coast duplicates "
          f"({q[dup >= 0].sum():,.0f} m3/s) {'set to use=0' if drop else 'KEPT (--keep-duplicates)'}")
    for lim in (1, 10, 100, 1000):
        print(f"   cells >= {lim:>5} m3/s: {int((q >= lim).sum()):>7,}  "
              f"({100 * q[q >= lim].sum() / max(q.sum(), 1e-9):5.1f} % of the water)")
    print(f"wrote {a.out}")
    if a.dead_ends:
        dead = glofas.dead_end_candidates(scan, min_q=a.dead_min_q, frac=a.dead_frac,
                                          max_sea_km=a.dead_max_sea_km)
        for r in dead:
            r["name"] = f"deadend_{r['lon']:.3f}_{r['lat']:.3f}"
            for k in ("lon", "lat", "sea_lon", "sea_lat"):
                r[k] = f"{r[k]:.3f}"
            r["mean_q"], r["sea_km"] = f"{r['mean_q']:.4g}", f"{r['sea_km']:.1f}"
            r["deliver_lon"], r["deliver_lat"] = r["lon"], r["lat"]
            r["use"] = 0
        _write_csv(a.dead_ends, dead, ["name", "lon", "lat", "deliver_lon", "deliver_lat",
                                       "mean_q", "channel_cells", "sea_km", "sea_lon",
                                       "sea_lat", "use"])
        print(f"{len(dead)} dead-end CANDIDATES (>= {a.dead_min_q:g} m3/s, within "
              f"{a.dead_max_sea_km:g} km of the sea) -> {a.dead_ends}. All use=0: triage by hand "
              f"(docs/METHODS.md), set use=1 for real dead ends, pass with --extra.")
    print(f"({time.time() - t0:.0f} s)")


# ---------------------------------------------------------------------- run
def _extra_sources(a, year, mapper=None):
    """--extra CSV: name, lon, lat[, deliver_lon, deliver_lat][, use] read from raw GloFAS."""
    from .readers import glofas
    if not a.extra:
        return None
    rows = [r for r in _read_csv(a.extra) if str(r.get("use", "1")).strip() not in ("0", "")]
    if not rows:
        return None
    lon = [float(r["lon"]) for r in rows]
    lat = [float(r["lat"]) for r in rows]
    dlon = [float(r.get("deliver_lon") or r["lon"]) for r in rows]
    dlat = [float(r.get("deliver_lat") or r["lat"]) for r in rows]
    return glofas.sources_from_raw(a.raw, year, lon, lat, name=[r["name"] for r in rows],
                                   kind=["extra"] * len(rows), deliver_lon=dlon,
                                   deliver_lat=dlat, var=a.var)


def _sources_for(a, year, mapper):
    from .readers import glofas, points, gridded
    if a.source == "glofas-raw":
        rows = [r for r in _read_csv(a.outlets) if str(r.get("use", "1")).strip() not in ("0", "")]
        lon = np.array([float(r["lon"]) for r in rows])
        lat = np.array([float(r["lat"]) for r in rows])
        sn = mapper.snap(lon, lat)                      # read only what can reach this grid
        near = (sn.status == 0) | (sn.snap_km <= 3 * mapper.max_snap_km)
        sel = np.flatnonzero(near)
        s = glofas.sources_from_raw(a.raw, year, lon[sel], lat[sel],
                                    name=[rows[k]["name"] for k in sel],
                                    kind=[rows[k].get("kind", "coastal") for k in sel], var=a.var)
    elif a.source == "glofas-runoff-map":
        s = glofas.runoff_map_sources(a.runoff_map.format(year=year), mapper, year=year,
                                      in_lat_order=a.in_lat_order, max_records=a.max_records,
                                      temp_path=a.temp_map.format(year=year) if a.temp_map else None)
    elif a.source == "table":
        s = points.from_csv(a.outlets, a.series, temp_csv=a.temp_table)
    elif a.source == "gridded":
        s = gridded.from_netcdf(a.field.format(year=year), a.var, lon=a.lon_var, lat=a.lat_var,
                                time=a.time_var, units=a.units)
    else:
        raise SystemExit(f"unknown --source {a.source}")
    ex = _extra_sources(a, year) if a.source.startswith("glofas") else None
    if ex is not None:
        s = s.concat(ex)
    if a.max_records and s.nt > a.max_records:
        s = s.head(a.max_records)
    return _with_temperature(a, s)


_CLIM = {}


def _with_temperature(a, s):
    """Attach river temperature (climatology sampled at the sources), then fill unknowns."""
    from .readers import temperature as rt
    if a.temp_clim:
        if a.temp_clim not in _CLIM:
            _CLIM[a.temp_clim] = rt.TemperatureClimatology.load(a.temp_clim)
        if s.has_temp:      # e.g. a table with some temperatures: keep them, sample the gaps
            t, _ = _CLIM[a.temp_clim].sample(s.lon, s.lat, s.time, max_km=a.temp_max_km)
            s = s.with_temp(np.where(np.isfinite(s.temp), s.temp, t))
        else:
            s = rt.add_temperature(s, _CLIM[a.temp_clim], max_km=a.temp_max_km)
    if s.has_temp and a.temp_fill_km > 0 and not np.isfinite(s.temp).all():
        s, donors = rt.fill_missing(s, max_km=a.temp_fill_km)
        print(f"   temperature: {s.attrs.get('temperature_fill')}")
    return s


def cmd_run(a):
    from .mapping import Mapper, Delivery
    from .writers import write_mitgcm_runoff, write_mitgcm_temperature, write_sparse_netcdf
    grid = _grid(a)
    print(grid.summary())
    mapper = Mapper(grid, a.sigma_km, max_snap_km=a.max_snap_km,
                    exclude_outside=a.exclude_outside_domain, domain_tol_km=a.domain_tol_km,
                    cutoff_sigma=a.cutoff_sigma)
    print(f"sigma {a.sigma_km:g} km = {a.sigma_km / grid.mean_cell_km():.1f} cells; max snap "
          f"{a.max_snap_km:g} km; exclude outside {a.exclude_outside_domain} "
          f"(tolerance {mapper.domain_tol_km:.2f} km)")
    yearly = any("{year}" in (x or "") for x in (a.raw, a.runoff_map, a.field))
    years = _years(a.years) if yearly else [None]
    prev = None
    for y in years:
        t0 = time.time()
        s = _sources_for(a, y, mapper)
        if (prev is not None and prev[0].n == s.n and np.array_equal(prev[0].lon, s.lon)
                and np.array_equal(prev[0].lat, s.lat)):
            d = Delivery(mapper, s, prev[1].snap, prev[1].W)       # same sources: reuse W
        else:
            d = mapper.map_sources(s)
        prev = (s, d)
        lab = str(y) if y is not None else "all"
        print(f"\n== {lab}: {s.n:,} sources, {s.nt} records, {d.n_point:,} receiving cells")
        info = d.report()
        if a.mitgcm_dir:
            os.makedirs(a.mitgcm_dir, exist_ok=True)
            p = os.path.join(a.mitgcm_dir, f"{a.prefix}{a.tag}_{lab}")
            n = write_mitgcm_runoff(p, d)
            print(f"wrote {p} ({n} records of {grid.shape}, >f4, m s-1)")
            if s.has_temp:
                p = os.path.join(a.mitgcm_dir, f"{a.temp_prefix}{a.tag}_{lab}")
                n = write_mitgcm_temperature(p, d)
                print(f"wrote {p} ({n} records, >f4, deg C, 0 where no water)")
        if a.sparse_dir:
            os.makedirs(a.sparse_dir, exist_ok=True)
            p = os.path.join(a.sparse_dir, f"{a.tag}_runoff_{lab}.nc")
            write_sparse_netcdf(p, d, title=f"river runoff {a.tag} {lab}",
                                with_discharge=not a.no_discharge,
                                attrs=dict(source_type=a.source))
            print(f"wrote {p} ({d.n_point:,} points of {grid.size:,} cells)")
        if a.tables_dir:
            os.makedirs(a.tables_dir, exist_ok=True)
            p = os.path.join(a.tables_dir, f"sources_{a.tag}_{lab}.csv")
            st = d.snap.status
            rows = [dict(r, status=int(st[k])) for k, r in enumerate(s.table())]
            _write_csv(p, rows, ["name", "kind", "lon", "lat", "mean_q", "status"]
                       + (["mean_temp"] if s.has_temp else []))
        print(f"({time.time() - t0:.0f} s; delivered/kept {info['accounted']:.12f})", flush=True)


# ---------------------------------------------------------------------- temp-clim
def cmd_temp_clim(a):
    from .readers import temperature as rt
    t0 = time.time()
    if a.bin:
        files = sorted(glob.glob(a.bin))
        if not files:
            raise SystemExit(f"no files match {a.bin}")
        lay = dict(rt.WATTMP)
        if a.grid_spec:
            nlon, nlat, lon0, latf, inc = [float(x) for x in a.grid_spec.split(",")]
            lay.update(nlon=int(nlon), nlat=int(nlat), lon0=lon0, lat_first=latf, inc=inc)
        if a.dtype:
            lay["dtype"] = a.dtype
        if a.south_first:
            lay["north_first"] = False
        if a.celsius:
            lay["kelvin"] = False
        print(f"{len(files)} binary files, layout {lay}")
        clim = rt.TemperatureClimatology.from_binary(files, lat_range=a.lat_range,
                                                     min_good=a.min_good, **lay)
    elif a.nc:
        files = sorted(glob.glob(a.nc))
        if not files:
            raise SystemExit(f"no files match {a.nc}")
        clim = rt.TemperatureClimatology.from_netcdf(files, a.var, lon=a.lon_var, lat=a.lat_var,
                                                     time=a.time_var, lat_range=a.lat_range,
                                                     min_good=a.min_good)
    else:
        raise SystemExit("temp-clim: give --bin or --nc")
    clim.save(a.out)
    print(f"wrote {a.out}: {clim.clim.shape} (slot, lat, lon), {int(clim.valid.sum()):,} usable "
          f"cells ({time.time() - t0:.0f} s)")


# ---------------------------------------------------------------------- dropouts
def cmd_dropouts(a):
    from .checks import scan_dropouts
    tables = {}
    for f in sorted(glob.glob(a.tables)):
        m = re.search(r"_(\d{4})\.csv$", f)
        if not m:
            continue
        rows = [r for r in _read_csv(f) if r.get("status", "0") == "0"]
        tables[int(m.group(1))] = (np.array([float(r["lon"]) for r in rows]),
                                   np.array([float(r["lat"]) for r in rows]),
                                   np.array([float(r["mean_q"]) for r in rows]))
    if len(tables) < 3:
        raise SystemExit(f"need >= 3 yearly tables matching {a.tables}")
    res = scan_dropouts(tables, merge_km=a.merge_km, min_q=a.min_q, gap_frac=a.gap_frac,
                        deficit_km=a.deficit_km)
    print(f"{len(tables)} years; {len(res)} rivers >= {a.min_q:g} m3/s with gap years")
    for r in res:
        print(f"  {r['lon']:9.3f} {r['lat']:8.3f} median {r['median_q']:8.0f}  lost "
              f"{r['lost_frac']:5.2f} {'KEEP' if r['keep'] else 'drop'}  "
              f"{' '.join(map(str, r['gap_years']))}")
    if a.out:
        _write_csv(a.out, [dict(r, gap_years=" ".join(map(str, r["gap_years"]))) for r in res],
                   ["lon", "lat", "median_q", "gap_years", "lost_frac", "keep"])
        print(f"wrote {a.out}")


# ---------------------------------------------------------------------- parser
def build_parser():
    p = argparse.ArgumentParser(prog="discharge2grid", description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = p.add_subparsers(dest="cmd", required=True)

    o = sub.add_parser("outlets", help="GloFAS outlets (LDD or coast rule) + dead-end candidates")
    o.add_argument("--ldd", default=None,
                   help="GloFAS drainage-direction map (netCDF, e.g. glofas_ldd.nc): exact outlets, "
                        "one per river (preferred). Without it the coast rule is used.")
    o.add_argument("--ldd-var", default=None)
    o.add_argument("--raw", default=None, help="raw GloFAS file pattern with {year}")
    o.add_argument("--years", default="", help="years to scan, e.g. 1990,2005,2015 or 2010-2012")
    o.add_argument("--var", default=None, help="discharge variable (default avg_dis/dis24/dis)")
    o.add_argument("--lat-range", type=float, nargs=2, default=None, metavar=("S", "N"))
    o.add_argument("--min-q", type=float, default=0.0, help="keep coastal cells with mean Q > this")
    o.add_argument("--dup-tol", type=float, default=0.1)
    o.add_argument("--dup-min-q", type=float, default=10.0)
    o.add_argument("--keep-duplicates", action="store_true",
                   help="keep along-coast duplicates (default: smaller twin gets use=0; keeping "
                        "them counts rivers like the Mackenzie twice)")
    o.add_argument("--out", required=True, help="outlets CSV")
    o.add_argument("--dead-ends", default=None, help="also write dead-end candidates to this CSV")
    o.add_argument("--dead-min-q", type=float, default=100.0)
    o.add_argument("--dead-frac", type=float, default=0.8)
    o.add_argument("--dead-max-sea-km", type=float, default=250.0)
    o.set_defaults(func=cmd_outlets)

    r = sub.add_parser("run", help="sources -> grid -> outputs")
    g = r.add_argument_group("grid")
    g.add_argument("--grid", help="grid netCDF")
    g.add_argument("--grid-mds", help="MITgcm grid directory (XC, YC, RAC, hFacC)")
    g.add_argument("--grid-lon"); g.add_argument("--grid-lat"); g.add_argument("--grid-area")
    g.add_argument("--grid-mask"); g.add_argument("--grid-depth")
    g.add_argument("--mask-level", type=int, default=0)
    s = r.add_argument_group("sources")
    s.add_argument("--source", required=True,
                   choices=["glofas-raw", "glofas-runoff-map", "table", "gridded"])
    s.add_argument("--years", default="", help="years for {year} patterns")
    s.add_argument("--raw", help="glofas-raw: raw file pattern with {year} (also for --extra)")
    s.add_argument("--var", default=None, help="variable (glofas-raw: avg_dis; gridded: required)")
    s.add_argument("--outlets", help="glofas-raw: outlets CSV (use=1 rows); table: outlets CSV")
    s.add_argument("--extra", help="glofas-*: extra rivers CSV (name,lon,lat[,deliver_lon,deliver_lat][,use])")
    s.add_argument("--runoff-map", help="glofas-runoff-map: legacy file pattern with {year}")
    s.add_argument("--in-lat-order", default="north-to-south",
                   choices=["north-to-south", "south-to-north"])
    s.add_argument("--series", help="table: wide time-series CSV")
    s.add_argument("--field", help="gridded: netCDF pattern (may contain {year})")
    s.add_argument("--units", help="gridded: override the units attribute")
    s.add_argument("--lon-var", default="lon"); s.add_argument("--lat-var", default="lat")
    s.add_argument("--time-var", default="time")
    s.add_argument("--max-records", type=int, default=None, help="testing: first N records only")
    t = r.add_argument_group("temperature (optional)")
    t.add_argument("--temp-clim", help="climatology .npz from `temp-clim`, sampled at the sources")
    t.add_argument("--temp-max-km", type=float, default=None,
                   help="no climatology value beyond this distance from a source (-> filled)")
    t.add_argument("--temp-map", help="glofas-runoff-map: legacy temperature map pattern with {year}")
    t.add_argument("--temp-table", help="table: wide temperature CSV (deg C), same layout as --series")
    t.add_argument("--temp-fill-km", type=float, default=150.0,
                   help="sources with unknown temperature take that of the largest known source "
                        "within this distance (0 = off)")
    t.add_argument("--temp-prefix", default="runoff_temp_",
                   help="MITgcm temperature file prefix (file = PREFIX TAG _YEAR)")
    m = r.add_argument_group("mapping")
    m.add_argument("--sigma-km", type=float, required=True)
    m.add_argument("--max-snap-km", type=float, default=25.0)
    m.add_argument("--exclude-outside-domain", action="store_true")
    m.add_argument("--domain-tol-km", type=float, default=None)
    m.add_argument("--cutoff-sigma", type=float, default=3.5)
    w = r.add_argument_group("output")
    w.add_argument("--tag", required=True)
    w.add_argument("--mitgcm-dir", help="write dense MITgcm runoff files here")
    w.add_argument("--prefix", default="runoff_", help="MITgcm file prefix (file = PREFIX TAG _YEAR)")
    w.add_argument("--sparse-dir", help="write sparse netCDF files here")
    w.add_argument("--no-discharge", action="store_true", help="sparse: omit discharge(time, source)")
    w.add_argument("--tables-dir", help="write per-year source tables (for `dropouts`)")
    r.set_defaults(func=cmd_run)

    c = sub.add_parser("temp-clim", help="daily river-temperature climatology from gridded files")
    c.add_argument("--bin", help="glob of flat binary daily files (default layout: wattmp, see "
                                 "readers/temperature.py)")
    c.add_argument("--grid-spec", help="binary layout NLON,NLAT,LON0,LAT_OF_FIRST_ROW,INC")
    c.add_argument("--dtype", help="binary dtype, e.g. '<f4' (little-endian) or '>f4'")
    c.add_argument("--south-first", action="store_true", help="binary rows run south to north")
    c.add_argument("--celsius", action="store_true", help="binary values are deg C (default Kelvin)")
    c.add_argument("--nc", help="glob of netCDF files with dated records (time, lat, lon)")
    c.add_argument("--var", help="netCDF temperature variable")
    c.add_argument("--lon-var", default="lon"); c.add_argument("--lat-var", default="lat")
    c.add_argument("--time-var", default="time")
    c.add_argument("--lat-range", type=float, nargs=2, default=None, metavar=("S", "N"),
                   help="keep these latitudes only (pad a few degrees beyond your grid)")
    c.add_argument("--min-good", type=int, default=30, help="good values a cell needs to be used")
    c.add_argument("--out", required=True, help="climatology .npz")
    c.set_defaults(func=cmd_temp_clim)

    d = sub.add_parser("dropouts", help="rivers missing in some years")
    d.add_argument("--tables", required=True, help="glob of per-year tables, e.g. 'tab/sources_X_*.csv'")
    d.add_argument("--merge-km", type=float, default=8.0)
    d.add_argument("--min-q", type=float, default=20.0)
    d.add_argument("--gap-frac", type=float, default=0.3)
    d.add_argument("--deficit-km", type=float, default=30.0)
    d.add_argument("--out", default=None)
    d.set_defaults(func=cmd_dropouts)
    return p


def main(argv=None):
    a = build_parser().parse_args(argv)
    if getattr(a, "cmd", None) == "run" and not (a.grid or a.grid_mds):
        raise SystemExit("run: give --grid or --grid-mds")
    a.func(a)


if __name__ == "__main__":
    sys.exit(main())
