"""GloFAS (v4, 0.05 deg, daily) -> SourceSet.

GloFAS gives river discharge (avg_dis, m3 s-1) on every LAND cell of a 0.05 deg grid; the
sea is NaN. A cell's value is the water flowing OUT of it, so along a river the value is
the river's discharge, and summing all land cells would count each river many times.
The ocean receives a river at its OUTLET: the land cell where the channel meets the sea.

Two ways to find outlets, both giving a FIXED list that is then read every year:

(A) From the GloFAS drainage map (LDD, a GloFAS static file) -- exact, preferred:
    out = outlets_from_ldd("glofas_ldd.nc", lat_range=(40, 90))
    Every river ends in exactly one pit cell; pits next to the sea are mouths, pits within
    max_sea_km of the sea are rivers that stop short of it (the Ob), farther pits are
    endorheic. No duplicates, no hand triage of dead ends. (Tested on synthetic maps only in
    v0.1: we did not have the LDD file; see METHODS 2.2.)

(B) From the discharge field alone (no LDD):
    scan = scan_raw(["glofas_original_2015.nc", ...], lat_range=(40, 90))
    out  = coastal_outlets(scan)                  # land cells with a sea neighbour
    dup  = flag_coastal_duplicates(out, scan)     # REMOVE these: same river counted again
    dead = dead_end_candidates(scan)              # rivers that stop inland: TRIAGE BY HAND
    On real GloFAS (2015, 40-90N) the coastal rule WITHOUT duplicate removal gives 38 % more
    water than our production input -- channels running along the coast (Mackenzie, Rhone,
    Po deltas) make 2-3 adjacent coastal cells carry the same river. With duplicates removed
    it is within ~5-10 % (METHODS V2/V4).

    src = sources_from_raw("glofas_original_{year}.nc", 2015, out_lon, out_lat, ...)

A fixed list is what keeps a river from vanishing in the years a year-by-year conversion
happens to miss it (METHODS: dropouts). Dead ends are rivers whose GloFAS channel ends
inland (the Ob ends 9 km from the Gulf of Ob); the coast rule misses them, so they are
found separately -- the automatic test flags real dead ends AND points along big rivers
whose discharge dips downstream, so its list must be checked by hand.

Legacy path (what made the pan_AO / WestAO forcing): read a pre-made "runoff map" (the
GloFAS grid with discharge kept only at outlet cells, m s-1, big-endian float32, one file
per year) with runoff_map_sources(). Kept so results can be reproduced; for new work use
the raw path.
"""
import math
import os

import numpy as np

from ..geometry import EARTH_R_M, gc_km, wrap_lon
from ..sources import SourceSet, daily_axis

# GloFAS v4 global grid (runoff-map files); lat0 is the centre of the SOUTHERNMOST row.
GLO_NLON, GLO_NLAT = 7200, 3000
GLO_LON0, GLO_LAT0, GLO_INC = -179.975, -59.975, 0.05


def glofas_cell_area(lat_deg, inc=GLO_INC):
    """Area (m2) of a GloFAS cell, the formula used to make the pan_AO/WestAO forcing."""
    d = inc * math.pi / 180.0 * EARTH_R_M
    return d * d * np.cos(np.radians(lat_deg))


# ====================================================================== raw netCDF files
def _pick(ds, names):
    for n in names:
        if n in ds.variables:
            return n
    raise KeyError(f"none of {names} in {list(ds.variables)}")


def _open(path, var=None):
    from netCDF4 import Dataset
    ds = Dataset(path)
    vn = var or _pick(ds, ["avg_dis", "dis24", "dis", "discharge", "river_discharge"])
    la = _pick(ds, ["latitude", "lat"])
    lo = _pick(ds, ["longitude", "lon"])
    tn = _pick(ds, ["valid_time", "time"])
    return ds, vn, la, lo, tn


def _times(ds, tn):
    from netCDF4 import num2date
    tv = ds[tn]
    t = num2date(tv[:], tv.units, getattr(tv, "calendar", "standard"),
                 only_use_cftime_datetimes=False, only_use_python_datetimes=True)
    return np.array([np.datetime64(x, "D") for x in t])


def chunk_points(v, js, is_):
    """Values of a (time, lat, lon) netCDF variable at points (js[k], is_[k]) for all times.

    Reads ONE hyperslab per storage chunk that holds any of the points. The obvious
    orthogonal read v[:, js, is_] decompresses a chunk per point and per day: on GloFAS
    files (chunks 29 x 231 x 554, deflate) that was 903 s per year for 50 outlets against
    9.6 s this way, with identical values. Returns a list of float64 arrays (NaN -> 0)."""
    ch = v.chunking()
    cj, ci = (ch[1], ch[2]) if isinstance(ch, (list, tuple)) else (v.shape[1], v.shape[2])
    groups = {}
    for k, (j, i) in enumerate(zip(js, is_)):
        groups.setdefault((j // cj, i // ci), []).append(k)
    out = [None] * len(js)
    for ks in groups.values():
        j0, j1 = min(js[k] for k in ks), max(js[k] for k in ks)
        i0, i1 = min(is_[k] for k in ks), max(is_[k] for k in ks)
        box = np.ma.filled(v[:, j0:j1 + 1, i0:i1 + 1], np.nan).astype(float)
        for k in ks:
            out[k] = np.nan_to_num(box[:, js[k] - j0, is_[k] - i0])
    return out


class RawScan:
    """Time-mean discharge and the sea mask of raw GloFAS file(s) over a latitude range.

    lat, lon : 1-D axes of the scanned window (file order)
    rows     : indices of `lat` in the files' latitude axis
    qmean    : mean discharge over all days scanned (m3 s-1), 0 at sea
    ocean    : True where the value was NaN on every day of every file (= sea)
    """

    def __init__(self, lat, lon, rows, qmean, ocean, ndays, files):
        self.lat, self.lon, self.rows = lat, lon, rows
        self.qmean, self.ocean, self.ndays, self.files = qmean, ocean, ndays, files

    @property
    def is_global_lon(self):
        return abs((self.lon.size * abs(float(np.median(np.diff(self.lon))))) - 360.0) < 1e-6


def scan_raw(paths, var=None, lat_range=None, pad_deg=0.1, verbose=True):
    """Mean discharge + sea mask over one or more raw yearly files, reading in
    chunk-aligned blocks (whole time chunks x whole row chunks x all longitudes)."""
    if isinstance(paths, str):
        paths = [paths]
    lat = lon = rows = None
    acc = finite = None
    ndays = 0
    for p in paths:
        ds, vn, la, lo, tn = _open(p, var)
        try:
            v = ds[vn]
            flat = np.asarray(ds[la][:], float)
            flon = wrap_lon(np.asarray(ds[lo][:], float))
            if lat is None:
                if lat_range is None:
                    sel = np.arange(flat.size)
                else:
                    lo_, hi_ = lat_range[0] - pad_deg, lat_range[1] + pad_deg
                    sel = np.flatnonzero((flat >= lo_) & (flat <= hi_))
                if sel.size == 0:
                    raise ValueError(f"no latitude rows in {lat_range}")
                rows = np.arange(sel.min(), sel.max() + 1)
                lat, lon = flat[rows], flon
                acc = np.zeros((rows.size, lon.size))
                finite = np.zeros((rows.size, lon.size), dtype=np.int32)
            elif flat.size <= rows.max() or not np.allclose(flat[rows], lat):
                raise ValueError(f"{p}: latitude axis differs from the first file")
            ch = v.chunking()
            ct, cj = (ch[0], ch[1]) if isinstance(ch, (list, tuple)) else (v.shape[0], v.shape[1])
            nt = v.shape[0]
            r0, r1 = int(rows[0]), int(rows[-1]) + 1
            b0 = (r0 // cj) * cj
            for t0 in range(0, nt, ct):
                t1 = min(t0 + ct, nt)
                for j0 in range(b0, r1, cj):
                    a, b = max(j0, r0), min(j0 + cj, r1)
                    blk = np.ma.filled(v[t0:t1, a:b, :], np.nan).astype(np.float64)
                    ok = np.isfinite(blk)
                    acc[a - r0:b - r0] += np.where(ok, blk, 0.0).sum(axis=0)
                    finite[a - r0:b - r0] += ok.sum(axis=0)
            ndays += nt
            if verbose:
                print(f"  scanned {os.path.basename(p)}: {nt} days, rows {r0}..{r1 - 1} "
                      f"(lat {lat.min():.2f}..{lat.max():.2f})", flush=True)
        finally:
            ds.close()
    qmean = acc / max(ndays, 1)
    return RawScan(lat, lon, rows, qmean, finite == 0, ndays, list(paths))


def _neighbour_any(mask, wrap):
    """8-neighbour 'any' of a boolean mask (lon axis wraps if global)."""
    from scipy.ndimage import maximum_filter
    if wrap:
        m = np.concatenate([mask[:, -1:], mask, mask[:, :1]], axis=1)
        f = maximum_filter(m.astype(np.uint8), size=3, mode="constant")[:, 1:-1]
    else:
        f = maximum_filter(mask.astype(np.uint8), size=3, mode="constant")
    return f > 0


def coastal_outlets(scan, min_q=0.0, lat_range=None):
    """Land cells with at least one sea neighbour (8-connectivity) and mean discharge > min_q.

    This is the "coast-only" outlet rule. It keeps every coastal cell, so small coastal
    catchments are included (they carry real water too). Returns a dict of arrays:
    lon, lat, row, col, mean_q, n_sea (sea neighbours), kind='coastal'."""
    wrap = scan.is_global_lon
    sea_nb = _neighbour_any(scan.ocean, wrap)
    land = ~scan.ocean
    cand = land & sea_nb & (scan.qmean > min_q)
    if lat_range is not None:
        cand &= ((scan.lat >= lat_range[0]) & (scan.lat <= lat_range[1]))[:, None]
    jj, ii = np.nonzero(cand)
    # number of sea neighbours (diagnostic)
    from scipy.ndimage import convolve
    k = np.ones((3, 3), dtype=np.int32)
    k[1, 1] = 0
    oc = scan.ocean.astype(np.int32)
    if wrap:
        oc2 = np.concatenate([oc[:, -1:], oc, oc[:, :1]], axis=1)
        nsea = convolve(oc2, k, mode="constant")[:, 1:-1]
    else:
        nsea = convolve(oc, k, mode="constant")
    order = np.argsort(-scan.qmean[jj, ii], kind="stable")
    jj, ii = jj[order], ii[order]
    return dict(lon=scan.lon[ii], lat=scan.lat[jj], row=scan.rows[jj], col=ii,
                mean_q=scan.qmean[jj, ii], n_sea=nsea[jj, ii],
                kind=np.full(jj.size, "coastal", dtype=object))


def flag_coastal_duplicates(outlets, scan, rel_tol=0.1, min_q=10.0):
    """Coastal cells that are probably the SAME river: an 8-neighbour coastal outlet whose
    mean discharge is within rel_tol and >= min_q (a channel running along the coast makes
    every cell of it 'coastal', and each carries the river's full discharge).
    Returns an int array: -1 = not a duplicate, else the index of the larger twin.
    Nothing is removed; decide after looking (METHODS: duplicates)."""
    key = {(int(r), int(c)): k for k, (r, c) in enumerate(zip(outlets["row"], outlets["col"]))}
    q = outlets["mean_q"]
    dup = np.full(q.size, -1, dtype=np.int64)
    nlon = scan.lon.size
    for k, (r, c) in enumerate(zip(outlets["row"], outlets["col"])):
        if q[k] < min_q:
            continue
        best = -1
        for dr in (-1, 0, 1):
            for dc in (-1, 0, 1):
                if dr == 0 and dc == 0:
                    continue
                cc = (c + dc) % nlon if scan.is_global_lon else c + dc
                j = key.get((int(r + dr), int(cc)))
                if j is None or q[j] < q[k]:
                    continue
                if abs(q[j] - q[k]) <= rel_tol * q[j] and (q[j] > q[k] or j < k):
                    if best < 0 or q[j] > q[best]:
                        best = j
        dup[k] = best
    return dup


def dead_end_candidates(scan, min_q=100.0, frac=0.8, max_sea_km=250.0, window=150):
    """Rivers whose channel never reaches the sea ("dead ends").

    For each local discharge maximum >= min_q on land (largest first), the channel is the
    connected set of land cells carrying >= frac of its discharge (within +-window cells).
    If that set touches no cell next to the sea, the river ends inland: its largest cell is
    reported with the distance to the nearest sea cell. Rivers ending farther than
    max_sea_km from the sea are endorheic or end in lakes and are not listed.

    THE LIST NEEDS HAND TRIAGE. It contains (a) real dead ends (the Ob, Anadyr, Kuskokwim,
    Yana, Thelon ...), (b) points inside deltas and estuaries whose discharge falls slightly
    downstream (the Yukon at -159.975/62.275 is NOT a dead end), (c) duplicates of the same
    river from neighbouring seeds. Returns a list of dicts sorted by discharge."""
    from scipy.ndimage import maximum_filter, label
    q, ocean = scan.qmean, scan.ocean
    wrap = scan.is_global_lon
    sea_nb = _neighbour_any(ocean, wrap)
    qf = np.where(ocean, -1.0, q)
    seeds = ~ocean & (q >= min_q) & (qf >= maximum_filter(qf, size=3, mode="nearest"))
    sj, si = np.nonzero(seeds)
    order = np.argsort(q[sj, si])[::-1]
    done = np.zeros(q.shape, bool)
    found, seen = [], set()
    for t in order:
        j, i = sj[t], si[t]
        if done[j, i]:
            continue
        j0, j1 = max(j - window, 0), min(j + window + 1, q.shape[0])
        i0, i1 = max(i - window, 0), min(i + window + 1, q.shape[1])
        sub = (q[j0:j1, i0:i1] >= frac * q[j, i]) & ~ocean[j0:j1, i0:i1]
        lab, _ = label(sub, structure=np.ones((3, 3)))
        comp = lab == lab[j - j0, i - i0]
        done[j0:j1, i0:i1] |= comp
        if (comp & sea_nb[j0:j1, i0:i1]).any():
            continue
        cq = np.where(comp, q[j0:j1, i0:i1], -1.0)
        tj, ti = np.unravel_index(int(np.argmax(cq)), cq.shape)
        tj, ti = tj + j0, ti + i0
        if (tj, ti) in seen:
            continue
        seen.add((tj, ti))
        w0, w1 = max(tj - window, 0), min(tj + window + 1, q.shape[0])
        v0, v1 = max(ti - window, 0), min(ti + window + 1, q.shape[1])
        oj, oi = np.nonzero(ocean[w0:w1, v0:v1])
        if not oj.size:
            continue
        d = gc_km(scan.lon[ti], scan.lat[tj], scan.lon[oi + v0], scan.lat[oj + w0])
        k = int(np.argmin(d))
        if d[k] > max_sea_km:
            continue
        found.append(dict(lon=float(scan.lon[ti]), lat=float(scan.lat[tj]), row=int(scan.rows[tj]),
                          col=int(ti), mean_q=float(q[tj, ti]), channel_cells=int(comp.sum()),
                          sea_km=float(d[k]), sea_lon=float(scan.lon[oi[k] + v0]),
                          sea_lat=float(scan.lat[oj[k] + w0])))
    found.sort(key=lambda r: -r["mean_q"])
    return found


# PCRaster / LISFLOOD local drain direction: numeric-keypad codes, 5 = pit.
#   7 8 9        NW  N  NE
#   4 5 6   =    W  pit E
#   1 2 3        SW  S  SE
_LDD_STEP = {1: (-1, -1), 2: (-1, 0), 3: (-1, 1), 4: (0, -1), 6: (0, 1),
             7: (1, -1), 8: (1, 0), 9: (1, 1)}          # (northward rows, eastward cols)


def outlets_from_ldd(path, lat_range=None, var=None, max_sea_km=250.0, pad_deg=0.1):
    """River mouths from a GloFAS/LISFLOOD drainage-direction map (LDD).

    An outlet is a land cell that drains to the sea: a pit (code 5), or a cell whose
    downstream neighbour is not land (missing in the LDD). Each river has exactly one, so
    there are no along-coast duplicates. Kinds:
        'mouth'       the outlet touches the sea (8-neighbour)
        'near_sea'    a pit within max_sea_km of the sea: a river that stops short of the
                      coast in GloFAS (the Ob) -- delivered at the pit; the snap to model
                      water moves it to the coast
        'inland'      a pit farther from the sea: endorheic / lake basin -> use = 0
    Returns a dict of arrays like coastal_outlets (lon, lat, row, col, kind, sea_km, use).
    Discharge is then read at these cells from the raw files (sources_from_raw).
    Orientation is taken from the file's latitude axis, so north-up and south-up maps both work.
    """
    from netCDF4 import Dataset
    from scipy.ndimage import distance_transform_edt
    with Dataset(path) as ds:
        vn = var or _pick(ds, ["ldd", "Band1", "LDD", "ldd_glofas"])
        la = _pick(ds, ["lat", "latitude", "y"])
        lo = _pick(ds, ["lon", "longitude", "x"])
        flat = np.asarray(ds[la][:], float)
        flon = wrap_lon(np.asarray(ds[lo][:], float))
        if lat_range is None:
            rows = np.arange(flat.size)
        else:
            sel = np.flatnonzero((flat >= lat_range[0] - pad_deg) & (flat <= lat_range[1] + pad_deg))
            rows = np.arange(sel.min(), sel.max() + 1)
        L = np.ma.filled(ds[vn][rows.min():rows.max() + 1, :], 0).astype(np.int16)
    lat = flat[rows]
    north = 1 if (lat.size > 1 and lat[1] > lat[0]) else -1      # +1 if row index grows north
    land = (L >= 1) & (L <= 9)
    ny, nx = L.shape
    glob_lon = abs(nx * abs(float(np.median(np.diff(flon)))) - 360.0) < 1e-6
    jj, ii = np.nonzero(land)
    code = L[jj, ii]
    dj = np.zeros(code.size, dtype=np.int64)
    di = np.zeros(code.size, dtype=np.int64)
    for c, (n, e) in _LDD_STEP.items():
        m = code == c
        dj[m], di[m] = n * north, e
    tj, ti = jj + dj, ii + di
    if glob_lon:
        ti %= nx
    inside = (tj >= 0) & (tj < ny) & (ti >= 0) & (ti < nx)
    tland = np.zeros(code.size, dtype=bool)
    tland[inside] = land[tj[inside], ti[inside]]
    pit = code == 5
    outlet = pit | (~pit & ~tland & inside)
    jj, ii = jj[outlet], ii[outlet]
    # distance to the sea (km, approximate: cell sizes at each row's latitude)
    sea = ~land
    dlat_km = abs(float(np.median(np.diff(lat)))) * 111.2 if lat.size > 1 else 5.5
    dlon_deg = abs(float(np.median(np.diff(flon))))
    _, (ej, ei) = distance_transform_edt(~sea, return_indices=True)
    dy = (ej[jj, ii] - jj) * dlat_km                 # approximate (index space, no lon wrap)
    dx = (ei[jj, ii] - ii) * dlon_deg * 111.2 * np.cos(np.radians(lat[jj]))
    sea_km = np.hypot(dx, dy)
    from scipy.ndimage import maximum_filter
    if glob_lon:
        s2 = np.concatenate([sea[:, -1:], sea, sea[:, :1]], axis=1)
        touch = maximum_filter(s2.astype(np.uint8), size=3, mode="constant")[:, 1:-1] > 0
    else:
        touch = maximum_filter(sea.astype(np.uint8), size=3, mode="constant") > 0
    kind = np.where(touch[jj, ii], "mouth", np.where(sea_km <= max_sea_km, "near_sea", "inland"))
    if lat_range is not None:
        keep = (lat[jj] >= lat_range[0]) & (lat[jj] <= lat_range[1])
        jj, ii, kind, sea_km = jj[keep], ii[keep], kind[keep], sea_km[keep]
    return dict(lon=flon[ii], lat=lat[jj], row=rows[jj], col=ii, kind=kind.astype(object),
                sea_km=sea_km, use=(kind != "inland").astype(int))


def read_outlet_series(path, lon, lat, var=None):
    """Daily discharge (nt, n) at the GloFAS cells nearest to (lon, lat), and the time axis."""
    ds, vn, la, lo, tn = _open(path, var)
    try:
        flat = np.asarray(ds[la][:], float)
        flon = wrap_lon(np.asarray(ds[lo][:], float))
        js = [int(np.argmin(np.abs(flat - y))) for y in np.asarray(lat, float)]
        is_ = [int(np.argmin(np.abs(flon - x))) for x in wrap_lon(np.asarray(lon, float))]
        cols = chunk_points(ds[vn], js, is_)
        time = _times(ds, tn)
    finally:
        ds.close()
    q = np.stack(cols, axis=1) if cols else np.zeros((time.size, 0))
    return q, time


def sources_from_raw(raw_pattern, year, lon, lat, name=None, kind=None,
                     deliver_lon=None, deliver_lat=None, var=None):
    """SourceSet for one year from a FIXED outlet list: discharge read at (lon, lat),
    delivered at (deliver_lon, deliver_lat) if given (e.g. a dead end delivered at the
    coast, or a river moved to the bay it really drains to), else at (lon, lat)."""
    q, time = read_outlet_series(raw_pattern.format(year=year), lon, lat, var=var)
    dl = np.asarray(lon if deliver_lon is None else deliver_lon, float)
    dt = np.asarray(lat if deliver_lat is None else deliver_lat, float)
    return SourceSet(dl, dt, q, time, name=name, kind=kind,
                     attrs=dict(dataset="GloFAS raw avg_dis at fixed outlets",
                                file=raw_pattern.format(year=year)))


# ====================================================================== legacy runoff map
def _band(grid_lat, margin_deg):
    lo = float(np.min(grid_lat)) - margin_deg
    hi = float(np.max(grid_lat)) + margin_deg
    j0 = int(math.floor((lo - GLO_LAT0) / GLO_INC))
    j1 = int(math.ceil((hi - GLO_LAT0) / GLO_INC))
    return max(0, min(j0, GLO_NLAT - 1)), max(0, min(j1, GLO_NLAT - 1))


def cell_area_from_file(path, lat, lon, var=None):
    """Cell areas (m2) at (lat, lon) from a gridded area file, e.g. GloFAS's own
    GLOFAS_pixarea_Global_03min.nc (nearest cell; orientation read from the file)."""
    from netCDF4 import Dataset
    with Dataset(path) as ds:
        vn = var or _pick(ds, ["Band1", "pixarea", "area", "cell_area"])
        flat = np.asarray(ds[_pick(ds, ["lat", "latitude", "y"])][:], float)
        flon = wrap_lon(np.asarray(ds[_pick(ds, ["lon", "longitude", "x"])][:], float))
        js = np.array([int(np.argmin(np.abs(flat - y))) for y in np.asarray(lat, float)])
        is_ = np.array([int(np.argmin(np.abs(flon - x))) for x in wrap_lon(np.asarray(lon, float))])
        if js.size == 0:
            return np.zeros(0)
        blk = np.ma.filled(ds[vn][js.min():js.max() + 1, :], np.nan).astype(float)
    return blk[js - js.min(), is_]


def runoff_map_sources(path, mapper, year=None, margin_deg=1.0, in_lat_order="north-to-south",
                       scale=1.0, halo_factor=3.0, max_records=None, area_file=None, temp_path=None,
                       verbose=True):
    """SourceSet from a legacy GloFAS "runoff map" year file (>f4, (nd, 3000, 7200), m s-1 per
    0.05 deg cell, discharge only at outlet cells; the input of the pan_AO/WestAO forcing).

    Every GloFAS cell in the grid's latitude band (+margin) is a candidate source, exactly as
    in glofas_to_mitgcm.py; only the cells that are kept by the mapper (or near enough to be
    reported as dropped) AND carry water on some day become sources, in the same order, so
    the runoff from this SourceSet reproduces that script's runoff files.

    m s-1 -> m3 s-1 uses the spherical cell area of glofas_to_mitgcm.py by default (exact
    reproduction). The maps were made by dividing discharge by GloFAS's own pixel areas
    (GLOFAS_pixarea_Global_03min.nc), which are ~0.7 % larger at these latitudes; pass
    area_file= that file to recover the raw discharge instead (METHODS V2).

    temp_path: the matching legacy temperature map (GloFas_temp_<year>, deg C on the same
    grid; made by make_glofas_river_temp_forcing.py). Its values at the source cells become
    SourceSet.temp, read exactly as glofas_to_mitgcm.py reads them (record count = the
    shorter of the two files)."""
    rec = GLO_NLON * GLO_NLAT * 4
    nb = os.path.getsize(path)
    if nb % rec:
        raise ValueError(f"{path}: {nb} bytes is not a whole number of 3000x7200 >f4 records")
    nd = nb // rec
    if max_records:
        nd = min(nd, int(max_records))
    mm = np.memmap(path, dtype=">f4", mode="r", shape=(nb // rec, GLO_NLAT, GLO_NLON))
    tm = None
    if temp_path is not None:
        nbt = os.path.getsize(temp_path)
        if nbt % rec:
            raise ValueError(f"{temp_path}: not a whole number of 3000x7200 >f4 records")
        tm = np.memmap(temp_path, dtype=">f4", mode="r", shape=(nbt // rec, GLO_NLAT, GLO_NLON))
        nd = min(nd, nbt // rec)
    j0, j1 = _band(mapper.grid.lat, margin_deg)
    slat = GLO_LAT0 + GLO_INC * np.arange(j0, j1 + 1)
    slon = GLO_LON0 + GLO_INC * np.arange(GLO_NLON)
    LON, LAT = np.meshgrid(slon, slat)
    snap = mapper.snap(LON.ravel(), LAT.ravel(), halo_factor=halo_factor)
    cand = np.flatnonzero((snap.status == 0) | (snap.snap_km <= mapper.max_snap_km * halo_factor))
    row, col = cand // GLO_NLON, cand % GLO_NLON
    nts = in_lat_order == "north-to-south"

    def band_day(k, src=None):
        src = mm if src is None else src
        if nts:
            b = np.asarray(src[k, GLO_NLAT - 1 - j1:GLO_NLAT - j0, :], dtype=np.float64)[::-1, :]
        else:
            b = np.asarray(src[k, j0:j1 + 1, :], dtype=np.float64)
        return b

    anyq = np.zeros(cand.size, dtype=bool)                   # pass 1: which ever carry water
    for k in range(nd):
        v = band_day(k)[row, col]
        anyq |= np.isfinite(v) & (v > 0)
    use = cand[anyq]
    r2, c2 = use // GLO_NLON, use % GLO_NLON
    area = (glofas_cell_area(LAT.ravel()[use]) if area_file is None
            else cell_area_from_file(area_file, LAT.ravel()[use], LON.ravel()[use]))
    q = np.zeros((nd, use.size))
    temp = None if tm is None else np.zeros((nd, use.size))
    for k in range(nd):                                      # pass 2: their daily discharge
        q[k] = band_day(k)[r2, c2] * area * scale            # m/s * m2 -> m3/s, as production
        if tm is not None:
            temp[k] = band_day(k, tm)[r2, c2]                # deg C, as read by production
    del mm, tm
    if verbose:
        print(f"  {os.path.basename(path)}: band rows {j0}..{j1}, {cand.size:,} candidate cells, "
              f"{use.size:,} carry water", flush=True)
    y = year if year is not None else int(os.path.basename(path)[-4:])
    return SourceSet(LON.ravel()[use], LAT.ravel()[use], q, daily_axis(y, nd),
                     kind=np.full(use.size, "runoff_map", dtype=object),
                     attrs=dict(dataset="GloFAS runoff map (legacy)", file=path), temp=temp)
