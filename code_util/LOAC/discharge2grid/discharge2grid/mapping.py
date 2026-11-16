"""Point sources -> model grid: snap to water, spread over a Gaussian footprint, conserve volume.

This is the core that is independent of the discharge dataset. It is a re-packaging of
the regridding in glofas_to_mitgcm.py (the code that made the pan_AO and WestAO forcing);
given the same sources it reproduces that code's runoff fields (see docs/METHODS.md, V1).

1. SNAP. A spherical KD-tree over the WET cells' 3-D unit vectors. One query gives each
   source its nearest wet cell and the great-circle distance to it.
     - distance > max_snap_km                    -> dropped ("too far"; reported)
     - exclude_outside and the nearest model cell, WET OR DRY, is farther than
       domain_tol_km                             -> dropped ("outside the domain")
   The second test is what lets max_snap_km be large (95 km for unresolved estuaries and
   fjords) without pulling in rivers from beyond a regional grid's edge: a river inside
   the domain is near SOME cell (land counts) however far it is from water.

2. SPREAD. Each kept source's water is spread over the wet cells within cutoff*sigma of its
   destination cell with weights g = exp(-d^2 / 2 sigma^2), d the great-circle distance
   from the destination cell (km). Weights are normalised so that sum(t * area) == 1,
   hence the runoff R = q * t (m/s) integrates to exactly q (m3/s): VOLUME IS CONSERVED
   exactly, even where the footprint is clipped by land or the domain edge (clipping
   concentrates the water instead of losing it; see the footprint-area diagnostic).

3. OPERATOR. All of it is one sparse matrix W (n_wet x n_source, units m-2):
        R_wet(t) = W @ q(t)          [m s-1]
   W is time-invariant, so a 37-year run builds it once. The rows of W that are not
   empty are the only cells that ever receive water -- the "points" of the sparse output.
"""
import math

import numpy as np

from .geometry import to_xyz, chord_to_km, km_to_chord

KEPT, TOO_FAR, OUTSIDE = 0, 1, 2
STATUS_NAMES = {KEPT: "kept", TOO_FAR: "too_far", OUTSIDE: "outside_domain"}


class Snap:
    """Result of Mapper.snap for n sources."""

    def __init__(self, dest, snap_km, status):
        self.dest = dest            # wet index of the destination cell, -1 if dropped
        self.snap_km = snap_km      # distance to the nearest wet cell (km), all sources
        self.status = status        # KEPT / TOO_FAR / OUTSIDE

    @property
    def kept(self):
        return self.status == KEPT


class Mapper:
    """Snap + Gaussian footprint on one grid.

    Parameters
    ----------
    grid           : ModelGrid
    sigma_km       : Gaussian footprint scale (km). ~ the grid spacing or more; see METHODS.
    max_snap_km    : drop sources farther than this from any wet cell
    exclude_outside: drop sources not ON the grid (nearest cell, wet or dry, > domain_tol_km)
    domain_tol_km  : default 1.5 x the mean wet-cell size
    cutoff_sigma   : footprint truncated at cutoff_sigma * sigma_km (3.5 in production)
    depth_pref_m   : optional preference for deep columns: weight *= min(1, depth/D)**power
                     (needs grid.depth; 0 = off, as in production)
    """

    def __init__(self, grid, sigma_km, max_snap_km=25.0, exclude_outside=False,
                 domain_tol_km=None, cutoff_sigma=3.5, depth_pref_m=0.0, depth_pref_power=2.0):
        from scipy.spatial import cKDTree
        if sigma_km <= 0:
            raise ValueError("sigma_km must be > 0")
        self.grid = grid
        self.sigma_km = float(sigma_km)
        self.max_snap_km = float(max_snap_km)
        self.exclude_outside = bool(exclude_outside)
        self.domain_tol_km = (1.5 * grid.mean_cell_km() if domain_tol_km is None
                              else float(domain_tol_km))
        self.cutoff_sigma = float(cutoff_sigma)
        self.depth_pref_m = float(depth_pref_m)
        self.depth_pref_power = float(depth_pref_power)
        self.wxyz = grid.wet_xyz()
        self.tree = cKDTree(self.wxyz)
        self.area_w = grid.wet_area()
        self.depth_w = grid.wet_depth()
        if self.depth_pref_m > 0 and self.depth_w is None:
            raise ValueError("depth_pref_m needs a grid with depth")
        self._gtree = None
        self._tpl = {}

    def describe(self):
        return dict(sigma_km=self.sigma_km, max_snap_km=self.max_snap_km,
                    exclude_outside=int(self.exclude_outside), domain_tol_km=self.domain_tol_km,
                    cutoff_sigma=self.cutoff_sigma, depth_pref_m=self.depth_pref_m)

    # ------------------------------------------------------------------ 1. snap
    def snap(self, lon, lat, halo_factor=3.0):
        """Nearest wet cell for each source, with the keep/drop decision."""
        from scipy.spatial import cKDTree
        sxyz = to_xyz(np.asarray(lon, float).ravel(), np.asarray(lat, float).ravel())
        n = sxyz.shape[0]
        if n == 0:
            return Snap(np.zeros(0, np.int64), np.zeros(0), np.zeros(0, np.int8))
        d, k = self.tree.query(sxyz, k=1, workers=-1)
        dkm = chord_to_km(d)
        in_dom = np.ones(n, dtype=bool)
        if self.exclude_outside:
            if self._gtree is None:
                self._gtree = cKDTree(to_xyz(self.grid.lon.ravel(), self.grid.lat.ravel()))
            cand = np.flatnonzero(dkm <= self.max_snap_km * halo_factor)
            if cand.size:
                dg, _ = self._gtree.query(sxyz[cand], k=1, workers=-1,
                                          distance_upper_bound=km_to_chord(self.domain_tol_km))
                in_dom[cand] = np.isfinite(dg)
        status = np.full(n, TOO_FAR, dtype=np.int8)
        status[(dkm <= self.max_snap_km) & in_dom] = KEPT
        status[~in_dom] = OUTSIDE
        dest = np.where(status == KEPT, k, -1).astype(np.int64)
        return Snap(dest, dkm, status)

    # ------------------------------------------------------------------ 2. spread
    def template(self, kw):
        """(wet indices, weights t) of the footprint centred on wet cell kw.
        sum(t * area) == 1. Cached per destination cell."""
        tpl = self._tpl.get(kw)
        if tpl is not None:
            return tpl
        idx = self.tree.query_ball_point(self.wxyz[kw],
                                         km_to_chord(self.cutoff_sigma * self.sigma_km))
        idx = np.asarray(idx, dtype=np.int64)
        dkm = chord_to_km(np.linalg.norm(self.wxyz[idx] - self.wxyz[kw], axis=1))
        g = np.exp(-(dkm * dkm) / (2.0 * self.sigma_km * self.sigma_km))
        if self.depth_pref_m > 0.0:
            g = g * np.clip(self.depth_w[idx] / self.depth_pref_m, 0.0, 1.0) ** self.depth_pref_power
        a = self.area_w[idx]
        denom = float((g * a).sum())
        if denom <= 0.0:
            raise RuntimeError(f"empty footprint at wet cell {kw}")
        t = g / denom
        self._tpl[kw] = (idx, t)
        return idx, t

    def footprint_area_m2(self, kw):
        """Participation-ratio area of a footprint (an unclipped Gaussian gives 4 pi sigma^2)."""
        idx, t = self.template(kw)
        return 1.0 / float((t * t * self.area_w[idx]).sum())

    # ------------------------------------------------------------------ 3. operator
    def operator(self, snap):
        """Sparse W (n_wet x n_source, m-2): R_wet = W @ q. Dropped sources get empty columns."""
        from scipy.sparse import csr_matrix
        rows, cols, vals = [], [], []
        for s in np.flatnonzero(snap.kept):
            idx, t = self.template(int(snap.dest[s]))
            rows.append(idx)
            cols.append(np.full(idx.size, s, dtype=np.int64))
            vals.append(t)
        n_src = snap.dest.size
        if rows:
            rows, cols, vals = np.concatenate(rows), np.concatenate(cols), np.concatenate(vals)
        else:
            rows = cols = np.zeros(0, np.int64)
            vals = np.zeros(0)
        W = csr_matrix((vals, (rows, cols)), shape=(self.grid.n_wet, n_src))
        W.sum_duplicates()
        W.sort_indices()       # column order inside each row = source order (deterministic sums)
        return W

    def map_sources(self, sources):
        """Snap + operator for a SourceSet -> Delivery."""
        snap = self.snap(sources.lon, sources.lat)
        return Delivery(self, sources, snap, self.operator(snap))


class Delivery:
    """A SourceSet mapped onto a grid. Produces runoff on the receiving cells ("points").

    points   : flat grid indices of the cells that receive water (fixed in time)
    Wp       : (n_point x n_source) operator restricted to those cells, m-2
    """

    def __init__(self, mapper, sources, snap, W):
        self.mapper, self.sources, self.snap, self.W = mapper, sources, snap, W
        g = mapper.grid
        nz = np.flatnonzero(np.diff(W.indptr) > 0)             # wet rows with any weight
        self.point_wet = nz                                     # wet indices
        self.points = g.wet_flat[nz]                            # flat grid indices
        self.Wp = W[nz]
        self.point_area = g.wet_area()[nz]

    @property
    def n_point(self):
        return int(self.points.size)

    def runoff_points(self, k):
        """Runoff (m s-1) on the points for record k, float64."""
        return self.Wp @ self.sources.q[k]

    def runoff_wet(self, k):
        """Runoff (m s-1) on all wet cells for record k, float64."""
        return self.W @ self.sources.q[k]

    def runoff_dense(self, k):
        """Runoff (m s-1) on the full grid for record k, float64, grid shape."""
        return self.mapper.grid.to_dense(self.runoff_wet(k))

    # ------------------------------------------------------------------ temperature
    def temperature_points(self, k, fill=np.nan, return_runoff=False):
        """River temperature (deg C) on the points for record k.

        Temperature is INTENSIVE, so it is not spread like water. Each cell gets the
        flux-weighted temperature of the water that reaches it:
            T_c = sum_i (q_i t_ic T_i) / sum_i (q_i t_ic)
        so heat is conserved: sum_c R_c T_c A_c = sum_i q_i T_i. Water with unknown
        temperature (NaN) is left out of both sums; a cell whose water is ALL of unknown
        temperature, or that gets no water on this record, gets `fill`.
        The products are formed as (q*t)*T and summed in source order, exactly as in
        glofas_to_mitgcm.py, so the legacy path reproduces its temperature files."""
        from scipy.sparse import csr_matrix
        if not self.sources.has_temp:
            raise ValueError("these sources have no temperature")
        W = self.Wp
        q = self.sources.q[k]
        T = self.sources.temp[k]
        col = W.indices
        qt = W.data * q[col]                                    # q_i t_ic
        known = np.isfinite(T)[col] & (qt > 0)
        a = np.where(known, qt, 0.0)
        b = np.where(known, qt * np.where(known, T[col], 0.0), 0.0)   # (q t) T
        one = np.ones(self.sources.n)
        shp = W.shape
        Rk = csr_matrix((a, col, W.indptr), shape=shp) @ one     # water of known temperature
        H = csr_matrix((b, col, W.indptr), shape=shp) @ one
        out = np.full(self.n_point, fill, dtype=float)
        has = Rk > 0
        out[has] = H[has] / Rk[has]
        if return_runoff:
            return out, csr_matrix((qt, col, W.indptr), shape=shp) @ one
        return out

    def temperature_dense(self, k, fill=0.0):
        """River temperature on the full grid for record k (fill where no water; 0 as in
        the production files), grid shape."""
        g = self.mapper.grid
        out = np.full(g.size, fill, dtype=float)
        tp = self.temperature_points(k, fill=fill)
        out[self.points] = tp
        return out.reshape(g.shape)

    # ------------------------------------------------------------------ diagnostics
    def report(self, max_list=8, long_snap_km=25.0, file=None):
        """Print what was kept, dropped and how the water lands. Returns a dict of numbers."""
        pr = (lambda *a, **k: print(*a, **k, file=file))
        s, sn, m = self.sources, self.snap, self.mapper
        mq = s.mean_q()
        kept = sn.kept
        q_kept = mq[kept].sum()
        q_far = mq[sn.status == TOO_FAR].sum()
        q_out = mq[sn.status == OUTSIDE].sum()
        # conservation over all records (exact by construction; this checks the code)
        tot_in = float(s.q[:, kept].sum())
        tot_out = 0.0
        for k in range(s.nt):
            tot_out += float((self.runoff_points(k) * self.point_area).sum())
        acct = tot_out / tot_in if tot_in > 0 else float("nan")
        dests = np.unique(sn.dest[kept])
        areas = np.array([m.footprint_area_m2(int(d)) for d in dests]) if dests.size else np.zeros(0)
        ideal = 4.0 * math.pi * (m.sigma_km * 1e3) ** 2
        pr(f"sources: {s.n:,} ({s.nt} records); kept {int(kept.sum()):,} carrying "
           f"{q_kept:,.1f} m3/s mean; too far (> {m.max_snap_km:g} km) "
           f"{int((sn.status == TOO_FAR).sum()):,} carrying {q_far:,.1f}; outside the domain "
           f"{int((sn.status == OUTSIDE).sum()):,} carrying {q_out:,.1f}")
        pr(f"conservation: delivered / kept = {acct:.12f}  "
           f"({self.n_point:,} receiving cells, {dests.size:,} distinct destination cells)")
        if areas.size:
            pr(f"footprint area: median {np.median(areas) / 1e6:,.0f} km2, min "
               f"{areas.min() / 1e6:,.0f} (unclipped 4 pi sigma^2 = {ideal / 1e6:,.0f} km2); "
               f"{int((areas < 0.1 * ideal).sum())} below 10 %")
        far = np.flatnonzero(kept & (sn.snap_km > long_snap_km) & (mq > 0))
        if far.size:
            far = far[np.argsort(mq[far])[::-1]]
            pr(f"snapped > {long_snap_km:g} km: {mq[far].sum():,.1f} m3/s from {far.size} "
               f"sources; largest (check each lands in the right basin):")
            gl, gt = m.grid.wet_lon(), m.grid.wet_lat()
            for i in far[:max_list]:
                d = sn.dest[i]
                pr(f"   {s.name[i]:<28} {s.lon[i]:9.3f} {s.lat[i]:8.3f} {mq[i]:10,.1f} m3/s "
                   f"snap {sn.snap_km[i]:6.1f} km -> {gl[d]:9.3f} {gt[d]:8.3f}")
        for st, lab in ((TOO_FAR, "dropped as too far"), (OUTSIDE, "dropped as outside the domain")):
            sel = np.flatnonzero((sn.status == st) & (mq > 0) & (sn.snap_km <= 3 * m.max_snap_km))
            if sel.size:
                sel = sel[np.argsort(mq[sel])[::-1]]
                pr(f"{lab} (within {3 * m.max_snap_km:g} km of water): {mq[sel].sum():,.1f} m3/s; largest:")
                for i in sel[:max_list]:
                    pr(f"   {s.name[i]:<28} {s.lon[i]:9.3f} {s.lat[i]:8.3f} {mq[i]:10,.1f} m3/s "
                       f"{sn.snap_km[i]:6.1f} km from water")
        out = dict(n_source=s.n, n_kept=int(kept.sum()), q_kept=float(q_kept), q_too_far=float(q_far),
                   q_outside=float(q_out), accounted=acct, n_point=self.n_point,
                   footprint_median_km2=float(np.median(areas) / 1e6) if areas.size else float("nan"))
        if s.has_temp:
            out.update(self._report_temperature(pr, kept))
        return out

    def _report_temperature(self, pr, kept):
        """Heat conservation and the temperature range on the receiving cells."""
        s = self.sources
        ok = np.isfinite(s.temp) & kept[None, :]
        heat_in = float((np.where(ok, s.q, 0.0) * np.where(ok, s.temp, 0.0)).sum())
        q_unknown = float(np.where(~np.isfinite(s.temp) & kept[None, :], s.q, 0.0).sum()) / max(s.nt, 1)
        heat_out, tmin, tmax, n_unknown, n_wet = 0.0, np.inf, -np.inf, 0, 0
        for k in range(s.nt):
            t, r = self.temperature_points(k, return_runoff=True)
            good = np.isfinite(t)
            heat_out += float((r[good] * t[good] * self.point_area[good]).sum())
            # cells with water but no known temperature
            n_unknown += int(((r > 0) & ~good).sum())
            n_wet += int((r > 0).sum())
            if good.any():
                tmin, tmax = min(tmin, float(t[good].min())), max(tmax, float(t[good].max()))
        # heat_out uses all water reaching a cell (r) times the known-water temperature; equal
        # to heat_in when every source has a temperature (cells with partly-unknown water get
        # the temperature of their known part)
        ratio = heat_out / heat_in if heat_in != 0 else float("nan")
        pr(f"temperature: {tmin:.2f} .. {tmax:.2f} C on receiving cells; heat out/in = {ratio:.12f}; "
           f"cell-records with water of unknown temperature: {n_unknown:,} of {n_wet:,}"
           + (f" (unknown-temperature water {q_unknown:,.1f} m3/s mean: fill it, see "
              f"readers.temperature.fill_missing)" if q_unknown > 0 else ""))
        return dict(temp_min=tmin, temp_max=tmax, heat_ratio=ratio, cells_unknown_temp=n_unknown,
                    q_unknown_temp=q_unknown)

    def residence_days(self, k, depth_m=None):
        """Worst (smallest) h / R over the receiving cells for record k, days. h = column depth
        if the grid has it, else depth_m (e.g. the top-level thickness)."""
        r = self.runoff_points(k)
        dw = self.mapper.grid.wet_depth()
        h = dw[self.point_wet] if dw is not None else np.full(r.size, float(depth_m or 1.0))
        with np.errstate(divide="ignore"):
            tt = np.where(r > 0, h / r, np.inf)
        j = int(np.argmin(tt))
        return float(tt[j] / 86400.0), int(self.points[j])
