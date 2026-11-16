"""The target ocean-model grid.

A grid is four arrays of the SAME shape -- cell-centre longitude and latitude, cell area
(m2) and a wet mask -- of ANY rank:

    regional MITgcm          (ny, nx)
    ECCO LLC (tiled)         (tile, j, i)        e.g. (13, 90, 90) for LLC90
    MITgcm compact LLC       (ny, nx)            e.g. (1170, 90)

Nothing downstream uses grid topology (neighbours, faces, rotations): sources are matched
to cells by great-circle distance only, so every layout works the same way. The array
shape and dimension names are kept only so outputs can be written back in the grid's own
layout (dense MITgcm binaries) or indexed by it (the sparse netCDF).
"""
import glob
import os

import numpy as np

from .geometry import to_xyz, regular_cell_area

# Variable names seen in the wild (MNC output, xmitgcm, ECCO grid files, hand-made files).
NC_NAMES = {
    "lon": ["XC", "xC", "x_c", "lon", "longitude", "nav_lon", "lon_rho", "geolon", "LONC", "tlon"],
    "lat": ["YC", "yC", "y_c", "lat", "latitude", "nav_lat", "lat_rho", "geolat", "LATC", "tlat"],
    "area": ["rA", "RAC", "rac", "area", "cell_area", "area_t", "tarea", "AREA"],
    "mask": ["hFacC", "HFacC", "hfacc", "maskC", "mask_c", "maskInC", "mask", "Depth", "depth",
             "bathy", "bathymetry", "H", "wet"],
    "depth": ["Depth", "depth", "bathy", "bathymetry", "H", "R_low"],
}
_DEPTH_LIKE = ("depth", "bathy", "bathymetry", "h", "r_low")


def _pick(names, kind, override):
    if override:
        if override not in names:
            raise ValueError(f"{kind} variable {override!r} not in file; present: {sorted(names)}")
        return override
    for n in NC_NAMES[kind]:
        if n in names:
            return n
    raise ValueError(f"no {kind} variable found (tried {NC_NAMES[kind]}); present: "
                     f"{sorted(names)}. Pass it explicitly.")


class ModelGrid:
    """Cell centres, areas and wet mask of a model grid, any shape.

    Attributes
    ----------
    lon, lat, area : float arrays, shape S
    wet            : bool array, shape S (surface level)
    depth          : optional float array, shape S (column depth, m, positive)
    dims           : names of the S dimensions (used by the sparse netCDF writer)
    wet_flat       : flat (C-order) indices of the wet cells; "wet index" k refers to
                     wet_flat[k] everywhere in this package
    """

    def __init__(self, lon, lat, area, wet, dims=None, depth=None, name=""):
        lon = np.asarray(lon, dtype=float)
        lat = np.asarray(lat, dtype=float)
        area = np.asarray(area, dtype=float)
        wet = np.asarray(wet, dtype=bool)
        if not (lon.shape == lat.shape == area.shape == wet.shape):
            raise ValueError(f"shapes disagree: lon{lon.shape} lat{lat.shape} "
                             f"area{area.shape} wet{wet.shape}")
        if not (np.isfinite(lon).all() and np.isfinite(lat).all()):
            raise ValueError("lon/lat contain non-finite values (masked fill values?)")
        if not wet.any():
            raise ValueError("no wet cells")
        if float(area[wet].min()) <= 0.0:
            raise ValueError("area has non-positive values at wet cells")
        if depth is not None:
            depth = np.asarray(depth, dtype=float)
            if depth.shape != wet.shape:
                depth = None
        self.lon, self.lat, self.area, self.wet, self.depth = lon, lat, area, wet, depth
        self.shape = wet.shape
        self.dims = tuple(dims) if dims is not None else tuple(f"dim{i}" for i in range(wet.ndim))
        if len(self.dims) != wet.ndim:
            raise ValueError(f"{len(self.dims)} dim names for a {wet.ndim}-D grid")
        self.name = name
        self.wet_flat = np.flatnonzero(wet)

    # ------------------------------------------------------------------ properties
    @property
    def n_wet(self):
        return int(self.wet_flat.size)

    @property
    def size(self):
        return int(self.wet.size)

    def wet_lon(self):
        return self.lon.ravel()[self.wet_flat]

    def wet_lat(self):
        return self.lat.ravel()[self.wet_flat]

    def wet_area(self):
        return self.area.ravel()[self.wet_flat]

    def wet_depth(self):
        return None if self.depth is None else self.depth.ravel()[self.wet_flat]

    def wet_xyz(self):
        return to_xyz(self.wet_lon(), self.wet_lat())

    def mean_cell_km(self):
        """sqrt of the mean wet-cell area, km."""
        return float(np.sqrt(self.wet_area().mean())) / 1e3

    def unravel(self, flat):
        """Flat cell indices -> tuple of index arrays, one per dimension."""
        return np.unravel_index(np.asarray(flat, dtype=np.int64), self.shape)

    def to_dense(self, values_wet, fill=0.0, dtype=float):
        """Values on the wet cells (length n_wet) -> full array of the grid's shape."""
        out = np.full(self.size, fill, dtype=dtype)
        out[self.wet_flat] = values_wet
        return out.reshape(self.shape)

    def summary(self):
        return (f"grid {self.name or ''} shape {self.shape} dims {self.dims}: "
                f"{self.n_wet:,} wet cells, mean cell {self.mean_cell_km():.2f} km, "
                f"lon {self.lon.min():.2f}..{self.lon.max():.2f}, "
                f"lat {self.lat.min():.2f}..{self.lat.max():.2f}")

    # ------------------------------------------------------------------ readers
    @classmethod
    def from_netcdf(cls, path, lon=None, lat=None, area=None, mask=None, level=0,
                    depth=None, name=None, verbose=True):
        """Read a grid from one netCDF file.

        lon/lat/area/mask : variable names (default: first match in NC_NAMES).
        level             : vertical index of the mask when it has one more dimension than
                            lon (hFacC(k,j,i), ECCO maskC(k,tile,j,i)); 0 = surface.
        Works for 2-D (j,i) grids, tiled LLC files (tile,j,i), and regular grids given as
        1-D lon/lat vectors (meshed; area computed if the file has none).
        A mask named like a depth/bathymetry is turned into wet = depth > 0 (or < 0 if the
        file uses negative-down values).
        """
        from netCDF4 import Dataset
        if not os.path.exists(path):
            raise FileNotFoundError(path)
        with Dataset(path) as ds:
            v = ds.variables
            nlon = _pick(v, "lon", lon)
            nlat = _pick(v, "lat", lat)
            nmask = _pick(v, "mask", mask)
            LON = np.asarray(v[nlon][:], dtype=float)
            LAT = np.asarray(v[nlat][:], dtype=float)
            dims = tuple(v[nlon].dimensions)
            regular = LON.ndim == 1 and LAT.ndim == 1 and v[nlon].dimensions != v[nlat].dimensions
            if regular:
                lon1, lat1 = LON, LAT
                LON, LAT = np.meshgrid(lon1, lat1)
                dims = (v[nlat].dimensions[0], v[nlon].dimensions[0])
            try:
                narea = _pick(v, "area", area)
                AREA = np.asarray(v[narea][:], dtype=float)
            except ValueError:
                if not regular:
                    raise
                narea = "(computed: regular lat-lon cells)"
                AREA = regular_cell_area(lon1, lat1)
            M = np.ma.filled(v[nmask][:], 0.0).astype(float)
            while M.ndim > LON.ndim + 1:          # leading time-like dims: first record
                M = M[0]
            took_level = M.ndim == LON.ndim + 1
            if took_level:                        # vertical dim: take the requested level
                M = M[level]
            if nmask.lower() in _DEPTH_LIKE:
                neg = float(np.nanmin(M)) < 0 and float(np.nanmax(M)) <= 0
                WET = (M < 0) if neg else (M > 0)
            else:
                WET = M > 0
            DEP = None
            try:
                nd = _pick(v, "depth", depth)
                D = np.asarray(np.ma.filled(v[nd][:], 0.0), dtype=float)
                while D.ndim > LON.ndim:
                    D = D[0]
                DEP = np.abs(D)
            except ValueError:
                pass
        if AREA.shape != LON.shape:
            raise ValueError(f"area shape {AREA.shape} != lon shape {LON.shape}")
        if verbose:
            print(f"grid {path}: lon={nlon} lat={nlat} area={narea} mask={nmask}"
                  f"{f' (level {level})' if took_level else ''} -> shape {LON.shape} dims {dims}")
        return cls(LON, LAT, AREA, WET, dims=dims, depth=DEP,
                   name=name or os.path.basename(path))

    @classmethod
    def from_mds(cls, grid_dir, name=None):
        """Read XC, YC, RAC, hFacC (and DRF) from an MITgcm grid dump (needs MITgcmutils)."""
        from MITgcmutils import mds
        if not glob.glob(os.path.join(grid_dir, "XC.*")):
            raise FileNotFoundError(f"no XC.* in {grid_dir}")
        rd = lambda n: np.asarray(mds.rdmds(os.path.join(grid_dir, n)), dtype=float)  # noqa
        XC, YC, RAC, hf = rd("XC"), rd("YC"), rd("RAC"), rd("hFacC")
        wet = (hf[0] > 0) if hf.ndim == XC.ndim + 1 else (hf > 0)
        dep = None
        try:
            drf = rd("DRF").ravel()
            if hf.ndim == XC.ndim + 1 and drf.size >= hf.shape[0]:
                dep = np.tensordot(drf[:hf.shape[0]], hf, axes=(0, 0))
        except Exception:
            pass
        dims = ("j", "i") if XC.ndim == 2 else tuple(f"dim{i}" for i in range(XC.ndim))
        return cls(XC, YC, RAC, wet, dims=dims, depth=dep, name=name or os.path.basename(grid_dir))
