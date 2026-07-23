"""
File handling module (translated from file.c)

Output is written in NetCDF by default (config.OUTPUT_FORMAT / CGEM_OUTPUT), with the
legacy tab-separated .dat path preserved for CGEM_OUTPUT=dat|both. hydwrite/transwrite/
Rates keep their original signatures and dispatch to whichever format(s) are enabled;
the analysis loaders in tools/make_*_pdf.py read whichever exists.
"""

from config import M, DELTI, repeatYear, nbday_ice, OUTPUT_FORMAT
from variables import DEPTH, B, D, U, ice_thickness, ice_frac
import numpy as np
from numpy import genfromtxt, interp, linspace

_WANT_DAT = OUTPUT_FORMAT in ("dat", "both")
_WANT_NC = OUTPUT_FORMAT in ("nc", "both")


class _NcWriter:
    """Append time slices (time x distance) to one NetCDF file per site. Variables are
    created lazily on first write, on a shared x-dimension of size M+1 (cell index).
    Species/rate fields occupy cells 1..M (cell 0 stays fill); depth/width use 0..M."""

    def __init__(self, path="output.nc"):
        """Set up the writer for one NetCDF file `path` (dataset created lazily on first write)."""
        self._path = path
        self._ds = None
        self._vars = {}
        self._t_index = {}

    def _ensure(self):
        """Create the NetCDF dataset and its time/x dimensions on first write (lazy)."""
        if self._ds is None:
            from netCDF4 import Dataset
            self._ds = Dataset(self._path, "w", format="NETCDF4")
            self._ds.createDimension("time", None)     # unlimited
            self._ds.createDimension("x", M + 1)
            tv = self._ds.createVariable("time", "f8", ("time",))
            tv.units = "s"
            tv.long_name = "seconds since model start"
            xv = self._ds.createVariable("x", "i4", ("x",))
            xv.long_name = "grid cell index (staggered; odd=concentration, even=velocity)"
            xv[:] = np.arange(M + 1)
            self._ds.description = ("NS-RAD output. Fields are time x cell; "
                                    "species/rates on cells 1..M, depth/width on 0..M.")

    def _tidx(self, t):
        """Return (creating if needed) the unlimited-time index for model time `t`, writing t into the time coordinate."""
        if t not in self._t_index:
            i = len(self._t_index)
            self._t_index[t] = i
            self._ds.variables["time"][i] = t
        return self._t_index[t]

    def write(self, name, arr, t, lo):
        """Write arr[lo:M+1] into variable `name` at the time index for t."""
        self._ensure()
        name = name[:-4] if name.endswith(".dat") else name
        if name not in self._vars:
            self._vars[name] = self._ds.createVariable(
                name, "f8", ("time", "x"), zlib=True, complevel=4,
                fill_value=np.nan)
        idx = self._tidx(t)
        self._vars[name][idx, lo:M + 1] = np.asarray(arr[lo:M + 1], dtype="f8")

    def close(self):
        """Flush and close the NetCDF file (idempotent)."""
        if self._ds is not None:
            self._ds.close()
            self._ds = None


_NC = _NcWriter() if _WANT_NC else None


def close_output():
    """Flush and close the NetCDF file. Call once at the end of the run."""
    if _NC is not None:
        _NC.close()


def hydwrite(t):
    """Write physical variables (depth, width, velocity)."""
    if _WANT_DAT:
        with open("depth.dat", "a") as fptr1, open("width.dat", "a") as fptr2, \
                open("velocity.dat", "a") as fptr3:
            fptr1.write(f"{t}\t")
            fptr2.write(f"{t}\t")
            fptr3.write(f"{t}\t")
            for i in range(M + 1):
                fptr1.write(f"{DEPTH[i]:f}\t")
                fptr2.write(f"{B[i]:f}\t")
                fptr3.write(f"{U[i]:f}\t")
            fptr1.write("\n")
            fptr2.write("\n")
            fptr3.write("\n")
    if _WANT_NC:
        _NC.write("depth", DEPTH, t, 0)
        _NC.write("width", B, t, 0)
        # velocity lives on the even (velocity) nodes of the staggered grid
        _NC.write("velocity", U, t, 0)


def icewrite(t):
    """Write the prognostic ice cover (thickness [m] and areal fraction [0-1])."""
    if _WANT_DAT:
        with open("ice_thickness.dat", "a") as f1, open("ice_frac.dat", "a") as f2:
            f1.write(f"{t}\t")
            f2.write(f"{t}\t")
            for i in range(1, M + 1):
                f1.write(f"{ice_thickness[i]:.6e}\t")
                f2.write(f"{ice_frac[i]:.6e}\t")
            f1.write("\n")
            f2.write("\n")
    if _WANT_NC:
        _NC.write("ice_thickness", ice_thickness, t, 1)
        _NC.write("ice_frac", ice_frac, t, 1)


def transwrite(co, s, t):
    """Write concentration values."""
    if _WANT_DAT:
        with open(s, "a") as fptr1:
            fptr1.write(f"{t}\t")
            for i in range(1, M + 1):
                fptr1.write(f"{co[i]:.15e}\t")
            fptr1.write("\n")
    if _WANT_NC:
        _NC.write(s, co, t, 1)


def Rates(co, s, t):
    """Write biogeochemical reaction rates."""
    if _WANT_DAT:
        with open(s, "a") as fptr1:
            fptr1.write(f"{t}\t")
            for i in range(1, M + 1):
                fptr1.write(f"{co[i]:.15e}\t")
            fptr1.write("\n")
    if _WANT_NC:
        _NC.write(s, co, t, 1)

# Parsed forcing series, keyed by filename. The originals were re-opened and
# re-parsed on EVERY call -- six calls per timestep across ~840k timesteps, i.e.
# millions of full CSV parses per run, all returning identical arrays. Caching is
# purely a speedup: the interpolation below still happens per timestep, so results
# are bit-for-bit unchanged.
_FORCING_CACHE = {}

# Fixed 365-point time axis the forcings are interpolated onto. A series of any
# other length raises ValueError in interp; see CLAUDE.md.
_TIME_AXIS = linspace(0, 365*86400, 365)


def _load(name):
    """Return (series, rolling_ice_sum) for a forcing file, parsing at most once."""
    cached = _FORCING_CACHE.get(name)
    if cached is None:
        with open(name, 'r', encoding='utf-8-sig') as f:
            data = genfromtxt(f, delimiter=',', dtype=float)
        # Rolling nbday_ice-day cumulative sum. Only meaningful for water
        # temperature, where main.py uses it to drive the ice gate.
        b = data.cumsum()
        b[nbday_ice:] = b[nbday_ice:] - b[:-nbday_ice]
        cached = (data, b)
        _FORCING_CACHE[name] = cached
    return cached


def exfread(name, t):
    """Interpolate a cached daily forcing series to model time `t` [s]. Returns (value,
    rolling_ice_sum); the second is meaningful only for the water-temperature call (the
    legacy previousdays gate). Series are parsed once and memoised in _FORCING_CACHE."""
    data, b = _load(name)
    if repeatYear == 1 and t > 31536000:
        t = t - 31536000
    y = interp(t, _TIME_AXIS, data)
    # for water temperature only
    y2 = interp(t, _TIME_AXIS, b)
    return y, y2