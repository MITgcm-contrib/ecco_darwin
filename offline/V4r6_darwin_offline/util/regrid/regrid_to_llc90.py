"""Regrid Steph's 1x1-degree offline IC/forcing binaries onto ECCOv4 LLC90.

Built and self-tested (with synthetic data, against the real local LLC90 grid
geometry) BEFORE Steph's actual files were available -- see
regrid_to_llc90_selftest() at the bottom. Swap in real source files and run
regrid_field() for real once they arrive; nothing else needs to change.

Deliberately built on the ECCO group's own reference tooling
(ecco_v4_py.llc_array_conversion, MITgcmutils.rdmds) rather than reimplementing
the LLC90 "compact" binary tile layout -- that layout is exactly the kind of
detail that silently produces plausible-but-wrong output if guessed at.
Round-trip-verified lossless (llc_tiles_to_compact(llc_compact_to_tiles(x)) ==
x to bit precision) against the real grid files in bathy/grid/ECCO_V4r5/.

Grid geometry note: bathy/grid/ECCO_V4r5/ is confirmed correct for ECCOv4r6
too -- its RF.data cumulative depths match r6's own "namelist/data" delR array
exactly, and r6's PARM05 bathyFile name matches the file present there
byte-for-byte (see ECCO/CLAUDE.md and the plan file for the full provenance
trail).
"""
import os

import numpy as np
from scipy.interpolate import RegularGridInterpolator
import MITgcmutils as mu
import ecco_v4_py as ecco

# Grid directories are machine-specific: override via the environment rather
# than editing this file.  LLC90_GRID_DIR must hold the ECCOv4 LLC90 grid
# (XC/YC/RF/hFacC .data+.meta, as published with the ECCOv4r5/r6 ancillary
# data); LLC270_GRID_DIR the LLC270 grid, needed only by
# regrid_llc270_to_llc90().
GRID_DIR = os.environ.get("LLC90_GRID_DIR", "./grid/ECCO_V4r5")

# Steph's 1x1-degree grid, from Steph's original 1-deg set-up (data) PARM04:
#   xgOrigin=0, delX=360*1.E0 ; ygOrigin=-80., delY=160*1.E0
SRC_NX, SRC_NY = 360, 160
SRC_LON0, SRC_DLON = 0.0, 1.0
SRC_LAT0, SRC_DLAT = -80.0, 1.0

# Steph's 23-level vertical grid (Steph's original 1-deg set-up (data) PARM04
# delZ), converted to cell-center depths (positive down) for vertical
# interpolation of 3D tracer fields.
_SRC_DELZ = np.array([
    10., 10., 15., 20., 20., 25., 35., 50., 75.,
    100., 150., 200., 275., 350., 415., 450., 500.,
    500., 500., 500., 500., 500., 500.,
])
SRC_RC = np.cumsum(_SRC_DELZ) - _SRC_DELZ / 2.0  # 23 values

assert _SRC_DELZ.size == 23

# ECCOv4r6's real 50-level delR, from "ECCOv4 Release 6"/namelist/data
# PARM04 (fetched) -- same array already used in
# input_darwin_offline/data.
TARGET_DELR = np.array([
    10.00, 10.00, 10.00, 10.00, 10.00, 10.00, 10.00, 10.01,
    10.03, 10.11, 10.32, 10.80, 11.76, 13.42, 16.04, 19.82, 24.85,
    31.10, 38.42, 46.50, 55.00, 63.50, 71.58, 78.90, 85.15, 90.18,
    93.96, 96.58, 98.25, 99.25, 100.01, 101.33, 104.56, 111.33, 122.83,
    139.09, 158.94, 180.83, 203.55, 226.50, 249.50, 272.50, 295.50, 318.50,
    341.50, 364.50, 387.50, 410.50, 433.50, 456.50,
])
assert TARGET_DELR.size == 50


def target_rc():
    """ECCOv4r6's 50-level cell-center depths (positive down)."""
    return np.cumsum(TARGET_DELR) - TARGET_DELR / 2.0


def _interp_vertical(horiz, trc):
    """Vectorized 1D linear interpolation of horiz (23,13,90,90) from
    SRC_RC onto the 50 target depths trc, applied independently at every
    horizontal point in a single broadcasted pass (no Python loop over the
    105300 horizontal points).

    Extrapolates with the nearest source level's value (constant), not a
    linear extrapolation, since the target grid's bottom (6134.5 m) is
    deeper than the source's deepest cell center and a linear extrapolation
    could produce unphysical (e.g. negative concentration) values at depth.
    """
    idx = np.searchsorted(SRC_RC, trc)              # (50,), in [0,23]
    idx_lo = np.clip(idx - 1, 0, SRC_RC.size - 1)    # (50,)
    idx_hi = np.clip(idx, 0, SRC_RC.size - 1)        # (50,)

    horiz_flat = horiz.reshape(horiz.shape[0], -1)   # (23, N)
    lo_vals = horiz_flat[idx_lo]                     # (50, N)
    hi_vals = horiz_flat[idx_hi]                     # (50, N)

    rc_lo, rc_hi = SRC_RC[idx_lo], SRC_RC[idx_hi]    # (50,)
    span = rc_hi - rc_lo
    # idx_lo == idx_hi happens exactly where trc is outside SRC_RC's range
    # (constant extrapolation) -- span is 0 there, so guard the divide and
    # force weight 0 (picks lo_vals, which equals hi_vals at that index).
    w_hi = np.where(span > 0, (trc - rc_lo) / np.where(span > 0, span, 1.0), 0.0)

    out_flat = lo_vals + w_hi[:, None] * (hi_vals - lo_vals)
    return out_flat.reshape((trc.size,) + horiz.shape[1:])


def src_lonlat():
    """Steph's 1x1-degree cell-center longitude/latitude axes."""
    lon = SRC_LON0 + SRC_DLON * (np.arange(SRC_NX) + 0.5)
    lat = SRC_LAT0 + SRC_DLAT * (np.arange(SRC_NY) + 0.5)
    return lon, lat


def target_grid(grid_dir=GRID_DIR):
    """LLC90 target coordinates and wet mask, tile-indexed (13,90,90)."""
    xc = ecco.llc_compact_to_tiles(mu.rdmds(grid_dir + "/XC"), less_output=True)
    yc = ecco.llc_compact_to_tiles(mu.rdmds(grid_dir + "/YC"), less_output=True)
    hfacc = ecco.llc_compact_to_tiles(mu.rdmds(grid_dir + "/hFacC"), less_output=True)
    return xc, yc, hfacc  # xc/yc: (13,90,90); hfacc: (50,13,90,90)


def _fill_land_nearest_wet_2d(tile2d, wet2d):
    """Fill dry (wet2d==False) points in ONE (90,90) tile with the nearest
    wet point's value, within that tile only.

    Avoids writing raw interpolated-through-land values (which can be
    unphysical, e.g. extrapolated from a distant open-ocean point across a
    strait) and avoids a bare zero (which creates an artificial gradient at
    the coastline that some non-hFac-weighted operators are not shielded
    from). Nearest-wet-neighbor fill is what real ECCO-Darwin IC pipelines
    do for exactly this reason. Must be per-tile: LLC90's 13 tiles are not
    laid out contiguously in physical space along the tile axis, so a
    nearest-neighbor search across tiles would pull values from a
    geometrically unrelated tile.
    """
    if wet2d.all():
        return tile2d
    if not wet2d.any():
        return tile2d  # nothing to fill from; leave as-is (e.g. all-land tile)
    from scipy.ndimage import distance_transform_edt
    _, (iy, ix) = distance_transform_edt(~wet2d, return_indices=True)
    return tile2d[iy, ix]


def _fill_land_nearest_wet(tiles, wet):
    """Apply _fill_land_nearest_wet_2d independently to each of the 13
    (90,90) tiles in a (13,90,90) array."""
    out = tiles.copy()
    for t in range(tiles.shape[0]):
        out[t] = _fill_land_nearest_wet_2d(tiles[t], wet[t])
    return out


def _fill_source_land_nearest_2d(src_field):
    """Fill Steph's source-grid land points (exact 0.0 -- confirmed the
    shared flag across every IC/forcing field checked: DIC, alkalinity,
    nitrate, phosphate, silicate, oxygen, the DARWIN_ECCO_run33_4 restart
    tracers, biomass, iron and wind all show the identical 36.1% masked
    fraction at the surface) with the nearest wet source-grid value, before
    horizontal interpolation.

    Without this, RegularGridInterpolator blends real ocean values with
    literal 0.0 land neighbors for every LLC90 point that lands near a
    source-grid coastline, biasing coastal concentrations low -- a
    plausible-but-wrong number of exactly the kind this pipeline is trying
    to avoid at the target-grid side already (_fill_land_nearest_wet_2d).
    This is the source-side counterpart of that same problem.

    Safe to do a whole-grid (not per-tile) nearest-neighbor search here,
    unlike the LLC90 target side: Steph's source is a single contiguous
    360x160 lat-lon grid, not disjoint cube-sphere tiles.
    """
    wet = src_field != 0.0
    if wet.all() or not wet.any():
        return src_field
    from scipy.ndimage import distance_transform_edt
    # Pad periodically in longitude so the nearest-wet search doesn't treat
    # the 0/360 seam as a hard edge.
    pad = 10
    field_ext = np.concatenate(
        [src_field[:, -pad:], src_field, src_field[:, :pad]], axis=1
    )
    wet_ext = np.concatenate([wet[:, -pad:], wet, wet[:, :pad]], axis=1)
    _, (iy, ix) = distance_transform_edt(~wet_ext, return_indices=True)
    filled_ext = field_ext[iy, ix]
    return filled_ext[:, pad:-pad]


def regrid_horizontal_2d(src_field, xc, yc, wet2d):
    """Regrid one 2D (lat,lon) field from the 1-degree source grid onto one
    LLC90 target level's (13,90,90) coordinates, with land-fill on both
    ends of the interpolation.

    src_field: (SRC_NY, SRC_NX)
    xc, yc, wet2d: (13,90,90)
    returns: (13,90,90)
    """
    src_filled = _fill_source_land_nearest_2d(src_field)

    lon, lat = src_lonlat()
    # Wrap source periodically in longitude so interpolation near the 0/360
    # seam doesn't fall outside the source grid's domain.
    lon_ext = np.concatenate(([lon[-1] - 360.0], lon, [lon[0] + 360.0]))
    src_ext = np.concatenate(
        [src_filled[:, -1:], src_filled, src_filled[:, :1]], axis=1
    )
    interp = RegularGridInterpolator(
        (lat, lon_ext), src_ext, bounds_error=False, fill_value=None
    )

    lon_q = np.mod(xc, 360.0)
    # Clamp query latitude into the source's actual domain (Steph's grid
    # only covers -80..80) before interpolating. fill_value=None makes
    # scipy do LINEAR EXTRAPOLATION past the grid edge, not nearest-edge as
    # a first version of this function assumed -- confirmed empirically:
    # regridding real nitrate/phosphate/silicate/oxygen fields produced
    # unphysical negative concentrations (nitrate as low as -152 mmol/m3)
    # exactly at LLC90's Arctic-cap points poleward of 79.5N, the one
    # region with no source-grid analog to interpolate between. Clamping
    # converts that extrapolation into the intended nearest-edge behavior.
    lat_q = np.clip(yc, lat.min(), lat.max())
    out = interp((lat_q, lon_q))
    return _fill_land_nearest_wet(out.reshape(13, 90, 90), wet2d.astype(bool))


def regrid_field(src_field, is_3d, grid_dir=GRID_DIR):
    """Regrid a full Steph 1-degree field onto LLC90, compact-ordered, ready
    to write to a pkg/offline / PTRACERS_initialFile binary.

    src_field: (SRC_NY, SRC_NX) for 2D, or (23, SRC_NY, SRC_NX) for 3D.
    returns: (1170,90) for 2D, or (50,1170,90) for 3D -- the compact layout
             MITgcm's mds reader (READ_REC_3D_RS / READ_REC_XY_RS) expects.
    """
    xc, yc, hfacc = target_grid(grid_dir)
    if not is_3d:
        wet2d = hfacc[0] > 0
        tiles = regrid_horizontal_2d(src_field, xc, yc, wet2d)
        return ecco.llc_tiles_to_compact(tiles, less_output=True)

    assert src_field.shape == (23, SRC_NY, SRC_NX), src_field.shape
    # Horizontal regrid at each of the 23 source levels first (cheap: 23
    # 2D interpolations), then vertical interpolation onto the 50 target
    # levels -- much cheaper than the reverse order, and the two steps
    # commute for a separable lat/lon-then-depth regridding like this.
    horiz = np.stack(
        [
            regrid_horizontal_2d(src_field[k], xc, yc, hfacc[0] > 0)
            for k in range(23)
        ],
        axis=0,
    )  # (23,13,90,90)

    out = _interp_vertical(horiz, target_rc())
    for k in range(50):
        out[k] = _fill_land_nearest_wet(out[k], hfacc[k] > 0)
    return ecco.llc_tiles_to_compact(out, less_output=True)


def regrid_monthly_2d(src_stack, grid_dir=GRID_DIR):
    """Regrid a 12-month climatology of 2D (lat,lon) fields (Steph's
    ironFile/windFile convention: one surface value per calendar month,
    read+time-interpolated at runtime by darwin_exf_load.F -- no vertical
    dimension, so this is 12 independent horizontal-only regrids, not a
    3D regrid_field call).

    src_stack: (12, SRC_NY, SRC_NX)
    returns: (12, 1170, 90) compact-ordered, ready for write_compact_binary.
    """
    assert src_stack.shape == (12, SRC_NY, SRC_NX), src_stack.shape
    xc, yc, hfacc = target_grid(grid_dir)
    wet2d = hfacc[0] > 0
    out = np.stack(
        [
            ecco.llc_tiles_to_compact(
                regrid_horizontal_2d(src_stack[m], xc, yc, wet2d), less_output=True
            )
            for m in range(12)
        ],
        axis=0,
    )
    return out  # (12, 1170, 90)


# ---------------------------------------------------------------------------
# LLC270 -> LLC90 (native-grid-to-native-grid, NOT lat-lon-to-LLC90). Added
# for the real LLC270-compact ALK/DIC IC files under steph_input/v06/ICs/ -- a genuinely different regridding problem from
# everything above: both grids are curvilinear cube-sphere tiles, so there is
# no rectangular source grid for RegularGridInterpolator to interpolate from.
#
# Method: per-level nearest-WET-neighbor, matching this group's own
# established LLC-to-LLC regridding convention (see
# carbon/m_files/grid/interpolate_CS510_to_LLC270/nninterp.m and
# carbon/m_files/grid/interpolate_pickup_manfredi/convert_input_cs510_llc270.m
# -- same nearest-neighbor-restricted-to-wet-source-cells idea, done per
# depth level since the wet/dry pattern changes with depth). Reimplemented
# here with a KDTree over each point's Cartesian unit-sphere embedding
# instead of that script's degree-space distance + bounding-box search +
# dateline-unwrap hack -- the Cartesian embedding handles the dateline and
# poles correctly with no special-casing, and is far faster (that MATLAB
# loop is O(N*M) per level; a KDTree query is O(M log N)).
# ---------------------------------------------------------------------------

LLC270_GRID_DIR = os.environ.get("LLC270_GRID_DIR", "./grid/LLC_270")


def _latlon_to_unit_sphere_xyz(lon, lat):
    """(..., ) lon/lat in degrees -> (..., 3) Cartesian unit-sphere coords.
    Distance in this embedding is a monotonic function of great-circle
    distance, so nearest-neighbor here matches nearest-neighbor on the
    sphere -- and unlike degree-space distance, needs no dateline unwrap or
    pole special-casing."""
    lon_r = np.radians(lon)
    lat_r = np.radians(lat)
    x = np.cos(lat_r) * np.cos(lon_r)
    y = np.cos(lat_r) * np.sin(lon_r)
    z = np.sin(lat_r)
    return np.stack([x, y, z], axis=-1)


def llc270_grid(grid_dir=LLC270_GRID_DIR):
    """LLC270 target coordinates and wet mask, tile-indexed (13,270,270) /
    (50,13,270,270) -- same convention as target_grid(), just at N=270.
    ecco_v4_py.llc_compact_to_tiles infers the tile size from the array
    shape, so this needs no changes from target_grid() beyond the path and
    grid size (confirmed: (3510,270) reshapes to (13,270,270) cleanly)."""
    xc = ecco.llc_compact_to_tiles(mu.rdmds(grid_dir + "/XC"), less_output=True)
    yc = ecco.llc_compact_to_tiles(mu.rdmds(grid_dir + "/YC"), less_output=True)
    hfacc = ecco.llc_compact_to_tiles(mu.rdmds(grid_dir + "/hFacC"), less_output=True)
    return xc, yc, hfacc


def regrid_llc270_to_llc90(field_270, grid_dir_270=LLC270_GRID_DIR, grid_dir_90=GRID_DIR):
    """Nearest-wet-neighbor regrid of a 50-level LLC270-compact field onto
    LLC90-compact.

    field_270: (50, 3510, 270) compact, e.g. np.fromfile(path, dtype='>f4')
               .reshape(50, 3510, 270) -- matches the real
               ALK.R4_270x3510x50.bin/DIC.R4_270x3510x50.bin layout (verified
               against the real LLC270 hFacC.data, which is the identical
               shape/byte count).
    returns: (50, 1170, 90) compact, ready for write_compact_binary.

    Land handling: unlike the lat-lon pipeline (which had to infer land from
    an empirical "exact 0.0" convention, since no real land mask was
    available for Steph's 1-degree source), LLC270's real hFacC is used
    directly to restrict candidate neighbors to wet source cells -- no
    source-side land-fill heuristic needed here. Target-side dry cells are
    NOT explicitly re-masked afterward: MITgcm's own ptracers_init_varia.F
    unconditionally zeros every tracer where maskC==0 at cold start
    regardless of what the initial file contains, so whatever a dry LLC90
    cell gets written here is irrelevant -- only wet-cell values matter, and
    for those the nearest-wet-neighbor search already keeps a wet LLC90 cell
    from ever reading a land value.
    """
    xc270, yc270, hfacc270 = llc270_grid(grid_dir_270)
    xc90, yc90, _ = target_grid(grid_dir_90)

    field_tiles = ecco.llc_compact_to_tiles(field_270, less_output=True)  # (50,13,270,270)

    xyz270 = _latlon_to_unit_sphere_xyz(xc270, yc270)  # (13,270,270,3)
    xyz90 = _latlon_to_unit_sphere_xyz(xc90, yc90)      # (13,90,90,3)
    xyz90_flat = xyz90.reshape(-1, 3)

    from scipy.spatial import cKDTree

    nz = field_tiles.shape[0]
    out = np.empty((nz, 13, 90, 90), dtype=np.float64)
    for k in range(nz):
        wet270 = hfacc270[k] > 0
        if not wet270.any():
            out[k] = 0.0
            continue
        tree = cKDTree(xyz270[wet270])
        _, idx = tree.query(xyz90_flat)
        out[k] = field_tiles[k][wet270][idx].reshape(13, 90, 90)

    return ecco.llc_tiles_to_compact(out, less_output=True)


def write_compact_binary(compact_field, out_path):
    """Write a compact-ordered field as the big-endian float32 MITgcm expects
    (readBinaryPrec=32 in every namelist here)."""
    compact_field.astype(">f4").tofile(out_path)


def _test_interp_vertical_matches_np_interp():
    """_interp_vertical is a vectorized replacement for a per-column
    np.interp loop -- check it's numerically identical, at a handful of
    random columns, to the direct per-column call it replaced."""
    rng = np.random.default_rng(0)
    horiz = rng.normal(size=(23, 13, 90, 90))
    trc = target_rc()
    out = _interp_vertical(horiz, trc)
    for _ in range(20):
        t, y, x = rng.integers(13), rng.integers(90), rng.integers(90)
        expected = np.interp(trc, SRC_RC, horiz[:, t, y, x])
        np.testing.assert_allclose(out[:, t, y, x], expected)
    print("_interp_vertical matches np.interp at 20 random columns: OK")


def _test_source_land_fill_removes_zeros_before_interp():
    """Inject a synthetic 'coastline' (a block of exact 0.0, Steph's real
    land flag) next to a smooth nonzero ocean field, and check that
    regrid_horizontal_2d's output near that block is close to the smooth
    field's true value -- not pulled toward 0 by blending with the raw
    land zeros. Also check the periodic pad doesn't leave a seam artifact.
    """
    lon, lat = src_lonlat()
    LON, LAT = np.meshgrid(lon, lat)
    smooth = 20.0 + 1e-6 * (LON + LAT)  # ~constant, tiny slope to break ties
    with_land = smooth.copy()
    with_land[60:70, 0:10] = 0.0  # a synthetic coastal block, incl. lon seam
    with_land[60:70, -10:] = 0.0

    filled = _fill_source_land_nearest_2d(with_land)
    assert not (filled == 0.0).any(), "land zeros leaked through source fill"
    np.testing.assert_allclose(filled, smooth, atol=1e-3)
    print("source land-fill OK: coastal block recovered, no 0.0 leakage")

    xc, yc, hfacc = target_grid()
    wet2d = hfacc[0] > 0
    out_unfilled = _fill_land_nearest_wet(
        RegularGridInterpolator(
            (lat, np.concatenate(([lon[-1] - 360.0], lon, [lon[0] + 360.0]))),
            np.concatenate([with_land[:, -1:], with_land, with_land[:, :1]], axis=1),
            bounds_error=False, fill_value=None,
        )((yc, np.mod(xc, 360.0))).reshape(13, 90, 90),
        wet2d.astype(bool),
    )
    out_filled = regrid_horizontal_2d(with_land, xc, yc, wet2d)
    # Wherever the unfilled path dipped noticeably below the smooth value
    # (i.e. got contaminated by the land zeros), the fixed path must be
    # closer to truth.
    contaminated = out_unfilled < smooth.mean() - 1.0
    assert contaminated.any(), "test didn't actually exercise contamination"
    assert np.abs(out_filled[contaminated] - smooth.mean()).max() < \
        np.abs(out_unfilled[contaminated] - smooth.mean()).max()
    print("regrid_horizontal_2d source-fill demonstrably reduces coastal bias")


def _test_polar_points_use_nearest_edge_not_linear_extrapolation():
    """Regression test for the real bug found regridding Steph's actual
    nitrate/phosphate/silicate/oxygen fields: LLC90's Arctic-cap points sit
    poleward of Steph's source domain (source stops at 79.5N; LLC90 goes to
    ~89.7N), and RegularGridInterpolator's fill_value=None does LINEAR
    EXTRAPOLATION there, not nearest-edge -- producing wild overshoot for
    any field with a nonzero gradient near the domain edge. Build a field
    with a strong, non-flat gradient right at the north edge and check the
    output at real Arctic-cap LLC90 points stays within the source field's
    own min/max (nearest-edge behavior), rather than overshooting it.
    """
    lon, lat = src_lonlat()
    LON, LAT = np.meshgrid(lon, lat)
    # Strong gradient concentrated at the poleward edge -- exactly the
    # shape that exposed the bug (steep drop-off near 79.5N, e.g. a nutrient
    # field going from ~30 to ~0 over the last few source latitudes).
    field = 30.0 * np.clip((LAT - 60.0) / 19.5, 0.0, 1.0)
    field = field.max() - field  # 30 at low lat, 0 at the poleward edge
    field[field == 0.0] = 1e-6  # avoid tripping the land-flag (exact 0.0)

    xc, yc, hfacc = target_grid()
    out = regrid_horizontal_2d(field, xc, yc, hfacc[0] > 0)
    polar = yc > lat.max()
    assert polar.any(), "test grid has no points poleward of the source domain"
    assert out[polar].min() >= -1e-3, \
        f"polar undershoot below source range: {out[polar].min()}"
    assert out[polar].max() <= field.max() + 1e-3, \
        f"polar overshoot above source range: {out[polar].max()}"
    print("polar clamp OK: no linear-extrapolation overshoot poleward of",
          lat.max(), "-- polar min/max", out[polar].min(), out[polar].max())


def _test_llc270_to_llc90_nearest_neighbor_matches_smooth_function():
    """Self-test against the REAL local LLC270 and LLC90 grids: build a
    synthetic field as a smooth function of lon/lat evaluated at LLC270's
    own grid points (not a rectangular grid -- genuinely exercises the
    tile-shaped nearest-neighbor path), regrid to LLC90, and check the
    result is close to the same smooth function evaluated directly at
    LLC90's own grid points. Nearest-neighbor isn't exact interpolation, so
    this allows a tolerance -- but LLC270 is 3x LLC90's resolution, so a
    correct implementation should track a smooth function closely.
    """
    xc270, yc270, hfacc270 = llc270_grid()
    xc90, yc90, hfacc90 = target_grid()

    def smooth_fn(lon, lat):
        return 10.0 + 5.0 * np.cos(np.radians(lat)) * np.sin(np.radians(lon))

    field_270_tiles = smooth_fn(xc270, yc270)[None, :, :, :] * np.ones((50, 1, 1, 1))
    field_270_compact = ecco.llc_tiles_to_compact(field_270_tiles, less_output=True)

    out_compact = regrid_llc270_to_llc90(field_270_compact.astype(np.float32))
    assert out_compact.shape == (50, 1170, 90), out_compact.shape
    assert np.isfinite(out_compact).all()

    out_tiles = ecco.llc_compact_to_tiles(out_compact, less_output=True)
    expected = smooth_fn(xc90, yc90)
    wet90 = hfacc90[0] > 0
    # Nearest-neighbor error bound: LLC270 cells are ~1/3 the size of LLC90
    # cells, so the smooth function changes little between a wet LLC90 point
    # and its nearest wet LLC270 neighbor -- 0.5 is generous versus the
    # field's own range (5..15).
    err = np.abs(out_tiles[0][wet90] - expected[wet90])
    assert err.max() < 0.5, f"nearest-neighbor error too large: {err.max()}"
    print("LLC270->LLC90 self-test OK: max error vs. smooth function",
          err.max(), "over", wet90.sum(), "wet LLC90 points")


def _test_regrid_monthly_2d_shape():
    rng = np.random.default_rng(1)
    stack = 10.0 + rng.normal(scale=0.1, size=(12, SRC_NY, SRC_NX))
    out = regrid_monthly_2d(stack.astype(np.float32))
    assert out.shape == (12, 1170, 90), out.shape
    assert np.isfinite(out).all()
    print("regrid_monthly_2d self-test OK: shape", out.shape)


def regrid_to_llc90_selftest():
    """Self-test with synthetic data against the REAL local grid geometry --
    Steph's actual files aren't available yet, so this validates the
    mechanics (shapes, longitude wraparound, land-fill, vertical
    extrapolation, compact round-trip) ahead of time.
    """
    _test_interp_vertical_matches_np_interp()
    _test_source_land_fill_removes_zeros_before_interp()
    _test_polar_points_use_nearest_edge_not_linear_extrapolation()
    _test_regrid_monthly_2d_shape()
    _test_llc270_to_llc90_nearest_neighbor_matches_smooth_function()

    lon, lat = src_lonlat()
    LON, LAT = np.meshgrid(lon, lat)

    # 2D synthetic field: smooth function of lon/lat, checkable by eye.
    synthetic_2d = 10.0 + 5.0 * np.cos(np.radians(LAT)) * np.sin(np.radians(LON))
    out_2d = regrid_field(synthetic_2d, is_3d=False)
    assert out_2d.shape == (1170, 90), out_2d.shape
    assert np.isfinite(out_2d).all(), "NaN/Inf in 2D regrid output"
    print("2D self-test OK: shape", out_2d.shape,
          "min/max", out_2d.min(), out_2d.max())

    # 3D synthetic field: same horizontal pattern, decaying with depth.
    synthetic_3d = synthetic_2d[None, :, :] * np.exp(-SRC_RC / 1000.0)[:, None, None]
    out_3d = regrid_field(synthetic_3d.astype(np.float32), is_3d=True)
    assert out_3d.shape == (50, 1170, 90), out_3d.shape
    assert np.isfinite(out_3d).all(), "NaN/Inf in 3D regrid output"
    print("3D self-test OK: shape", out_3d.shape,
          "min/max", out_3d.min(), out_3d.max())

    # Land points should now be finite and drawn from a nearby wet value,
    # not left at whatever raw interpolation-through-land produced.
    _, _, hfacc = target_grid()
    surf_tiles = ecco.llc_compact_to_tiles(out_2d, less_output=True)
    dry = hfacc[0] == 0
    print("dry-point fraction at surface:", dry.mean(),
          "; dry-point value range:", surf_tiles[dry].min(), surf_tiles[dry].max())


if __name__ == "__main__":
    regrid_to_llc90_selftest()
