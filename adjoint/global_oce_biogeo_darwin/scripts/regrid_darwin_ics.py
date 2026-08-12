"""Regrid Darwin ecosystem tracer ICs from LLC90-compact onto the
128x64 lat-lon grid used by global_oce_biogeo_bling/_SOCAT/_darwin.

Source: real Darwin IC fields (a prior Darwin ecosystem run, LLC270
regridded to LLC90), 50-level LLC90-compact format
((50, 1170, 90), i.e. (nz, 13*90, 90)).

Method: per-level nearest-WET-neighbor via a KDTree over each point's
Cartesian unit-sphere embedding (handles the dateline/poles with no
special-casing, restricts candidates to wet source cells so dry LLC90
land values never leak into a wet target cell). Horizontal regrid first
(50 independent 2D nearest-neighbor lookups against LLC90's own
hFacC-derived wet mask), THEN vertical interpolation from LLC90's 50
levels onto the target's 15 (the two steps commute for a separable
regridding, and horizontal-first is far cheaper).

Land handling on the output side: deliberately NOT re-masked against the
target's own bathy.bin here -- MITgcm's ptracers_init_varia.F
unconditionally zeros every tracer where maskC==0 at cold start
regardless of what the initial file contains, so only wet target-grid
values need to be correct, and the per-level nearest-wet-neighbor search
already guarantees that.

Output: global (15, 64, 128) big-endian float32 binaries, matching the
MDS convention already used for the other *_init.bin files in this
project (NOT the 4-tile-decomposed xx_*/adxx_* control-vector format).

Requires: MITgcmutils (ships with MITgcm under utils/python/MITgcmutils
-- add that directory to PYTHONPATH) and ecco_v4_py.

Configure via environment variables (see defaults below) rather than
editing this file:
  LLC90_GRID_DIR  -- directory with LLC90 XC/YC/hFacC (ECCOv4-style)
  DARWIN_IC_SRC_DIR -- directory with the source LLC90-compact IC fields
  DARWIN_INPUT_AD_DIR -- output directory (this experiment's input_ad/)
"""
import os
import numpy as np
import MITgcmutils as mu
import ecco_v4_py as ecco
from scipy.spatial import cKDTree

GRID_DIR = os.environ.get("LLC90_GRID_DIR", "./ECCO_V4r5_grid")
SRC_DIR = os.environ.get("DARWIN_IC_SRC_DIR", "./regridded_llc90")
OUT_DIR = os.environ.get("DARWIN_INPUT_AD_DIR", "../input_ad")

# Target grid: 128x64, 2.8125deg, matching global_oce_biogeo_bling's bathy.bin.
TGT_NX, TGT_NY, TGT_NZ = 128, 64, 15
TGT_DELX = TGT_DELY = 2.8125

# LLC90's real 50-level delR (ECCOv4r6), same array already validated in
# the reference regrid_to_llc90.py's TARGET_DELR.
SRC_DELR = np.array([
    10.00, 10.00, 10.00, 10.00, 10.00, 10.00, 10.00, 10.01,
    10.03, 10.11, 10.32, 10.80, 11.76, 13.42, 16.04, 19.82, 24.85,
    31.10, 38.42, 46.50, 55.00, 63.50, 71.58, 78.90, 85.15, 90.18,
    93.96, 96.58, 98.25, 99.25, 100.01, 101.33, 104.56, 111.33, 122.83,
    139.09, 158.94, 180.83, 203.55, 226.50, 249.50, 272.50, 295.50, 318.50,
    341.50, 364.50, 387.50, 410.50, 433.50, 456.50,
])
assert SRC_DELR.size == 50
SRC_RC = np.cumsum(SRC_DELR) - SRC_DELR / 2.0  # 50 cell-center depths

# Target 15-level delZ (global_oce_biogeo_bling's data PARM04), converted
# to cell-center depths the same way.
TGT_DELZ = np.array([
    50., 70., 100., 140., 190.,
    240., 290., 340., 390., 440.,
    490., 540., 590., 640., 690.,
])
assert TGT_DELZ.size == TGT_NZ
TGT_RC = np.cumsum(TGT_DELZ) - TGT_DELZ / 2.0  # 15 cell-center depths

# The 21 fields this regrid needs, mapping source filename -> local output
# name (matches data.ptracers' PTRACERS_initialFile entries 1:1).
FIELDS = {
    "DIC_llc270_to_llc90.bin": "dic_darwin_init.bin",
    "lev01_nitrate_ann-3d_llc90.bin": "no3_darwin_init.bin",
    "ptr03_run33_4_28800_llc90.bin": "no2_darwin_init.bin",
    "ptr04_run33_4_28800_llc90.bin": "nh4_darwin_init.bin",
    "lev01_phosphate_ann-3d_llc90.bin": "po4_darwin_init.bin",
    "ptr06_run33_4_28800_llc90.bin": "fet_darwin_init.bin",
    "lev01_silicate_ann-3d_llc90.bin": "sio2_darwin_init.bin",
    "ptr08_run33_4_28800_llc90.bin": "doc_darwin_init.bin",
    "ptr09_run33_4_28800_llc90.bin": "don_darwin_init.bin",
    "ptr10_run33_4_28800_llc90.bin": "dop_darwin_init.bin",
    "ptr11_run33_4_28800_llc90.bin": "dofe_darwin_init.bin",
    "ptr12_run33_4_28800_llc90.bin": "poc_darwin_init.bin",
    "ptr13_run33_4_28800_llc90.bin": "pon_darwin_init.bin",
    "ptr14_run33_4_28800_llc90.bin": "pop_darwin_init.bin",
    "ptr15_run33_4_28800_llc90.bin": "pofe_darwin_init.bin",
    "ptr16_run33_4_28800_llc90.bin": "posi_darwin_init.bin",
    "ptr17_run33_4_28800_llc90.bin": "pic_darwin_init.bin",
    "ALK_llc270_to_llc90.bin": "alk_darwin_init.bin",
    "lev01_oxygen_micromolar_ann-3d_llc90.bin": "o2_darwin_init.bin",
    "ptr20_run33_4_28800_llc90.bin": "cdom_darwin_init.bin",
    "biomass_run33_4_28800_llc90.bin": "biomass_darwin_init.bin",
}


def _latlon_to_unit_sphere_xyz(lon, lat):
    lon_r = np.radians(lon)
    lat_r = np.radians(lat)
    x = np.cos(lat_r) * np.cos(lon_r)
    y = np.cos(lat_r) * np.sin(lon_r)
    z = np.sin(lat_r)
    return np.stack([x, y, z], axis=-1)


def target_lonlat():
    """Target grid's cell-center lon [0,360) / lat [-90,90] axes."""
    lon = (np.arange(TGT_NX) + 0.5) * TGT_DELX
    lat = -90.0 + (np.arange(TGT_NY) + 0.5) * TGT_DELY
    return lon, lat


def llc90_grid():
    xc = ecco.llc_compact_to_tiles(mu.rdmds(GRID_DIR + "/XC"), less_output=True)
    yc = ecco.llc_compact_to_tiles(mu.rdmds(GRID_DIR + "/YC"), less_output=True)
    hfacc = ecco.llc_compact_to_tiles(mu.rdmds(GRID_DIR + "/hFacC"), less_output=True)
    return xc, yc, hfacc  # xc/yc: (13,90,90); hfacc: (50,13,90,90)


def _interp_vertical(horiz):
    """horiz: (50, TGT_NY, TGT_NX) on SRC_RC depths -> (15, TGT_NY, TGT_NX)
    on TGT_RC depths. Constant extrapolation past the source's depth range
    (the target's deepest cell center is well within LLC90's 50-level
    range, so this only matters at the shallow end in practice, but kept
    for safety/symmetry with the reference interpolator)."""
    idx = np.searchsorted(SRC_RC, TGT_RC)
    idx_lo = np.clip(idx - 1, 0, SRC_RC.size - 1)
    idx_hi = np.clip(idx, 0, SRC_RC.size - 1)

    horiz_flat = horiz.reshape(horiz.shape[0], -1)
    lo_vals, hi_vals = horiz_flat[idx_lo], horiz_flat[idx_hi]
    rc_lo, rc_hi = SRC_RC[idx_lo], SRC_RC[idx_hi]
    span = rc_hi - rc_lo
    w_hi = np.where(span > 0, (TGT_RC - rc_lo) / np.where(span > 0, span, 1.0), 0.0)

    out_flat = lo_vals + w_hi[:, None] * (hi_vals - lo_vals)
    return out_flat.reshape((TGT_NZ,) + horiz.shape[1:])


def regrid_field(src_path, xc90, yc90, hfacc90, xyz90_flat_by_level):
    """Regrid one (50,1170,90) LLC90-compact field onto the target grid."""
    field_compact = np.fromfile(src_path, dtype=">f4").reshape(50, 1170, 90)
    field_tiles = ecco.llc_compact_to_tiles(field_compact, less_output=True)  # (50,13,90,90)

    lon_t, lat_t = target_lonlat()
    LON, LAT = np.meshgrid(lon_t, lat_t)  # (64,128) each
    xyz_t = _latlon_to_unit_sphere_xyz(LON, LAT).reshape(-1, 3)  # (64*128, 3)

    horiz = np.empty((50, TGT_NY, TGT_NX), dtype=np.float64)
    for k in range(50):
        wet90 = hfacc90[k] > 0
        if not wet90.any():
            horiz[k] = 0.0
            continue
        tree = cKDTree(xyz90_flat_by_level[k]) if xyz90_flat_by_level[k] is not None \
            else cKDTree(_latlon_to_unit_sphere_xyz(xc90, yc90)[wet90])
        _, idx = tree.query(xyz_t)
        horiz[k] = field_tiles[k][wet90][idx].reshape(TGT_NY, TGT_NX)

    return _interp_vertical(horiz)


def write_binary(field, out_path):
    """field: (15,64,128) -> big-endian float32, matching the MDS
    convention already used for the other init.bin/wt_*.bin files."""
    field.astype(">f4").tofile(out_path)


def main():
    xc90, yc90, hfacc90 = llc90_grid()
    xyz90_by_level = [None] * 50  # regrid_field builds each level's tree lazily

    for src_name, out_name in FIELDS.items():
        src_path = f"{SRC_DIR}/{src_name}"
        out_path = f"{OUT_DIR}/{out_name}"
        out = regrid_field(src_path, xc90, yc90, hfacc90, xyz90_by_level)
        assert out.shape == (TGT_NZ, TGT_NY, TGT_NX), out.shape
        assert np.isfinite(out).all(), f"NaN/Inf in regridded {src_name}"
        write_binary(out, out_path)
        print(f"{src_name:45s} -> {out_name:25s}  "
              f"surface min/mean/max = {out[0].min():.4g}/{out[0].mean():.4g}/{out[0].max():.4g}")


if __name__ == "__main__":
    main()
