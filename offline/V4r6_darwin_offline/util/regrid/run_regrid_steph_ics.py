"""Run regrid_to_llc90.py's real pipeline on Steph's actual delivered
1-degree files (steph_input/ICs/, steph_input/forcing/).

As of, three fields are deliberately NOT regridded by this
script even though Steph's own 1-degree versions arrived and could be:

- Iron: data.darwin now sources it from ecco_darwin's copy of the same
  Hamilton et al. 2020 product instead, read via pkg/darwin's own runtime
  horizontal interpolation (iron_lon0/lat0/nlon/nlat) rather than a
  pre-regridded file -- see data.darwin's own comments.
- icefile/SIarea comes from this project's own V4r6 archive, not from
  Steph at all -- see readme.txt.
- DIC, ALK: superseded by real LLC270-native fields from
  steph_input/v06/ICs/, regridded by run_regrid_v06_llc270_ics.py instead
  (nearest-wet-neighbor on the native cube-sphere grid -- a fundamentally
  different regrid than this script's lat-lon RegularGridInterpolator
  path). Only DIC and ALK come from that LLC270 source:
  NO3/PO4/FeT/SiO2/O2 were briefly switched to it as well, and are
  reverted here to their original 1-degree sources (lev01_nitrate/
  lev01_phosphate/ptr06/lev01_silicate/lev01_oxygen_micromolar), and the
  corresponding *_llc270_to_llc90.bin outputs were deleted as no longer
  referenced.

OASIM irradiance (switched to pkg/oasim, not file-based at all) and the
wavelength-indexed optics tables (still pending from Steph) are also not
touched here.

Output lands in regridded_llc90/, named to match the TODO_regrid_* placeholders
in input_darwin_offline/{data.ptracers,data.darwin} minus the "TODO_regrid_"
prefix -- readme.txt Part 7 symlinks this directory's contents in whole.
"""
import os
import numpy as np
import regrid_to_llc90 as r

SRC_DIR_ICS = "steph_input/ICs"
SRC_DIR_FORCING = "steph_input/forcing"
OUT_DIR = "regridded_llc90"

# (source filename, output filename, is_3d) for the 19 PTRACERS_initialFile
# fields sourced from Steph's 1-degree files -- mapping cross-checked
# against Steph's original 1-deg set-up (data.ptracers) PTRACERS_initialFile
# entries. DIC/ALK removed (see module docstring -- superseded by
# run_regrid_v06_llc270_ics.py); NO3/PO4/FeT/SiO2/O2 reinstated here after
# being briefly switched to LLC270 sources and then reverted.
PTRACER_FILES = [
    ("lev01_nitrate_ann-3d.bin", "lev01_nitrate_ann-3d_llc90.bin"),
    ("ptr03_run33_4_28800.bin", "ptr03_run33_4_28800_llc90.bin"),
    ("ptr04_run33_4_28800.bin", "ptr04_run33_4_28800_llc90.bin"),
    ("lev01_phosphate_ann-3d.bin", "lev01_phosphate_ann-3d_llc90.bin"),
    ("ptr06_run33_4_28800.bin", "ptr06_run33_4_28800_llc90.bin"),
    ("lev01_silicate_ann-3d.bin", "lev01_silicate_ann-3d_llc90.bin"),
    ("ptr08_run33_4_28800.bin", "ptr08_run33_4_28800_llc90.bin"),
    ("ptr09_run33_4_28800.bin", "ptr09_run33_4_28800_llc90.bin"),
    ("ptr10_run33_4_28800.bin", "ptr10_run33_4_28800_llc90.bin"),
    ("ptr11_run33_4_28800.bin", "ptr11_run33_4_28800_llc90.bin"),
    ("ptr12_run33_4_28800.bin", "ptr12_run33_4_28800_llc90.bin"),
    ("ptr13_run33_4_28800.bin", "ptr13_run33_4_28800_llc90.bin"),
    ("ptr14_run33_4_28800.bin", "ptr14_run33_4_28800_llc90.bin"),
    ("ptr15_run33_4_28800.bin", "ptr15_run33_4_28800_llc90.bin"),
    ("ptr16_run33_4_28800.bin", "ptr16_run33_4_28800_llc90.bin"),
    ("ptr17_run33_4_28800.bin", "ptr17_run33_4_28800_llc90.bin"),
    ("lev01_oxygen_micromolar_ann-3d.bin", "lev01_oxygen_micromolar_ann-3d_llc90.bin"),
    ("ptr20_run33_4_28800.bin", "ptr20_run33_4_28800_llc90.bin"),
    ("biomass_run33_4_28800.bin", "biomass_run33_4_28800_llc90.bin"),
]

# (source filename, output filename) for the monthly-climatology 2D forcing
# fields, cross-checked against Steph's original 1-deg set-up (data.darwin.)
# hamilton20_monthmeanSFe.bin (iron) deliberately NOT included -- see the
# module docstring: data.darwin now sources iron from ecco_darwin instead.
FORCING_FILES = [
    ("tren_speed_mth-2d.bin", "tren_speed_mth-2d_llc90.bin"),
]


def regrid_ptracer_file(name, out_name):
    src_path = os.path.join(SRC_DIR_ICS, name)
    raw = np.fromfile(src_path, dtype=">f4")
    expected = r.SRC_NX * r.SRC_NY * 23
    assert raw.size == expected, f"{name}: got {raw.size} elements, expected {expected}"
    field = raw.reshape(23, r.SRC_NY, r.SRC_NX)
    out = r.regrid_field(field, is_3d=True)
    assert out.shape == (50, 1170, 90), out.shape
    assert np.isfinite(out).all(), f"{name}: non-finite values after regrid"
    out_path = os.path.join(OUT_DIR, out_name)
    r.write_compact_binary(out, out_path)
    print(f"{name} -> {out_name}: shape {out.shape}, "
          f"min/max {out.min():.6g}/{out.max():.6g}")


def regrid_forcing_file(name, out_name):
    src_path = os.path.join(SRC_DIR_FORCING, name)
    raw = np.fromfile(src_path, dtype=">f4")
    expected = r.SRC_NX * r.SRC_NY * 12
    assert raw.size == expected, f"{name}: got {raw.size} elements, expected {expected}"
    stack = raw.reshape(12, r.SRC_NY, r.SRC_NX)
    out = r.regrid_monthly_2d(stack)
    assert out.shape == (12, 1170, 90), out.shape
    assert np.isfinite(out).all(), f"{name}: non-finite values after regrid"
    out_path = os.path.join(OUT_DIR, out_name)
    r.write_compact_binary(out, out_path)
    print(f"{name} -> {out_name}: shape {out.shape}, "
          f"min/max {out.min():.6g}/{out.max():.6g}")


if __name__ == "__main__":
    os.makedirs(OUT_DIR, exist_ok=True)
    for src, out in PTRACER_FILES:
        regrid_ptracer_file(src, out)
    for src, out in FORCING_FILES:
        regrid_forcing_file(src, out)
    print(f"\nDone: {len(PTRACER_FILES)} ptracer ICs + {len(FORCING_FILES)} "
          f"forcing fields regridded to {OUT_DIR}/")
