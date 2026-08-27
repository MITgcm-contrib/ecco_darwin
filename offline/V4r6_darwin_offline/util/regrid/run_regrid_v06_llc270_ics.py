"""Regrid Steph's real LLC270-native IC files onto LLC90, via
regrid_to_llc90.regrid_llc270_to_llc90 (nearest-wet-neighbor, NOT the
lat-lon pipeline -- these are native LLC270-compact binaries, confirmed by
byte count against the real LLC270 hFacC.data: 270*3510*50*4 bytes).

DIC and ALK use these real LLC270-native fields instead of the 1-degree
GLODAP climatology. Only those two: NO3/PO4/FeT/SiO2/O2 were briefly also
switched to this source (steph_input/v06/ICs/ has all 7), then reverted. NO3/PO4/SiO2/O2 (GLODAP/Levitus)
and FeT (Steph's 1-degree DARWIN_ECCO_run33_4 restart, previously ptr06)
are back on run_regrid_steph_ics.py's 1-degree pipeline. The other 14
PTRACERS_initialFile fields (NO2, NH4, DOC, DON, DOP, DOFe, POC, PON, POP,
POFe, POSi, PIC, CDOM, plus c07-c10 biomass) never had an LLC270 version
delivered and are unaffected either way.

Output lands in regridded_llc90/ alongside the 1-degree-regridded files,
named <TRACER>_llc270_to_llc90.bin so the two regridding paths are never
ambiguous from the filename alone.
"""
import os
import numpy as np
import regrid_to_llc90 as r

SRC_DIR = "steph_input/v06/ICs"
OUT_DIR = "regridded_llc90"

# (source filename, output filename, PTRACERS_initialFile index, PTRACERS name)
LLC270_FILES = [
    ("DIC.R4_270x3510x50.bin", "DIC_llc270_to_llc90.bin", 1, "DIC"),
    ("ALK.R4_270x3510x50.bin", "ALK_llc270_to_llc90.bin", 18, "ALK"),
]

SRC_SHAPE = (50, 3510, 270)  # (nz, ny_compact=13*270, nx_compact=270)


def regrid_llc270_file(name, out_name, tracer_name):
    src_path = os.path.join(SRC_DIR, name)
    raw = np.fromfile(src_path, dtype=">f4")
    expected = SRC_SHAPE[0] * SRC_SHAPE[1] * SRC_SHAPE[2]
    assert raw.size == expected, f"{name}: got {raw.size} elements, expected {expected}"
    field = raw.reshape(SRC_SHAPE)
    out = r.regrid_llc270_to_llc90(field)
    assert out.shape == (50, 1170, 90), out.shape
    assert np.isfinite(out).all(), f"{name}: non-finite values after regrid"
    out_path = os.path.join(OUT_DIR, out_name)
    r.write_compact_binary(out, out_path)
    print(f"{tracer_name}: {name} -> {out_name}: shape {out.shape}, "
          f"min/max {out.min():.6g}/{out.max():.6g}")


if __name__ == "__main__":
    os.makedirs(OUT_DIR, exist_ok=True)
    for src, out, idx, name in LLC270_FILES:
        regrid_llc270_file(src, out, name)
    print(f"\nDone: {len(LLC270_FILES)} LLC270-native IC fields regridded to {OUT_DIR}/")
