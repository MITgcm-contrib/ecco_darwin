#!/usr/bin/env python3
"""Replace Darwin's independently-regridded (LLC270-derived) BGC initial
conditions with BLING's own IC fields, unit-converted (mol/m3 -> mmol/m3,
x1000). Both experiments share the same 128x64x15 grid, so no regridding
is needed -- straight read/convert/write.

Rationale: an earlier ~2x mismatch between Darwin's forward-FD and
BLING's adjoint sensitivity for DIC traced back to Darwin's DIC/ALK ICs
coming from a totally different, never-cross-validated source (a prior
Darwin LLC270 ecosystem run) than BLING's ICs (bundled with the MITgcm
BLING tutorial) -- a ~7% ALK / ~2% DIC mismatch already present at t=0
that barely grows over a 30-day integration (see
darwin_bling_comparison/readme.txt). Using BLING's own ICs for Darwin
eliminates that confound so the two models' sensitivities are compared
starting from an identical carbon-chemistry state.

Configure via environment variables:
  BLING_INPUT_AD_DIR  -- BLING_SOCAT's input_ad/ (source of *_init.bin)
  DARWIN_INPUT_AD_DIR -- Darwin's input_ad/ (destination)
"""
import os
import numpy as np

BLING_DIR = os.environ.get("BLING_INPUT_AD_DIR", "../../global_oce_biogeo_bling_SOCAT/input_ad")
DARWIN_DIR = os.environ.get("DARWIN_INPUT_AD_DIR", "../input_ad")
NZ, NY, NX = 15, 64, 128

# BLING source file -> Darwin target file (BLING units mol/m3 or eq/m3 -> Darwin mmol/m3)
FIELDS = {
    "dic_init.bin": "dic_darwin_init.bin",
    "alk_init.bin": "alk_darwin_init.bin",
    "o2_init.bin":  "o2_darwin_init.bin",
    "no3_init.bin": "no3_darwin_init.bin",
    "po4_init.bin": "po4_darwin_init.bin",
    "fe_init.bin":  "fet_darwin_init.bin",
}

for src_name, dst_name in FIELDS.items():
    src = np.fromfile(f"{BLING_DIR}/{src_name}", dtype=">f4").reshape(NZ, NY, NX)
    assert np.isfinite(src).all(), f"NaN/Inf in {src_name}"
    dst = (src.astype(np.float64) * 1000.0).astype(">f4")
    dst.tofile(f"{DARWIN_DIR}/{dst_name}")
    print(f"{src_name:14s} -> {dst_name:22s}  "
          f"surface min/mean/max: BLING={src[0].min():.4g}/{src[0].mean():.4g}/{src[0].max():.4g} "
          f"-> Darwin={dst[0].min():.4g}/{dst[0].mean():.4g}/{dst[0].max():.4g}")
