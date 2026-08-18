#!/usr/bin/env python3
"""Same byte-patch as patch_ecco_cost.py, defaulting to run_bling_1yr.
See patch_ecco_cost.py for the full file-layout documentation.
Usage: patch_ecco_cost_1yr.py <cycle> <darwin_cost_value>
Configure via environment variable: BLING_RUN_DIR (default ../run_bling_1yr)
"""
import os
import struct
import sys

RUN_BLING = os.environ.get("BLING_RUN_DIR", "../run_bling_1yr")
YCTRLID = "MIT_CE_000"


def main(cycle, darwin_cost):
    path = f"{RUN_BLING}/ecco_cost_{YCTRLID}.opt{cycle:04d}"
    with open(path, "r+b") as f:
        f.seek(54)
        marker_lo = struct.unpack(">i", f.read(4))[0]
        f.seek(66)
        marker_hi = struct.unpack(">i", f.read(4))[0]
        if marker_lo != 8 or marker_hi != 8:
            raise RuntimeError(
                f"{path}: unexpected record markers ({marker_lo}, {marker_hi}), "
                "expected (8, 8) -- file layout may have changed, refusing to patch"
            )
        bling_fc = struct.unpack(">d", open(path, "rb").read()[58:66])[0]
        f.seek(58)
        f.write(struct.pack(">d", darwin_cost))
    print(f"{path}: patched fc {bling_fc:.6e} (BLING) -> {darwin_cost:.6e} (Darwin), "
          f"gradient payload (bytes 70+) untouched")


if __name__ == "__main__":
    main(int(sys.argv[1]), float(sys.argv[2]))
