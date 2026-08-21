#!/usr/bin/env python3
"""Binary-patch the scalar cost value inside BLING's ecco_cost_<yctrlid>.opt<cycle>
file, replacing it with Darwin's independently-computed cost, while
leaving the packed adjoint-gradient payload (bytes 70 onward) untouched.
This is the surrogate-gradient injection point: optim.x subsequently
reads (Darwin cost, BLING gradient) as if they were consistent.

File layout (Fortran unformatted sequential, big-endian, verified
byte-for-byte against a real 10-cycle SOCAT campaign's
ecco_cost_MIT_CE_000.opt0000..opt0010):
  bytes [54:58)  record-length marker, must read as int32 == 8
  bytes [58:66)  fc (real*8, big-endian) -- THIS is what we overwrite
  bytes [66:70)  record-length marker, must read as int32 == 8
  bytes [70:...) packed gradient/control payload -- NOT touched

Usage: patch_ecco_cost.py <cycle> <darwin_cost_value>
Configure via environment variable: BLING_RUN_DIR (default ../run_bling)
"""
import os
import struct
import sys

RUN_BLING = os.environ.get("BLING_RUN_DIR", "../run_bling")
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
