#!/usr/bin/env python3
"""Same byte-patch as hybrid_darwin_bling/scripts/patch_ecco_cost.py,
pointed at this project's own run_dic/. See that script / README.md for
the full binary-file-layout documentation (verified byte-for-byte against
real M1QN3 campaign files): bytes [54:58) and [66:70) are int32 Fortran
record markers (must both == 8), bytes [58:66) is fc as big-endian
real*8/float64, bytes [70:...) is the packed gradient/control payload
(untouched by this patch -- pkg/dic's own gradient stays in place, only
its cost value fc is overwritten with Darwin's own).

Usage: patch_ecco_cost_dic.py <cycle> <darwin_cost_value>
"""
import struct
import sys

RUN_DIC = "/Users/carrolld/Documents/research/adjoint/MITgcm/verification/hybrid_darwin_dic/run_dic"
YCTRLID = "MIT_CE_000"


def main(cycle, darwin_cost):
    path = f"{RUN_DIC}/ecco_cost_{YCTRLID}.opt{cycle:04d}"
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
        dic_fc = struct.unpack(">d", open(path, "rb").read()[58:66])[0]
        f.seek(58)
        f.write(struct.pack(">d", darwin_cost))
    print(f"{path}: patched fc {dic_fc:.6e} (pkg/dic) -> {darwin_cost:.6e} (Darwin), "
          f"gradient payload (bytes 70+) untouched")


if __name__ == "__main__":
    main(int(sys.argv[1]), float(sys.argv[2]))
