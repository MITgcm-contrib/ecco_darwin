#!/usr/bin/env python3
"""
gen_derived_tracer_pickup.py -- build an initial-condition pickup for a
pan_AO tracer that the parent LLC270 v5r1 solution does not carry, by
scaling an ALREADY-DOWNSCALED pan_AO pickup.

WHY IT WORKS THIS WAY
---------------------
pan_AO runs 33 tracers (Bertin/AO_ecosystem); the parent ran 31. The two
extras, `rDOC` (pan_AO index 9) and `CDOM` (index 21), have no parent
field, so `gen_pickups.py` could not produce them -- see DECISIONS.md
2026-09-30.

User decision (2026-09-30): seed both as a FRACTION OF PARENT DOC, and use
the SAME fractions for the initial condition and for the OBCS boundary
values so the two stay internally consistent.

This script derives from `forcings/pickups/pickup_pTr08.<iter>.data` --
the parent's DOC, ALREADY interpolated onto the pan_AO 1680x1680x90 grid
and already verified by `scan_pickups.py` (NaN-free, zero-set == land-set).
It deliberately does NOT re-run the Delaunay/Hinterp pipeline from the
parent grid, because:
  - that pipeline takes ~54 min and needs a compute node;
  - re-deriving would re-expose the 135,784 out-of-hull cells and the
    sigmaG=15 Gaussian, i.e. a fresh chance to introduce NaN;
  - scaling the verified field guarantees the extras land on EXACTLY the
    same grid with EXACTLY the same land mask as the other 31 tracers.
A multiplication cannot move a land/ocean boundary. That is the point.

FRACTIONS -- provenance, not invention
  --frac for CDOM defaults to 0.02, taken from this configuration's OWN
  `data.darwin` &DARWIN_CDOM_PARAMS `fracCDOM = 0.02`. That parameter is
  formally the fraction of remineralised POP contributing to CDOM, i.e. a
  process rate rather than an inventory ratio, so it is an ANCHOR and not
  a derivation -- but it is the only CDOM:organic-matter ratio anywhere in
  the setup, and it puts CDOM at the right order of magnitude.
  --frac for rDOC defaults to 0.50, i.e. refractory DOC is about half the
  DOC inventory. This is an order-of-magnitude choice, NOT a tuned value.
Both are single numbers on the command line precisely so they are cheap to
revise. Revise them here and re-run; do not hand-edit the outputs.

CLIPPING
  The parent-derived tracers carry small negative values from the
  interpolation (DOC reaches -0.0898). The other 31 pickups were not
  clipped, so for consistency this script does not clip by default --
  EXCEPT that a negative CDOM would feed a negative absorption through
  `CDOMcoeff`, so `--clip-min 0` is recommended for CDOM and the number of
  clipped cells is always reported.

Processes level by level via memmap: the field is 1.89 GB and there is no
reason to hold it resident.
"""

import argparse
import os
import re
import shutil
import sys

import numpy as np

DIMLIST_RE = re.compile(r"dimList\s*=\s*\[(.*?)\]", re.S)
PREC_RE = re.compile(r"dataprec\s*=\s*\[\s*'([^']+)'\s*\]")
NREC_RE = re.compile(r"nrecords\s*=\s*\[\s*(\d+)\s*\]")


def parse_meta(path):
    with open(path) as fh:
        txt = fh.read()
    m = DIMLIST_RE.search(txt)
    if not m:
        sys.exit(f"FATAL: no dimList in {path}")
    nums = [int(x) for x in re.findall(r"-?\d+", m.group(1))]
    if len(nums) % 3 != 0:
        sys.exit(f"FATAL: malformed dimList in {path}")
    # triples are (global, start, end), x fastest-varying first
    dims = [nums[i] for i in range(0, len(nums), 3)]
    p = PREC_RE.search(txt)
    prec = p.group(1) if p else "float64"
    n = NREC_RE.search(txt)
    nrec = int(n.group(1)) if n else 1
    return dims, prec, nrec, txt


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--pickup-dir",
                    default="/nobackup/dcarrol2/downscaling/pan_AO/forcings/pickups")
    ap.add_argument("--source", default="pickup_pTr08",
                    help="source pickup basename (default parent DOC = pTr08)")
    ap.add_argument("--name", required=True, help="output tracer name, e.g. rDOC")
    ap.add_argument("--frac", type=float, required=True,
                    help="multiply source by this fraction")
    ap.add_argument("--iter", type=int, default=736344)
    ap.add_argument("--clip-min", type=float, default=None,
                    help="clip output below this value (recommended 0 for CDOM)")
    ap.add_argument("--force", action="store_true",
                    help="overwrite an existing output")
    args = ap.parse_args()

    stem = f"{args.iter:010d}"
    src_data = os.path.join(args.pickup_dir, f"{args.source}.{stem}.data")
    src_meta = os.path.join(args.pickup_dir, f"{args.source}.{stem}.meta")
    out_data = os.path.join(args.pickup_dir, f"pickup_{args.name}.{stem}.data")
    out_meta = os.path.join(args.pickup_dir, f"pickup_{args.name}.{stem}.meta")

    for f in (src_data, src_meta):
        if not os.path.exists(f):
            sys.exit(f"FATAL: source missing: {f}")
    if os.path.exists(out_data) and not args.force:
        sys.exit(f"FATAL: {out_data} exists; pass --force to overwrite")
    if not (0.0 < args.frac <= 1.0):
        sys.exit(f"FATAL: --frac {args.frac} outside (0, 1]; refusing")

    dims, prec, nrec, meta_txt = parse_meta(src_meta)
    if len(dims) != 3:
        sys.exit(f"FATAL: expected 3 dims, got {dims}")
    nx, ny, nz = dims
    dt = {"float64": ">f8", "float32": ">f4"}.get(prec)
    if dt is None:
        sys.exit(f"FATAL: unhandled dataprec {prec}")
    itemsize = np.dtype(dt).itemsize
    expect = nx * ny * nz * nrec * itemsize
    actual = os.path.getsize(src_data)
    if actual != expect:
        sys.exit(f"FATAL: {src_data} is {actual} B, expected "
                 f"{nx}x{ny}x{nz}x{nrec}x{itemsize} = {expect} B")

    print(f"source : {args.source}  {nx}x{ny}x{nz}  {prec}  ({actual/2**30:.2f} GiB)")
    print(f"output : pickup_{args.name}  = {args.frac} x source"
          + (f", clipped at >= {args.clip_min}" if args.clip_min is not None else ""))

    src = np.memmap(src_data, dtype=dt, mode="r", shape=(nz, ny, nx))
    out = np.memmap(out_data, dtype=dt, mode="w+", shape=(nz, ny, nx))

    n_nan = 0
    n_clip = 0
    vmin, vmax = np.inf, -np.inf
    for k in range(nz):
        lev = np.asarray(src[k], dtype=np.float64) * args.frac
        n_nan += int(np.count_nonzero(~np.isfinite(lev)))
        if args.clip_min is not None:
            mask = lev < args.clip_min
            n_clip += int(np.count_nonzero(mask))
            lev[mask] = args.clip_min
        out[k] = lev
        vmin = min(vmin, float(np.nanmin(lev)))
        vmax = max(vmax, float(np.nanmax(lev)))
    out.flush()
    del out, src

    # meta: copy the source's verbatim -- same dims, prec, nrecords, iter
    shutil.copyfile(src_meta, out_meta)

    print(f"  range  : {vmin:.6g} .. {vmax:.6g}")
    print(f"  NaN    : {n_nan}  (must be 0)")
    if args.clip_min is not None:
        print(f"  clipped: {n_clip} cells raised to {args.clip_min}")
    if n_nan:
        sys.exit("FATAL: NaN in output -- source was supposed to be NaN-free")
    print(f"  wrote  : {out_data}")
    print(f"  wrote  : {out_meta}")
    print("DERIVED_PICKUP_OK")
    return 0


if __name__ == "__main__":
    sys.exit(main())
