#!/usr/bin/env python3
"""
gen_derived_tracer_obcs.py -- build OBCS boundary files for the two pan_AO
tracers the parent LLC270 v5r1 solution does not carry, by scaling the
ALREADY-GENERATED parent DOC boundary files.

This is the OBCS counterpart of `gen_derived_tracer_pickup.py`, and it
exists for the same reason: pan_AO runs 33 tracers (Bertin/AO_ecosystem),
the parent ran 31, and the two extras are INSERTED mid-list --

    rDOC  at pan_AO index  9   (refractory DOC, terrestrial/riverine)
    CDOM  at pan_AO index 21   (dominates Arctic shelf light attenuation)

-- so there is no parent field to interpolate from. See DECISIONS.md
2026-09-30.

User decision (2026-09-30): seed both as a FRACTION OF PARENT DOC, and use
the SAME fractions for the initial condition and for the OBCS boundary
values so the two stay internally consistent. `gen_derived_tracer_pickup.py`
already did the IC half with rDOC=0.50, CDOM=0.02; this script does the
boundary half with the same defaults. **If you change a fraction, change it
in BOTH places** -- an IC and a boundary that disagree will drive a spurious
gradient in from the rim for the whole run.

WHY IT SCALES THE OBCS FILE AND NOT THE dv FILE
-----------------------------------------------
The source is `forcings/OBCS/TRAC<nn>_<bnd>.bin`, i.e. parent DOC ALREADY
projected onto the pan_AO disc rim by `gen_obcs_fast.py -disc on` and
already passed by `scan_obcs.py` (finite, nothing outside the disc open
boundary, varying in time). Deriving from the finished boundary file rather
than re-running the projection from `parent_run/run/dv/`:

  - guarantees the extras land on EXACTLY the same rim cells, with exactly
    the same zero/nonzero footprint, as the other 31 tracers. A
    multiplication cannot move a boundary. That is the point;
  - avoids a second trip through the disc gather, i.e. a fresh chance to
    reintroduce the `_MASK_FACE` class of bug;
  - costs ~1 min of streaming instead of a compute-node job.

WHICH FILE IS DOC -- RESOLVED BY NAME, NEVER BY INDEX
------------------------------------------------------
`TRAC<nn>_<bnd>.bin` is named for the PARENT's ptracer slot nn (
`gen_obcs_fast.py:_dv_read_name` maps TRACnn -> the dv file PTRACEnn, and
the parent's dv masks are emitted in parent ptracer order). pan_AO's own
index for DOC is 8 as well, but ONLY because DOC sits before the first
insertion -- that coincidence does not survive to any tracer after index 8,
and hardcoding `TRAC08` here would be the same latent error
`build_ptracer_map.py` exists to prevent. So the DOC slot is looked up by
name in the PARENT's `data.ptracers`, and the script aborts if the name is
absent or ambiguous.

CLIPPING -- deliberately OFF, to match the initial condition
------------------------------------------------------------
The parent-derived fields carry small negative values from the
interpolation. `gen_derived_tracer_pickup.py`'s docstring RECOMMENDS
`--clip-min 0` for CDOM, but the ICs that were actually built on 2026-10-01
were NOT clipped -- `pickup_CDOM` has range -0.00180642 .. 1.23824 (see
pickups_out.txt). So clipping here would make the boundary condition
disagree with the interior it is supposed to match, which is the one thing
the "same fractions for IC and OBCS" decision was meant to prevent.

Default is therefore NO clipping, matching both the IC and the other 31
unclipped OBCS tracer files (parent DOC itself reaches -0.0898). The
negative count and minimum are ALWAYS reported so the choice stays visible.

If negative CDOM turns out to matter -- it feeds a negative absorption
through `CDOMcoeff` -- clip BOTH sides: re-run `gen_derived_tracer_pickup.py
--clip-min 0` for the IC and this script with `--clip-min-cdom 0`. Clipping
only one side is worse than clipping neither.

Streams record-block by record-block: each file is 4.95 GB and there is no
reason to hold one resident.
"""

import argparse
import os
import re
import sys

import numpy as np

# Must match scan_obcs.py / gen_obcs_fast.py. NR is the pan_AO CHILD grid's
# 90 levels, not the parent's 50 -- see the note in scan_obcs.py.
NX = 1680
NR = 90
NREC = 8784
DTYPE = '>f4'
CHUNK = 256          # records per block; 256*90*1680*4 B = 155 MB

BOUNDARIES = ('east', 'west', 'north', 'south')

NAME_RE = re.compile(r"PTRACERS_names\(\s*(\d+)\s*\)\s*=\s*'([^']*)'")

# Defaults deliberately identical to the ICs that were actually built
# (pickups_out.txt: rDOC = 0.5 x source, CDOM = 0.02 x source, NEITHER
# clipped). Do not "fix" clip_min here without also rebuilding the IC.
DERIVED = {
    'rDOC': dict(frac=0.50, clip_min=None),
    'CDOM': dict(frac=0.02, clip_min=None),
}

SOURCE_SPECIES = 'DOC'


def parent_slot(path, species):
    """Return the parent's 1-based ptracer index for `species`, by NAME."""
    with open(path) as fh:
        txt = fh.read()
    names = {}
    for idx, nm in NAME_RE.findall(txt):
        nm = nm.strip()
        if nm:
            names[int(idx)] = nm
    if not names:
        sys.exit(f"FATAL: no PTRACERS_names entries in {path}")
    hits = sorted(i for i, nm in names.items() if nm == species)
    if len(hits) != 1:
        sys.exit(f"FATAL: parent {path} has {len(hits)} tracers named "
                 f"'{species}' (expected exactly 1): {hits}. Refusing to "
                 f"guess which boundary file holds it.")
    return hits[0], names


def derive(src, dst, frac, clip_min, force, nx=NX, nrec=NREC):
    if os.path.exists(dst) and not force:
        sys.exit(f"FATAL: {dst} exists. Re-running would silently re-scale an "
                 f"already-derived file. Delete it or pass --force.")
    nbytes = os.path.getsize(src)
    want = nx * NR * nrec * 4
    if nbytes != want:
        sys.exit(f"FATAL: {src} is {nbytes} B, expected {want} "
                 f"({nrec} rec x {NR} lev x {nx}). Wrong grid or truncated "
                 f"file -- not scaling it.")

    a = np.memmap(src, dtype=DTYPE, mode='r').reshape(nrec, NR, nx)
    nclip = 0
    nneg = 0
    nz_src = 0
    nz_out = 0
    vmin, vmax = np.inf, -np.inf

    tmp = dst + '.partial'
    with open(tmp, 'wb') as fh:
        for i in range(0, nrec, CHUNK):
            blk = np.asarray(a[i:i + CHUNK], dtype='f8')
            if not np.isfinite(blk).all():
                sys.exit(f"FATAL: non-finite value in source {src} near "
                         f"record {i}. scan_obcs.py passed, so this means "
                         f"the file changed under us -- aborting.")
            nz_src += int((blk != 0).sum())
            out = blk * frac
            if clip_min is not None:
                below = out < clip_min
                nclip += int(below.sum())
                out = np.where(below, clip_min, out)
            nneg += int((out < 0).sum())
            nz = out != 0
            nz_out += int(nz.sum())
            if nz.any():
                vmin = min(vmin, float(out[nz].min()))
                vmax = max(vmax, float(out[nz].max()))
            out.astype(DTYPE).tofile(fh)

    got = os.path.getsize(tmp)
    if got != want:
        os.unlink(tmp)
        sys.exit(f"FATAL: wrote {got} B to {tmp}, expected {want}")
    os.replace(tmp, dst)
    return dict(nclip=nclip, nneg=nneg, nz_src=nz_src, nz_out=nz_out,
                vmin=vmin, vmax=vmax)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('-d', '--config-dir',
                    default='/nobackup/dcarrol2/downscaling/pan_AO')
    ap.add_argument('--parent-ptracers',
                    default='/nobackupp19/dcarrol2/v05_V5r1/darwin3/run3/data.ptracers',
                    help='parent data.ptracers, used to locate DOC BY NAME')
    ap.add_argument('--frac-rdoc', type=float, default=DERIVED['rDOC']['frac'])
    ap.add_argument('--frac-cdom', type=float, default=DERIVED['CDOM']['frac'])
    ap.add_argument('--clip-min-cdom', type=float, default=None,
                    help='clip CDOM at this minimum (e.g. 0). OFF by default '
                         'to match the UNCLIPPED pickup_CDOM initial '
                         'condition -- if you set it here, rebuild the IC '
                         'with gen_derived_tracer_pickup.py --clip-min too')
    ap.add_argument('--clip-min-rdoc', type=float, default=None,
                    help='clip rDOC at this minimum. OFF by default, matching '
                         'the unclipped pickup_rDOC initial condition')
    ap.add_argument('--boundaries', nargs='+', default=list(BOUNDARIES))
    ap.add_argument('--force', action='store_true')
    # Edge lengths / record count default to pan_AO (square 1680 edges, 8784
    # hourly records). Other children (e.g. WestAO: E/W 680, N/S 860, 4344
    # records) must pass their own.
    ap.add_argument('--nx-ew', type=int, default=NX, help='points on east/west edges')
    ap.add_argument('--nx-ns', type=int, default=NX, help='points on north/south edges')
    ap.add_argument('--nrec', type=int, default=NREC)
    args = ap.parse_args()

    out_dir = os.path.join(args.config_dir, 'forcings/OBCS')
    if not os.path.isdir(out_dir):
        sys.exit(f"FATAL: no such directory {out_dir}")

    slot, names = parent_slot(args.parent_ptracers, SOURCE_SPECIES)
    print(f"parent {SOURCE_SPECIES} resolved BY NAME to parent ptracer slot "
          f"{slot} -> TRAC{slot:02d}_<bnd>.bin")
    print(f"  (from {args.parent_ptracers}, {len(names)} tracers)")

    spec = {
        'rDOC': dict(frac=args.frac_rdoc, clip_min=args.clip_min_rdoc),
        'CDOM': dict(frac=args.frac_cdom, clip_min=args.clip_min_cdom),
    }

    print(f"\nfractions: " + ", ".join(
        f"{k}={v['frac']}" + ("" if v['clip_min'] is None
                              else f" (clip>={v['clip_min']})")
        for k, v in spec.items()))
    print("NOTE: these must match gen_derived_tracer_pickup.py's --frac "
          "values used for the ICs.\n")

    print(f"{'output':22s} {'src nonzero':>12s} {'out nonzero':>12s} "
          f"{'clipped':>9s} {'negative':>10s} {'min':>13s} {'max':>13s}")
    print('-' * 100)

    footprint = {}
    for nm, s in spec.items():
        for bnd in args.boundaries:
            src = os.path.join(out_dir, f'TRAC{slot:02d}_{bnd}.bin')
            dst = os.path.join(out_dir, f'{nm}_{bnd}.bin')
            if not os.path.exists(src):
                sys.exit(f"FATAL: missing source {src}")
            nx = args.nx_ew if bnd in ('east', 'west') else args.nx_ns
            r = derive(src, dst, s['frac'], s['clip_min'], args.force,
                       nx=nx, nrec=args.nrec)
            print(f"{nm + '_' + bnd + '.bin':22s} {r['nz_src']:12d} "
                  f"{r['nz_out']:12d} {r['nclip']:9d} {r['nneg']:10d} "
                  f"{r['vmin']:13.6g} {r['vmax']:13.6g}")
            footprint.setdefault(bnd, []).append((nm, r))

    # A multiplication by a nonzero constant cannot move a land/ocean
    # boundary, so with clipping off the derived footprint MUST equal the
    # source's exactly. Asserted rather than eyeballed -- this is the cheap
    # check that catches a truncated or partially-written output.
    print()
    bad = 0
    for bnd, rows in sorted(footprint.items()):
        for nm, r in rows:
            clipped = any(s['clip_min'] is not None
                          for k, s in spec.items() if k == nm)
            if not clipped and r['nz_out'] != r['nz_src']:
                print(f"FAIL: {nm}_{bnd}.bin footprint {r['nz_out']} != "
                      f"source {r['nz_src']} with clipping OFF")
                bad += 1
    if bad:
        sys.exit(f"GEN_DERIVED_TRACER_OBCS_FAIL: {bad} footprint mismatch(es)")
    print("footprint check: every unclipped output matches its source "
          "nonzero count exactly")

    print("\nGEN_DERIVED_TRACER_OBCS_DONE")


if __name__ == '__main__':
    main()
