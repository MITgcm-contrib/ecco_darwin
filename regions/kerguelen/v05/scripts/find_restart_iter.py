#!/usr/bin/env python3
"""Resolve the furthest verified-complete MITgcm checkpoint in a run
directory and, if it's a rolling (ckptA/ckptB) checkpoint rather than a
numbered (pChkptFreq) one, rename all 5 filesets to the numbered
convention so MITgcm's restart reader (which literally string-formats
nIter0 into 'pickup.<iter>.data', never reads .meta timeStepNumber --
confirmed against read_pickup.F/ptracers_read_pickup.F/etc in the
darwin3 source) can pick it up via nIter0 alone.

Prints the resolved restart iteration as the last line of stdout, or
prints nothing and exits 1 if no verified-complete checkpoint exists.
"""
import argparse
import os
import re
import sys

# fileset base name -> expected .data byte size (from packages_write_pickup.F
# field inventory, cross-checked against real files on disk, 2026-08-09).
FILESETS = {
    "pickup": 6_821_418_240,
    "pickup_ptracers": 26_323_315_200,
    "pickup_seaice": 56_609_280,
    "pickup_ggl90": 849_139_200,
    "pickup_darwin": 849_139_200,
}


def read_timestep_number(meta_path):
    try:
        with open(meta_path) as f:
            text = f.read()
    except OSError:
        return None
    m = re.search(r"timeStepNumber\s*=\s*\[\s*(\d+)\s*\]", text)
    return int(m.group(1)) if m else None


def candidate_ok(run_dir, suffix, expect_iter=None):
    """Check all 5 filesets exist under `suffix`, agree on timeStepNumber,
    and have the exact expected byte size. Returns the shared iter or None."""
    iters = set()
    for base, expect_size in FILESETS.items():
        data_path = os.path.join(run_dir, f"{base}.{suffix}.data")
        meta_path = os.path.join(run_dir, f"{base}.{suffix}.meta")
        if not (os.path.isfile(data_path) and os.path.isfile(meta_path)):
            return None
        if os.path.getsize(data_path) != expect_size:
            return None
        it = read_timestep_number(meta_path)
        if it is None:
            return None
        iters.add(it)
    if len(iters) != 1:
        return None
    it = iters.pop()
    if expect_iter is not None and it != expect_iter:
        return None
    return it


def rename_to_numbered(run_dir, suffix, iter_num):
    numbered = f"{iter_num:010d}"
    for base in FILESETS:
        for ext in ("data", "meta"):
            src = os.path.join(run_dir, f"{base}.{suffix}.{ext}")
            dst = os.path.join(run_dir, f"{base}.{numbered}.{ext}")
            os.rename(src, dst)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--run", required=True, help="run directory")
    args = ap.parse_args()
    run_dir = args.run

    candidates = []  # (iter, suffix)

    # numbered (permanent, pChkptFreq) checkpoints
    for fn in os.listdir(run_dir):
        m = re.match(r"^pickup\.(\d{10})\.data$", fn)
        if m:
            it = int(m.group(1))
            candidates.append((it, m.group(1)))

    # rolling (chkptFreq) checkpoints
    for suffix in ("ckptA", "ckptB"):
        meta_path = os.path.join(run_dir, f"pickup.{suffix}.meta")
        it = read_timestep_number(meta_path)
        if it is not None:
            candidates.append((it, suffix))

    candidates.sort(key=lambda c: c[0], reverse=True)

    for it, suffix in candidates:
        verified = candidate_ok(run_dir, suffix, expect_iter=it)
        if verified is None:
            print(f"REJECTED suffix={suffix} claimed_iter={it} "
                  f"(incomplete/inconsistent/wrong-size fileset)",
                  file=sys.stderr)
            continue
        if suffix in ("ckptA", "ckptB"):
            print(f"RENAMING suffix={suffix} iter={it} -> numbered convention",
                  file=sys.stderr)
            rename_to_numbered(run_dir, suffix, it)
        print(it)
        return 0

    print("NO_VALID_CHECKPOINT_FOUND", file=sys.stderr)
    return 1


if __name__ == "__main__":
    sys.exit(main())
