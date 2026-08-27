"""Convert a V4r6 forward run's own MDS diagnostic dumps into the
pkg/offline per-record binary naming convention.

Unlike regrid_to_llc90.py, this needs NO horizontal/vertical regridding --
the diagnostic output already comes straight off the LLC90 grid, in the
same "compact" byte layout pkg/offline's READ_REC_3D_RS/READ_REC_XY_RS
expect. The only real work is renaming: MITgcm's diagnostics package names
each dump by ITERATION NUMBER (<prefix>.<iter>.data), but
offline_fields_load.F expects <prefix>.<record_index*Ifprd + offlineIter0>
(see OFFLINE.h / offline_fields_load.F in code_offline_ggl90/), and the
data may need a precision cast to match offlineLoadPrec=32.

Data doesn't exist yet (blocked on the actual V4r6 forward run happening
on Pleiades -- see the plan file), so this is self-tested below with
synthetic MDS-like files rather than real ones. The renaming/precision
logic doesn't change once real files exist; only MDS_DIR/PREFIX do.
"""
import glob
import os
import re

import numpy as np
import MITgcmutils as mu


def find_dump_iterations(mds_dir, mds_prefix):
    """Sorted list of iteration numbers for which <mds_prefix>.<iter>.data
    exists in mds_dir, and a sanity check that dumps are evenly spaced
    (irregular spacing would mean either a run restart/gap or the wrong
    prefix was picked -- worth erroring on rather than silently producing
    a mistimed offline archive).
    """
    pattern = os.path.join(mds_dir, f"{mds_prefix}.*.data")
    iters = sorted(
        int(re.search(r"\.(\d{10})\.data$", f).group(1))
        for f in glob.glob(pattern)
    )
    if len(iters) < 2:
        return iters
    stride = iters[1] - iters[0]
    bad = [
        (a, b) for a, b in zip(iters, iters[1:]) if b - a != stride
    ]
    if bad:
        raise ValueError(
            f"{mds_prefix}: dump iterations are not evenly spaced "
            f"(expected stride {stride}); first mismatch at {bad[0]}. "
            "This usually means a run restart/gap, or the wrong prefix."
        )
    return iters


def convert_archive_field(
    mds_dir, mds_prefix, out_dir, out_prefix, offline_ifprd,
    offline_iter0=0, prec=32,
):
    """Read every <mds_prefix>.<iter>.data dump in mds_dir (sorted by
    iteration -> record index 0,1,2,...) and write each out as
    <out_dir>/<out_prefix>.<record_index*offline_ifprd + offline_iter0>,
    at the given precision, big-endian (MITgcm's own convention).

    Returns the list of (record_index, mds_iter, out_path) written, so a
    caller can sanity-check e.g. that all fields for one day map to a
    consistent set of record indices before trusting the archive.
    """
    dtype = ">f4" if prec == 32 else ">f8"
    iters = find_dump_iterations(mds_dir, mds_prefix)
    written = []
    for record_index, it in enumerate(iters):
        field = mu.rdmds(os.path.join(mds_dir, mds_prefix), it)
        suffix = record_index * offline_ifprd + offline_iter0
        out_path = os.path.join(out_dir, f"{out_prefix}.{suffix:010d}")
        field.astype(dtype).tofile(out_path)
        written.append((record_index, it, out_path))
    return written


def _write_synthetic_mds(mds_dir, prefix, iters, shape):
    """Write fake <prefix>.<iter>.data/.meta pairs, mimicking MITgcm's real
    MDS layout, for self-testing without real V4r6 data.

    shape is the numpy array shape rdmds would return: (1170,90) for a 2D
    LLC90-compact field, or (Nr,1170,90) for a 3D one -- matching the real
    XC.meta/hFacC.meta convention (checked against the actual local files):
    dimList lists one (size,1,size) triple per dimension in FILE order
    (x,y,[z]), i.e. reversed from numpy's (...,z,y,x) shape order.
    """
    rng = np.random.default_rng(0)
    file_dims = tuple(reversed(shape))  # (nx,ny[,nz]) file order
    ndims = len(shape)
    dim_list = "\n".join(f"   {n},    1, {n}," for n in file_dims)
    for it in iters:
        field = rng.normal(size=shape).astype(">f8")
        data_path = os.path.join(mds_dir, f"{prefix}.{it:010d}.data")
        field.tofile(data_path)
        with open(os.path.join(mds_dir, f"{prefix}.{it:010d}.meta"), "w") as f:
            f.write(f" nDims = [ {ndims} ];\n")
            f.write(f" dimList = [\n{dim_list}\n ];\n")
            f.write(" dataprec = [ 'float64' ];\n")
            f.write(" nrecords = [ 1 ];\n")


def mds_to_offline_selftest():
    """Self-test with synthetic MDS files -- validates the renaming
    (record-index derivation, filename format) and precision-cast logic
    ahead of the real V4r6 archive existing."""
    import shutil
    import tempfile

    mds_dir = tempfile.mkdtemp(prefix="mds_selftest_")
    out_dir = tempfile.mkdtemp(prefix="offline_selftest_")
    try:
        shape = (1170, 90)  # LLC90-compact 2D field (matches SFLUX: 1 level)
        iters = [24, 48, 72, 96]  # daily dumps, hourly clock -> stride 24
        _write_synthetic_mds(mds_dir, "archive_sflux", iters, shape)

        written = convert_archive_field(
            mds_dir, "archive_sflux", out_dir, "SFLUX",
            offline_ifprd=24, offline_iter0=0, prec=32,
        )
        assert [w[0] for w in written] == [0, 1, 2, 3]
        assert [w[1] for w in written] == iters
        expected_names = [
            "SFLUX.0000000000", "SFLUX.0000000024",
            "SFLUX.0000000048", "SFLUX.0000000072",
        ]
        got_names = sorted(os.path.basename(p) for _, _, p in written)
        assert got_names == expected_names, got_names
        print("filename/record-index mapping OK:", got_names)

        # precision cast: written as float32 even though source was float64
        raw = np.fromfile(written[0][2], dtype=">f4")
        assert raw.size == shape[0] * shape[1]
        print("precision-cast OK: wrote", raw.size, "float32 values")

        # uneven-spacing guard
        _write_synthetic_mds(mds_dir, "archive_bad", [24, 48, 96], shape)
        try:
            find_dump_iterations(mds_dir, "archive_bad")
            raise AssertionError("expected uneven-spacing ValueError")
        except ValueError as e:
            print("uneven-spacing guard OK:", e)
    finally:
        shutil.rmtree(mds_dir)
        shutil.rmtree(out_dir)


if __name__ == "__main__":
    mds_to_offline_selftest()
