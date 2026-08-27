"""Real Part-5 conversion run, against the finished ECCOv4r6 forward-run
archive (see readme.txt Part 5) -- not the synthetic self-test in
mds_to_offline.py.

Forward run confirmed complete: NORMAL END, exit 0, 12418 daily
records in diags/archive_theta/ at a clean 24-iteration stride, first at
iteration 24, last at 298032 (= nIter0(1) + nTimeSteps(298031)).

Run this on Pleiades via the accompanying PBS job (job_convert_archive.pbs),
not interactively on pfe -- it reads/writes on the order of 2.4 TB across
the 5 3-D velocity/tracer fields, 3 3-D GM diffusivity fields, 1 3-D GGL90
diffusivity field, and 2 2-D fields (SFLUX, SIarea), all at 12418 records.
"""
import os
import sys

sys.path.insert(0, os.path.dirname(__file__))
from mds_to_offline import convert_archive_field, find_dump_iterations

MDS_DIR = os.environ.get("V4R6_RUN_DIR", "./MITgcm/ECCOV4/release6")
ARCHIVE_ROOT = os.environ.get("ARCHIVE_ROOT", "./archive")

# offline_ifprd = offlineForcingPeriod / deltaToffline = 86400 / 3600 = 24
# (see input_darwin_offline/data.off) -- must match for offline_fields_load.F
# to find each record at the iteration number it expects.
OFFLINE_IFPRD = 24
EXPECTED_RECORDS = 12418

FIELDS = [
    ("run/diags/archive_theta/archive_theta", "DDtheta"),
    ("run/diags/archive_salt/archive_salt", "DDsalt"),
    ("run/diags/archive_uvel/archive_uvel", "DDuvel"),
    ("run/diags/archive_vvel/archive_vvel", "DDvvel"),
    ("run/diags/archive_wvel/archive_wvel", "DDwvel"),
    ("run/diags/archive_gmkwx/archive_gmkwx", "DGMkwx"),
    ("run/diags/archive_gmkwy/archive_gmkwy", "DGMkwy"),
    ("run/diags/archive_gmkwz/archive_gmkwz", "DGMkwz"),
    ("run/diags/archive_ggl90kr/archive_ggl90kr", "DGGL90kr"),
    ("run/diags/archive_sflux/archive_sflux", "DFOSflux"),
    ("run/diags/archive_siarea/archive_siarea", "DSIarea"),
]

# DSIarea is consumed differently from the other 10 fields: it's read via
# data.darwin's icefile / data.radtrans's RT_icefile (both point at the
# same path), which go through pkg/exf's EXF_INTERP_READ -- a
# direct-access reader that opens ONE file and seeks to record `irecord`
# within it (same convention as ironFile/pCO2File/OASIM's fields). That's
# a different file format from the other 10 fields, which are read by
# pkg/offline's own MDS_READ_FIELD (one file PER record, iteration-number
# suffix) -- the format convert_archive_field() above produces. Found on
# the first-ever run attempt: EXF_INTERP_READ crashed with "File does not
# exist" on the bare path (no suffix), since only per-record files existed.
FIELDS_NEEDING_SINGLE_FILE = ["DSIarea"]


def concatenate_daily_records(out_dir, out_prefix, n_records, offline_ifprd):
    """Build the single multi-record file EXF_INTERP_READ expects
    (bare <out_dir>/<out_prefix>, no suffix) by concatenating the
    per-record files convert_archive_field() already wrote, in order.
    Cheap: those files are already correctly formatted and ordered, so
    this is a pure byte concatenation, not a re-read of the original MDS
    dumps.
    """
    combined_path = os.path.join(out_dir, out_prefix)
    with open(combined_path, "wb") as out_f:
        for record_index in range(n_records):
            suffix = record_index * offline_ifprd
            record_path = os.path.join(out_dir, f"{out_prefix}.{suffix:010d}")
            with open(record_path, "rb") as in_f:
                out_f.write(in_f.read())
    return combined_path


def main():
    for mds_prefix, out_prefix in FIELDS:
        n_dumps = len(find_dump_iterations(MDS_DIR, mds_prefix))
        if n_dumps != EXPECTED_RECORDS:
            raise ValueError(
                f"{out_prefix}: found {n_dumps} dumps, expected "
                f"{EXPECTED_RECORDS} (mismatched vs. archive_theta's own "
                "count -- one field's diagnostic stream is missing records)"
            )

        out_dir = os.path.join(ARCHIVE_ROOT, out_prefix)
        os.makedirs(out_dir, exist_ok=True)

        print(f"converting {out_prefix}: {n_dumps} records -> {out_dir}", flush=True)
        written = convert_archive_field(
            mds_dir=MDS_DIR, mds_prefix=mds_prefix,
            out_dir=out_dir, out_prefix=out_prefix,
            offline_ifprd=OFFLINE_IFPRD, offline_iter0=0, prec=32,
        )
        assert len(written) == EXPECTED_RECORDS
        print(f"  done: {len(written)} records written", flush=True)

        if out_prefix in FIELDS_NEEDING_SINGLE_FILE:
            combined_path = concatenate_daily_records(
                out_dir, out_prefix, EXPECTED_RECORDS, OFFLINE_IFPRD,
            )
            print(f"  combined into single file for EXF_INTERP_READ: "
                  f"{combined_path}", flush=True)

    print("all 11 fields converted OK", flush=True)


if __name__ == "__main__":
    main()
