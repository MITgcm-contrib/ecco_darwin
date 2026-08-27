"""Build real MDS .data/.meta pairs for the 10 BGC river runoff fields
(data.darwin's *runofffile parameters), found missing on the first-ever
run attempt of the offline Darwin build (readme.txt Part 7).

Root cause: unlike ironFile/pCO2File (which have a real _lon0/_lat0/
_nlon/_nlat interpolation block in data.darwin, routing them through
pkg/exf's EXF_INTERP_READ -- a direct-access reader needing only a raw
concatenated binary, no .meta), the runoff fields have NO interpolation
block (deliberately -- the source data is already native LLC90, per
readme.txt's own note). darwin_exf_load.F passes each
field's _lon0/_nlon/etc. (all zero/unset for runoff) to a shared loader
that, lacking those, falls back to the plain MDS_READ_FIELD path -- which
requires a genuine <prefix>.data/<prefix>.meta pair, not just a raw
concatenated file. EXF_INTERP_READ's own interpolation mechanism
couldn't be used instead even if we wanted to: it assumes a regular
lat-lon source grid, and LLC90-compact isn't one -- a real .meta pair is
the only correct fix for a genuinely native-grid periodic field like this.

Source data: /nobackup/rsavelli/LOAC/ECCO_V4r5/bgc_runoff/, one raw
(headerless) binary file per field per calendar year
(<prefix>_ECCO_V4r5_<YYYY>), each a sequence of daily 2-D LLC90-compact
(1170x90 compact, i.e. 90x1170 file-order) float32 records -- confirmed
directly from file sizes (365/366 records x 421200 bytes/record, matching
real/leap years exactly).
"""
import glob
import os
import re

MDS_DIR = "/nobackup/rsavelli/LOAC/ECCO_V4r5/bgc_runoff"
OUT_DIR = os.environ.get("OFFLINE_RUN_DIR", ".")

RECORD_BYTES = 1170 * 90 * 4  # one 2-D LLC90-compact float32 record

# (source prefix, data.darwin's *runofffile value -- confirmed against
# the real filenames on Pleiades; note POPrunofffile/PONrunofffile use
# the short "PP"/"PN" names, not "POP"/"PON").
FIELDS = [
    "DOC_ECCO_V4r5",
    "DON_ECCO_V4r5",
    "DOP_ECCO_V4r5",
    "DIN_ECCO_V4r5",
    "DIP_ECCO_V4r5",
    "DSi_ECCO_V4r5",
    "POC_ECCO_V4r5",
    "PP_ECCO_V4r5",
    "PN_ECCO_V4r5",
    "DIC_ECCO_V4r5",
]


def find_yearly_files(prefix):
    """Sorted (by year, not lexicographically) list of yearly source
    files for one field prefix."""
    pattern = os.path.join(MDS_DIR, f"{prefix}_*")
    files = glob.glob(pattern)
    with_year = []
    for f in files:
        m = re.search(r"_(\d{4})$", f)
        if not m:
            raise ValueError(f"{f}: doesn't end in a 4-digit year")
        with_year.append((int(m.group(1)), f))
    with_year.sort()
    return with_year


def build_one(prefix):
    with_year = find_yearly_files(prefix)
    if not with_year:
        raise ValueError(f"{prefix}: no yearly source files found")
    years = [y for y, _ in with_year]
    if years != list(range(years[0], years[-1] + 1)):
        raise ValueError(f"{prefix}: year gap in {years}")

    out_data = os.path.join(OUT_DIR, f"{prefix}.data")
    total_records = 0
    with open(out_data, "wb") as out_f:
        for year, path in with_year:
            size = os.path.getsize(path)
            if size % RECORD_BYTES != 0:
                raise ValueError(
                    f"{path}: size {size} not a multiple of "
                    f"{RECORD_BYTES} (one 2-D LLC90 record)"
                )
            n_records = size // RECORD_BYTES
            expected = 366 if _is_leap(year) else 365
            if n_records != expected:
                raise ValueError(
                    f"{path}: {n_records} records, expected {expected} "
                    f"for {'leap' if _is_leap(year) else 'non-leap'} "
                    f"year {year}"
                )
            with open(path, "rb") as in_f:
                out_f.write(in_f.read())
            total_records += n_records

    out_meta = os.path.join(OUT_DIR, f"{prefix}.meta")
    with open(out_meta, "w") as f:
        f.write(" nDims = [   2 ];\n")
        f.write(" dimList = [\n")
        f.write("    90,    1,   90,\n")
        f.write("  1170,    1, 1170\n")
        f.write(" ];\n")
        f.write(" dataprec = [ 'float32' ];\n")
        f.write(f" nrecords = [ {total_records} ];\n")

    return out_data, out_meta, total_records, years[0], years[-1]


def _is_leap(year):
    return year % 4 == 0 and (year % 100 != 0 or year % 400 == 0)


def main():
    for prefix in FIELDS:
        out_data, out_meta, n, y0, y1 = build_one(prefix)
        print(f"{prefix}: {n} daily records ({y0}-{y1}) -> "
              f"{out_data} + {out_meta}", flush=True)
    print("all 10 runoff fields built OK", flush=True)


if __name__ == "__main__":
    main()
