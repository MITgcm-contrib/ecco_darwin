"""Build combined per-field (and, for the 3 spectral aerosol fields,
per-waveband) OASIM forcing files, found missing on the first-ever run
attempt of the offline Darwin build (readme.txt Part 7).

Root cause: data.oasim's field paths (e.g. aertauFile=
'oasim/C61extrap/R4/taua') are read via pkg/exf's EXF_INTERP_READ --
a direct-access reader expecting ONE file with all records (same
convention as apCO2/DSIarea, see build_runoff_mds.py's docstring for the
general pattern). But the real source at /nobackup/ojahn/forcing/oasim/
is organized as one raw (headerless) binary file PER CALENDAR YEAR (12
monthly records each) for the 8 non-spectral fields, and PER WAVEBAND
PER CALENDAR YEAR for the 3 spectral aerosol fields (aertau/asymp/ssalb
-- confirmed via OASIM_SIZE.h: nltaer=33 wavebands). Crashed with
"EXF_INTERP_READ: filename: oasim/C61extrap/R4/taua001 / File does not
exist" -- the code appends the 3-digit waveband index but not a year,
so it wants ALL years already concatenated into ONE file per waveband.

Since the real oasim/ directory is read-only (owned by another user) and
was previously symlinked wholesale, this also restructures our own
oasim/ into a real local directory: the 3 .dat optical lookup tables
(oasim_atmoFile/waterFile/slingoFile) are symlinked individually (no
concatenation needed, not time-varying), but the periodic forcing fields
are real files we build here.

Grid: confirmed 360x180 (global 1x1-degree, matching data.oasim's own
_nlon=360/_nlat=180 interpolation blocks) x 12 monthly records/year for
every field, spectral or not (file sizes divide out exactly).
"""
import glob
import os
import re

SRC_DIR = os.environ.get("OASIM_SRC_DIR", "/nobackup/ojahn/forcing/oasim/C61extrap/R4")
BCS_SRC_DIR = os.environ.get("OASIM_BCS_DIR", "/nobackup/ojahn/forcing/oasim/bcs")
OUT_ROOT = os.environ.get("OFFLINE_RUN_DIR", ".") + "/oasim"
OUT_DIR = os.path.join(OUT_ROOT, "C61extrap", "R4")
BCS_OUT_DIR = os.path.join(OUT_ROOT, "bcs")

RECORD_BYTES = 360 * 180 * 4  # one global 1x1-degree float32 monthly record
RECORDS_PER_YEAR = 12
YEAR_BYTES = RECORD_BYTES * RECORDS_PER_YEAR

# (data.oasim field, real source prefix, is_spectral)
FIELDS = [
    ("cldcovFile", "ccov", False),
    ("cldlwpFile", "rlwp", False),
    ("cldreFile", "cdre", False),
    ("presFile", "slp", False),
    ("OAwindFile", "wsm", False),
    ("relhumFile", "rh", False),
    ("ozoneFile", "oz", False),
    ("wvFile", "wv", False),
    ("aertauFile", "taua", True),
    ("asympFile", "asymp", True),
    ("ssalbFile", "ssalb", True),
]
N_WAVEBANDS = 33  # OASIM_SIZE.h: nltaer=33


def find_yearly_files(prefix):
    pattern = os.path.join(SRC_DIR, f"{prefix}_*")
    files = glob.glob(pattern)
    with_year = []
    for f in files:
        m = re.search(r"_(\d{4})$", f)
        if not m:
            raise ValueError(f"{f}: doesn't end in a 4-digit year")
        with_year.append((int(m.group(1)), f))
    with_year.sort()
    return with_year


def concatenate_years(prefix, out_path):
    with_year = find_yearly_files(prefix)
    if not with_year:
        raise ValueError(f"{prefix}: no yearly source files found")
    years = [y for y, _ in with_year]
    if years != list(range(years[0], years[-1] + 1)):
        raise ValueError(f"{prefix}: year gap in {years}")

    with open(out_path, "wb") as out_f:
        for year, path in with_year:
            size = os.path.getsize(path)
            if size != YEAR_BYTES:
                raise ValueError(
                    f"{path}: size {size}, expected {YEAR_BYTES} "
                    f"(12 monthly {RECORD_BYTES}-byte records)"
                )
            with open(path, "rb") as in_f:
                out_f.write(in_f.read())
    return years[0], years[-1], len(with_year)


def main():
    os.makedirs(OUT_DIR, exist_ok=True)
    os.makedirs(BCS_OUT_DIR, exist_ok=True)

    for fname in ["atmo25b.dat", "abw25b.dat", "slingo.dat"]:
        src = os.path.join(BCS_SRC_DIR, fname)
        dst = os.path.join(BCS_OUT_DIR, fname)
        if not os.path.lexists(dst):
            os.symlink(src, dst)
        print(f"linked {fname}", flush=True)

    for field, prefix, is_spectral in FIELDS:
        if not is_spectral:
            out_path = os.path.join(OUT_DIR, prefix)
            y0, y1, n = concatenate_years(prefix, out_path)
            print(f"{field} ({prefix}): {n} years ({y0}-{y1}) -> {out_path}",
                  flush=True)
        else:
            for band in range(1, N_WAVEBANDS + 1):
                band_prefix = f"{prefix}{band:03d}"
                out_path = os.path.join(OUT_DIR, band_prefix)
                y0, y1, n = concatenate_years(band_prefix, out_path)
                print(f"{field} band {band:03d} ({band_prefix}): "
                      f"{n} years ({y0}-{y1}) -> {out_path}", flush=True)

    print("all OASIM fields built OK", flush=True)


if __name__ == "__main__":
    main()
