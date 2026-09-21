"""
Build INTERANNUAL discharge and riverine DOC (TOC cub) forcings for NS-RAD's four North
Slope rivers, from Michael's PWBM (Pan-Arctic Water Balance Model) output.

Unlike every other forcing in forcing/ (a single 365-day climatology that repeats every
year -- see config.repeatYear), these are GENUINE MULTI-YEAR daily series, 1980-2023,
with real year-to-year variability. They drive the *_interannual site variants
(sites/<name>_interannual.py) via CGEM_SITE=<name>_interannual -- see CLAUDE.md ->
"Interannual forcing" and tools/run_interannual.sh. The four rivers' regular sites
(colville.py etc.) and runs/definitive/ are completely unaffected: this script only
writes forcing/*_interannual_*.csv, it does not touch colville_river_discharge_2022_m3sec.csv
or any existing BOUNDARIES/TOC constant.

Source (not in this repo, read directly at build time) -- Michael's per-river DAILY
output, one row per calendar day 1980-01-01 to 2023-12-31, replacing the earlier
monthly-aggregate PWBM extraction this script used:
    /Users/rsavelli/Documents/FORTE/Michael/Colville/derived_daily_annual_csv/Colville_daily_1980_2023.csv
    /Users/rsavelli/Documents/FORTE/Michael/Colville/derived_daily_annual_csv/Canning_daily_1980_2023_fromMichael.csv
    /Users/rsavelli/Documents/FORTE/Michael/Kuparuk_daily_1980_2023.csv
    /Users/rsavelli/Documents/FORTE/Michael/Sagavanirktok_daily_1980_2023.csv
Each has columns `date,freshwater_m3_d,DOC_g_d[,...per-sub-basin breakdown]`. Colville and
Sagavanirktok also carry a per-sub-basin/per-branch breakdown (mainstem/basin134/Kupigruak;
EastBranch/WestBranch/mouth) whose sum equals the total `freshwater_m3_d`/`DOC_g_d` columns
(checked directly against the files) -- only the totals are used here, same "one number per
river" granularity the model consumes.

Unit conversion
----------------
`freshwater_m3_d` is already a daily volume [m3/day] and `DOC_g_d` a daily DOC export
[g/day], both PWBM whole-basin sums, so no cell-area/days-in-month bookkeeping is needed
(that was only required to convert the OLD monthly runoff-depth/DOC-mass product; this
data is already in daily flux units):

    Q(day)     [m3/s]  = freshwater_m3_d / 86400
    C_DOC(day) [mg/L]  = DOC_g_d / freshwater_m3_d      (g/m3 == mg/L numerically, so the
                                                          m3<->L conversion cancels exactly
                                                          the same way the old 1 km2/cell
                                                          factor did)
    TOC_cub    [mmol/m3] = C_DOC[mg/L] * 83.3     (same DOC->TOC factor every sites/*.py
                                                    site already uses: 1000/12.011)

2022 annual-mean sanity check (this script's own printed summary should reproduce these,
and does, to within ~0.5 m3/s -- confirms this daily product and the retired monthly one
are the same underlying PWBM run, just at different output resolution): Colville ~394,
Kuparuk ~72, Sagavanirktok ~198, Canning ~68 m3/s.

DOC-only, not POC -- unchanged from the retired monthly version. NS-RAD's `TOC` species is
purely dissolved labile DOC in this model's actual kinetics (Monod oxidation/denitrification,
env=1 transport identical to salinity, zero coupling to SPM/settling/burial anywhere in
biogeo_module.py or sed_module.py -- see CLAUDE.md). No POC/TSS product exists for these
four rivers, so PWBM DOC is fed directly into TOC_cub with no particulate addition.

Winter gap-fill -- same technique as before, now at daily instead of monthly resolution.
Zero-flow days (frozen, `freshwater_m3_d == 0`) leave `C_DOC` undefined (0/0). Discharge is
genuinely ~0 those days so the exact concentration barely matters physically, but it must
still be finite for the model to read: each river's open-water DOC days are linearly
interpolated (np.interp, the same gap-fill technique tools/build_river_temp.py uses) across
each winter's zero-flow gap.

Calendar. The source is real calendar days including Feb 29 in leap years (16071 rows,
1980-2023). NS-RAD's own daily forcings all treat every year as exactly 365 days (see
config.repeatYear's annual wrap in file_module.exfread, and tools/fetch_discharge.py) so
that a genuinely multi-year series like this one stays in phase, year over year, with the
single repeating 2022-typical seasonal cycle every OTHER forcing (met, tides, temperature)
still uses. Feb 29 is therefore dropped -- its one row per leap year is discarded, not
interpolated -- giving a compressed 44*365 = 16060-day output series, exactly as the retired
monthly-interpolation version produced.

Output: forcing/<site>_river_discharge_interannual_1980-2023_m3sec.csv (CRLF, 2 dp, to
match every other *_river_discharge_*.csv) and
forcing/<site>_toc_interannual_1980-2023_mmolC_m3.csv (LF, 6 dp, to match the other
BOUNDARY_FORCING-consumed series in forcing/idealized_*_cub_river.csv). Both: one value
per day, no header, no trailing newline -- file_module.exfread's genfromtxt reader.

Usage:  python tools/build_interannual_forcings.py [--plot]
"""
import csv
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent.parent
FORCINGS = ROOT / "forcing"

# NS-RAD site key -> Michael's per-river daily CSV
SRC = {
    "colville": Path("/Users/rsavelli/Documents/FORTE/Michael/Colville/derived_daily_annual_csv/Colville_daily_1980_2023.csv"),
    "kuparuk": Path("/Users/rsavelli/Documents/FORTE/Michael/Kuparuk_daily_1980_2023.csv"),
    "sagavanirktok": Path("/Users/rsavelli/Documents/FORTE/Michael/Sagavanirktok_daily_1980_2023.csv"),
    "canning": Path("/Users/rsavelli/Documents/FORTE/Michael/Colville/derived_daily_annual_csv/Canning_daily_1980_2023_fromMichael.csv"),
}

YEAR_START, YEAR_END = 1980, 2023  # inclusive, matches the full record on disk
N_DAYS_COMPRESSED = (YEAR_END - YEAR_START + 1) * 365  # = 16060, Feb 29 dropped

MGL_TO_MMOLM3 = 83.3  # DOC mg/L -> TOC mmol C/m3 = 1000/12.011, rounded as every
                       # sites/<name>.py DOC->TOC conversion already is


def _read_daily(path):
    """Return (q_m3d, doc_g_d) as float arrays, one value per row of `path`, in file
    order (already chronological, 1980-01-01..2023-12-31). Only the total
    `freshwater_m3_d`/`DOC_g_d` columns are read; a per-sub-basin breakdown some files
    carry is not (see module docstring)."""
    q, doc, dates = [], [], []
    with open(path, newline="") as f:
        for row in csv.DictReader(f):
            dates.append(row["date"])
            q.append(float(row["freshwater_m3_d"]))
            doc.append(float(row["DOC_g_d"]))
    return dates, np.array(q), np.array(doc)


def _compress_calendar(dates, values):
    """Drop the Feb-29 row of every leap year (see module docstring's Calendar section).
    `dates` are 'YYYY-MM-DD' strings in chronological order."""
    keep = np.array([not d.endswith("02-29") for d in dates])
    return values[keep]


def _write(path, values, decimals, newline):
    body = newline.join(f"{v:.{decimals}f}" for v in values)
    path.write_text(body, encoding="utf-8", newline="")


def main():
    plot = "--plot" in sys.argv
    FORCINGS.mkdir(exist_ok=True)

    if plot:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        fig, axes = plt.subplots(len(SRC), 2, figsize=(11, 9), sharex=True)

    for i, (site, path) in enumerate(SRC.items()):
        dates, q_m3d, doc_gd = _read_daily(path)

        q_daily = _compress_calendar(dates, q_m3d / 86400.0)          # m3/s
        doc_daily = _compress_calendar(dates, doc_gd)                  # g/d, for gap-fill below
        q_raw_daily = _compress_calendar(dates, q_m3d)                 # m3/d, same axis as doc_daily

        assert q_daily.size == N_DAYS_COMPRESSED, \
            f"{site}: expected {N_DAYS_COMPRESSED} compressed days, got {q_daily.size}"

        conc_daily = np.full(N_DAYS_COMPRESSED, np.nan)  # mg/L
        flowing = q_raw_daily > 0
        conc_daily[flowing] = doc_daily[flowing] / q_raw_daily[flowing]

        # Gap-fill undefined (zero-flow) days by linear interpolation between the
        # nearest defined open-water days on either side.
        ok = ~np.isnan(conc_daily)
        if not ok.all():
            conc_daily[~ok] = np.interp(np.flatnonzero(~ok), np.flatnonzero(ok), conc_daily[ok])
        toc_daily = conc_daily * MGL_TO_MMOLM3

        q_path = FORCINGS / f"{site}_river_discharge_interannual_{YEAR_START}-{YEAR_END}_m3sec.csv"
        toc_path = FORCINGS / f"{site}_toc_interannual_{YEAR_START}-{YEAR_END}_mmolC_m3.csv"
        _write(q_path, q_daily, 2, "\r\n")
        _write(toc_path, toc_daily, 6, "\n")

        # 2022 annual-mean check, printed against the values already verified against
        # the retired monthly product (see module docstring).
        idx2022 = slice((2022 - YEAR_START) * 365, (2022 - YEAR_START + 1) * 365)
        q2022_mean = q_daily[idx2022].mean()
        print(f"{site:14s} 2022 annual-mean Q = {q2022_mean:7.1f} m3/s   "
              f"TOC range {toc_daily.min():6.1f}-{toc_daily.max():6.1f} mmol/m3   "
              f"(DOC {conc_daily.min():5.2f}-{conc_daily.max():5.2f} mg/L)   "
              f"n_gapfilled_days={int((~ok).sum())}")

        if plot:
            years = YEAR_START + np.arange(N_DAYS_COMPRESSED) / 365.0
            axes[i, 0].plot(years, q_daily)
            axes[i, 0].set_ylabel(f"{site}\nQ [m3/s]")
            axes[i, 1].plot(years, toc_daily)
            axes[i, 1].set_ylabel("TOC [mmol/m3]")
    if plot:
        axes[0, 0].set_title("Daily discharge")
        axes[0, 1].set_title("Daily TOC (from daily DOC)")
        axes[-1, 0].set_xlabel("year")
        axes[-1, 1].set_xlabel("year")
        fig.tight_layout()
        out = FORCINGS.parent / "docs" / "interannual_forcings_preview.png"
        fig.savefig(out, dpi=130)
        print(f"wrote {out}")


if __name__ == "__main__":
    main()
