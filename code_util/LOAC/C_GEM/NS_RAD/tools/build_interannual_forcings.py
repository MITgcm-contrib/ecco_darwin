"""
Build INTERANNUAL discharge and riverine DOC (TOC cub) forcings for NS-RAD's four North
Slope rivers, from Michael's PWBM (Pan-Arctic Water Balance Model) monthly output.

Unlike every other forcing in forcing/ (a single 365-day climatology that repeats every
year -- see config.repeatYear), these are GENUINE MULTI-YEAR daily series, 1980-2023,
with real year-to-year variability. They drive the new *_interannual site variants
(sites/<name>_interannual.py) via CGEM_SITE=<name>_interannual -- see CLAUDE.md ->
"Interannual forcing" and tools/run_interannual.sh. The four rivers' regular sites
(colville.py etc.) and runs/definitive/ are completely unaffected: this script only adds
new files, it does not touch colville_river_discharge_2022_m3sec.csv or any existing
BOUNDARIES/TOC constant.

Source (not in this repo, read directly at build time):
    /Users/rsavelli/Documents/FORTE/Michael/river_timeseries_output/
        runoff_monthly_timeseries_by_river.csv   -- per-river monthly SUM of per-cell
                                                     PWBM total runoff depth [mm]
        DOC_monthly_timeseries_by_river.csv      -- per-river monthly SUM of per-cell
                                                     PWBM total DOC export [kg C]
        watershed_cell_counts.csv                -- n_cells per river watershed

Unit conversion (validated this session -- see CLAUDE.md for the full derivation and
cross-checks against USGS gauges / basin areas)
-----------------------------------------------
PWBM grid cells are 1 km^2: n_cells * 1 km^2 reproduces each river's known basin area
(e.g. Colville 59,647 cells matches the "~60,000 km2" already cited in sites/colville.py).
With that:

    Q(month)     [m3/s]     = runoff_sum[mm] * 1000 / (days_in_month * 86400)
    C_DOC(month) [mg/L]     = DOC_sum[kg] / runoff_sum[mm]        (the 1 km2 factor cancels)
    TOC_cub      [mmol/m3]  = C_DOC[mg/L] * 83.3     (same DOC->TOC factor every sites/*.py
                                                       site already uses: 1000/12.011)

2022 sanity check (this script's own printed summary should reproduce these): annual mean
discharge Colville ~395, Kuparuk ~72.6, Sagavanirktok ~197, Canning ~68 m3/s -- all higher
than or close to the existing gauge-based 2022 means (238.8 / 63.8 / 47.6 / 45.4 m3/s),
consistent with the documented gauge caveats (Colville/Sagavanirktok gauged well upstream
of the delta and known to understate mouth flow; Kuparuk near-tidewater and best-constrained
-- PWBM lands closest to ITS gauge too; Canning has no gauge at all today, only an ad-hoc
Hulahula-proxy reconstruction that PWBM's real whole-basin runoff replaces).

DOC-only, not POC. NS-RAD's `TOC` species is purely dissolved labile DOC in this model's
actual kinetics (Monod oxidation/denitrification, env=1 transport identical to salinity,
zero coupling to SPM/settling/burial anywhere in biogeo_module.py or sed_module.py -- see
CLAUDE.md). PWBM's own extraction on this machine is DOC-only; no POC/TSS product exists
for these four rivers anywhere. So this script feeds PWBM DOC directly into TOC_cub with
no particulate addition -- documented as a known limitation, the same way the model
already documents its missing prognostic sediment-OM pool (docs/arctic_biogeochemistry.md).

Winter gap-fill. Many winter months have runoff_sum == 0 (frozen, no PWBM signal), which
leaves C_DOC undefined (0/0). Since discharge is genuinely ~0 those months the exact
concentration barely matters physically (negligible advective mass flux), but it must
still be finite and smooth for daily interpolation: each river's open-water DOC values are
linearly interpolated (np.interp, the same gap-fill technique tools/build_river_temp.py
uses) across the zero-flow gap.

Calendar. Monthly values are placed at each REAL calendar month's midpoint (using true
days-in-month, including real Feb 29 in leap years, for the unit conversion itself), then
interpolated onto NS-RAD's own simplified daily calendar, which -- like every other
forcing and like config.repeatYear's annual wrap -- treats every year as exactly 365 days
(Feb 29 is dropped from the output axis, same convention tools/fetch_discharge.py uses).

Output: forcing/<site>_river_discharge_interannual_1980-2023_m3sec.csv (CRLF, 2 dp, to
match every other *_river_discharge_*.csv) and
forcing/<site>_toc_interannual_1980-2023_mmolC_m3.csv (LF, 6 dp, to match the other
BOUNDARY_FORCING-consumed series in forcing/idealized_*_cub_river.csv). Both: one value
per day, no header, no trailing newline -- file_module.exfread's genfromtxt reader.

Usage:  python tools/build_interannual_forcings.py [--plot]
"""
import calendar
import datetime as dt
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent.parent
FORCINGS = ROOT / "forcing"
SRC = Path("/Users/rsavelli/Documents/FORTE/Michael/river_timeseries_output")

YEAR_START, YEAR_END = 1980, 2023  # inclusive, matches the full PWBM record on disk

# PWBM column name -> NS-RAD site key
RIVERS = {
    "Colville": "colville",
    "Kuparuk": "kuparuk",
    "Sagavanirktok": "sagavanirktok",
    "Canning": "canning",
}

MGL_TO_MMOLM3 = 83.3  # DOC mg/L -> TOC mmol C/m3 = 1000/12.011, rounded as every
                       # sites/<name>.py DOC->TOC conversion already is


def _read_monthly(path, river):
    """Return {(year, month): value} for one river column of a *_by_river.csv."""
    out = {}
    with open(path) as f:
        header = f.readline().strip().split(",")
        col = header.index(river)
        for line in f:
            parts = line.strip().split(",")
            if not parts or not parts[0]:
                continue
            y, m = int(parts[0][:4]), int(parts[0][5:7])
            out[(y, m)] = float(parts[col])
    return out


def _months():
    """Every (year, month) from YEAR_START-01 to YEAR_END-12, in order."""
    return [(y, m) for y in range(YEAR_START, YEAR_END + 1) for m in range(1, 13)]


def _model_daily_axis(months):
    """Map each (year, month) to the index of that month's midpoint in NS-RAD's
    simplified daily calendar (365 days/year, Feb 29 dropped -- see module docstring),
    and return that calendar's total length. Real days-in-month (leap-aware) are used
    for the midpoint itself, since that only positions the interpolation node."""
    day_idx = 0
    centers = []
    for (y, m) in months:
        ndays = calendar.monthrange(y, m)[1]
        # Position of this month's midpoint within the compressed (no-Feb-29) calendar.
        # A leap Feb's midpoint is unaffected; later months in a leap year shift back by
        # one day relative to the real calendar, which is the intended compression.
        centers.append(day_idx + ndays / 2.0)
        skip_feb29 = 1 if (m == 2 and calendar.isleap(y)) else 0
        day_idx += ndays - skip_feb29
    return np.array(centers), day_idx  # day_idx is now the total compressed length


def _write(path, values, decimals, newline):
    body = newline.join(f"{v:.{decimals}f}" for v in values)
    path.write_text(body, encoding="utf-8", newline="")


def main():
    plot = "--plot" in sys.argv
    FORCINGS.mkdir(exist_ok=True)

    months = _months()
    month_centers, n_days = _model_daily_axis(months)
    day_grid = np.arange(n_days, dtype=float)

    runoff_by_river = {r: _read_monthly(SRC / "runoff_monthly_timeseries_by_river.csv", r)
                        for r in RIVERS}
    doc_by_river = {r: _read_monthly(SRC / "DOC_monthly_timeseries_by_river.csv", r)
                     for r in RIVERS}

    if plot:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        fig, axes = plt.subplots(len(RIVERS), 2, figsize=(11, 9), sharex=True)

    for i, (pwbm_name, site) in enumerate(RIVERS.items()):
        runoff = runoff_by_river[pwbm_name]
        doc = doc_by_river[pwbm_name]

        q_month = np.empty(len(months))
        conc_month = np.full(len(months), np.nan)
        for j, (y, m) in enumerate(months):
            ndays = calendar.monthrange(y, m)[1]
            runoff_mm = runoff[(y, m)]
            q_month[j] = runoff_mm * 1000.0 / (ndays * 86400.0)
            if runoff_mm > 0:
                conc_month[j] = doc[(y, m)] / runoff_mm  # mg/L; the 1 km2/cell factor cancels

        # Gap-fill undefined (zero-flow) months by linear interpolation between the
        # nearest defined open-water months on either side.
        ok = ~np.isnan(conc_month)
        if not ok.all():
            conc_month[~ok] = np.interp(np.flatnonzero(~ok), np.flatnonzero(ok), conc_month[ok])
        toc_month = conc_month * MGL_TO_MMOLM3

        q_daily = np.interp(day_grid, month_centers, q_month)
        toc_daily = np.interp(day_grid, month_centers, toc_month)

        q_path = FORCINGS / f"{site}_river_discharge_interannual_{YEAR_START}-{YEAR_END}_m3sec.csv"
        toc_path = FORCINGS / f"{site}_toc_interannual_{YEAR_START}-{YEAR_END}_mmolC_m3.csv"
        _write(q_path, q_daily, 2, "\r\n")
        _write(toc_path, toc_daily, 6, "\n")

        # 2022 annual-mean check, printed against the values already verified this
        # session (see module docstring).
        idx2022 = [j for j, (y, m) in enumerate(months) if y == 2022]
        q2022_mean = np.average(q_month[idx2022],
                                 weights=[calendar.monthrange(y, m)[1] for (y, m) in
                                          [months[j] for j in idx2022]])
        print(f"{site:14s} 2022 annual-mean Q = {q2022_mean:7.1f} m3/s   "
              f"TOC range {toc_daily.min():6.1f}-{toc_daily.max():6.1f} mmol/m3   "
              f"(DOC {conc_month.min():5.1f}-{conc_month.max():5.1f} mg/L)   "
              f"n_gapfilled_months={int((~ok).sum())}")

        if plot:
            years = np.array([y + (m - 1) / 12.0 for (y, m) in months])
            axes[i, 0].plot(years, q_month)
            axes[i, 0].set_ylabel(f"{site}\nQ [m3/s]")
            axes[i, 1].plot(years, toc_month)
            axes[i, 1].set_ylabel("TOC [mmol/m3]")
    if plot:
        axes[0, 0].set_title("Monthly discharge")
        axes[0, 1].set_title("Monthly TOC (from PWBM DOC)")
        axes[-1, 0].set_xlabel("year")
        axes[-1, 1].set_xlabel("year")
        fig.tight_layout()
        out = FORCINGS.parent / "docs" / "interannual_forcings_preview.png"
        fig.savefig(out, dpi=130)
        print(f"wrote {out}")


if __name__ == "__main__":
    main()
