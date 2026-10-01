"""
Forcing-only report: every upstream (riverine) and downstream (marine) boundary
forcing used by the 2005-2023 full-forcing interannual variant
(sites/<name>_interannual_full.py), as time series, plus a summary table of each
variable's source and real-data coverage. No model output/state variables here --
this is exclusively what goes IN to the model, not what comes out.

Reads the forcing/*_interannual_2005-2023_*.csv files built by
tools/build_interannual_{met,humidity,solar,surge,marine,river_temp}.py and the
existing discharge/TOC slices, plus forcing/tidal_constituents.json (tides need no
forcing file at all -- computed analytically, same as fun_module.Tide).

Usage:  python tools/make_forcing_report.py
Writes: docs/ns_rad_forcing_report.pdf
"""
import json
import math
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

import nsrad_style as S

ROOT = Path(__file__).resolve().parent.parent
FORC = ROOT / "forcing"
DOCS = ROOT / "docs"

YEAR_START, YEAR_END = 2005, 2023
N_YEARS = YEAR_END - YEAR_START + 1
N_DAYS = N_YEARS * 365
RIVERS = ["colville", "kuparuk", "sagavanirktok"]

S.apply()
S.install_autoscale(1.15)

# Calendar axis: 2005-01-01 .. 2023-12-31, Feb 29 dropped every leap year, matching
# every interannual forcing file's own convention.
def _year_dates_no_leap(y):
    start = np.datetime64(f"{y}-01-01")
    end = np.datetime64(f"{y+1}-01-01")
    real_days = np.arange(start, end, dtype="datetime64[D]")  # 365 or 366
    is_leap = (y % 4 == 0 and y % 100 != 0) or (y % 400 == 0)
    if is_leap:
        real_days = np.delete(real_days, 59)  # Feb 29
    return real_days


DATES = np.concatenate([_year_dates_no_leap(y) for y in range(YEAR_START, YEAR_END + 1)])
assert len(DATES) == N_DAYS, (len(DATES), N_DAYS)


def load(name):
    return np.genfromtxt(FORC / name, delimiter=",")


# ---------------------------------------------------------------------------
# Variable registry: (label, category, unit, source, coverage_note, pct_real, loader)
# pct_real values are from each build script's own printed summary (this run).
# ---------------------------------------------------------------------------
SHARED_VARS = [
    dict(key="wind", label="Wind speed", unit="m/s", category="Shared (atmospheric)",
         source="NDBC PRDA2, Prudhoe Bay", coverage=f"{YEAR_START}-{YEAR_END} (station start)",
         native_freq="Hourly (native)", hourly=True,
         pct_real=89.0, file="windspeed_hourly_interannual_2005-2023_msec.csv"),
    dict(key="airtemp", label="Air temperature", unit="degC", category="Shared (atmospheric)",
         source="NDBC PRDA2, Prudhoe Bay", coverage=f"{YEAR_START}-{YEAR_END} (station start)",
         native_freq="Hourly (native)", hourly=True,
         pct_real=92.8, file="airtemp_hourly_interannual_2005-2023_degC.csv"),
    dict(key="relhum", label="Relative humidity", unit="fraction", category="Shared (atmospheric)",
         source="Deadhorse Airport (NOAA ISD)", coverage=f"{YEAR_START}-{YEAR_END} (real from 1980)",
         native_freq="Hourly (native)", hourly=True,
         pct_real=94.8, file="relhum_hourly_interannual_2005-2023_frac.csv"),
    dict(key="solar", label="Solar radiation", unit="W/m2", category="Shared (atmospheric)",
         source="NOAA GML Barrow (Utqiagvik) BSRN", coverage=f"{YEAR_START}-{YEAR_END} (real since 1998)",
         native_freq="1-minute (native) -> hourly mean", hourly=True,
         pct_real=98.3, file="solar_hourly_interannual_2005-2023_Wm2.csv"),
    dict(key="surge", label="Storm surge", unit="m", category="Shared (tidal/surge)",
         source="NOAA CO-OPS 9497645, Prudhoe Bay", coverage=f"{YEAR_START}-{YEAR_END} (real since 1995)",
         native_freq="Hourly -> daily mean (tide must cancel out)",
         pct_real=99.4, file="surge_prudhoe_interannual_2005-2023_m.csv"),
    dict(key="river_temp", label="Riverine temperature", unit="degC", category="Upstream (riverine)",
         source="Derived (air->water formula)", coverage=f"{YEAR_START}-{YEAR_END}",
         native_freq="Derived at daily resolution (from daily air temp)",
         pct_real=None, file="river_temp_interannual_2005-2023_degC.csv"),
]

PER_RIVER_UPSTREAM = [
    dict(key="discharge", label="River discharge", unit="m3/s",
         source="PWBM", coverage="1980-2023 (sliced)",
         native_freq="Daily (native -- PWBM's own daily output)",
         pct_real=100.0,
         file_tmpl="{river}_river_discharge_interannual_2005-2023_m3sec.csv"),
    dict(key="toc_cub", label="Upstream DOC (TOC, cub)", unit="mmol C/m3",
         source="PWBM", coverage="1980-2023 (sliced)",
         native_freq="Daily (native -- PWBM's own daily output)",
         pct_real=100.0,
         file_tmpl="{river}_toc_interannual_2005-2023_mmolC_m3.csv"),
]

PER_RIVER_MARINE = [
    dict(key="s", label="Salinity", unit="PSU", var="SALTanom"),
    dict(key="t", label="Sea temperature", unit="degC", var="SST"),
    dict(key="dic", label="DIC", unit="mmol C/m3", var="DIC"),
    dict(key="alk", label="Alkalinity", unit="mmol/m3", var="ALK"),
    dict(key="no3", label="Nitrate (NO3)", unit="mmol N/m3", var="NO3"),
    dict(key="nh4", label="Ammonium (NH4)", unit="mmol N/m3", var="NH4"),
    dict(key="po4", label="Phosphate (PO4)", unit="mmol P/m3", var="PO4"),
    dict(key="dsi", label="Dissolved silica (dSi)", unit="mmol Si/m3", var="SiO2"),
    dict(key="o2", label="Dissolved oxygen", unit="mmol O2/m3", var="O2"),
    dict(key="toc_clb", label="Marine DOC (TOC, clb)", unit="mmol C/m3", var="DOC"),
]
_MARINE_FILE_KEY = {"toc_clb": "toc"}  # species.lower() in the filename vs. this table's key
for d in PER_RIVER_MARINE:
    d["source"] = "ECCO-Darwin v5 (extended)"
    d["coverage"] = "1995-2023 (sliced)" if d["key"] != "t" else "1995-2023 (was PRDA2)"
    d["native_freq"] = "Monthly (native) -> interpolated to daily"
    d["pct_real"] = 100.0
    file_key = _MARINE_FILE_KEY.get(d["key"], d["key"])
    d["file_tmpl"] = "{river}_" + file_key + "_marine_interannual_2005-2023.csv"


def _daily_mean_of_hourly(y):
    """Collapse an hourly (N_DAYS*24-row) series to one daily mean per day, for the
    full-record overview plot -- a true 19-year hourly plot would be an unreadable
    smear at this figure width, same reasoning as the tide plot's daily/zoom split."""
    return y.reshape(N_DAYS, 24).mean(axis=1)


def plot_shared(pdf, spec):
    fig, ax = plt.subplots(figsize=(11, 4.2))
    y = load(spec["file"])
    if spec.get("hourly"):
        y = _daily_mean_of_hourly(y)
    ax.plot(DATES, y, color=S.ACCENT, linewidth=0.6)
    ax.set_ylabel(f"{spec['label']} [{spec['unit']}]")
    title_suffix = " (daily mean shown; native hourly)" if spec.get("hourly") else ""
    ax.set_title(f"{spec['label']} -- {spec['category']}{title_suffix}", loc="left", fontsize=13)
    S.tidy(ax)
    S.brand(fig, extra=f"{spec['source']}  ·  native frequency: {spec['native_freq']}")
    pdf.savefig(fig)
    plt.close(fig)


def plot_per_river(pdf, spec, title_suffix):
    fig, ax = plt.subplots(figsize=(11, 4.2))
    for river in RIVERS:
        fname = spec["file_tmpl"].format(river=river)
        y = load(fname)
        ax.plot(DATES, y, color=S.RIVC[river], linewidth=0.6, label=S.LABEL[river])
    ax.set_ylabel(f"{spec['label']} [{spec['unit']}]")
    ax.set_title(f"{spec['label']} -- {title_suffix}", loc="left", fontsize=13)
    ax.legend(loc="upper right", ncol=3, fontsize=9)
    S.tidy(ax)
    S.brand(fig, extra=f"{spec['source']}  ·  native frequency: {spec['native_freq']}")
    pdf.savefig(fig)
    plt.close(fig)


def tide_elevation(river, days):
    """Harmonic tide + storm surge at the mouth, same formula as fun_module.Tide,
    evaluated directly from the constituents JSON + this river's surge series (shared
    Prudhoe Bay proxy for all three, same as the model itself uses)."""
    consts = json.load(open(FORC / "tidal_constituents.json"))[river]["constituents"]
    t_hours = days * 24.0
    eta = np.zeros_like(days, dtype=float)
    for c in consts:
        eta += c["amp_m"] * np.cos(np.radians(c["speed_deg_hr"] * t_hours - c["phase_deg"]))
    surge = load("surge_prudhoe_interannual_2005-2023_m.csv")
    # surge is DAILY; broadcast each day's value across that day's sub-daily samples
    surge_per_sample = np.repeat(surge, len(days) // len(surge)) if len(days) > len(surge) else surge
    return eta + surge_per_sample[: len(days)]


def plot_tide_surge(pdf):
    # Full-record DAILY view (tide averages out at this resolution -- shows the surge
    # + spring-neap ENVELOPE, not the sub-daily oscillation).
    fig, ax = plt.subplots(figsize=(11, 4.2))
    days_daily = np.arange(N_DAYS, dtype=float)
    for river in RIVERS:
        eta = tide_elevation(river, days_daily)
        ax.plot(DATES, eta, color=S.RIVC[river], linewidth=0.5, label=S.LABEL[river])
    ax.set_ylabel("Sea-surface elevation [m]")
    ax.set_title("Mouth water level (tide + storm surge) -- daily, full record", loc="left", fontsize=13)
    ax.legend(loc="upper right", ncol=3, fontsize=9)
    S.tidy(ax)
    S.brand(fig, extra="Harmonic tide (continuous, analytic) + Prudhoe Bay surge proxy "
                       "(hourly obs -> daily mean)")
    pdf.savefig(fig)
    plt.close(fig)

    # Zoomed 90-day window at hourly resolution -- the full record would be an
    # unreadable smear at the semi-diurnal tidal period.
    fig, ax = plt.subplots(figsize=(11, 4.2))
    hours = np.arange(0, 90 * 24, 1.0)
    days_frac = hours / 24.0
    zoom_start = 5 * 365  # start the zoom in year 6 (2010), away from record edges
    for river in RIVERS:
        consts = json.load(open(FORC / "tidal_constituents.json"))[river]["constituents"]
        eta = np.zeros_like(hours)
        for c in consts:
            eta += c["amp_m"] * np.cos(np.radians(c["speed_deg_hr"] * (hours + zoom_start * 24) - c["phase_deg"]))
        surge = load("surge_prudhoe_interannual_2005-2023_m.csv")
        surge_win = surge[zoom_start: zoom_start + 90]
        surge_hourly = np.repeat(surge_win, 24)[: len(hours)]
        ax.plot(DATES[zoom_start] + np.timedelta64(1, "h") * hours.astype(int),
                eta + surge_hourly, color=S.RIVC[river], linewidth=0.7, label=S.LABEL[river])
    ax.set_ylabel("Sea-surface elevation [m]")
    ax.set_title("Mouth water level -- 90-day zoom (resolves the semi-diurnal tide + spring-neap beat)",
                 loc="left", fontsize=12)
    ax.legend(loc="upper right", ncol=3, fontsize=9)
    S.tidy(ax)
    S.brand(fig)
    pdf.savefig(fig)
    plt.close(fig)


def plot_diurnal_zoom(pdf):
    """A short (7-day, hourly-resolution) window showing the actual diurnal cycle in
    air temperature, solar radiation, wind, and humidity -- the whole point of
    building these at hourly resolution, which the daily-mean overview plots above
    cannot show at all. Picked in early July of a middle record year so the
    continuous-daylight/low-sun-angle Arctic diurnal shape is visible."""
    # Mid-April, NOT solstice: early July at 70N is continuous daylight (no true
    # night), so it would show no real diurnal solar cycle -- mid-April still has a
    # genuine sunrise/sunset each day, which is what this plot is meant to show.
    zoom_start_day = 12 * 365 + 100  # ~April 11 of year 13 (2017), mid-record
    start_hour = zoom_start_day * 24
    n_hours = 7 * 24
    hours = np.arange(n_hours)
    dates_hourly = DATES[zoom_start_day] + np.timedelta64(1, "h") * hours.astype(int)

    fig, axes = plt.subplots(4, 1, figsize=(11, 9), sharex=True)
    panels = [
        ("airtemp", "Air temperature", "degC", S.ACCENT),
        ("solar", "Solar radiation", "W/m2", "#c77f00"),
        ("wind", "Wind speed", "m/s", "#5c5c5c"),
        ("relhum", "Relative humidity", "fraction", "#2a9d6b"),
    ]
    spec_by_key = {s["key"]: s for s in SHARED_VARS}
    for ax, (key, label, unit, color) in zip(axes, panels):
        spec = spec_by_key[key]
        y = load(spec["file"])
        seg = y[start_hour:start_hour + n_hours]
        ax.plot(dates_hourly, seg, color=color, linewidth=1.0)
        ax.set_ylabel(f"{label}\n[{unit}]", fontsize=9)
        S.tidy(ax)
    axes[0].set_title("Diurnal cycle, 7-day zoom (hourly resolution) -- April 2017",
                       loc="left", fontsize=13)
    S.brand(fig, extra="Shows the day/night structure the daily-mean plots above cannot")
    pdf.savefig(fig)
    plt.close(fig)


def build_table_rows():
    rows = []
    for spec in SHARED_VARS:
        rows.append((spec["label"], spec["category"], spec["source"], spec["native_freq"],
                     spec["coverage"],
                     "n/a (derived)" if spec["pct_real"] is None else f"{spec['pct_real']:.0f}%",
                     spec["unit"]))
    for spec in PER_RIVER_UPSTREAM:
        rows.append((spec["label"], "Upstream (riverine)", spec["source"], spec["native_freq"],
                     spec["coverage"], f"{spec['pct_real']:.0f}%", spec["unit"]))
    for spec in PER_RIVER_MARINE:
        rows.append((spec["label"], "Downstream (marine)", spec["source"], spec["native_freq"],
                     spec["coverage"], f"{spec['pct_real']:.0f}%", spec["unit"]))
    rows.append(("Tides (harmonic)", "Downstream (marine)", "NOAA CO-OPS (per-river station)",
                 "Continuous (analytic formula)", "unlimited (analytic)", "100%", "m"))
    return rows


def plot_table(pdf):
    rows = build_table_rows()
    fig, ax = plt.subplots(figsize=(14, 0.34 * len(rows) + 1.0))
    fig.subplots_adjust(top=0.90, bottom=0.04, left=0.02, right=0.98)
    ax.axis("off")
    col_labels = ["Variable", "Category", "Source", "Native frequency", "Period coverage",
                  "% real", "Units"]
    tbl = ax.table(cellText=rows, colLabels=col_labels, loc="upper center", cellLoc="left",
                    colWidths=[0.14, 0.12, 0.23, 0.22, 0.15, 0.07, 0.07])
    tbl.auto_set_font_size(False)
    tbl.set_fontsize(8.5)
    tbl.scale(1, 1.6)
    for (r, c), cell in tbl.get_celld().items():
        cell.set_edgecolor(S.GRID)
        cell.PAD = 0.01
        if r == 0:
            cell.set_facecolor(S.ACCENT)
            cell.set_text_props(color="white", fontweight="bold")
    ax.set_title("Forcing summary: source and real-data coverage (2005-2023 window)",
                 loc="left", fontsize=13, pad=14)
    S.brand(fig)
    pdf.savefig(fig)
    plt.close(fig)


def plot_title_page(pdf):
    fig = plt.figure(figsize=(11, 8.5))
    fig.text(0.1, 0.65, "NS-RAD Forcing Report", fontsize=26, fontweight="bold", color=S.INK)
    fig.text(0.1, 0.58, "Upstream (riverine) and downstream (marine) boundary forcings,\n"
                        f"full-forcing interannual variant, {YEAR_START}-{YEAR_END}",
             fontsize=13, color=S.INK2)
    fig.text(0.1, 0.45, "Forcings only -- no model output/state variables are shown here.\n"
                        "See docs/ns_rad_interannual.pdf and docs/publication_figures/ for "
                        "model output built from these inputs.",
             fontsize=10, color=S.MUTED)
    S.brand(fig)
    pdf.savefig(fig)
    plt.close(fig)


def main():
    out = DOCS / "ns_rad_forcing_report.pdf"
    with PdfPages(out) as pdf:
        plot_title_page(pdf)
        plot_table(pdf)

        for spec in SHARED_VARS:
            if spec["key"] == "river_temp":
                continue  # plotted in the upstream section below
            plot_shared(pdf, spec)
        plot_diurnal_zoom(pdf)

        # Upstream (riverine)
        for spec in PER_RIVER_UPSTREAM:
            plot_per_river(pdf, spec, "Upstream (riverine)")
        river_temp_spec = next(s for s in SHARED_VARS if s["key"] == "river_temp")
        plot_shared(pdf, river_temp_spec)

        # Downstream (marine)
        for spec in PER_RIVER_MARINE:
            plot_per_river(pdf, spec, "Downstream (marine boundary)")
        plot_tide_surge(pdf)

    print(f"wrote {out.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
