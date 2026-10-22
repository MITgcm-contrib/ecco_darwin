"""
OSSE report: the FULL FORTE Expedition observation-value framework, end to end --
(1) the nature run + 5 synthetic sampling operators (tools/run_osse_nature.sh,
tools/build_osse_observations.py), (2) the Green's-function sensitivity engine
(tools/build_green_perturbations.py, tools/build_green_functions.py), and (3) the
Fisher-information/degrees-of-freedom-for-signal channel ranking
(tools/rank_observation_value.py) -- which channel most reduces uncertainty in the
model's boundary/forcing inputs, the actual research question this whole framework
answers.

Reads:
  - forcing/<site>_river_discharge_interannual_2005-2023_m3sec.csv (flow-regime
    classification of the 3 analog years)
  - runs/osse_nature/<site>/output.nc, WIDTH PROFILE ONLY (a single time-slice lazy
    read) -- for the SWOT channel-width finding
  - docs/osse/<site>_<channel>.csv (synthetic observations, tools/build_osse_
    observations.py)
  - runs/osse_green/green_functions_<site>.npz (sensitivity fields, tools/
    build_green_functions.py) -- for the example Green's-function figure
  - docs/osse/observation_value_ranking.csv and
    docs/osse/observation_value_ranking_combined.csv (tools/rank_observation_
    value.py's output -- every ranking figure below is built FROM this output, not
    re-derived)

Usage:  python tools/make_osse_report.py
Writes: docs/osse/osse_sampling_report.pdf
"""
import csv
import sys
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

import nsrad_style as S

sys.path.insert(0, str(Path(__file__).resolve().parent))
import rank_observation_value as RANK   # noqa: E402 (CONTROLS/VAR_MAP/CHANNEL_SOURCES/etc.)

ROOT = Path(__file__).resolve().parent.parent
FORC = ROOT / "forcing"
NATURE = ROOT / "runs" / "osse_nature"
GREEN = ROOT / "runs" / "osse_green"
OSSE = ROOT / "docs" / "osse"
SITES = ["colville", "kuparuk", "sagavanirktok"]
STUDY_YEARS = [2007, 2008, 2009]

CHANNEL_LABELS = {
    "grab": "Grab sampling", "mooring_inshore": "Mooring (inshore)",
    "mooring_offshore": "Mooring (offshore)", "drifter": "Drifters",
    "airborne": "Airborne (AVIRIS-3+drone)", "swot": "SWOT (notional)",
    "pace": "PACE (notional)",
}
CHANNEL_ORDER = ["grab", "mooring_inshore", "mooring_offshore", "drifter",
                 "airborne", "swot", "pace"]

S.apply()
S.install_autoscale(1.15)


def _read_csv(site, channel):
    p = OSSE / f"{site}_{channel}.csv"
    if not p.exists():
        return []
    with open(p) as f:
        return list(csv.DictReader(f))


def _discharge(site):
    p = FORC / f"{site}_river_discharge_interannual_2005-2023_m3sec.csv"
    v = np.genfromtxt(open(p, encoding="utf-8-sig"), delimiter=",")
    years = 2005 + np.arange(v.size) / 365.0
    return years, v


def _width_profile(site):
    from netCDF4 import Dataset
    d = Dataset(NATURE / site / "output.nc")
    w = np.asarray(d.variables["width"][0, 1:])   # one time-slice, lazy read
    return w


# =====================================================================================
# Title page
# =====================================================================================
def plot_title_page(pdf):
    fig = plt.figure(figsize=(11, 8.5))
    fig.text(0.08, 0.72, "OSSE Observation-Value Report", fontsize=26, fontweight="bold",
             color=S.INK)
    fig.text(0.08, 0.645, "Nature run, synthetic sampling, Green's-function sensitivity,\n"
                          "and channel ranking for the FORTE Expedition field-campaign plan",
             fontsize=13, color=S.INK2)
    fig.text(0.08, 0.52,
             "Four parts, end to end:",
             fontsize=11, color=S.INK, weight="bold")
    fig.text(0.08, 0.46,
             "1. Nature run + 5 synthetic observation channels matching the FORTE\n"
             "   campaign's real instrument plan as closely as a 1-D model allows.\n"
             "2. Green's-function sensitivity: how the model state responds to a\n"
             "   perturbation in each of 11 boundary/forcing inputs.\n"
             "3. The observation operator connecting (1) and (2) -- which model\n"
             "   variables and controls each channel actually informs.\n"
             "4. The ranking itself: a Fisher-information/degrees-of-freedom-for-\n"
             "   signal score per channel, answering the actual research question --\n"
             "   which observation TYPE most reduces uncertainty in the model's\n"
             "   inputs.",
             fontsize=10.5, color=S.INK, va="top")
    fig.text(0.08, 0.14,
             "Built by tools/run_osse_nature.sh, build_osse_observations.py,\n"
             "build_green_perturbations.py, build_green_functions.py, and\n"
             "rank_observation_value.py -- each documented in full in its own\n"
             "module docstring.",
             fontsize=9, color=S.MUTED, va="top")
    S.brand(fig)
    pdf.savefig(fig)
    plt.close(fig)


# =====================================================================================
# Methodology: campaign phase table
# =====================================================================================
PHASE_ROWS = [
    ("Phase 0", "Late April (TBD)", "--", "Under-ice",
     "NOT SYNTHESIZED -- source doc marks infeasible (negligible under-ice flow)"),
    ("Phase 1a", "Jun 1-5", "152-156", "Grab sampling",
     "Helicopter, pre-breakup; DOC/nutrients/DIC/ALK/Chl-a"),
    ("Phase 1b", "Jun 6-20", "157-171", "Grab sampling + ice drone",
     "Helicopter+boats; ice-elevation drone work NOT synthesized (different product)"),
    ("Phase 2", "Jul 1-12", "182-193", "AVIRIS-3 + moorings deployed",
     "Boats+R/V Ukpik; moorings go in, stay through Phase 3"),
    ("Phase 3", "Jul 20-Aug 7", "201-219", "AVIRIS-3+drone, offshore moorings",
     "Boats+R/V Ukpik+OceanXplorer"),
]


def plot_campaign_table(pdf):
    fig, ax = plt.subplots(figsize=(14, 4.2))
    fig.subplots_adjust(top=0.82, bottom=0.04, left=0.02, right=0.98)
    ax.axis("off")
    col_labels = ["Phase", "Real dates (2027)", "Day-of-year", "Instruments", "Notes"]
    tbl = ax.table(cellText=PHASE_ROWS, colLabels=col_labels, loc="upper center",
                    cellLoc="left", colWidths=[0.07, 0.11, 0.08, 0.18, 0.50])
    tbl.auto_set_font_size(False)
    tbl.set_fontsize(8.5)
    tbl.scale(1, 2.0)
    for (r, c), cell in tbl.get_celld().items():
        cell.set_edgecolor(S.GRID)
        if r == 0:
            cell.set_facecolor(S.ACCENT)
            cell.set_text_props(color="white", fontweight="bold")
    ax.set_title("FORTE campaign phases -> OSSE day-of-year windows",
                  loc="left", fontsize=13, pad=14)
    fig.text(0.02, 0.11,
             "Real campaign is 2027 (Phase 0 \"maybe 2028\"); NS-RAD's forcing stops at\n"
             "2023 and the model has no restart/checkpoint capability, so these windows\n"
             "are mapped onto three ANALOG years below, not literally 2027.",
             fontsize=8, color=S.MUTED, style="italic", va="top")
    S.brand(fig)
    pdf.savefig(fig)
    plt.close(fig)


# =====================================================================================
# Methodology + results: nature-run analog years (flow regime)
# =====================================================================================
def plot_flow_regime(pdf):
    fig, axes = plt.subplots(3, 1, figsize=(11, 8.5), sharex=True)
    fig.suptitle("Nature-run analog years: 2007 (low) / 2008 (average) / 2009 (high)",
                 x=0.05, ha="left", fontsize=13, weight="bold")
    fig.text(0.05, 0.945,
             "Chosen as the cleanest consecutive low/average/high flow-regime triplet "
             "in the 2005-2023 record for all three rivers simultaneously -- the real "
             "2027/2028 campaign has no NS-RAD forcing to simulate directly.",
             color=S.INK2, fontsize=8.0)
    for ax, site in zip(axes, SITES):
        years, q = _discharge(site)
        mask = (years >= 2006.5) & (years < 2010)
        ax.plot(years[mask], q[mask], color=S.MUTED, lw=0.6)
        means = []
        for y, lbl in zip(STUDY_YEARS, ["low", "average", "high"]):
            ym = (years >= y) & (years < y + 1)
            m = q[ym].mean()
            means.append(m)
            ax.axvspan(y, y + 1, color=S.RIVC[site], alpha=0.12)
            ax.text(y + 0.5, ax.get_ylim()[1] * 0.85 if ax.get_ylim()[1] > 0 else 1,
                    f"{lbl}\nmean={m:.0f}", ha="center", fontsize=7.5, color=S.RIVC[site])
        ax.set_ylabel(f"{S.LABEL[site]}\ndischarge [m3/s]", fontsize=8.5, color=S.RIVC[site])
        S.tidy(ax)
    axes[-1].set_xlabel("year")
    S.brand(fig)
    pdf.savefig(fig)
    plt.close(fig)


# =====================================================================================
# Methodology: channel design table
# =====================================================================================
CHANNEL_ROWS = [
    ("Grab sampling", "5 stations/river\n(2-cell-in from\neach boundary +\n3 interior)",
     "4 phases x 3 days/phase\n(first/mid/last)",
     "DOC,CDOM(proxy),NH4,\nNO3,PO4,dSi,DIC,ALK,\nChl-a(DIA proxy)",
     "2-15%\n(species-dependent)"),
    ("Moorings", "1 interior (inshore)\n+ boundary forcing\n(offshore)",
     "Continuous, Jul1-Sep30\n(deployed/recovered)",
     "S, T, O2, current speed\n(inshore); S,T,O2 (offshore)",
     "T+-0.2C, S+-0.1psu,\nO2+-5%"),
    ("Drifters", "Lagrangian\n(advected by model\nvelocity field)",
     "30/year, released in\n4 batches (10/10/5/5)\nby campaign phase",
     "S, T along trajectory +\narrival time at mouth",
     "T+-0.3C, S+-0.2psu"),
    ("Airborne\n(AVIRIS-3+drone)", "Every grid cell\n(full longitudinal\nsnapshot)",
     "Every 8h, Phase 2+3\nonly (not Phase 1)",
     "CDOM proxy, Chl-a proxy,\nSPM (optical only)",
     "15-20%"),
    ("Satellite\n(SWOT+PACE,\nnotional)", "SWOT: width>=100m\ncells; PACE: outer\n1-2 cells only",
     "SWOT ~21-day repeat;\nPACE ~daily, 50% cloud\ndropout",
     "SWOT: depth anomaly;\nPACE: CDOM/Chl-a proxy",
     "SWOT 3-10cm (width-\ndep.); PACE 25%"),
]


def plot_channel_table(pdf):
    fig, ax = plt.subplots(figsize=(14, 8.5))
    fig.subplots_adjust(top=0.90, bottom=0.18, left=0.02, right=0.98)
    ax.axis("off")
    col_labels = ["Channel", "Spatial pattern", "Temporal pattern", "Variables", "Error"]
    tbl = ax.table(cellText=CHANNEL_ROWS, colLabels=col_labels, loc="upper center",
                    cellLoc="left", colWidths=[0.13, 0.19, 0.22, 0.28, 0.18])
    tbl.auto_set_font_size(False)
    tbl.set_fontsize(7.8)
    tbl.scale(1, 3.2)
    for (r, c), cell in tbl.get_celld().items():
        cell.set_edgecolor(S.GRID)
        cell.PAD = 0.02
        if r == 0:
            cell.set_facecolor(S.ACCENT)
            cell.set_text_props(color="white", fontweight="bold")
    ax.set_title("Five observation channels: design summary", loc="left", fontsize=13, pad=14)
    fig.text(0.02, 0.14,
             "NOT synthesized (documented gaps, not oversights): POC/PON (no defensible\n"
             "NS-RAD state-variable proxy), flow cytometry/productivity incubations\n"
             "(process/community measurements, no model analog), Phase 1b ice-elevation\n"
             "drone work (a different product). Satellite channel gated by PHYSICAL\n"
             "limits only (ice, cloud, channel width), not real launch dates -- see\n"
             "module docstring for why.",
             fontsize=7.5, color=S.MUTED, style="italic", va="top")
    S.brand(fig)
    pdf.savefig(fig)
    plt.close(fig)


# =====================================================================================
# Results: grab sampling
# =====================================================================================
def plot_grab_results(pdf):
    fig, axes = plt.subplots(2, 1, figsize=(11, 8.5))
    fig.subplots_adjust(hspace=0.45, top=0.90)
    fig.suptitle("Grab sampling -- station placement and example time series",
                 x=0.05, ha="left", fontsize=13, weight="bold")
    ax0, ax1 = axes

    for site in SITES:
        rows = _read_csv(site, "grab")
        if not rows:
            continue
        stations = sorted(set(float(r["distance_km"]) for r in rows))
        ax0.scatter(stations, [S.LABEL[site]] * len(stations), color=S.RIVC[site], s=40, zorder=3)
    ax0.set_xlabel("distance from mouth [km]")
    ax0.set_title("Station placement (5/river: 2-cell-in from each boundary + 3 interior)",
                   loc="left", fontsize=9.5)
    S.tidy(ax0, grid="x")

    site = "colville"
    rows = _read_csv(site, "grab")
    first_station = min(int(r["station"]) for r in rows) if rows else None
    dic = [(float(r["t_days"]), float(r["value"])) for r in rows
           if r["variable"] == "DIC" and int(r["station"]) == first_station]
    if dic:
        dic.sort()
        tt, vv = zip(*dic)
        ax1.plot(tt, vv, "o-", color=S.RIVC[site], ms=4, lw=0.8, label="DIC [mmol C/m3]")
    chl = [(float(r["t_days"]), float(r["value"])) for r in rows
           if r["variable"] == "Chl_a_mg_m3" and int(r["station"]) == first_station]
    if chl:
        chl.sort()
        ax1b = ax1.twinx()
        tt2, vv2 = zip(*chl)
        ax1b.plot(tt2, vv2, "s--", color=S.WARN, ms=4, lw=0.8, label="Chl-a [mg/m3]")
        ax1b.set_ylabel("Chl-a [mg/m3]", color=S.WARN)
    ax1.set_xlabel("t [days from run start]")
    ax1.set_ylabel("DIC [mmol C/m3]", color=S.RIVC[site])
    ax1.set_title(f"Example: {S.LABEL[site]} mouth-adjacent station, DIC and Chl-a "
                  "across all sampling days (3 years x 4 phases x 3 days/phase)",
                  loc="left", fontsize=9.5)
    S.tidy(ax1)
    S.brand(fig)
    pdf.savefig(fig)
    plt.close(fig)


# =====================================================================================
# Results: moorings
# =====================================================================================
def plot_moorings_results(pdf):
    site = "colville"
    rows_in = _read_csv(site, "moorings")
    year = STUDY_YEARS[1]
    fig, axes = plt.subplots(2, 1, figsize=(11, 7.0), sharex=True)
    fig.suptitle(f"Moorings -- {S.LABEL[site]}, inshore vs. offshore, {year} (average year)",
                 x=0.05, ha="left", fontsize=13, weight="bold")
    for ax, var, ylab in zip(axes, ["temperature", "salinity"], ["T [degC]", "S [psu]"]):
        in_pts = sorted((float(r["t_days"]), float(r["value"])) for r in rows_in
                         if r["channel"] == "mooring_inshore" and r["variable"] == var
                         and r["year"] == str(year))
        off_pts = sorted((float(r["t_days"]), float(r["value"])) for r in rows_in
                          if r["channel"] == "mooring_offshore" and r["variable"] == var
                          and r["year"] == str(year))
        if in_pts:
            tt, vv = zip(*in_pts)
            ax.plot(tt, vv, lw=0.5, color=S.RIVC[site], label="inshore (interior cell)")
        if off_pts:
            tt2, vv2 = zip(*off_pts)
            ax.plot(tt2, vv2, "o", ms=3, color=S.WARN, label="offshore (= marine boundary forcing)")
        ax.set_ylabel(ylab)
        S.tidy(ax)
        ax.legend(loc="best", fontsize=7.5)
    axes[-1].set_xlabel("t [days from run start]")
    fig.subplots_adjust(bottom=0.22)
    fig.text(0.05, 0.11,
             "Offshore samples the marine BOUNDARY FORCING file directly (daily native\n"
             "cadence) -- it validates the assumed boundary input, not interior model\n"
             "skill, a conceptually different product from the inshore mooring despite\n"
             "the shared instrument class.",
             fontsize=7.5, color=S.MUTED, style="italic", va="top")
    S.brand(fig)
    pdf.savefig(fig)
    plt.close(fig)


# =====================================================================================
# Results: drifters
# =====================================================================================
def plot_drifters_results(pdf):
    fig, axes = plt.subplots(1, 2, figsize=(11, 5.5))
    ax0, ax1 = axes
    fig.suptitle("microSWIFT drifters -- example trajectories and arrival-time distribution",
                 x=0.05, ha="left", fontsize=13, weight="bold")

    site = "colville"
    rows = _read_csv(site, "drifters")
    year = STUDY_YEARS[0]
    sample_ids = sorted(set(r["drifter_id"] for r in rows
                             if r["year"] == str(year) and r["phase"] == "phase2"))[:6]
    for did in sample_ids:
        pts = sorted((float(r["t_days"]), float(r["distance_km"])) for r in rows
                      if r["drifter_id"] == did and r["variable"] == "salinity")
        if pts:
            t0 = pts[0][0]
            tt = [p[0] - t0 for p in pts]
            xx = [p[1] for p in pts]
            ax0.plot(tt, xx, lw=1.0, alpha=0.8)
    ax0.set_xlabel("days since release")
    ax0.set_ylabel("distance from mouth [km]")
    ax0.set_title(f"{S.LABEL[site]}, Phase 2 {year}, 6 example drifters", loc="left", fontsize=9.5)
    S.tidy(ax0)

    for site in SITES:
        rows = _read_csv(site, "drifters")
        arrivals = [float(r["value"]) * 24.0 for r in rows
                    if r["variable"] == "arrival_time_at_mouth_days"]
        if arrivals:
            ax1.hist(arrivals, bins=30, histtype="step", color=S.RIVC[site], lw=1.3,
                     label=f"{S.LABEL[site]} (n={len(arrivals)})")
    ax1.set_xlabel("arrival time at mouth [hours]")
    ax1.set_ylabel("count")
    ax1.set_title("Arrival-time distribution, all years/phases", loc="left", fontsize=9.5)
    ax1.legend(fontsize=7.5)
    S.tidy(ax1)
    S.brand(fig)
    pdf.savefig(fig)
    plt.close(fig)


# =====================================================================================
# Results: airborne (tidal-phase coverage)
# =====================================================================================
def plot_airborne_results(pdf):
    site = "colville"
    rows = _read_csv(site, "airborne")
    year = STUDY_YEARS[1]
    fig, ax = plt.subplots(figsize=(11, 5.5))
    fig.subplots_adjust(top=0.82, bottom=0.30)
    fig.suptitle(f"Airborne (AVIRIS-3+drone) -- {S.LABEL[site]} {year} flight times",
                 x=0.05, ha="left", fontsize=13, weight="bold")
    fig.text(0.05, 0.895, "Every 8h within each phase window -- 8h does not divide\n"
             "evenly into any dominant tidal period (M2/S2/O1/K1), so sampled tidal\n"
             "phase keeps advancing rather than re-locking.",
             fontsize=8.0, color=S.INK2, va="top")
    for phase, color in zip(["phase2", "phase3"], [S.RIVC[site], S.WARN]):
        tt = sorted(set(float(r["t_days"]) for r in rows
                         if r["phase"] == phase and r["year"] == str(year)))
        if tt:
            t0 = tt[0]
            ax.scatter([t - t0 for t in tt], [phase] * len(tt), color=color, s=14)
    ax.set_xlabel("days since phase start")
    S.tidy(ax, grid="x")
    fig.text(0.05, 0.12,
             "Phase 2 (12 days) and Phase 3 (19 days) both approach or exceed the\n"
             "~14.8-day M2/S2 spring-neap beat period, so spring-neap diversity falls\n"
             "out of the real campaign windows without special-case scheduling.",
             fontsize=7.5, color=S.MUTED, style="italic", va="top")
    S.brand(fig)
    pdf.savefig(fig)
    plt.close(fig)


# =====================================================================================
# Results: satellite (width-resolvability finding)
# =====================================================================================
def plot_satellite_results(pdf):
    fig, axes = plt.subplots(1, 2, figsize=(11, 5.5))
    ax0, ax1 = axes
    fig.suptitle("Satellite (notional): SWOT width-resolvability and coverage",
                 x=0.05, ha="left", fontsize=13, weight="bold")

    for site in SITES:
        w = _width_profile(site)
        delxi_km = 0.1
        dist = np.arange(len(w)) * delxi_km
        ax0.plot(dist, w, color=S.RIVC[site], lw=1.3, label=S.LABEL[site])
    ax0.axhline(100, color=S.GRID, lw=1.0, ls="--")
    ax0.text(0.1, 105, "SWOT 100 m floor", fontsize=7.5, color=S.INK2)
    ax0.set_xlabel("distance from mouth [km]")
    ax0.set_ylabel("channel width [m]")
    ax0.set_yscale("log")
    ax0.set_title("All three rivers flare from a narrow upstream reach to a wide\n"
                  "mouth -- SWOT's blind spot is upstream-only, not whole-river",
                  loc="left", fontsize=9)
    ax0.legend(fontsize=7.5)
    S.tidy(ax0)

    site = "kuparuk"
    rows = _read_csv(site, "satellite")
    swot = [(float(r["t_days"]), float(r["distance_km"])) for r in rows if r["channel"] == "swot"]
    pace = [(float(r["t_days"]), float(r["distance_km"])) for r in rows if r["channel"] == "pace"]
    # Fixed colors (not S.RIVC[site]): the per-river brand color can coincidentally
    # clash with S.WARN (seen directly for kuparuk, whose own color is also
    # orange-red) -- SWOT/PACE need two colors distinct from EACH OTHER here, not
    # colors tied to river identity, which this single-river panel doesn't need.
    if swot:
        tt, xx = zip(*swot)
        ax1.scatter(tt, xx, s=4, color=S.ACCENT, label=f"SWOT (n={len(swot)})")
    if pace:
        tt2, xx2 = zip(*pace)
        ax1.scatter(tt2, xx2, s=4, color=S.WARN, label=f"PACE (n={len(pace)})")
    ax1.set_xlabel("t [days from run start]")
    ax1.set_ylabel("distance from mouth [km]")
    ax1.set_title(f"{S.LABEL[site]}: coverage gaps = ice season (+ PACE cloud dropout)",
                  loc="left", fontsize=9)
    ax1.legend(fontsize=7.5, loc="upper center")
    S.tidy(ax1)
    S.brand(fig)
    pdf.savefig(fig)
    plt.close(fig)


# =====================================================================================
# Summary table (row counts)
# =====================================================================================
def plot_summary_table(pdf):
    channels = ["grab", "moorings", "drifters", "airborne", "satellite"]
    rows = []
    for site in SITES:
        row = [S.LABEL[site]]
        for ch in channels:
            n = len(_read_csv(site, ch))
            row.append(f"{n:,}")
        rows.append(row)
    fig, ax = plt.subplots(figsize=(9, 2.6))
    fig.subplots_adjust(top=0.78, bottom=0.08, left=0.04, right=0.96)
    ax.axis("off")
    col_labels = ["River"] + [c.capitalize() for c in channels]
    tbl = ax.table(cellText=rows, colLabels=col_labels, loc="upper center", cellLoc="center")
    tbl.auto_set_font_size(False)
    tbl.set_fontsize(9)
    tbl.scale(1, 2.2)
    for (r, c), cell in tbl.get_celld().items():
        cell.set_edgecolor(S.GRID)
        if r == 0:
            cell.set_facecolor(S.ACCENT)
            cell.set_text_props(color="white", fontweight="bold")
    ax.set_title("Synthetic observation row counts, by river x channel", loc="left",
                  fontsize=12, pad=10)
    S.brand(fig)
    pdf.savefig(fig)
    plt.close(fig)


# =====================================================================================
# PART 2 -- Green's-function sensitivity methodology
# =====================================================================================
def plot_green_methodology(pdf):
    fig, ax = plt.subplots(figsize=(14, 8.5))
    fig.subplots_adjust(top=0.90, bottom=0.04, left=0.03, right=0.97)
    ax.axis("off")
    ax.set_title("Green's-function sensitivity: methodology", loc="left",
                  fontsize=14, weight="bold", pad=16)

    ctrl_rows = [["discharge", "River discharge forcing"],
                 ["clb_T", "Marine boundary: sea temperature"]]
    for sp in RANK.MARINE_SPECIES:
        ctrl_rows.append([f"clb_{sp}", f"Marine boundary: {sp}"])
    tbl = ax.table(cellText=ctrl_rows, colLabels=["Control", "What it perturbs"],
                    loc="upper left", bbox=[0.0, 0.22, 0.30, 0.68], cellLoc="left")
    tbl.auto_set_font_size(False)
    tbl.set_fontsize(7.5)
    for (r, c), cell in tbl.get_celld().items():
        cell.set_edgecolor(S.GRID)
        if r == 0:
            cell.set_facecolor(S.ACCENT)
            cell.set_text_props(color="white", fontweight="bold")

    fig.text(0.34, 0.86,
              f"{RANK.N_CTRL} CONTROLS per river (table, left): discharge and the\n"
              "marine boundary for every genuinely time-varying solute in the\n"
              "full-forcing variant. Upstream (cub) boundaries are NOT perturbed --\n"
              "every solute except TOC holds its upstream boundary constant even in\n"
              "the full-forcing variant, so there is nothing time-varying there to\n"
              "be uncertain about.",
              fontsize=9.0, color=S.INK, va="top")

    fig.text(0.34, 0.62,
              "PERTURBATION: a single one-sided +10% relative step applied to the\n"
              "WHOLE forcing series (not a delta-impulse) -- this answers \"if this\n"
              "boundary were systematically biased 10% high, how would the model\n"
              "respond,\" the right question for ranking which OBSERVATIONS best\n"
              "catch a persistent forcing misspecification.\n\n"
              "NS-RAD has no adjoint (a ported legacy Fortran model, not auto-\n"
              "differentiable), so sensitivity is FINITE-DIFFERENCE: one baseline\n"
              "and one perturbed run per control, differenced and divided by the\n"
              "10% step:\n\n"
              "    G_c(doy, cell) = [X_perturbed(doy, cell) - X_baseline(doy, cell)] / 0.10",
              fontsize=9.0, color=S.INK, va="top")

    fig.text(0.34, 0.30,
              "COST: each perturbation run is MAXT_DAYS=730 (365 warmup + 365 post-\n"
              "warmup reference year, 2006) at DAILY cadence -- cheap (~minutes),\n"
              "since Green's functions only need the seasonal SHAPE of the model's\n"
              "sensitivity, not tidal-scale dynamics. 2006 is used because NS-RAD has\n"
              "no restart/checkpoint capability -- it is the first full year\n"
              f"reachable at all. {len(RANK.SITES)} rivers x {RANK.N_CTRL} controls = "
              f"{len(RANK.SITES)*RANK.N_CTRL} perturbation runs, plus 3 baselines = "
              f"{len(RANK.SITES)*RANK.N_CTRL + 3} total.",
              fontsize=9.0, color=S.INK, va="top")

    fig.text(0.03, 0.195,
              "KEY CAVEAT: sensitivity is assumed STATIONARY IN DAY-OF-YEAR across\n"
              "different calendar years -- i.e. \"how a +10% DIC bias in May\n"
              "propagates\" is assumed similar in 2006, 2007, 2008, and 2009, even\n"
              "though the actual flow regime differs across those three analog study\n"
              "years (that is the whole reason they were chosen as low/avg/high). A\n"
              "fully rigorous version would compute a SEPARATE reference run per\n"
              "analog year -- not done here, due to the same no-restart cost problem\n"
              "that made reaching 2007/2008/2009 themselves expensive in the first\n"
              "place. This is the single most important caveat for interpreting the\n"
              "ranking results in Part 4.",
              fontsize=8.5, color=S.MUTED, style="italic", va="top")
    S.brand(fig)
    pdf.savefig(fig)
    plt.close(fig)


def plot_green_example(pdf):
    site = "colville"
    green = RANK._load_green(site)
    fig, axes = plt.subplots(2, 1, figsize=(11, 9.0))
    fig.subplots_adjust(hspace=0.55, top=0.78, bottom=0.17)
    fig.suptitle(f"Example sensitivity fields -- {S.LABEL[site]}, reference year 2006",
                 x=0.05, ha="left", fontsize=13, weight="bold")
    fig.text(0.05, 0.925,
             "Finite-difference sensitivity dX/d(fractional control), by day-of-year\n"
             "and distance from mouth -- both panels reveal the SAME real, physically\n"
             "coherent pattern: marine/discharge influence collapses when the estuary\n"
             "is fully fresh.", color=S.INK2, fontsize=8.0, va="top")

    ax0 = axes[0]
    g = green["S__discharge"]
    delxi_km = 0.1
    cells = [0, 5, 10, 20]
    for ci in cells:
        ax0.plot(np.arange(365), g[:, ci], lw=1.1, label=f"{ci*delxi_km:.1f} km from mouth")
    ax0.axhline(0, color=S.GRID, lw=0.8)
    ax0.set_ylabel("dS / d(fractional\ndischarge)")
    ax0.set_title("Salinity sensitivity to +10% discharge", loc="left", fontsize=9.5)
    ax0.legend(fontsize=7, ncols=4, loc="lower center")
    S.tidy(ax0)

    ax1 = axes[1]
    g2 = green["DIC__clb_DIC"]
    for ci in cells:
        ax1.plot(np.arange(365), g2[:, ci], lw=1.1, label=f"{ci*delxi_km:.1f} km from mouth")
    ax1.axhline(0, color=S.GRID, lw=0.8)
    ax1.set_ylabel("dDIC / d(fractional\nclb_DIC)")
    ax1.set_xlabel("day of year")
    ax1.set_title("DIC sensitivity to +10% marine DIC boundary", loc="left", fontsize=9.5)
    S.tidy(ax1)

    fig.text(0.05, 0.095,
             "Both sensitivities are non-zero only while the estuary is salt-intruded\n"
             "(winter/low-flow); during the spring freshet the domain flushes\n"
             "completely fresh and baseline==perturbed exactly (verified directly\n"
             "against the raw netCDF, not a plotting artifact) -- the marine boundary\n"
             "genuinely cannot influence a cell with zero residual marine water in it.",
             fontsize=7.5, color=S.MUTED, style="italic", va="top")
    S.brand(fig)
    pdf.savefig(fig)
    plt.close(fig)


# =====================================================================================
# PART 3 -- observation-to-control mapping (what each channel can actually inform)
# =====================================================================================
def plot_observation_mapping(pdf):
    rows = []
    for file_key, channel_label, channel_filter in RANK.CHANNEL_SOURCES:
        csv_rows = _read_csv("colville", file_key)
        if channel_filter:
            csv_rows = [r for r in csv_rows if r["channel"] == channel_filter]
        obs_varnames = sorted({r["variable"] for r in csv_rows} - RANK.EXCLUDED_VARS)
        obs_varnames = [v for v in obs_varnames if v in RANK.VAR_MAP]

        if channel_label == "mooring_offshore":
            model_vars = ["S", "T", "O2"]
            ctrls = ["clb_S", "clb_T", "clb_O2"]
            note = ("reads the marine boundary FORCING file directly (not a\n"
                     "propagated model field) -- validates the assumed input itself")
        else:
            model_vars = sorted({RANK.VAR_MAP[v][0] for v in obs_varnames})
            ctrls_set = {"discharge"}
            for mv in model_vars:
                if mv == "T":
                    ctrls_set.add("clb_T")
                elif mv in RANK.MARINE_SPECIES:
                    ctrls_set.add(f"clb_{mv}")
            ctrls = sorted(ctrls_set)
            note = ""
        rows.append([CHANNEL_LABELS[channel_label], ", ".join(model_vars) or "(n/a)",
                     ", ".join(ctrls), note])

    fig, ax = plt.subplots(figsize=(14, 7.5))
    fig.subplots_adjust(top=0.88, bottom=0.22, left=0.02, right=0.98)
    ax.axis("off")
    col_labels = ["Channel", "Model variables observed", "Controls it can inform", "Note"]
    tbl = ax.table(cellText=rows, colLabels=col_labels, loc="upper center",
                    cellLoc="left", colWidths=[0.16, 0.26, 0.32, 0.26])
    tbl.auto_set_font_size(False)
    tbl.set_fontsize(8.0)
    tbl.scale(1, 2.6)
    for (r, c), cell in tbl.get_celld().items():
        cell.set_edgecolor(S.GRID)
        cell.PAD = 0.02
        if r == 0:
            cell.set_facecolor(S.ACCENT)
            cell.set_text_props(color="white", fontweight="bold")
    ax.set_title("Observation-to-control mapping: what each channel can inform",
                  loc="left", fontsize=13, pad=14)
    fig.text(0.02, 0.17,
             "\"Controls it can inform\" = model variables that channel samples AND\n"
             "that have a corresponding Green's-function control (CDOM/Chl-a proxies\n"
             "map to TOC/DIA, per tools/rank_observation_value.py's VAR_MAP). A channel\n"
             "reaching a control does not guarantee it meaningfully CONSTRAINS it --\n"
             "see Part 4 for the actual information content.",
             fontsize=7.5, color=S.MUTED, style="italic", va="top")
    S.brand(fig)
    pdf.savefig(fig)
    plt.close(fig)


# =====================================================================================
# PART 4 -- observation-value ranking results
# =====================================================================================
def _read_ranking_combined():
    p = OSSE / "observation_value_ranking_combined.csv"
    if not p.exists():
        return []
    with open(p) as f:
        return list(csv.DictReader(f))


def _read_ranking_per_site():
    p = OSSE / "observation_value_ranking.csv"
    if not p.exists():
        return []
    with open(p) as f:
        return list(csv.DictReader(f))


def plot_ranking_bars(pdf):
    rows = _read_ranking_combined()
    rows_by_ch = {r["channel"]: r for r in rows}
    channels = [c for c in CHANNEL_ORDER if c in rows_by_ch]
    dfs = [float(rows_by_ch[c]["dfs"]) for c in channels]
    dfs_u = [float(rows_by_ch[c]["dfs_uncorrected"]) for c in channels]
    labels = [CHANNEL_LABELS[c] for c in channels]

    order = np.argsort(dfs)[::-1]
    channels = [channels[i] for i in order]
    dfs = [dfs[i] for i in order]
    dfs_u = [dfs_u[i] for i in order]
    labels = [labels[i] for i in order]

    fig, ax = plt.subplots(figsize=(11, 7.5))
    fig.subplots_adjust(top=0.82, bottom=0.26, left=0.22, right=0.95)
    y = np.arange(len(channels))
    ax.barh(y + 0.18, dfs_u, height=0.32, color=S.MUTED, label="uncorrected (raw obs count)")
    ax.barh(y - 0.18, dfs, height=0.32, color=S.ACCENT, label="effective-sample-size corrected")
    ax.set_yticks(y)
    ax.set_yticklabels(labels, fontsize=9.5)
    ax.invert_yaxis()
    ax.set_xlabel(f"Degrees of freedom for signal (DFS), summed across 3 rivers "
                  f"(max possible = {RANK.N_CTRL * len(SITES)})")
    fig.suptitle("Observation-value ranking: which channel most reduces\n"
                  "boundary/forcing uncertainty", x=0.05, ha="left", fontsize=12.5,
                  weight="bold")
    ax.legend(fontsize=8.5, loc="lower right")
    S.tidy(ax, grid="x")
    fig.text(0.05, 0.12,
             "Mooring (inshore) drops from #2 to #4 once same-day/sub-daily repeat\n"
             "samples (which share an IDENTICAL daily-resolution Green's-function\n"
             "row) are no longer double-counted as independent evidence -- the\n"
             "single biggest qualitative change this correction makes.",
             fontsize=7.5, color=S.MUTED, style="italic", va="top")
    S.brand(fig)
    pdf.savefig(fig)
    plt.close(fig)


def plot_ranking_correction_table(pdf):
    rows = _read_ranking_per_site()
    agg = {}
    for r in rows:
        if r["channel"] == "ALL_COMBINED" or r["n_obs_raw"] in ("", "-1"):
            continue
        ch = r["channel"]
        agg.setdefault(ch, {"n_raw": 0, "n_eff": 0.0})
        agg[ch]["n_raw"] += int(r["n_obs_raw"])
        agg[ch]["n_eff"] += float(r["n_obs_effective"])

    table_rows = []
    for ch in CHANNEL_ORDER:
        if ch not in agg:
            continue
        n_raw, n_eff = agg[ch]["n_raw"], agg[ch]["n_eff"]
        ratio = n_raw / n_eff if n_eff else float("nan")
        table_rows.append([CHANNEL_LABELS[ch], f"{n_raw:,}", f"{n_eff:,.0f}",
                            f"{ratio:,.1f}x"])

    fig, ax = plt.subplots(figsize=(11, 5.8))
    fig.subplots_adjust(top=0.84, bottom=0.26, left=0.04, right=0.96)
    ax.axis("off")
    col_labels = ["Channel", "Raw observations\n(summed, 3 rivers)",
                  "Effective observations\n(summed, 3 rivers)", "Overcounting factor"]
    tbl = ax.table(cellText=table_rows, colLabels=col_labels, loc="upper center",
                    cellLoc="center")
    tbl.auto_set_font_size(False)
    tbl.set_fontsize(9)
    tbl.scale(1, 2.4)
    for (r, c), cell in tbl.get_celld().items():
        cell.set_edgecolor(S.GRID)
        if r == 0:
            cell.set_facecolor(S.ACCENT)
            cell.set_text_props(color="white", fontweight="bold")
    ax.set_title("Effective-sample-size correction, by channel", loc="left",
                  fontsize=12.5, pad=12)
    fig.text(0.02, 0.19,
             "Two caps (tools/rank_observation_value.py): (A) MECHANICAL, every\n"
             "channel -- the Green's-function arrays are daily-resolution only, so\n"
             "same-day repeat samples of the same (variable, cell) share an identical\n"
             "sensitivity row; (B) PHYSICAL, mooring (offshore) only -- it reads the\n"
             "marine boundary forcing directly, whose real update rate is ~monthly\n"
             "(ECCO-Darwin v5), so daily-interpolated samples within the same ~30-day\n"
             "window add no new information.",
             fontsize=7.5, color=S.MUTED, style="italic", va="top")
    S.brand(fig)
    pdf.savefig(fig)
    plt.close(fig)


def plot_top_constrained_table(pdf):
    rows = _read_ranking_combined()
    table_rows = []
    for ch in CHANNEL_ORDER:
        match = next((r for r in rows if r["channel"] == ch), None)
        if match is None:
            continue
        reductions = [(c, float(match.get(f"var_reduction__{c}", 0.0))) for c in RANK.CONTROLS]
        reductions.sort(key=lambda x: -x[1])
        top3 = [f"{c} ({v:.4f})" for c, v in reductions[:3] if v > 1e-5]
        table_rows.append([CHANNEL_LABELS[ch], ", ".join(top3) or "(none)"])

    fig, ax = plt.subplots(figsize=(13, 4.5))
    fig.subplots_adjust(top=0.86, bottom=0.08, left=0.02, right=0.98)
    ax.axis("off")
    col_labels = ["Channel", "Top 3 constrained controls (variance reduction, "
                             f"prior = {RANK.PRIOR_SIGMA_FRAC*100:.0f}% rel. uncertainty)"]
    tbl = ax.table(cellText=table_rows, colLabels=col_labels, loc="upper center",
                    cellLoc="left", colWidths=[0.18, 0.78])
    tbl.auto_set_font_size(False)
    tbl.set_fontsize(8.5)
    tbl.scale(1, 2.6)
    for (r, c), cell in tbl.get_celld().items():
        cell.set_edgecolor(S.GRID)
        cell.PAD = 0.015
        if r == 0:
            cell.set_facecolor(S.ACCENT)
            cell.set_text_props(color="white", fontweight="bold")
    ax.set_title("Which controls each channel actually constrains (corrected, "
                  "summed across 3 rivers)", loc="left", fontsize=12.5, pad=12)
    S.brand(fig)
    pdf.savefig(fig)
    plt.close(fig)


def plot_ranking_summary_text(pdf):
    fig = plt.figure(figsize=(11, 8.5))
    fig.text(0.07, 0.90, "Summary: observation-value ranking", fontsize=15,
             weight="bold", color=S.INK)
    rows = _read_ranking_combined()
    by_ch = {r["channel"]: float(r["dfs"]) for r in rows}
    ordered = sorted(by_ch, key=lambda c: -by_ch[c])
    lines = []
    for i, ch in enumerate(ordered, 1):
        lines.append(f"{i}. {CHANNEL_LABELS[ch]} -- DFS = {by_ch[ch]:.2f}")
    fig.text(0.07, 0.80, "\n".join(lines), fontsize=12, color=S.INK, va="top",
             family="monospace")
    fig.text(0.07, 0.50,
             "Key findings:\n\n"
             "- Grab sampling is the clear top channel -- sparse but spans many\n"
             "  species, stations, and phases, reaching essentially every control.\n\n"
             "- Mooring (offshore) and PACE are surprisingly competitive: both are\n"
             "  INFORMATION-EFFICIENT -- they reach full value from relatively few\n"
             "  truly-independent samples, unlike continuous in-situ monitoring.\n\n"
             "- Mooring (inshore) ranks lower than its raw observation count would\n"
             "  suggest once duplicate same-day samples are no longer double-\n"
             "  counted (Part 4's effective-sample-size correction).\n\n"
             "- SWOT (notional) ranks lowest -- its depth-anomaly signal reaches\n"
             "  only the discharge control, not any biogeochemical boundary.",
             fontsize=10, color=S.INK, va="top")
    fig.text(0.07, 0.11,
             "All results carry the Part 2 stationarity caveat (sensitivity built from\n"
             "ONE reference year, 2006, reused for 2007-2009) and the single, equal\n"
             "20% relative-uncertainty prior on every control -- see those sections\n"
             "before treating small DFS differences as decisive.",
             fontsize=8.5, color=S.MUTED, style="italic", va="top")
    S.brand(fig)
    pdf.savefig(fig)
    plt.close(fig)


# =====================================================================================
def main():
    out = OSSE / "osse_sampling_report.pdf"
    OSSE.mkdir(parents=True, exist_ok=True)
    with PdfPages(out) as pdf:
        plot_title_page(pdf)
        plot_campaign_table(pdf)
        plot_flow_regime(pdf)
        plot_channel_table(pdf)
        plot_grab_results(pdf)
        plot_moorings_results(pdf)
        plot_drifters_results(pdf)
        plot_airborne_results(pdf)
        plot_satellite_results(pdf)
        plot_summary_table(pdf)

        plot_green_methodology(pdf)
        plot_green_example(pdf)

        plot_observation_mapping(pdf)

        plot_ranking_bars(pdf)
        plot_ranking_correction_table(pdf)
        plot_top_constrained_table(pdf)
        plot_ranking_summary_text(pdf)
    print(f"wrote {out.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
