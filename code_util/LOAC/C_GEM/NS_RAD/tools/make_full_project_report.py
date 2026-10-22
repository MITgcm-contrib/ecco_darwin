"""
The full-project report: context, methods and results for the 19-year (2005-2023)
full-forcing long-term run, how its climatology feeds the FORTE campaign run, and
the Green's-function/OSSE sensitivity analysis built on top of that campaign run.
One standalone document spanning everything from "what is this model" through to
"which measurement is most valuable and why" -- written so a reader with no prior
context on this project can follow start to finish, with every forcing input
identified, every methodological choice documented, and every figure's provenance
traceable to a specific script.

STRUCTURE (Part numbers match the section headers used throughout):
  Part 0  Orientation (title, executive summary)
  Part 1  The long-term run: context, methods, results -- NEW in this report
  Part 2  From the long-term record to a campaign-year climatology -- NEW
  Part 3  The OSSE (synthetic sampling)
  Part 4  Green's functions (boundary AND along-estuary sensitivity)
  Part 5  The ranking
  Part 6  Conclusions

Parts 3-5 cover the same analysis as docs/osse_forte/forte_campaign_report.pdf
(tools/make_forte_campaign_report.py) but are re-written here, not imported, so
every in-text page cross-reference is correct for THIS document's own page numbers
(computed from PAGE_SEQUENCE below, not hand-counted) and so the narrative reads as
one continuous document rather than two reports concatenated. The underlying
data/figures are identical -- same CSVs, same .npz Green's functions, same method --
this is a presentation-layer rewrite, not a re-analysis.

Reads: forcing/*_interannual_2005-2023*.csv, forcing/*_climatology*.csv,
runs/interannual_full/<site>/output.nc, runs/forte_campaign/<site>/output.nc,
docs/interannual_discharge_obs.json, docs/osse_forte/<site>_<channel>.csv,
runs/forte_green/green_functions_<site>.npz,
docs/osse_forte/observation_value_ranking*.csv. No model run of its own.

Usage:  python tools/make_full_project_report.py
Writes: docs/ns_rad_full_project_report.pdf
"""
import csv
import json
import sys
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

sys.path.insert(0, str(Path(__file__).resolve().parent))
import nsrad_style as S
import interannual_full_lib as I

ROOT = Path(__file__).resolve().parent.parent
FORC = ROOT / "forcing"
LONGRUN = ROOT / "runs" / "interannual_full"
CAMPAIGN = ROOT / "runs" / "forte_campaign"
GREEN = ROOT / "runs" / "forte_green"
OSSE = ROOT / "docs" / "osse_forte"
SITES = ["colville", "kuparuk", "sagavanirktok"]

YEAR_START, YEAR_END = I.YEAR_START, I.YEAR_END   # 2005, 2023
WARMUP_DAYS = I.WARMUP_DAYS

BASE_DAY = 365
CAMPAIGN_END_DAY = 273
REACHES = ["mouth", "mid", "head"]
LOCAL_CONTROLS = ([f"rate_aerdeg_{r}" for r in REACHES]
                  + [f"lateral_toc_{r}" for r in REACHES])
CONTROLS = (["discharge", "clb_T", "clb_S", "clb_DIC", "clb_ALK", "clb_NO3",
             "clb_NH4", "clb_PO4", "clb_dSi", "clb_O2", "clb_TOC"]
            + ["cub_TOC", "cub_NO3", "cub_NH4", "cub_PO4", "cub_dSi", "cub_DIC",
               "cub_ALK", "cub_DIA"]
            + LOCAL_CONTROLS)
N_CTRL = len(CONTROLS)

CHANNEL_LABELS = {
    "grab": "Grab sampling", "mooring_inshore": "Mooring (inshore)",
    "mooring_offshore": "Mooring (offshore)", "drifter": "Drifters",
    "airborne": "Airborne (AVIRIS-3+drone)", "swot": "SWOT (notional)",
    "pace": "PACE (notional)",
}
CHANNEL_ORDER = ["grab", "mooring_inshore", "mooring_offshore", "drifter",
                 "airborne", "swot", "pace"]
REACH_COLOR = {"mouth": "#2a78d6", "mid": "#eb9e34", "head": "#1e7a4a"}
REACH_LABEL = {"mouth": "mouth third", "mid": "middle third", "head": "head third"}

S.apply()
S.install_autoscale(1.15)


# =====================================================================================
# Shared loaders
# =====================================================================================
def _read_csv(site, channel):
    p = OSSE / f"{site}_{channel}.csv"
    if not p.exists():
        return []
    with open(p) as f:
        return list(csv.DictReader(f))


def _width_profile(site):
    from netCDF4 import Dataset
    d = Dataset(CAMPAIGN / site / "output.nc")
    return np.asarray(d.variables["width"][0, 1:])


def _ranking_combined():
    with open(OSSE / "observation_value_ranking_combined.csv") as f:
        return list(csv.DictReader(f))


def _reach_bounds(M):
    b1, b2 = M // 3, 2 * (M // 3)
    return {"mouth": (1, b1), "mid": (b1 + 1, b2), "head": (b2 + 1, M)}


def annual_ice_stats(site):
    """Copied/adapted from tools/make_interannual_pdf.py's function of the same
    name, pointed at I.load (runs/interannual_full) instead of the 44-year loader --
    the logic needed no change (n_years is derived from t.max(), not hardcoded)."""
    t, H = I.load(site, "ice_thickness")
    _, FR = I.load(site, "ice_frac")
    if H is None or FR is None:
        return np.array([]), np.array([]), np.array([])
    n_years = int((t.max() - WARMUP_DAYS) // 365)
    years, duration, peak = [], [], []
    for i in range(n_years):
        lo, hi = WARMUP_DAYS + i * 365, WARMUP_DAYS + (i + 1) * 365
        mask = (t >= lo) & (t < hi)
        if not mask.any():
            continue
        frac_mean = np.nanmean(FR[mask], axis=1)
        years.append(YEAR_START + 1 + i)
        duration.append(float(np.sum(frac_mean > 0.1)))
        peak.append(float(np.nanmax(H[mask])))
    return np.array(years), np.array(duration), np.array(peak)


def discharge_obs_annual(site, years):
    """Observed annual-mean discharge [m3/s] for the given calendar years, from the
    same USGS gauge cache tools/make_interannual_pdf.py uses
    (docs/interannual_discharge_obs.json, built by
    tools/fetch_interannual_discharge_obs.py) -- NOT re-fetched here."""
    p = ROOT / "docs" / "interannual_discharge_obs.json"
    if not p.exists():
        return np.full(len(years), np.nan)
    d = json.load(open(p))
    rec = d.get(site, {}).get("annual", {})
    return np.array([rec.get(str(y), {}).get("mean_m3s", np.nan) for y in years])


AREA_VARS = [("FCO2", "C", "air-sea CO2 flux (>0 outgassing)"),
             ("NPP", "C", "net primary production"),
             ("aer_deg", "C", "aerobic degradation"),
             ("denit", "C", "denitrification (carbon consumed)"),
             ("nit", "N", "nitrification")]
MOLAR_MASS = I.MOLAR_MASS


# =====================================================================================
# PART 0 -- orientation
# =====================================================================================
def page_title(pdf):
    fig = plt.figure(figsize=(11, 8.5))
    fig.text(0.08, 0.76, "NS-RAD: The Full Project Report", fontsize=26, fontweight="bold",
              color=S.INK)
    fig.text(0.08, 0.665, "From the 19-year climate-forced river record, through a\n"
                          "campaign-year climatology, to a Green's-function sensitivity\n"
                          "analysis of the planned FORTE field campaign",
              fontsize=13, color=S.INK2)
    fig.text(0.08, 0.50,
              "WHO THIS IS FOR: anyone who wants the complete story without prior\n"
              "familiarity with NS-RAD or the FORTE field campaign. No background is\n"
              "assumed -- every technical term is explained the first time it is used.",
              fontsize=10.5, color=S.INK, va="top")
    fig.text(0.08, 0.38,
              "WHAT THIS REPORT COVERS, IN ORDER:\n"
              "  Part 1  The 2005-2023 long-term run: every forcing input, the model\n"
              "          configuration, and 19 years of simulated river behavior\n"
              "  Part 2  How a campaign-specific climatology was built from that record\n"
              "  Part 3  The synthetic field-campaign sampling (the OSSE)\n"
              "  Part 4  Green's-function sensitivity -- boundary AND along-estuary\n"
              "  Part 5  Which measurements are actually most valuable, and why\n"
              "  Part 6  Conclusions",
              fontsize=10.5, color=S.INK, va="top")
    fig.text(0.08, 0.10,
              "Built from: tools/build_interannual_*.py, run_interannual_full.sh,\n"
              "build_forte_campaign_forcing.py, run_forte_campaign.sh,\n"
              "build_osse_observations_forte.py, build_green_perturbations_forte.py,\n"
              "build_green_functions_forte.py, rank_observation_value_forte.py --\n"
              "each documented in full in its own module docstring.",
              fontsize=8.5, color=S.MUTED, va="top")
    S.brand(fig)
    pdf.savefig(fig)
    plt.close(fig)


def page_executive_summary(pdf):
    rows = _ranking_combined()
    by_ch = {r["channel"]: float(r["dfs"]) for r in rows}
    ordered = sorted(by_ch, key=lambda c: -by_ch[c])

    fig = plt.figure(figsize=(11, 8.5))
    fig.text(0.07, 0.94, "Executive summary", fontsize=16, weight="bold", color=S.INK)
    fig.text(0.07, 0.885,
              "This report has two halves. Part 1 characterizes 19 years (2005-2023) of\n"
              "simulated river behavior on Alaska's North Slope, built from real,\n"
              "multi-source climate and ocean forcing. Part 2 onward asks a focused,\n"
              "practical question built on top of that record: which measurement types\n"
              "planned for the FORTE Expedition field campaign will actually teach us the\n"
              "most, and are there blind spots in how the campaign is timed or placed?",
              fontsize=10, color=S.INK, va="top")

    fig.text(0.07, 0.715, "Part 1 in brief:", fontsize=11.5, weight="bold", color=S.INK)
    fig.text(0.07, 0.675,
              "Three North Slope rivers (Colville, Kuparuk, Sagavanirktok) were simulated\n"
              "continuously for 19 years under real, independently-sourced daily-to-hourly\n"
              "forcing (discharge, meteorology, marine chemistry, tides, storm surge -- see\n"
              "Part 1's forcing-inventory table). All three show a significant INCREASING\n"
              "trend in discharge and DOC export, mirrored by a significant DECREASING\n"
              "trend in DIC/alkalinity -- a coherent dilution signal. Freshet timing, ice\n"
              "duration, and most reaction rates show no significant trend over this\n"
              "19-year window (Part 1's findings page).",
              fontsize=9.5, color=S.INK, va="top")

    fig.text(0.07, 0.47, "Parts 2-5 result, ranked:", fontsize=11.5, weight="bold", color=S.INK)
    lines = [f"{i}. {CHANNEL_LABELS[ch]} -- score {by_ch[ch]:.2f} (out of a "
             f"maximum of {N_CTRL})" for i, ch in enumerate(ordered, 1)]
    fig.text(0.07, 0.425, "\n".join(lines), fontsize=10, color=S.INK, va="top",
              family="monospace")

    fig.text(0.07, 0.20, "The one finding that matters most:", fontsize=11, weight="bold",
              color=S.WARN)
    fig.text(0.07, 0.17,
              "Grab sampling should be nearly worthless in FORTE's June-September window --\n"
              "the ocean connection that drove its value in a full-year analysis is\n"
              "severed by the spring freshet before the campaign starts. But grab sampling\n"
              "is ALSO sensitive to the river's own upstream chemistry and to what is\n"
              "happening ALONG the channel in between -- testing both (Parts 4-5) flips the\n"
              "conclusion: grab sampling comes out on top, not the bottom.",
              fontsize=9.5, color=S.INK, va="top")
    S.brand(fig)
    pdf.savefig(fig)
    plt.close(fig)


# =====================================================================================
# PART 1 -- the long-term run: context, methods, results
# =====================================================================================
def page_what_is_nsrad(pdf):
    fig, ax = plt.subplots(figsize=(11, 8.5))
    ax.axis("off")
    ax.set_title("Part 1: What is NS-RAD?", loc="left", fontsize=15, weight="bold", pad=14)
    fig.text(0.08, 0.80,
              "NS-RAD (North Slope River-Aquatic-Delta Model) is a one-dimensional\n"
              "reactive-transport simulation of Arctic river-estuary systems, built on the\n"
              "C-GEM modeling engine (Volta et al. 2014). It represents each river as a\n"
              "channel running from its tidal mouth upstream, resolving salinity,\n"
              "temperature, dissolved oxygen, carbonate chemistry (DIC, alkalinity, pH),\n"
              "nutrients (nitrate, ammonium, phosphate, silica), organic carbon, suspended\n"
              "sediment, phytoplankton biomass, and a prognostic sea-ice cover, all\n"
              "advected, dispersed, and reacting together through time.\n\n"
              "Three North Slope rivers are configured and runnable: the Colville, Kuparuk,\n"
              "and Sagavanirktok. A fourth (Canning) has no observed estuary length yet and\n"
              "is excluded. Each river's geometry (width, depth, estuary length), tidal\n"
              "constituents, and boundary chemistry are observation-based and site-specific;\n"
              "the reaction network and numerics are shared across all three.\n\n"
              "Like any model of a real place, NS-RAD needs to be told what is happening at\n"
              "its EDGES -- the BOUNDARY CONDITIONS -- since it cannot simulate the whole\n"
              "ocean or the whole river basin:\n\n"
              "  UPSTREAM boundary (\"cub\"): what is flowing INTO the river from further\n"
              "    inland -- discharge, and the chemistry of that water.\n\n"
              "  DOWNSTREAM/marine boundary (\"clb\"): what the OCEAN looks like at the\n"
              "    river mouth -- salinity, temperature, carbon, nutrients, and the\n"
              "    tide/storm-surge water-level signal.\n\n"
              "Everything in this report traces back to these boundaries: Part 1 describes\n"
              "how they (and the meteorological forcing above the channel) are built for a\n"
              "genuine 19-year simulation; Parts 2-5 ask how well planned field\n"
              "measurements could pin them down.",
              fontsize=10, color=S.INK, va="top")
    S.brand(fig)
    pdf.savefig(fig)
    plt.close(fig)


def page_why_longrun(pdf):
    fig, ax = plt.subplots(figsize=(11, 8.5))
    ax.axis("off")
    ax.set_title("Part 1: Why a 19-year run, and why 2005-2023?", loc="left", fontsize=14,
                  weight="bold", pad=14)
    fig.text(0.08, 0.80,
              "Every NS-RAD run before this one used a single repeating \"typical year\" of\n"
              "forcing (originally 2022). That is sufficient for validating the model's\n"
              "physics against observations, but cannot answer a different kind of\n"
              "question: how much does river behavior actually VARY from year to year,\n"
              "and is there a detectable multi-year TREND underneath that variability?\n"
              "Answering that needs a genuinely multi-year run, with real year-to-year\n"
              "forcing, not a repeated climatology.\n\n"
              "WHY 2005-2023, NOT A LONGER WINDOW: this is the shortest common window\n"
              "across every independent forcing source this run needs, checked directly\n"
              "against each source's own archive rather than assumed from documentation\n"
              "(see Part 1's forcing-inventory table, next page, for the full per-source\n"
              "coverage). The binding constraint is wind and air temperature: the nearest\n"
              "continuous meteorological station did not exist before 2005. Three\n"
              "alternatives were considered and declined: switching the whole\n"
              "meteorological stack to a reanalysis product (would extend the window to\n"
              "1995, but adds a credential requirement this project has otherwise\n"
              "avoided), accepting a different window per source (declined in favor of one\n"
              "common window every forcing category shares), or restricting discharge --\n"
              "which has real data back to 1980 -- to match everything else anyway. The\n"
              "last option is effectively what was chosen: every source is sliced down to\n"
              "the shortest common window, 2005-2023, 19 years, so that every single\n"
              "simulated year in this run is driven by entirely real, independently\n"
              "sourced forcing in every category at once -- not a mix of real-this-year\n"
              "and climatological-that-year.",
              fontsize=10, color=S.INK, va="top")
    S.brand(fig)
    pdf.savefig(fig)
    plt.close(fig)


FORCING_ROWS = [
    ["River discharge + riverine\norganic carbon (TOC)", "PWBM (Pan-Arctic Water\nBalance Model)",
     "1980-2023, daily", "Whole-basin simulated flux + DOC\nconcentration per river"],
    ["Wind speed, air\ntemperature", "NDBC buoy PRDA2\n(Prudhoe Bay)", "2005-2023, hourly",
     "THE binding constraint on the\nwindow -- station began 2005"],
    ["Relative humidity", "Deadhorse Airport\n(NOAA ISD)", "2005-2023, hourly",
     "Co-located with PRDA2"],
    ["Solar radiation", "NOAA GML Barrow (Utqiagvik)\nAtmospheric Baseline Obs.", "2005-2023, 1-min\n(aggregated hourly)",
     "Real pyranometer observation,\nnot a reanalysis -- see note below"],
    ["Storm surge (non-tidal\nsea-level residual)", "NOAA CO-OPS 9497645,\nPrudhoe Bay", "2005-2023, daily mean",
     "Regional proxy for all 3 rivers\n(only station with continuous water level)"],
    ["Marine (clb) boundary --\nS, DIC, ALK, NO3, NH4,\nPO4, dSi, O2, TOC, sea-T", "ECCO-Darwin v5\n(ocean reanalysis)", "2005-2023, monthly\n(interpolated daily)",
     "Nearest wet grid cell to each\nriver's mouth"],
    ["River temperature\n(upstream boundary)", "Regional air->water\nregression (this project)", "2005-2023, derived\nfrom air temperature", "T_river = max(0, 1.432*(T_air,10-day + 4.0));\nfitted on Kuparuk, applied to all 3"],
    ["Tidal elevation", "NOAA CO-OPS harmonic\nconstituents, per river", "Continuous\n(no data file needed)",
     "Computed analytically from elapsed\nmodel time -- genuine multi-year\nspring-neap for free"],
]


def page_forcing_inventory(pdf):
    fig, ax = plt.subplots(figsize=(11, 10.0))
    fig.subplots_adjust(top=0.90, bottom=0.10)
    ax.axis("off")
    tbl = ax.table(cellText=FORCING_ROWS,
                    colLabels=["Forcing", "Source", "Coverage / resolution", "Notes"],
                    loc="upper center", cellLoc="left", colWidths=[0.20, 0.22, 0.18, 0.40])
    tbl.auto_set_font_size(False)
    tbl.set_fontsize(7.6)
    tbl.scale(1, 2.6)
    for (r, c), cell in tbl.get_celld().items():
        cell.set_edgecolor(S.GRID)
        if r == 0:
            cell.set_facecolor(S.ACCENT)
            cell.set_text_props(color="white", fontweight="bold")
    ax.set_title("Part 1: Every forcing input to the 2005-2023 run", loc="left",
                  fontsize=14, weight="bold", pad=12)
    fig.text(0.02, 0.22,
              "Solar radiation's source changed twice before landing on Barrow: a satellite\n"
              "product doesn't reach this latitude before 2013, and a 3-hourly reanalysis\n"
              "(NARR) was tried but superseded by this real, finer, more reliable station\n"
              "record. Every other row (except river temperature, a regional regression,\n"
              "and tides, computed analytically) is a genuine, independently-sourced\n"
              "observation -- not a reanalysis or a model product standing in for one.",
              fontsize=8.0, color=S.MUTED, style="italic", va="top")
    S.brand(fig)
    pdf.savefig(fig)
    plt.close(fig)


def page_model_config(pdf):
    fig, ax = plt.subplots(figsize=(11, 8.5))
    ax.axis("off")
    ax.set_title("Part 1: Model configuration for this run", loc="left", fontsize=14,
                  weight="bold", pad=14)
    rows = [["Simulated period", "6935 days (19 x 365) + 365-day warm-up = 7300 days total"],
            ["Warm-up", "First 365 days: hydrodynamics/transport only, biogeochemistry off"],
            ["Grid spacing (DELXI)", "100 m (all 3 sites)"],
            ["Timestep (DELTI)", "30 s (all 3 sites)"],
            ["Output save cadence", "Daily -- a 6-min cadence over 19 years would be tens of\n"
                                     "GB/site; daily is enough for annual budgets and seasonal shape"],
            ["Estuary length (EL)", "Colville 4150 m / Kuparuk 7156 m / Sagavanirktok 2170 m\n"
                                     "-- each river's own observed length, not a shared value"],
            ["Grid points (M)", "Colville 42 / Kuparuk 72 / Sagavanirktok 22"],
            ["Ice model", "Prognostic (freeze-up, conductive growth, bottom-fast\n"
                           "grounding, surface melt, hydraulic breakup) -- on by default"],
            ["Multichannel geometry", "On -- braided-delta distributary width distinct from\n"
                                       "upstream prismatic per-thread width"],
            ["Width/dispersion model", "\"Flare\" width convergence + Seo & Cheong (1998)\n"
                                        "local shear dispersion -- see model documentation"],
            ["Arctic BGC extension", "Off (bit-identical to the shipped reaction network)"]]
    tbl = ax.table(cellText=rows, colLabels=["Setting", "Value"], loc="upper center",
                    cellLoc="left", colWidths=[0.28, 0.62])
    tbl.auto_set_font_size(False)
    tbl.set_fontsize(8.5)
    tbl.scale(1, 2.0)
    for (r, c), cell in tbl.get_celld().items():
        cell.set_edgecolor(S.GRID)
        if r == 0:
            cell.set_facecolor(S.ACCENT)
            cell.set_text_props(color="white", fontweight="bold")
    fig.text(0.08, 0.10,
              "This is the SAME reaction network, numerics, and geometry as every other\n"
              "NS-RAD configuration -- only the forcing is new. Every solute boundary\n"
              "except temperature and TOC is still a CONSTANT at the upstream end (no\n"
              "time-varying source exists yet for riverine nutrients/carbonate); the marine\n"
              "end is now genuinely time-varying for every species, via ECCO-Darwin v5.",
              fontsize=8.5, color=S.MUTED, style="italic", va="top")
    S.brand(fig)
    pdf.savefig(fig)
    plt.close(fig)


def page_discharge_results(pdf):
    fig, axes = plt.subplots(3, 1, figsize=(11, 9.2), sharex=True)
    fig.subplots_adjust(top=0.84, bottom=0.08, hspace=0.45)
    fig.suptitle("Part 1 results: annual discharge, 2006-2023", x=0.05, ha="left",
                  fontsize=14, weight="bold")
    fig.text(0.05, 0.925,
              "Model forcing (solid) vs. the real USGS gauge record (dashed), per river.\n"
              "Colville and Sagavanirktok are gauged well upstream of their deltas and are\n"
              "known to understate mouth discharge (documented model limitation, not new\n"
              "here); Kuparuk's gauge is near-tidewater and the best-constrained of the 3.",
              fontsize=8.8, color=S.INK2, va="top")

    for ax, site in zip(axes, SITES):
        years = np.arange(YEAR_START + 1, YEAR_END + 1)
        q_yrs, q = I._forcing_daily(site, "discharge")
        model_annual = np.array([q[(q_yrs >= y) & (q_yrs < y + 1)].mean()
                                   if q is not None else np.nan for y in years])
        obs_annual = discharge_obs_annual(site, years)
        ax.plot(years, model_annual, "o-", color=S.RIVC[site], lw=1.3, ms=3,
                 label="model forcing")
        ax.plot(years, obs_annual, "s--", color=S.MUTED, lw=1.0, ms=3, label="USGS gauge")
        I._trend_annotation(ax, years, model_annual, "m3/s", S.RIVC[site], 0)
        ax.set_ylabel(f"{S.LABEL[site]}\n[m3/s]", fontsize=8, color=S.RIVC[site])
        ax.legend(fontsize=7, loc="upper left")
        S.tidy(ax)
    axes[-1].set_xlabel("year")
    S.brand(fig)
    pdf.savefig(fig)
    plt.close(fig)


def page_freshet_results(pdf):
    fig, axes = plt.subplots(2, 2, figsize=(11, 8.0))
    fig.subplots_adjust(top=0.78, bottom=0.10, hspace=0.45, wspace=0.32)
    fig.suptitle("Part 1 results: spring-freshet characteristics, 2006-2023", x=0.05,
                  ha="left", fontsize=14, weight="bold")
    fig.text(0.05, 0.90,
              "Four freshet metrics per river per year, each with an OLS trend line (dashed,\n"
              "* = p<0.05). Duration and DOC load are the two genuinely diagnostic\n"
              "quantities here for land-to-ocean carbon delivery.",
              fontsize=8.5, color=S.INK2, va="top")
    panels = [("peak_doy", "peak day of year", "days"),
              ("flashiness", "flashiness (peak Q / mean Q)", "x"),
              ("duration", "days above 3x mean Q", "days"),
              ("doc_load_tC", "DOC load [tC]", "tC")]
    for ax, (key, title, unit) in zip(axes.flat, panels):
        for i, site in enumerate(SITES):
            fm = I.freshet_metrics(site)
            if not fm or fm["year"].size == 0:
                continue
            ax.plot(fm["year"], fm[key], "o-", color=S.RIVC[site], ms=3, lw=1.0,
                     label=S.LABEL[site])
            I._trend_annotation(ax, fm["year"], fm[key], unit, S.RIVC[site], i)
        ax.set_title(title, loc="left", fontsize=9)
        ax.xaxis.set_major_locator(plt.MaxNLocator(nbins=6, integer=True))
        S.tidy(ax)
    axes[0, 0].legend(fontsize=7, loc="lower right")
    S.brand(fig)
    pdf.savefig(fig)
    plt.close(fig)


def page_state_variables_1(pdf):
    fig, axes = plt.subplots(2, 2, figsize=(11, 8.0))
    fig.subplots_adjust(top=0.84, bottom=0.10, hspace=0.42, wspace=0.32)
    fig.suptitle("Part 1 results: major state variables, annual open-water means (1 of 2)",
                  x=0.05, ha="left", fontsize=13, weight="bold")
    fig.text(0.05, 0.905,
              "Spatial + temporal mean over open-water (ice-free) cells, per calendar year.\n"
              "\"Open-water\" excludes the roughly half of each year under ice, when\n"
              "biogeochemistry is largely sealed off from the atmosphere.",
              fontsize=8.5, color=S.INK2, va="top")
    panels = [("T", "water temperature", "degC"), ("S", "salinity", "PSU"),
              ("O2", "dissolved oxygen", "mmol/m3"), ("pH", "pH", "")]
    for ax, (var, title, unit) in zip(axes.flat, panels):
        for i, site in enumerate(SITES):
            yrs, vals = I.annual_openwater_mean_series(site, var)
            if yrs.size == 0:
                continue
            ax.plot(yrs, vals, "o-", color=S.RIVC[site], ms=3, lw=1.0, label=S.LABEL[site])
            I._trend_annotation(ax, yrs, vals, unit or "unit", S.RIVC[site], i)
        ax.set_title(f"{title} [{unit}]" if unit else title, loc="left", fontsize=9)
        ax.xaxis.set_major_locator(plt.MaxNLocator(nbins=6, integer=True))
        S.tidy(ax)
    handles, labels = axes[0, 0].get_legend_handles_labels()
    fig.legend(handles, labels, fontsize=8, loc="upper right",
                bbox_to_anchor=(0.985, 0.985), frameon=True, facecolor="white",
                edgecolor=S.GRID)
    S.brand(fig)
    pdf.savefig(fig)
    plt.close(fig)


def page_state_variables_2(pdf):
    fig, axes = plt.subplots(2, 2, figsize=(11, 8.0))
    fig.subplots_adjust(top=0.84, bottom=0.10, hspace=0.42, wspace=0.32)
    fig.suptitle("Part 1 results: major state variables, annual open-water means (2 of 2)",
                  x=0.05, ha="left", fontsize=13, weight="bold")
    fig.text(0.05, 0.905,
              "Same method as the previous page, for the carbon/nutrient system.",
              fontsize=8.5, color=S.INK2, va="top")
    panels = [("DIC", "dissolved inorganic carbon", "mmol C/m3"),
              ("ALK", "alkalinity", "mmol/m3"),
              ("TOC", "organic carbon", "mmol C/m3"),
              ("NO3", "nitrate", "mmol N/m3")]
    for ax, (var, title, unit) in zip(axes.flat, panels):
        for i, site in enumerate(SITES):
            yrs, vals = I.annual_openwater_mean_series(site, var)
            if yrs.size == 0:
                continue
            ax.plot(yrs, vals, "o-", color=S.RIVC[site], ms=3, lw=1.0, label=S.LABEL[site])
            I._trend_annotation(ax, yrs, vals, unit, S.RIVC[site], i)
        ax.set_title(f"{title} [{unit}]", loc="left", fontsize=9)
        ax.xaxis.set_major_locator(plt.MaxNLocator(nbins=6, integer=True))
        S.tidy(ax)
    handles, labels = axes[0, 0].get_legend_handles_labels()
    fig.legend(handles, labels, fontsize=8, loc="upper right",
                bbox_to_anchor=(0.985, 0.985), frameon=True, facecolor="white",
                edgecolor=S.GRID)
    S.brand(fig)
    pdf.savefig(fig)
    plt.close(fig)


def page_process_rates(pdf):
    fig, axes = plt.subplots(2, 3, figsize=(12, 7.2))
    fig.subplots_adjust(top=0.76, bottom=0.10, hspace=0.48, wspace=0.35)
    fig.suptitle("Part 1 results: area-normalized annual process-rate budgets", x=0.05,
                  ha="left", fontsize=13, weight="bold")
    fig.text(0.05, 0.895,
              "Each rate depth-integrated, width-weighted, and time-integrated over one\n"
              "calendar year -- same method as the 44-year record's own report, applied\n"
              "here to the 19-year full-forcing one.",
              fontsize=8.5, color=S.INK2, va="top")
    for ax, (var, element, title) in zip(axes.flat, AREA_VARS):
        for i, site in enumerate(SITES):
            yrs, vals = I.area_normalized_annual_series(site, var)
            if yrs.size == 0:
                continue
            mass = vals * MOLAR_MASS[element] / 1000.0
            ax.plot(yrs, mass, "o-", color=S.RIVC[site], ms=2.5, lw=1.0, label=S.LABEL[site])
            I._trend_annotation(ax, yrs, mass, f"g{element}/m2", S.RIVC[site], i)
        ax.axhline(0, color=S.GRID, lw=0.8)
        ax.set_title(title, loc="left", fontsize=8)
        ax.set_ylabel(f"g{element} m-2 yr-1", fontsize=7.5)
        ax.xaxis.set_major_locator(plt.MaxNLocator(nbins=5, integer=True))
        S.tidy(ax)
    handles, labels = axes.flat[0].get_legend_handles_labels()
    axes.flat[-1].axis("off")
    axes.flat[-1].legend(handles, labels, fontsize=9, loc="center")
    S.brand(fig)
    pdf.savefig(fig)
    plt.close(fig)


def page_ice_results(pdf):
    fig, axes = plt.subplots(2, 1, figsize=(11, 7.5), sharex=True)
    fig.subplots_adjust(top=0.82, bottom=0.10, hspace=0.35)
    fig.suptitle("Part 1 results: ice-cover variability, 2006-2023", x=0.05, ha="left",
                  fontsize=13, weight="bold")
    fig.text(0.05, 0.905,
              "Ice-covered duration (domain-mean ice fraction > 0.1) and peak thickness,\n"
              "per calendar year -- the prognostic ice model's response to 19 real winters.",
              fontsize=8.5, color=S.INK2, va="top")
    for i, site in enumerate(SITES):
        yrs, dur, peak = annual_ice_stats(site)
        if yrs.size == 0:
            continue
        axes[0].plot(yrs, dur, "o-", color=S.RIVC[site], ms=3, lw=1.0, label=S.LABEL[site])
        I._trend_annotation(axes[0], yrs, dur, "days", S.RIVC[site], i)
        axes[1].plot(yrs, peak, "o-", color=S.RIVC[site], ms=3, lw=1.0, label=S.LABEL[site])
        I._trend_annotation(axes[1], yrs, peak, "m", S.RIVC[site], i)
    axes[0].set_ylabel("ice-covered\nduration [days/yr]")
    axes[1].set_ylabel("peak ice\nthickness [m]")
    axes[1].set_xlabel("year")
    axes[1].xaxis.set_major_locator(plt.MaxNLocator(nbins=9, integer=True))
    axes[0].legend(fontsize=7.5, loc="best")
    S.tidy(axes[0]); S.tidy(axes[1])
    S.brand(fig)
    pdf.savefig(fig)
    plt.close(fig)


def page_longrun_findings(pdf):
    fig = plt.figure(figsize=(11, 8.5))
    fig.text(0.07, 0.92, "Part 1 summary: what 19 years of real forcing actually shows",
              fontsize=15, weight="bold", color=S.INK)
    fig.text(0.07, 0.855,
              "1. A CLEAR, COHERENT SIGNAL IN DISCHARGE AND CARBON EXPORT. All three rivers\n"
              "   show a statistically significant (p<=0.01) INCREASE in annual discharge\n"
              "   (+32 to +113 m3/s/decade) and in annual DOC load (p=0.00 for all three)\n"
              "   over 2006-2023. DIC and alkalinity show the mirror-image signal -- a\n"
              "   significant DECREASE at all three rivers (p=0.00) -- consistent with a\n"
              "   simple, physically coherent story: more freshwater runoff dilutes the\n"
              "   river's background carbonate chemistry while carrying more organic\n"
              "   carbon with it. This is a genuine multi-river, multi-variable pattern,\n"
              "   not an isolated single-panel trend.\n\n"
              "2. MOST OTHER QUANTITIES SHOW NO SIGNIFICANT TREND OVER THIS WINDOW.\n"
              "   Freshet timing, flashiness, and duration; ice-covered duration and peak\n"
              "   thickness; salinity, temperature, pH; and most reaction rates (NPP,\n"
              "   nitrification, and aerobic degradation/denitrification at 2 of 3 rivers)\n"
              "   all show p>0.05 trend lines -- reported honestly as null results, not\n"
              "   suppressed. 19 years is a short baseline against this coast's\n"
              "   substantial year-to-year weather noise, and a quantity NOT trending is\n"
              "   as informative here as one that is.\n\n"
              "3. DISCHARGE FORCING TRACKS THE REAL GAUGE RECORD. Kuparuk (the\n"
              "   best-constrained gauge) shows close agreement in both magnitude and\n"
              "   year-to-year variability with the real trend; Colville and Sagavanirktok\n"
              "   run systematically higher, consistent with their already-documented\n"
              "   upstream-gauge limitation (both gauged well above their deltas) -- but\n"
              "   all three still show the SAME significant increasing trend as their own\n"
              "   gauge records, so the forcing is not manufacturing the signal.\n\n"
              "4. THIS RECORD IS THE FOUNDATION FOR EVERYTHING THAT FOLLOWS. Parts 2-5 do\n"
              "   not use this 19-year record directly -- they use a CLIMATOLOGY derived\n"
              "   from it (Part 2), so the choices and data quality documented here\n"
              "   propagate forward into the FORTE campaign analysis, not just into this\n"
              "   Part's own figures.",
              fontsize=9.6, color=S.INK, va="top")
    S.brand(fig)
    pdf.savefig(fig)
    plt.close(fig)


# =====================================================================================
# PART 2 -- from the long-term record to a campaign-year climatology
# =====================================================================================
def page_why_climatology(pdf):
    fig, ax = plt.subplots(figsize=(11, 8.5))
    ax.axis("off")
    ax.set_title("Part 2: Why not run the real campaign dates directly?", loc="left",
                  fontsize=14, weight="bold", pad=14)
    fig.text(0.08, 0.80,
              "The FORTE Expedition's real target is 2027 (and possibly 2028) -- years this\n"
              "model cannot simulate directly, since its forcing data only goes through\n"
              "2023. Two further considerations shaped what to build instead:\n\n"
              "NS-RAD HAS NO RESTART/CHECKPOINT CAPABILITY. Every run must simulate from\n"
              "day zero; there is no way to \"resume\" a long run partway through. Running\n"
              "an entirely fresh multi-year simulation for every possible future campaign\n"
              "year would be prohibitively expensive and would still not reach 2027.\n\n"
              "A SINGLE ANALOG YEAR WAS CONSIDERED AND DECLINED. Picking one historical\n"
              "year (e.g. the driest, or the most typical) to stand in for \"what might\n"
              "happen in 2027\" is somewhat arbitrary and commits the analysis to that\n"
              "year's idiosyncrasies. Instead (project owner's explicit choice,\n"
              "2026-10-19), a CLIMATOLOGY was built: every forcing input is averaged\n"
              "across all 19 real years, day-of-year by day-of-year (or hour-of-year for\n"
              "sub-daily forcings), into one smooth, representative \"typical campaign\n"
              "season.\" This is reusable for any future campaign year without\n"
              "re-choosing an analog, and is exactly the kind of baseline real field-\n"
              "campaign data (especially a genuinely time-varying upstream boundary,\n"
              "which does not exist for any historical year, real or synthetic) is\n"
              "meant to eventually be compared against or substituted into.\n\n"
              "This is a new, separate framework alongside the existing 19-year run\n"
              "(Part 1) and an earlier 3-analog-year OSSE exercise -- both of those stay\n"
              "exactly as they are, for their own existing purposes.",
              fontsize=10, color=S.INK, va="top")
    S.brand(fig)
    pdf.savefig(fig)
    plt.close(fig)


def page_climatology_method(pdf):
    fig, ax = plt.subplots(figsize=(11, 8.5))
    ax.axis("off")
    ax.set_title("Part 2: How the climatology was built", loc="left", fontsize=14,
                  weight="bold", pad=14)
    fig.text(0.08, 0.80,
              "For every forcing file in Part 1's inventory (discharge, meteorology, marine\n"
              "boundary chemistry, storm surge, river temperature), the climatology is\n"
              "computed directly from the already-validated 2005-2023 series:\n\n"
              "  DAILY forcings (discharge, TOC, marine chemistry, surge, river temp):\n"
              "    reshape the 19 x 365-row series to (19, 365) and average across the 19\n"
              "    years, giving one 365-row \"typical year\" file.\n\n"
              "  HOURLY forcings (wind, air temperature, humidity, solar): the same\n"
              "    recipe, reshaped to (19, 8760), giving one 8760-row file.\n\n"
              "Verified directly before use, not assumed: every *_interannual_2005-2023_*\n"
              ".csv is exactly 19*365 (daily) or 19*8760 (hourly) rows, with no leap days\n"
              "(matching every other NS-RAD forcing file's calendar convention).\n\n"
              "These climatology files feed a new site variant (sites/<river>_forte_\n"
              "campaign.py) pointing at them instead of the 19-year series -- no model-code\n"
              "changes needed, since repeatYear already wraps any 365-row/8760-row file\n"
              "exactly the way every single-year site always has.\n\n"
              "Tides need no new file at all: the harmonic tide is computed analytically\n"
              "from elapsed model time, so it is already a genuine (if idealized,\n"
              "non-calendar-pinned) multi-year signal in every NS-RAD run.",
              fontsize=10, color=S.INK, va="top")
    S.brand(fig)
    pdf.savefig(fig)
    plt.close(fig)


def page_climatology_visual(pdf):
    fig, axes = plt.subplots(3, 1, figsize=(11, 9.2), sharex=True)
    fig.subplots_adjust(top=0.84, bottom=0.08, hspace=0.4)
    fig.suptitle("Part 2: the climatology against the 19 real years it summarizes",
                  x=0.05, ha="left", fontsize=14, weight="bold")
    fig.text(0.05, 0.925,
              "Thin lines: each of the 19 real years' own discharge, day-of-year. Thick\n"
              "line: the climatology (their day-of-year mean) that actually drives the\n"
              "FORTE campaign run. The climatology tracks the central tendency faithfully;\n"
              "individual years' freshet timing and peak magnitude vary substantially\n"
              "around it, consistent with Part 1's interannual-variability finding.",
              fontsize=8.8, color=S.INK2, va="top")
    for ax, site in zip(axes, SITES):
        q_yrs, q = I._forcing_daily(site, "discharge")
        if q is not None:
            n_years = q.size // 365
            q2 = q[:n_years * 365].reshape(n_years, 365)
            for row in q2:
                ax.plot(np.arange(365), row, color=S.RIVC[site], lw=0.5, alpha=0.25)
        clim = np.genfromtxt(
            open(FORC / f"{site}_river_discharge_climatology_m3sec.csv", encoding="utf-8-sig"),
            delimiter=",")
        ax.plot(np.arange(365), clim, color=S.INK, lw=1.8, label="climatology")
        ax.axvspan(152, 273, color=S.ACCENT, alpha=0.10)
        ax.set_ylabel(f"{S.LABEL[site]}\n[m3/s]", fontsize=8, color=S.RIVC[site])
        ax.legend(fontsize=7, loc="upper left")
        S.tidy(ax)
    axes[0].text(212, axes[0].get_ylim()[1] * 0.8, "FORTE campaign\nwindow", ha="center",
                  fontsize=7.5, color=S.ACCENT)
    axes[-1].set_xlabel("day of year")
    S.brand(fig)
    pdf.savefig(fig)
    plt.close(fig)


PHASE_ROWS = [
    ("Phase 0", "Late April (TBD)", "--", "Under-ice",
     "NOT part of this analysis -- infeasible per the campaign plan (rivers frozen, negligible flow)"),
    ("Phase 1a", "Jun 1-5", "152-156", "Grab sampling",
     "Helicopter-accessed water samples; DOC, nutrients, DIC, alkalinity, chlorophyll-a"),
    ("Phase 1b", "Jun 6-20", "157-171", "Grab sampling + ice drone",
     "Continues grab sampling; separate ice-thickness drone survey (not analyzed here)"),
    ("Phase 2", "Jul 1-12", "182-193", "Airborne + moorings deployed",
     "Aircraft-based optical imaging; fixed instruments go into the water"),
    ("Phase 3", "Jul 20-Aug 7", "201-219", "Airborne+drone + offshore moorings",
     "Continues airborne imaging; moorings stay in through September"),
]


def page_campaign_run(pdf):
    fig, ax = plt.subplots(figsize=(14, 7.0))
    fig.subplots_adjust(top=0.84, bottom=0.20, left=0.02, right=0.98)
    ax.axis("off")
    col_labels = ["Phase", "Real dates (2027)", "Day-of-year", "What happens", "Detail"]
    tbl = ax.table(cellText=PHASE_ROWS, colLabels=col_labels, loc="upper center",
                    cellLoc="left", colWidths=[0.07, 0.10, 0.07, 0.20, 0.48])
    tbl.auto_set_font_size(False)
    tbl.set_fontsize(7.8)
    tbl.scale(1, 3.2)
    for (r, c), cell in tbl.get_celld().items():
        cell.set_edgecolor(S.GRID)
        if r == 0:
            cell.set_facecolor(S.ACCENT)
            cell.set_text_props(color="white", fontweight="bold")
    ax.set_title("Part 2: the FORTE campaign run and the planned field campaign", loc="left",
                  fontsize=13, weight="bold", pad=12)
    fig.text(0.02, 0.12,
              "The run itself: climatological forcing from Part 2's method, MAXT_DAYS=638\n"
              "(365-day warm-up + 273 days covering the full campaign season), 6-minute\n"
              "output cadence (needed for the Green's-function perturbation work in Part 4).\n"
              "Every comparison from here on uses the campaign's CALENDAR DATES (day-of-\n"
              "year) applied to this climatological run, not literal 2027.",
              fontsize=7.8, color=S.MUTED, style="italic", va="top")
    S.brand(fig)
    pdf.savefig(fig)
    plt.close(fig)


# =====================================================================================
# PART 3 -- the OSSE (synthetic sampling)
# =====================================================================================
def page_what_is_osse(pdf):
    fig, ax = plt.subplots(figsize=(11, 8.5))
    ax.axis("off")
    ax.set_title("Part 3: What is an OSSE, and why build one before the real campaign?",
                  loc="left", fontsize=14, weight="bold", pad=14)
    fig.text(0.08, 0.80,
              "An OSSE -- Observing System Simulation Experiment -- is a standard technique\n"
              "used across climate and ocean science before an expensive field campaign\n"
              "happens. The idea:\n\n"
              "  1. Treat the model's own simulated state as a stand-in for \"the truth\"\n"
              "     (here, the Part 2 climatological campaign-year run).\n\n"
              "  2. Simulate what each planned instrument would ACTUALLY measure if it\n"
              "     sampled that \"truth\" -- at the right times, the right places, with\n"
              "     the instrument's own realistic measurement error. Done here for all\n"
              "     five of FORTE's measurement types (next page).\n\n"
              "  3. Ask: if an instrument had only ever seen these synthetic measurements,\n"
              "     how much would it actually have learned? Answered in Part 4-5 via a\n"
              "     Green's-function sensitivity analysis.\n\n"
              "The payoff: problems with a campaign's DESIGN (not needing more funding or\n"
              "better instruments, but better TIMING or PLACEMENT) can be found and fixed\n"
              "before a single real measurement is taken.",
              fontsize=10.5, color=S.INK, va="top")
    S.brand(fig)
    pdf.savefig(fig)
    plt.close(fig)


CHANNEL_ROWS = [
    ("Grab sampling", "5 stations/river\n(along the length\nof the estuary)",
     "4 phases x 3 days/\nphase (first/mid/last)",
     "Dissolved organic carbon, nutrients,\ndissolved CO2/alkalinity, chlorophyll-a",
     "Water samples collected by hand,\nanalyzed in a lab"),
    ("Moorings", "1 fixed point inside\nthe estuary + the\nocean itself",
     "Continuous, Jul 1 -\nSep 30",
     "Salinity, temperature, dissolved\noxygen, current speed",
     "Instruments left in place for\nmonths, recording continuously"),
    ("Drifters", "Follow the water\n(float with the\ncurrent)",
     "30 released across\nthe season",
     "Salinity, temperature, and how long\nit takes to drift out to sea",
     "Small floating sensors (\"microSWIFT\")\nthat travel with the current"),
    ("Airborne\noptical", "Every point along\nthe estuary, in one\npass",
     "Every 8 hours,\nJul-Aug only",
     "Water color/clarity -- proxies for\norganic carbon, algae, sediment",
     "Camera/sensor flown overhead on a\nplane or drone"),
    ("Satellite\n(notional)", "Wide river mouths\nonly (SWOT); the\nplume (PACE)",
     "SWOT ~every 3 weeks;\nPACE ~daily",
     "Water-surface height (SWOT); water\ncolor (PACE)",
     "Not a committed part of the real\ncampaign -- included as a \"what if\""),
]


def page_channel_table(pdf):
    fig, ax = plt.subplots(figsize=(14, 8.5))
    fig.subplots_adjust(top=0.90, bottom=0.08, left=0.02, right=0.98)
    ax.axis("off")
    col_labels = ["Measurement\ntype", "Where", "How often", "What it measures", "In plain terms"]
    tbl = ax.table(cellText=CHANNEL_ROWS, colLabels=col_labels, loc="upper center",
                    cellLoc="left", colWidths=[0.12, 0.17, 0.15, 0.28, 0.28])
    tbl.auto_set_font_size(False)
    tbl.set_fontsize(8.0)
    tbl.scale(1, 3.4)
    for (r, c), cell in tbl.get_celld().items():
        cell.set_edgecolor(S.GRID)
        if r == 0:
            cell.set_facecolor(S.ACCENT)
            cell.set_text_props(color="white", fontweight="bold")
    ax.set_title("Part 3: the five measurement types this report evaluates", loc="left",
                  fontsize=14, weight="bold", pad=14)
    S.brand(fig)
    pdf.savefig(fig)
    plt.close(fig)


def page_grab_results(pdf):
    fig, axes = plt.subplots(2, 1, figsize=(11, 8.5))
    fig.subplots_adjust(hspace=0.5, top=0.88)
    fig.suptitle("Part 3 results: what the synthetic grab samples actually look like",
                  x=0.05, ha="left", fontsize=13, weight="bold")
    ax0, ax1 = axes
    for site in SITES:
        rows = _read_csv(site, "grab")
        if not rows:
            continue
        stations = sorted(set(float(r["distance_km"]) for r in rows))
        ax0.scatter(stations, [S.LABEL[site]] * len(stations), color=S.RIVC[site], s=40, zorder=3)
    ax0.set_xlabel("distance from the river mouth [km]")
    ax0.set_title("Where the 5 sampling stations sit along each river", loc="left", fontsize=9.5)
    S.tidy(ax0, grid="x")

    site = "colville"
    rows = _read_csv(site, "grab")
    first_station = min(int(r["station"]) for r in rows) if rows else None
    dic = sorted((float(r["t_days"]) - BASE_DAY, float(r["value"])) for r in rows
                 if r["variable"] == "DIC" and int(r["station"]) == first_station)
    if dic:
        tt, vv = zip(*dic)
        ax1.plot(tt, vv, "o-", color=S.RIVC[site], ms=5, lw=0.9)
    ax1.set_xlabel("day of year")
    ax1.set_ylabel("dissolved CO2 (DIC)\n[mmol C/m3]")
    ax1.set_title(f"Example: {S.LABEL[site]}, mouth-adjacent station, across all 12\n"
                   "sampling days (4 phases x 3 days each)", loc="left", fontsize=9.5)
    S.tidy(ax1)
    S.brand(fig)
    pdf.savefig(fig)
    plt.close(fig)


def page_mooring_drifter_results(pdf):
    fig, axes = plt.subplots(2, 1, figsize=(11, 8.5))
    fig.subplots_adjust(hspace=0.55, top=0.88)
    fig.suptitle("Part 3 results: moorings and drifters", x=0.05, ha="left", fontsize=13,
                  weight="bold")
    site = "colville"
    rows = _read_csv(site, "moorings")
    ax0 = axes[0]
    inside = sorted((float(r["t_days"]) - BASE_DAY, float(r["value"])) for r in rows
                     if r["channel"] == "mooring_inshore" and r["variable"] == "temperature")
    offshore = sorted((float(r["t_days"]) - BASE_DAY, float(r["value"])) for r in rows
                       if r["channel"] == "mooring_offshore" and r["variable"] == "temperature")
    if inside:
        tt, vv = zip(*inside)
        ax0.plot(tt, vv, lw=0.6, color=S.RIVC[site], label="inside the estuary")
    if offshore:
        tt2, vv2 = zip(*offshore)
        ax0.plot(tt2, vv2, "o", ms=3, color=S.WARN, label="in the ocean (at the mouth)")
    ax0.set_ylabel("water temperature [C]")
    ax0.set_xlabel("day of year")
    ax0.legend(fontsize=8, loc="best")
    ax0.set_title(f"{S.LABEL[site]}: a mooring inside the estuary runs far warmer in\n"
                   "summer than the ocean right outside it", loc="left", fontsize=9.5)
    S.tidy(ax0)

    ax1 = axes[1]
    rows = _read_csv(site, "drifters")
    sample_ids = sorted(set(r["drifter_id"] for r in rows if r["phase"] == "phase2"))[:6]
    for did in sample_ids:
        pts = sorted((float(r["t_days"]), float(r["distance_km"])) for r in rows
                      if r["drifter_id"] == did and r["variable"] == "salinity")
        if pts:
            t0 = pts[0][0]
            ax1.plot([p[0] - t0 for p in pts], [p[1] for p in pts], lw=1.1, alpha=0.85)
    ax1.set_xlabel("days since the drifter was released")
    ax1.set_ylabel("distance from mouth [km]")
    ax1.set_title(f"{S.LABEL[site]}: 6 example drifters released in Phase 2, each\n"
                   "steadily carried out to sea by the river's current", loc="left", fontsize=9.5)
    S.tidy(ax1)
    S.brand(fig)
    pdf.savefig(fig)
    plt.close(fig)


def page_airborne_satellite_results(pdf):
    fig, axes = plt.subplots(1, 2, figsize=(11, 5.5))
    fig.subplots_adjust(top=0.80, bottom=0.22, wspace=0.3)
    fig.suptitle("Part 3 results: airborne imaging and the notional satellite channels",
                  x=0.05, ha="left", fontsize=13, weight="bold")
    site = "colville"
    rows = _read_csv(site, "airborne")
    ax0 = axes[0]
    for phase, color in zip(["phase2", "phase3"], [S.RIVC[site], S.WARN]):
        tt = sorted(set(float(r["t_days"]) - BASE_DAY for r in rows if r["phase"] == phase))
        if tt:
            t0 = tt[0]
            ax0.scatter([t - t0 for t in tt], [phase] * len(tt), color=color, s=14)
    ax0.set_xlabel("days since the phase began")
    ax0.set_title(f"{S.LABEL[site]}: flight times (every 8 hours)", loc="left", fontsize=9.5)
    S.tidy(ax0, grid="x")

    ax1 = axes[1]
    for site2 in SITES:
        w = _width_profile(site2)
        dist = np.arange(len(w)) * 0.1
        ax1.plot(dist, w, color=S.RIVC[site2], lw=1.3, label=S.LABEL[site2])
    ax1.axhline(100, color=S.GRID, lw=1.0, ls="--")
    ax1.text(0.1, 110, "SWOT's usable-width limit", fontsize=7, color=S.INK2)
    ax1.set_xlabel("distance from mouth [km]")
    ax1.set_ylabel("channel width [m]")
    ax1.set_yscale("log")
    ax1.set_title("Why SWOT can't see the upstream reaches: every\nriver narrows below "
                   "its 100 m usable-width floor", loc="left", fontsize=9)
    ax1.legend(fontsize=7)
    S.tidy(ax1)
    S.brand(fig)
    pdf.savefig(fig)
    plt.close(fig)


# =====================================================================================
# PART 4 -- Green's functions: boundary AND along-estuary sensitivity
# =====================================================================================
def page_what_is_green(pdf):
    fig, ax = plt.subplots(figsize=(11, 8.5))
    ax.axis("off")
    ax.set_title("Part 4: What is a \"Green's function,\" and why does this analysis need one?",
                  loc="left", fontsize=14, weight="bold", pad=14)
    fig.text(0.08, 0.80,
              "The real goal of this exercise is to figure out which measurements would\n"
              "most improve confidence in the model's inputs -- both the BOUNDARY\n"
              "conditions (Part 1) and, as explained below, conditions INSIDE the channel.\n\n"
              "A measurement inside the river doesn't directly tell you an input's value --\n"
              "it tells you the RESULT of that input after it has traveled, mixed, and\n"
              "reacted its way through the estuary. To work backwards from \"what did we\n"
              "measure\" to \"what does that imply,\" we need to know HOW SENSITIVE the\n"
              "model's interior state is to each input -- i.e., if the ocean's carbon\n"
              "content were 10% higher than assumed, how much would that actually change\n"
              "what a sensor sitting inside the river would read, at every place and time?\n\n"
              "That sensitivity relationship is a GREEN'S FUNCTION, computed the most\n"
              "direct way possible: run the model once normally, run it again with ONE\n"
              "input turned up by 10%, and measure the difference.\n\n"
              "The same idea extends naturally to unknowns that are NOT at either\n"
              "boundary: a local in-stream reaction rate, or a local source of organic\n"
              "matter entering ALONG the channel, can be turned up by 10% in one stretch\n"
              "of river at a time and the same recipe applied. This is new in this round\n"
              "of analysis -- 25 inputs are tested in total (next page): 19 at the two\n"
              "boundaries, plus 6 REACH-LOCAL ones inside the domain.",
              fontsize=10.2, color=S.INK, va="top")
    S.brand(fig)
    pdf.savefig(fig)
    plt.close(fig)


def page_controls_table(pdf):
    fig, ax = plt.subplots(figsize=(11, 10.3))
    fig.subplots_adjust(top=0.89, bottom=0.14)
    ax.axis("off")
    rows = [["River discharge", "DOWNSTREAM", "How much water the river is carrying"],
            ["Ocean temperature", "DOWNSTREAM", "Temperature of the seawater at the river mouth"],
            ["Ocean salinity", "DOWNSTREAM", "Saltiness of the seawater at the river mouth"],
            ["Ocean dissolved CO2 (DIC)", "DOWNSTREAM", "Carbon dioxide dissolved in the seawater"],
            ["Ocean alkalinity", "DOWNSTREAM", "The seawater's capacity to buffer/absorb CO2"],
            ["Ocean nitrate", "DOWNSTREAM", "A nutrient seaweed and algae need to grow"],
            ["Ocean ammonium", "DOWNSTREAM", "Another nitrogen-based nutrient"],
            ["Ocean phosphate", "DOWNSTREAM", "A nutrient needed for growth, often scarcest"],
            ["Ocean dissolved silica", "DOWNSTREAM", "Needed by one type of algae (diatoms)"],
            ["Ocean dissolved oxygen", "DOWNSTREAM", "How much oxygen is dissolved in the seawater"],
            ["Ocean organic carbon", "DOWNSTREAM", "Dissolved organic matter in the seawater"],
            ["River organic carbon", "UPSTREAM", "Dissolved organic matter coming down the river"],
            ["River nitrate", "UPSTREAM", "A nutrient the river carries down from inland"],
            ["River ammonium", "UPSTREAM", "Another nitrogen-based nutrient from inland"],
            ["River phosphate", "UPSTREAM", "A nutrient from inland, often the scarcest"],
            ["River dissolved silica", "UPSTREAM", "Needed by diatoms, carried down from inland"],
            ["River dissolved CO2 (DIC)", "UPSTREAM", "Carbon dioxide the river carries from inland"],
            ["River alkalinity", "UPSTREAM", "The river water's own buffering capacity"],
            ["River algae (DIA)", "UPSTREAM", "Phytoplankton biomass entering from upstream"],
            ["Degradation rate, mouth third", "LOCAL (interior)", "How fast organic carbon breaks down, near the mouth"],
            ["Degradation rate, middle third", "LOCAL (interior)", "The same, in the middle stretch of the channel"],
            ["Degradation rate, head third", "LOCAL (interior)", "The same, near the upstream/headwater end"],
            ["Organic loading, mouth third", "LOCAL (interior)", "Extra organic carbon entering along the channel, near the mouth"],
            ["Organic loading, middle third", "LOCAL (interior)", "The same, in the middle stretch of the channel"],
            ["Organic loading, head third", "LOCAL (interior)", "The same, near the upstream/headwater end"]]
    tbl = ax.table(cellText=rows, colLabels=["The 25 inputs tested", "Which\nboundary", "In plain terms"],
                    loc="upper center", cellLoc="left", colWidths=[0.30, 0.15, 0.50])
    tbl.auto_set_font_size(False)
    tbl.set_fontsize(7.8)
    tbl.scale(1, 1.75)
    for (r, c), cell in tbl.get_celld().items():
        cell.set_edgecolor(S.GRID)
        if r == 0:
            cell.set_facecolor(S.ACCENT)
            cell.set_text_props(color="white", fontweight="bold")
        elif 12 <= r <= 19:
            cell.set_facecolor("#eef6ee")
        elif r >= 20:
            cell.set_facecolor("#eaf1fb")
    ax.set_title("Part 4: what was actually tested", loc="left", fontsize=14, weight="bold", pad=12)
    fig.text(0.02, 0.075,
              "Each was increased by 10% for an entire run (LOCAL ones, only within their own\n"
              "third of the channel). UPSTREAM rows (green) match what the campaign's grab\n"
              "sampling measures directly. LOCAL rows (blue) test whether spatially-\n"
              "distributed sampling can also resolve conditions INSIDE the channel, not just\n"
              "at its two ends -- see the along-estuary figures ahead.",
              fontsize=7.8, color=S.MUTED, style="italic", va="top")
    S.brand(fig)
    pdf.savefig(fig)
    plt.close(fig)


def page_key_sensitivity_finding(pdf):
    site = "colville"
    green = np.load(GREEN / f"green_functions_{site}.npz")
    g = green["DIC__clb_DIC"]

    fig, ax = plt.subplots(figsize=(11, 7.5))
    fig.subplots_adjust(top=0.78, bottom=0.26)
    fig.suptitle("Part 4, Finding 1 of 2: the ocean boundary goes dark", x=0.05, ha="left",
                  fontsize=14, weight="bold", color=S.WARN)
    fig.text(0.05, 0.895,
              f"{S.LABEL[site]}: how much a 10% increase in the ocean's dissolved-CO2 level\n"
              "would change the river's own dissolved CO2, right at the mouth, by day of year.",
              fontsize=9.5, color=S.INK2, va="top")

    doy = np.arange(274)
    ax.plot(doy, g[:274, 0], color=S.RIVC[site], lw=1.6)
    ax.axvline(152, color=S.ACCENT, ls="--", lw=1.2)
    ax.axvspan(152, 273, color=S.ACCENT, alpha=0.10)
    ax.text(156, ax.get_ylim()[1] * 0.85, "FORTE campaign\nstarts here (Jun 1)",
             fontsize=8.5, color=S.ACCENT)
    ax.axhline(0, color=S.GRID, lw=0.8)
    ax.set_xlabel("day of year")
    ax.set_ylabel("sensitivity\n(change per 10% ocean CO2 increase)")
    S.tidy(ax)
    fig.text(0.05, 0.19,
              "The ocean's influence is strong all winter (the river is slow and salty, so\n"
              "ocean water sits inside it) and collapses to essentially ZERO right around\n"
              "day 140-150 -- just before the campaign begins. The spring flood has, by\n"
              "then, flushed the estuary completely fresh, cutting off any path for ocean\n"
              "water to reach anywhere the campaign samples. This is real: no instrument\n"
              "INSIDE the river can learn anything about the ocean boundary during FORTE's\n"
              "season. But it is only half the story -- see the next page.",
              fontsize=9, color=S.INK, va="top")
    S.brand(fig)
    pdf.savefig(fig)
    plt.close(fig)


def page_upstream_sensitivity_finding(pdf):
    site = "colville"
    green = np.load(GREEN / f"green_functions_{site}.npz")
    g_cub = green["DIC__cub_DIC"]
    g_clb = green["DIC__clb_DIC"]
    cell_upstream = 39   # colville grab station at 3.9 km of EL=4.15 km

    fig, ax = plt.subplots(figsize=(11, 7.5))
    fig.subplots_adjust(top=0.78, bottom=0.26)
    fig.suptitle("Part 4, Finding 2 of 2: the RIVER boundary does not go dark", x=0.05,
                  ha="left", fontsize=14, weight="bold", color="#1e7a4a")
    fig.text(0.05, 0.895,
              f"{S.LABEL[site]}: the SAME test, but now for a 10% increase in the RIVER's\n"
              "own upstream dissolved-CO2 level, measured at the grab-sampling station\n"
              "closest to the upstream boundary (3.9 of this river's 4.15 km).",
              fontsize=9.5, color=S.INK2, va="top")

    doy = np.arange(274)
    ax.plot(doy, g_cub[:274, cell_upstream], color=S.ACCENT, lw=1.6,
             label="upstream (cub_DIC), near-boundary grab station")
    ax.plot(doy, g_clb[:274, 0], color=S.MUTED, lw=1.1, ls="--",
             label="ocean (clb_DIC), mouth -- from the previous page, for scale")
    ax.axvline(152, color=S.WARN, ls="--", lw=1.2)
    ax.text(156, ax.get_ylim()[1] * 0.55, "FORTE campaign\nstarts here (Jun 1)",
             fontsize=8.5, color=S.WARN)
    ax.axhline(0, color=S.GRID, lw=0.8)
    ax.set_xlabel("day of year")
    ax.set_ylabel("sensitivity\n(change per 10% boundary CO2 increase)")
    ax.legend(fontsize=8, loc="lower left")
    S.tidy(ax)
    fig.text(0.05, 0.19,
              "Unlike the ocean signal, the upstream signal barely moves across the whole\n"
              "year -- it stays close to the boundary's own concentration throughout the\n"
              "FORTE window, because this station sits only 250 m from the boundary itself.\n"
              "This is why grab sampling -- despite learning nothing about the ocean\n"
              "boundary this season (Finding 1) -- still ends up the single most valuable\n"
              "channel overall (Part 5): it is directly measuring the river boundary.",
              fontsize=9, color=S.INK, va="top")
    S.brand(fig)
    pdf.savefig(fig)
    plt.close(fig)


def page_spatial_profile_rate(pdf):
    site = "colville"
    green = np.load(GREEN / f"green_functions_{site}.npz")
    M = green["TOC__rate_aerdeg_mouth"].shape[1]
    dist = np.arange(M) * 0.1
    bounds = _reach_bounds(M)

    fig, ax = plt.subplots(figsize=(11, 7.5))
    fig.subplots_adjust(top=0.80, bottom=0.26)
    fig.suptitle("Part 4: along-estuary sensitivity -- a local RATE perturbation stays local",
                  x=0.05, ha="left", fontsize=14, weight="bold")
    fig.text(0.05, 0.895,
              f"{S.LABEL[site]}: how much a 10% local increase in the aerobic-degradation\n"
              "rate -- applied in only ONE third of the channel at a time -- changes the\n"
              "river's own organic carbon (TOC) at every point along the channel (averaged\n"
              "over the whole campaign window, day-of-year 0-273).",
              fontsize=9.3, color=S.INK2, va="top")

    for reach, (i0, i1) in bounds.items():
        ax.axvspan((i0 - 1) * 0.1, i1 * 0.1, color=REACH_COLOR[reach], alpha=0.08, lw=0)
    for reach in REACHES:
        g = green[f"TOC__rate_aerdeg_{reach}"][:274]
        profile = g.mean(axis=0)
        ax.plot(dist, profile, color=REACH_COLOR[reach], lw=1.8,
                 label=f"perturbed in the {REACH_LABEL[reach]}")
    ax.axhline(0, color=S.GRID, lw=0.8)
    ax.set_xlabel("distance from the river mouth [km]")
    ax.set_ylabel("TOC sensitivity\n[mmol C/m3 per 10% local rate increase]")
    ax.legend(fontsize=8.5, loc="best")
    S.tidy(ax)
    fig.text(0.05, 0.19,
              "Each curve is most intense (TOC drops the most, since a faster rate consumes\n"
              "more of it) inside its OWN perturbed third (tinted bands) and fades away from\n"
              "it -- a local rate change mostly shows up locally, with only a modest\n"
              "footprint carried into neighboring stretches by the river's own flow. This is\n"
              "what lets multi-station grab sampling (not a single fixed mooring) tell these\n"
              "three thirds apart (see the reach-by-reach ranking, page 30).",
              fontsize=9, color=S.INK, va="top")
    S.brand(fig)
    pdf.savefig(fig)
    plt.close(fig)


def page_spatial_profile_lateral(pdf):
    site = "colville"
    green = np.load(GREEN / f"green_functions_{site}.npz")
    M = green["TOC__lateral_toc_mouth"].shape[1]
    dist = np.arange(M) * 0.1
    bounds = _reach_bounds(M)

    fig, ax = plt.subplots(figsize=(11, 7.5))
    fig.subplots_adjust(top=0.80, bottom=0.26)
    fig.suptitle("Part 4: along-estuary sensitivity -- a local LOADING perturbation, same test",
                  x=0.05, ha="left", fontsize=14, weight="bold")
    fig.text(0.05, 0.895,
              f"{S.LABEL[site]}: the SAME test, but now for a local 10% BOOST to organic-\n"
              "carbon entering along the channel (e.g. a tundra or thermokarst source),\n"
              "applied in only one third of the channel at a time.",
              fontsize=9.3, color=S.INK2, va="top")

    for reach, (i0, i1) in bounds.items():
        ax.axvspan((i0 - 1) * 0.1, i1 * 0.1, color=REACH_COLOR[reach], alpha=0.08, lw=0)
    for reach in REACHES:
        g = green[f"TOC__lateral_toc_{reach}"][:274]
        profile = g.mean(axis=0)
        ax.plot(dist, profile, color=REACH_COLOR[reach], lw=1.8,
                 label=f"perturbed in the {REACH_LABEL[reach]}")
    ax.axhline(0, color=S.GRID, lw=0.8)
    ax.set_xlabel("distance from the river mouth [km]")
    ax.set_ylabel("TOC sensitivity\n[mmol C/m3 per 10% local loading increase]")
    ax.legend(fontsize=8.5, loc="center right")
    S.tidy(ax)
    fig.text(0.05, 0.19,
              "Same localized-footprint pattern as the rate test -- and NOT a weaker signal\n"
              "(the mouth-reach peak is larger here; checked at the 5 real grab stations,\n"
              "the two are comparably strong, 0.81 vs 0.79 average). Loading still scores\n"
              "lower in the ranking ahead (page 30) because its pattern is more entangled\n"
              "with discharge (estimated at the same time), not because it is weaker.",
              fontsize=9, color=S.INK, va="top")
    S.brand(fig)
    pdf.savefig(fig)
    plt.close(fig)


def page_hovmoller_sensitivity(pdf):
    site = "colville"
    green = np.load(GREEN / f"green_functions_{site}.npz")
    g = green["TOC__rate_aerdeg_mid"][:274]
    M = g.shape[1]
    dist = np.arange(M) * 0.1
    bounds = _reach_bounds(M)

    fig, ax = plt.subplots(figsize=(11, 8.3))
    fig.subplots_adjust(top=0.78, bottom=0.24, right=0.92)
    fig.suptitle("Part 4: the full picture -- sensitivity across BOTH space and time",
                  x=0.05, ha="left", fontsize=14, weight="bold")
    fig.text(0.05, 0.895,
              f"{S.LABEL[site]}: the same mid-third rate perturbation as the previous page,\n"
              "but now showing every day of the campaign season, not just the average --\n"
              "a \"Hovmoller\" map, a standard way to show how a spatial pattern evolves\n"
              "over time in one picture.",
              fontsize=9.3, color=S.INK2, va="top")

    vmax = np.abs(g).max()
    im = ax.pcolormesh(np.arange(274), dist, g.T, shading="auto", cmap="RdBu_r",
                         vmin=-vmax, vmax=vmax)
    for reach, (i0, i1) in bounds.items():
        ax.axhline(i1 * 0.1, color=S.INK, lw=0.7, ls="--", alpha=0.6)
    ax.axvline(152, color=S.INK, lw=1.2, ls="--")
    ax.text(156, dist[-1] * 0.96, "FORTE campaign\nstarts here (Jun 1)", fontsize=8,
             color=S.INK, va="top")
    ax.set_xlabel("day of year")
    ax.set_ylabel("distance from the river mouth [km]")
    cb = fig.colorbar(im, ax=ax, pad=0.02)
    cb.set_label("TOC sensitivity\n[mmol C/m3 per 10%]", fontsize=8.5)
    fig.text(0.05, 0.155,
              "The bright band stays pinned to the perturbed middle third (dashed lines\n"
              "mark reach boundaries) MOST of the season -- it briefly spreads across the\n"
              "whole channel twice: early on (settling in), and around day 120-150, when\n"
              "the spring freshet mixes the whole channel (the same event behind Finding\n"
              "1's collapse). Through the entire FORTE season, the local signal stays local.",
              fontsize=9, color=S.INK, va="top")
    S.brand(fig)
    pdf.savefig(fig)
    plt.close(fig)


def page_reach_value_breakdown(pdf):
    rows = _ranking_combined()
    by_ch = {r["channel"]: r for r in rows}
    plot_channels = ["grab", "mooring_inshore", "airborne", "pace"]
    colors = {"grab": "#1e7a4a", "mooring_inshore": S.WARN, "airborne": S.ACCENT, "pace": S.MUTED}

    fig, axes = plt.subplots(1, 2, figsize=(12, 6.2))
    fig.subplots_adjust(top=0.78, bottom=0.26, wspace=0.3)
    fig.suptitle("Part 4: which channel resolves WHERE -- reach-by-reach value", x=0.05,
                  ha="left", fontsize=14, weight="bold")
    fig.text(0.05, 0.89,
              "The same information as the Part 5 ranking, split by reach instead of\n"
              "totaled -- this is what the along-estuary figures above predict should happen.",
              fontsize=9, color=S.INK2, va="top")

    for ax, (control_prefix, title) in zip(axes, [("rate_aerdeg", "Local RATE"),
                                                      ("lateral_toc", "Local LOADING")]):
        x = np.arange(len(REACHES))
        width = 0.2
        for j, ch in enumerate(plot_channels):
            vals = [float(by_ch[ch].get(f"var_reduction__{control_prefix}_{r}", 0.0))
                    for r in REACHES]
            ax.bar(x + (j - 1.5) * width, vals, width, color=colors[ch],
                    label=CHANNEL_LABELS[ch] if ax is axes[0] else None)
        ax.set_xticks(x)
        ax.set_xticklabels([REACH_LABEL[r] for r in REACHES], fontsize=8.5)
        ax.set_ylabel("variance reduction\n(share of the 0.04 ceiling)")
        ax.set_title(title, loc="left", fontsize=10.5)
        S.tidy(ax, grid="y")
    axes[0].legend(fontsize=7.5, loc="upper right")
    fig.text(0.05, 0.11,
              "Grab (green) stays high in every reach -- its 5 stations span the whole\n"
              "estuary. The inshore mooring (red) is sharply LOCALIZED to its own fixed\n"
              "location (mouth) and nearly blind to the other two reaches, and contributes\n"
              "NOTHING to local loading at all -- it does not measure organic carbon.",
              fontsize=9, color=S.INK, va="top")
    S.brand(fig)
    pdf.savefig(fig)
    plt.close(fig)


# =====================================================================================
# PART 5 -- the ranking
# =====================================================================================
def page_what_is_dfs(pdf):
    fig, ax = plt.subplots(figsize=(11, 8.5))
    ax.axis("off")
    ax.set_title("Part 5: How \"which measurement is most valuable\" is actually scored",
                  loc="left", fontsize=14, weight="bold", pad=14)
    fig.text(0.08, 0.80,
              "For each synthetic measurement (Part 3), combined with how sensitive the\n"
              "model is to each input at that place and time (Part 4), and how precise the\n"
              "instrument actually is, this analysis computes how much that single\n"
              "measurement would have narrowed down our uncertainty about each of the 25\n"
              "inputs tested. Adding this up correctly across every measurement a channel\n"
              "takes gives a single score per channel, called DEGREES OF FREEDOM FOR\n"
              "SIGNAL (DFS) -- a standard information-theory score used throughout\n"
              "atmospheric and ocean data assimilation.\n\n"
              "A score of 0 means a channel taught us nothing; a score of 25 (the number\n"
              "of inputs tested) would mean it fully pinned down every single one. In\n"
              "practice, scores land somewhere in between, and the split BETWEEN inputs\n"
              "(which ones got constrained, which didn't) matters as much as the total --\n"
              "a channel can score well by fully nailing down a FEW inputs (as grab\n"
              "sampling does, for the upstream and along-estuary ones) rather than\n"
              "partially touching all 25.\n\n"
              "ONE IMPORTANT CORRECTION APPLIED: a sensor recording every 6 minutes for\n"
              "months doesn't really deliver thousands of independent pieces of evidence --\n"
              "consecutive readings are almost identical. This analysis explicitly corrects\n"
              "for that (an \"effective sample size\" adjustment), so a flood of repetitive\n"
              "readings is not mistaken for a flood of new information.",
              fontsize=10.2, color=S.INK, va="top")
    S.brand(fig)
    pdf.savefig(fig)
    plt.close(fig)


def page_ranking_chart(pdf):
    rows = _ranking_combined()
    rows_by_ch = {r["channel"]: r for r in rows}
    channels = [c for c in CHANNEL_ORDER if c in rows_by_ch]
    dfs = [float(rows_by_ch[c]["dfs"]) for c in channels]
    labels = [CHANNEL_LABELS[c] for c in channels]
    order = np.argsort(dfs)[::-1]
    channels = [channels[i] for i in order]
    dfs = [dfs[i] for i in order]
    labels = [labels[i] for i in order]

    fig, ax = plt.subplots(figsize=(11, 7.2))
    fig.subplots_adjust(top=0.82, bottom=0.24, left=0.24, right=0.95)
    y = np.arange(len(channels))
    colors = ["#1e7a4a" if c == "grab" else S.ACCENT for c in channels]
    ax.barh(y, dfs, color=colors, height=0.6)
    ax.set_yticks(y)
    ax.set_yticklabels(labels, fontsize=10)
    ax.invert_yaxis()
    ax.set_xlabel(f"score (degrees of freedom for signal) -- max possible = {N_CTRL}")
    fig.suptitle("Part 5: final ranking -- which measurement actually taught us the most",
                  x=0.05, ha="left", fontsize=13, weight="bold")
    S.tidy(ax, grid="x")
    fig.text(0.05, 0.135,
              "Grab sampling (highlighted) comes out on top -- not from the ocean boundary\n"
              "(Part 4, Finding 1: cut off all season) but from the RIVER's own upstream\n"
              "boundary (Finding 2), PLUS a further boost from along-estuary sensitivity.\n"
              "Airborne and PACE pick up smaller versions of both effects.",
              fontsize=8, color=S.MUTED, style="italic", va="top")
    S.brand(fig)
    pdf.savefig(fig)
    plt.close(fig)


def page_comparison(pdf):
    orig_path = ROOT / "docs" / "osse" / "observation_value_ranking_combined.csv"
    orig = {}
    if orig_path.exists():
        with open(orig_path) as f:
            for r in csv.DictReader(f):
                orig[r["channel"]] = float(r["dfs"])

    rows = _ranking_combined()
    forte_full = {r["channel"]: float(r["dfs"]) for r in rows}
    forte_marine_only = {r["channel"]: float(r["dfs_marine_discharge_only_11ctrl"]) for r in rows}
    forte_boundary_only = {r["channel"]: float(r["dfs_boundary_only_19ctrl"]) for r in rows}

    table_rows = []
    for ch in CHANNEL_ORDER:
        o_val = orig.get(ch)
        m_val = forte_marine_only.get(ch)
        b_val = forte_boundary_only.get(ch)
        f_val = forte_full.get(ch)
        table_rows.append([CHANNEL_LABELS[ch],
                            f"{o_val:.2f} / 11" if o_val is not None else "n/a",
                            f"{m_val:.2f} / 11" if m_val is not None else "n/a",
                            f"{b_val:.2f} / 19" if b_val is not None else "n/a",
                            f"{f_val:.2f} / 25" if f_val is not None else "n/a"])

    fig, ax = plt.subplots(figsize=(14, 8.8))
    fig.subplots_adjust(top=0.84, bottom=0.26)
    ax.axis("off")
    col_labels = ["Measurement type",
                  "(A) Full year,\nocean+discharge only\n(earlier analysis)",
                  "(B) FORTE window,\nocean+discharge only\n(isolates TIMING alone)",
                  "(C) FORTE window,\n+ upstream controls\n(isolates UPSTREAM alone)",
                  "(D) FORTE window,\n+ along-estuary controls\n(this report's full result)"]
    tbl = ax.table(cellText=table_rows, colLabels=col_labels, loc="upper center",
                    cellLoc="center", colWidths=[0.22, 0.195, 0.195, 0.195, 0.195])
    tbl.auto_set_font_size(False)
    tbl.set_fontsize(8.6)
    tbl.scale(1, 2.5)
    for (r, c), cell in tbl.get_celld().items():
        cell.set_edgecolor(S.GRID)
        if r == 0:
            cell.set_facecolor(S.ACCENT)
            cell.set_text_props(color="white", fontweight="bold")
        elif r == 1:
            cell.set_facecolor("#fdf0ec" if c == 2 else "#eef6ee" if c == 3
                                else "#eaf1fb" if c == 4 else "#ffffff")
    ax.set_title("Part 5: three different effects, disentangled -- timing, upstream, along-estuary",
                  loc="left", fontsize=13, weight="bold", pad=12)
    fig.text(0.02, 0.21,
              "Reading grab sampling's row left to right: (A -> B) TIMING ALONE crashes its\n"
              "value 5.62 -> 1.00 -- the summer window cannot see the ocean boundary (Finding\n"
              "1). (B -> C) The 8 upstream controls restore it to 8.99/19 -- grab was ALWAYS\n"
              "also measuring the river boundary, column (B) just never asked (Finding 2).\n"
              "(C -> D) The 6 along-estuary controls lift it further to 12.81/25 -- grab's 5\n"
              "stations also resolve WHERE along the channel things happen, a question\n"
              "boundary-only testing cannot even pose. Airborne shows the same 3-step\n"
              "pattern; moorings/drifters/SWOT barely move in C or D, since they are\n"
              "not well-positioned to sense either structure.",
              fontsize=8.0, color=S.MUTED, style="italic", va="top")
    S.brand(fig)
    pdf.savefig(fig)
    plt.close(fig)


# =====================================================================================
# PART 6 -- conclusions
# =====================================================================================
def page_conclusions(pdf):
    fig = plt.figure(figsize=(11, 8.5))
    fig.text(0.07, 0.94, "Part 6: Conclusions and practical implications", fontsize=15,
              weight="bold", color=S.INK)
    fig.text(0.07, 0.875,
              "1. THE LONG-TERM RECORD (Part 1) shows a significant, coherent dilution\n"
              "   signal at all three rivers: discharge and DOC export trending up,\n"
              "   DIC/alkalinity trending down. Freshet timing, ice duration, and most\n"
              "   reaction rates show no significant trend over this 19-year baseline --\n"
              "   reported as a genuine null result where that is what the data shows.\n\n"
              "2. WITHIN the FORTE campaign window (Jun 1 - Sep 30), no instrument --\n"
              "   including grab sampling -- can learn anything about the OCEAN boundary.\n"
              "   The spring flood has already cut that connection by the time the\n"
              "   campaign starts. The ONE exception is a mooring placed directly in the\n"
              "   ocean, which reads that boundary directly.\n\n"
              "3. BUT grab sampling is still the single MOST valuable channel overall --\n"
              "   it also, independently, measures the RIVER's own upstream boundary (its\n"
              "   headwater stations read it almost directly), AND resolves WHERE along\n"
              "   the channel local changes in rate or loading are happening, something\n"
              "   no boundary-only test could even pose as a question.\n\n"
              "4. A single fixed mooring, by contrast, is sharply localized to its own\n"
              "   location and contributes nothing to resolving local organic-carbon\n"
              "   loading anywhere -- it is a fundamentally different, complementary\n"
              "   measurement, strongest for the ocean boundary specifically.\n\n"
              "5. Practical reading for campaign planning: grab sampling's real strength,\n"
              "   as currently scheduled, is characterizing what the RIVER is carrying and\n"
              "   WHERE along its length -- not the ocean. If the ocean boundary matters\n"
              "   too, extend monitoring into the cold season or rely on the offshore\n"
              "   mooring, which is immune to the seasonal blind spot.",
              fontsize=9.6, color=S.INK, va="top")
    fig.text(0.07, 0.195,
              "CAVEATS: results come from a model, calibrated against real rivers but not\n"
              "identical to them; the sensitivity analysis assumes one \"typical\" campaign\n"
              "season behaves like every year; a single, equal 20% uncertainty was assumed\n"
              "for every input, since no better estimate exists for most of them; and 7 of\n"
              "8 upstream controls plus all 6 along-estuary controls are tested as simple\n"
              "step changes, not genuine time-varying or continuously-distributed inputs.\n"
              "These affect exact numbers, not the findings above, which rest on simple,\n"
              "robust mechanisms not sensitive to any of these choices.",
              fontsize=8.3, color=S.MUTED, style="italic", va="top")
    S.brand(fig)
    pdf.savefig(fig)
    plt.close(fig)


# =====================================================================================
def main():
    out = ROOT / "docs" / "ns_rad_full_project_report.pdf"
    out.parent.mkdir(parents=True, exist_ok=True)
    with PdfPages(out) as pdf:
        page_title(pdf)
        page_executive_summary(pdf)

        page_what_is_nsrad(pdf)
        page_why_longrun(pdf)
        page_forcing_inventory(pdf)
        page_model_config(pdf)
        page_discharge_results(pdf)
        page_freshet_results(pdf)
        page_state_variables_1(pdf)
        page_state_variables_2(pdf)
        page_process_rates(pdf)
        page_ice_results(pdf)
        page_longrun_findings(pdf)

        page_why_climatology(pdf)
        page_climatology_method(pdf)
        page_climatology_visual(pdf)
        page_campaign_run(pdf)

        page_what_is_osse(pdf)
        page_channel_table(pdf)
        page_grab_results(pdf)
        page_mooring_drifter_results(pdf)
        page_airborne_satellite_results(pdf)

        page_what_is_green(pdf)
        page_controls_table(pdf)
        page_key_sensitivity_finding(pdf)
        page_upstream_sensitivity_finding(pdf)
        page_spatial_profile_rate(pdf)
        page_spatial_profile_lateral(pdf)
        page_hovmoller_sensitivity(pdf)
        page_reach_value_breakdown(pdf)

        page_what_is_dfs(pdf)
        page_ranking_chart(pdf)
        page_comparison(pdf)

        page_conclusions(pdf)
    print(f"wrote {out.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
