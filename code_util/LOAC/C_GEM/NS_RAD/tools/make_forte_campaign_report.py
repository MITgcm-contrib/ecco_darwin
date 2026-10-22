"""
Standalone, from-first-principles report on the FORTE campaign framework: the
climatological "campaign year" run, the synthetic observation channels (the
OSSE), the Green's-function sensitivity analysis, and the final observation-value
ranking. Written for a reader with NO prior context on this project -- every
technical term (OSSE, Green's function, degrees of freedom for signal, climatology,
upstream/downstream boundary) is explained in plain language before it is used.

Reads: forcing/*_climatology*.csv, runs/forte_campaign/<site>/output.nc (width
profile only), docs/osse_forte/<site>_<channel>.csv, runs/forte_green/
green_functions_<site>.npz, docs/osse_forte/observation_value_ranking*.csv.
No model run of its own. Companion to (but completely independent of)
docs/osse/osse_sampling_report.pdf, the earlier analog-year version of this
framework -- see tools/make_osse_report.py.

Usage:  python tools/make_forte_campaign_report.py
Writes: docs/osse_forte/forte_campaign_report.pdf
"""
import csv
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

import nsrad_style as S

ROOT = Path(__file__).resolve().parent.parent
FORC = ROOT / "forcing"
NATURE = ROOT / "runs" / "forte_campaign"
GREEN = ROOT / "runs" / "forte_green"
OUT = ROOT / "docs" / "osse_forte"
SITES = ["colville", "kuparuk", "sagavanirktok"]

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
N_MARINE_CTRL = 11   # discharge + clb_T + 9 marine species -- the ORIGINAL FORTE
                      # control set, before upstream (cub_*) controls were added
N_BOUNDARY_CTRL = 19  # + the 8 upstream controls, before the 6 reach-local
                      # (interior) controls were added -- see page_comparison

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
    p = OUT / f"{site}_{channel}.csv"
    if not p.exists():
        return []
    with open(p) as f:
        return list(csv.DictReader(f))


def _width_profile(site):
    from netCDF4 import Dataset
    d = Dataset(NATURE / site / "output.nc")
    return np.asarray(d.variables["width"][0, 1:])


def _ranking_combined():
    p = OUT / "observation_value_ranking_combined.csv"
    with open(p) as f:
        return list(csv.DictReader(f))


def _reach_bounds(M):
    """Same coarse-thirds convention as build_green_perturbations_forte.py's
    _reach_bounds -- cell 1 is the seaward/mouth end, cell M the upstream/
    headwater end."""
    b1, b2 = M // 3, 2 * (M // 3)
    return {"mouth": (1, b1), "mid": (b1 + 1, b2), "head": (b2 + 1, M)}


# =====================================================================================
# PART 0 -- orientation
# =====================================================================================
def page_title(pdf):
    fig = plt.figure(figsize=(11, 8.5))
    fig.text(0.08, 0.76, "The FORTE Campaign Framework", fontsize=27, fontweight="bold",
              color=S.INK)
    fig.text(0.08, 0.675, "A plain-language guide to the climatological model run,\n"
                          "the simulated field campaign, the sensitivity analysis,\n"
                          "and what it all means for planning real measurements",
              fontsize=13, color=S.INK2)
    fig.text(0.08, 0.50,
              "WHO THIS IS FOR: anyone who wants to understand this analysis without\n"
              "prior familiarity with the NS-RAD model or the FORTE field campaign. No\n"
              "background is assumed -- every technical term is explained the first\n"
              "time it is used.",
              fontsize=10.5, color=S.INK, va="top")
    fig.text(0.08, 0.36,
              "THE ONE-SENTENCE QUESTION THIS REPORT ANSWERS:\n"
              "If the FORTE field campaign collects data the way it is currently\n"
              "planned, which of its five measurement types will actually teach us the\n"
              "most about how these rivers work -- and does the campaign's timing\n"
              "itself limit what any of them can find out?",
              fontsize=11, color=S.INK, va="top", weight="bold")
    fig.text(0.08, 0.14,
              "Built from: tools/build_forte_campaign_forcing.py, run_forte_campaign.sh,\n"
              "build_osse_observations_forte.py, build_green_perturbations_forte.py,\n"
              "build_green_functions_forte.py, rank_observation_value_forte.py --\n"
              "each documented in full in its own module docstring.",
              fontsize=9, color=S.MUTED, va="top")
    S.brand(fig)
    pdf.savefig(fig)
    plt.close(fig)


def page_executive_summary(pdf):
    rows = _ranking_combined()
    by_ch = {r["channel"]: float(r["dfs"]) for r in rows}
    ordered = sorted(by_ch, key=lambda c: -by_ch[c])

    fig = plt.figure(figsize=(11, 8.5))
    fig.text(0.07, 0.92, "Executive summary", fontsize=16, weight="bold", color=S.INK)
    fig.text(0.07, 0.855,
              "This report builds and analyzes a computer simulation of three North\n"
              "Slope Alaska rivers (Colville, Kuparuk, Sagavanirktok) and the planned\n"
              "FORTE Expedition field campaign that will sample them, to answer one\n"
              "practical question: which of the campaign's measurement types is most\n"
              "valuable, and are there any blind spots in how the campaign is timed?",
              fontsize=10, color=S.INK, va="top")

    fig.text(0.07, 0.68, "The result, ranked:", fontsize=11.5, weight="bold", color=S.INK)
    lines = [f"{i}. {CHANNEL_LABELS[ch]} -- score {by_ch[ch]:.2f} (out of a "
             f"maximum of {N_CTRL})" for i, ch in enumerate(ordered, 1)]
    fig.text(0.07, 0.63, "\n".join(lines), fontsize=10.5, color=S.INK, va="top",
              family="monospace")

    fig.text(0.07, 0.40, "The one finding that matters most:", fontsize=11.5,
              weight="bold", color=S.WARN)
    fig.text(0.07, 0.365,
              "Grab sampling should be nearly worthless in FORTE's June-September window:\n"
              "these rivers' connection to the OCEAN only matters in the cold, low-flow\n"
              "season, and by June 1 the spring flood has already flushed the river fresh,\n"
              "cutting off any path for ocean water to reach the campaign. That part is\n"
              "true, and every measurement type loses value from it equally.\n\n"
              "But that is only half the story. Grab sampling is ALSO sensitive to the\n"
              "RIVER's own upstream chemistry, AND to what is happening ALONG the channel\n"
              "in between its two ends -- local degradation rates, local organic-carbon\n"
              "inputs. Testing both (not done in an earlier version of this report) flips\n"
              "the conclusion: grab sampling comes out on TOP, teaching us about the river\n"
              "itself -- end to end, not just its edges. See Parts 3-4 for the mechanism.",
              fontsize=9.5, color=S.INK, va="top")

    fig.text(0.07, 0.045,
              "The rest of this report explains, from scratch, how the numbers above\n"
              "were derived -- so this conclusion can be checked, not just taken on faith.",
              fontsize=8.5, color=S.MUTED, style="italic", va="top")
    S.brand(fig)
    pdf.savefig(fig)
    plt.close(fig)


# =====================================================================================
# PART 1 -- the model and the campaign
# =====================================================================================
def page_what_is_nsrad(pdf):
    fig, ax = plt.subplots(figsize=(11, 8.5))
    ax.axis("off")
    ax.set_title("What is the model being used here?", loc="left", fontsize=15,
                  weight="bold", pad=14)
    fig.text(0.08, 0.80,
              "NS-RAD is a computer simulation of three Arctic river-estuary systems on\n"
              "Alaska's North Slope: the Colville, Kuparuk, and Sagavanirktok rivers. It\n"
              "represents each river as a one-dimensional channel from its tidal mouth\n"
              "upstream, tracking water flow, salinity, temperature, dissolved oxygen,\n"
              "carbon, and nutrients through time.\n\n"
              "Like any model of a real place, it needs to be told what is happening at\n"
              "its EDGES, since it cannot simulate the whole ocean or the whole river\n"
              "basin. These edge inputs are called BOUNDARY CONDITIONS, and there are\n"
              "two of them per river:\n\n"
              "  UPSTREAM boundary (\"cub\"): what is flowing INTO the river from\n"
              "    further inland -- how much water (discharge), and the chemistry of\n"
              "    that water (carbon, nutrients, etc.).\n\n"
              "  DOWNSTREAM / marine boundary (\"clb\"): what the OCEAN looks like right\n"
              "    at the river mouth -- its salinity, temperature, carbon, and\n"
              "    nutrient content, and the tide/storm-surge water-level signal.\n\n"
              "Everything this report investigates comes down to one question: HOW\n"
              "CONFIDENT can we be about these boundary inputs, and which real-world\n"
              "measurements would most improve that confidence?",
              fontsize=10.5, color=S.INK, va="top")
    S.brand(fig)
    pdf.savefig(fig)
    plt.close(fig)


PHASE_ROWS = [
    ("Phase 0", "Late April (TBD)", "--", "Under-ice",
     "NOT part of this analysis -- the source campaign plan itself marks this\ninfeasible (the rivers are frozen solid with negligible flow)"),
    ("Phase 1a", "Jun 1-5", "152-156", "Grab sampling",
     "Helicopter-accessed water samples; DOC, nutrients, DIC, alkalinity,\nchlorophyll-a"),
    ("Phase 1b", "Jun 6-20", "157-171", "Grab sampling + ice drone",
     "Continues grab sampling; separate ice-thickness drone survey (not\nanalyzed here)"),
    ("Phase 2", "Jul 1-12", "182-193", "Airborne + moorings deployed",
     "Aircraft-based optical imaging; fixed instruments go into the water"),
    ("Phase 3", "Jul 20-Aug 7", "201-219", "Airborne+drone + offshore moorings",
     "Continues airborne imaging; moorings stay in through September"),
]


def page_campaign_plan(pdf):
    fig, ax = plt.subplots(figsize=(14, 6.3))
    fig.subplots_adjust(top=0.84, bottom=0.22, left=0.02, right=0.98)
    ax.axis("off")
    col_labels = ["Phase", "Real dates\n(2027)", "Day-of-\nyear", "What happens",
                  "Detail"]
    tbl = ax.table(cellText=PHASE_ROWS, colLabels=col_labels, loc="upper center",
                    cellLoc="left", colWidths=[0.07, 0.10, 0.07, 0.20, 0.48])
    tbl.auto_set_font_size(False)
    tbl.set_fontsize(8.0)
    tbl.scale(1, 3.4)
    for (r, c), cell in tbl.get_celld().items():
        cell.set_edgecolor(S.GRID)
        if r == 0:
            cell.set_facecolor(S.ACCENT)
            cell.set_text_props(color="white", fontweight="bold")
    ax.set_title("The planned FORTE field campaign", loc="left", fontsize=14,
                  weight="bold", pad=14)
    fig.text(0.02, 0.15,
              "The campaign's real target is 2027 (and possibly 2028). This matters for\n"
              "everything that follows: the model cannot simulate 2027 itself (its weather/\n"
              "ocean input data only goes up to 2023), so every comparison in this report\n"
              "uses the campaign's CALENDAR DATES (day-of-year) applied to a model year,\n"
              "not literal 2027.",
              fontsize=7.8, color=S.MUTED, style="italic", va="top")
    S.brand(fig)
    pdf.savefig(fig)
    plt.close(fig)


def page_why_climatology(pdf):
    fig, axes = plt.subplots(3, 1, figsize=(11, 9.2), sharex=True)
    fig.subplots_adjust(top=0.68, bottom=0.09, hspace=0.4)
    fig.suptitle("Why a \"typical year\" instead of picking one real year?",
                  x=0.05, ha="left", fontsize=14, weight="bold")
    fig.text(0.05, 0.925,
              "The model has 19 years (2005-2023) of real river/ocean/weather data on file.\n"
              "Rather than picking one specific year to stand in for \"what will happen in\n"
              "2027\" (somewhat arbitrary), this analysis builds a CLIMATOLOGY instead: for\n"
              "every input the model needs, it averages all 19 real years together, day-of-\n"
              "year by day-of-year, into one smooth, representative \"typical year\" -- reusable\n"
              "for any future campaign without re-deciding which real year to use.",
              fontsize=8.8, color=S.INK2, va="top")

    for ax, site in zip(axes, SITES):
        q = np.genfromtxt(open(FORC / f"{site}_river_discharge_climatology_m3sec.csv"),
                            delimiter=",")
        ax.plot(np.arange(365), q, color=S.RIVC[site], lw=1.2)
        ax.axvspan(152, 273, color=S.ACCENT, alpha=0.12)
        ax.set_ylabel(f"{S.LABEL[site]}\ndischarge\n[m3/s]", fontsize=8, color=S.RIVC[site])
        S.tidy(ax)
    axes[0].text(212, axes[0].get_ylim()[1] * 0.8, "FORTE campaign window\n(Jun 1 - Sep 30)",
                  ha="center", fontsize=7.5, color=S.ACCENT)
    axes[-1].set_xlabel("day of year")
    S.brand(fig)
    pdf.savefig(fig)
    plt.close(fig)


# =====================================================================================
# PART 2 -- the OSSE (synthetic sampling)
# =====================================================================================
def page_what_is_osse(pdf):
    fig, ax = plt.subplots(figsize=(11, 8.5))
    ax.axis("off")
    ax.set_title("What is an OSSE, and why build one before the real campaign?",
                  loc="left", fontsize=14, weight="bold", pad=14)
    fig.text(0.08, 0.80,
              "An OSSE -- Observing System Simulation Experiment -- is a standard\n"
              "technique used across climate and ocean science before an expensive field\n"
              "campaign happens. The idea:\n\n"
              "  1. Treat the computer model's own simulated state as a stand-in for\n"
              "     \"the truth\" (here, the climatological FORTE campaign-year run).\n\n"
              "  2. Simulate what each planned instrument would ACTUALLY measure if it\n"
              "     sampled that \"truth\" -- at the right times, the right places, with\n"
              "     the instrument's own realistic measurement error. This is done for\n"
              "     all five of FORTE's measurement types below.\n\n"
              "  3. Ask: if an instrument had only ever seen these synthetic\n"
              "     measurements, how much would it actually have learned? This is what\n"
              "     the Green's-function/ranking analysis later in this report answers.\n\n"
              "The payoff: problems with a campaign's DESIGN (not needing more funding\n"
              "or better instruments, but better TIMING or PLACEMENT) can be found and\n"
              "fixed before a single real measurement is taken.",
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
     "Water color/clarity --proxies for\norganic carbon, algae, sediment",
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
    col_labels = ["Measurement\ntype", "Where", "How often", "What it measures",
                  "In plain terms"]
    tbl = ax.table(cellText=CHANNEL_ROWS, colLabels=col_labels, loc="upper center",
                    cellLoc="left", colWidths=[0.12, 0.17, 0.15, 0.28, 0.28])
    tbl.auto_set_font_size(False)
    tbl.set_fontsize(8.0)
    tbl.scale(1, 3.4)
    for (r, c), cell in tbl.get_celld().items():
        cell.set_edgecolor(S.GRID)
        cell.PAD = 0.02
        if r == 0:
            cell.set_facecolor(S.ACCENT)
            cell.set_text_props(color="white", fontweight="bold")
    ax.set_title("The five measurement types this report evaluates", loc="left",
                  fontsize=14, weight="bold", pad=14)
    S.brand(fig)
    pdf.savefig(fig)
    plt.close(fig)


def page_grab_results(pdf):
    fig, axes = plt.subplots(2, 1, figsize=(11, 8.5))
    fig.subplots_adjust(hspace=0.5, top=0.88)
    fig.suptitle("What the synthetic grab samples actually look like",
                  x=0.05, ha="left", fontsize=13, weight="bold")
    ax0, ax1 = axes
    for site in SITES:
        rows = _read_csv(site, "grab")
        if not rows:
            continue
        stations = sorted(set(float(r["distance_km"]) for r in rows))
        ax0.scatter(stations, [S.LABEL[site]] * len(stations), color=S.RIVC[site],
                     s=40, zorder=3)
    ax0.set_xlabel("distance from the river mouth [km]")
    ax0.set_title("Where the 5 sampling stations sit along each river", loc="left",
                   fontsize=9.5)
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
    fig.suptitle("Moorings and drifters: what continuous and moving sensors see",
                  x=0.05, ha="left", fontsize=13, weight="bold")
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
                   "steadily carried out to sea by the river's current", loc="left",
                   fontsize=9.5)
    S.tidy(ax1)
    S.brand(fig)
    pdf.savefig(fig)
    plt.close(fig)


def page_airborne_satellite_results(pdf):
    fig, axes = plt.subplots(1, 2, figsize=(11, 5.5))
    fig.subplots_adjust(top=0.80, bottom=0.22, wspace=0.3)
    fig.suptitle("Airborne imaging and the notional satellite channels", x=0.05,
                  ha="left", fontsize=13, weight="bold")

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
# PART 3 -- Green's functions (sensitivity), explained from scratch
# =====================================================================================
def page_what_is_green(pdf):
    fig, ax = plt.subplots(figsize=(11, 8.5))
    ax.axis("off")
    ax.set_title("What is a \"Green's function,\" and why does this analysis need one?",
                  loc="left", fontsize=14, weight="bold", pad=14)
    fig.text(0.08, 0.80,
              "The real goal of this whole exercise is to figure out which measurements\n"
              "would most improve our confidence in the model's BOUNDARY INPUTS (page 3)\n"
              "-- the river discharge and ocean chemistry the model can't predict on its\n"
              "own, only be told about.\n\n"
              "But a measurement inside the river doesn't directly tell you the boundary\n"
              "input's value -- it tells you the RESULT of that input after it has\n"
              "traveled, mixed, and reacted its way through the estuary. To work\n"
              "backwards from \"what did we measure\" to \"what does that imply about the\n"
              "boundary,\" we need to know HOW SENSITIVE the model's interior state is to\n"
              "each boundary input -- i.e., if the ocean's carbon content were 10% higher\n"
              "than assumed, how much would that actually change what a sensor sitting\n"
              "inside the river would read, at every place and every time of year?\n\n"
              "That sensitivity relationship is called a GREEN'S FUNCTION. It is computed\n"
              "here the most direct way possible: run the model once normally, run it\n"
              "again with ONE input turned up by 10%, and measure the difference. Doing\n"
              "this once for each of 19 boundary inputs -- both DOWNSTREAM/ocean (clb)\n"
              "and UPSTREAM/riverine (cub), see page 3 -- for each of the 3 rivers, gives\n"
              "a complete map of \"what matters, where, and when.\"\n\n"
              "The same idea extends naturally to unknowns that are NOT at either boundary:\n"
              "a local in-stream reaction rate, or a local source of organic matter entering\n"
              "ALONG the channel (not just at its two ends), can be turned up by 10% in one\n"
              "stretch of river at a time and the same difference-and-measure recipe applied.\n"
              "This adds 6 more REACH-LOCAL inputs (page 12) -- the same method, just moved\n"
              "inside the domain instead of staying at its edges.",
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
    tbl = ax.table(cellText=rows, colLabels=["The 25 inputs tested", "Which\nboundary",
                                               "In plain terms"],
                    loc="upper center", cellLoc="left", colWidths=[0.30, 0.15, 0.50])
    tbl.auto_set_font_size(False)
    tbl.set_fontsize(7.8)
    tbl.scale(1, 1.75)
    for (r, c), cell in tbl.get_celld().items():
        cell.set_edgecolor(S.GRID)
        if r == 0:
            cell.set_facecolor(S.ACCENT)
            cell.set_text_props(color="white", fontweight="bold")
        elif 12 <= r <= 19:   # the 8 upstream rows, added to test the real campaign's
                               # own grab-sampling measurements (see Part 3, Finding 2)
            cell.set_facecolor("#eef6ee")
        elif r >= 20:         # the 6 reach-local (interior) rows, added to test
                               # ALONG-ESTUARY sensitivity (see page 15 onward)
            cell.set_facecolor("#eaf1fb")
    ax.set_title("What was actually tested", loc="left", fontsize=14, weight="bold", pad=12)
    fig.text(0.02, 0.075,
              "Each was increased by 10% for an entire run (LOCAL ones, only within their own\n"
              "third of the channel). UPSTREAM rows (green) match what the campaign's grab\n"
              "sampling measures directly. LOCAL rows (blue, new here) test whether spatially-\n"
              "distributed sampling can also resolve conditions INSIDE the channel, not just\n"
              "at its two ends -- see Part 3's along-estuary figures.",
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
    fig.suptitle("Finding 1 of 2: the ocean boundary goes dark", x=0.05,
                  ha="left", fontsize=14, weight="bold", color=S.WARN)
    fig.text(0.05, 0.895,
              f"{S.LABEL[site]}: how much a 10% increase in the ocean's dissolved-CO2\n"
              "level would change the river's own dissolved CO2, right at the mouth,\n"
              "by day of year.", fontsize=9.5, color=S.INK2, va="top")

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
              "The ocean's influence is strong all winter (the river is slow and salty,\n"
              "so ocean water sits inside it) and collapses to essentially ZERO right\n"
              "around day 140-150 -- just before the campaign begins. The spring flood\n"
              "(\"freshet\") has, by then, flushed the estuary completely fresh, cutting off\n"
              "any path for ocean water -- or anything dissolved in it -- to reach anywhere\n"
              "the campaign samples. This is real: no instrument INSIDE the river can learn\n"
              "anything about the ocean boundary during FORTE's season. But it is only half\n"
              "the story -- see the next page for the boundary that does NOT go dark.",
              fontsize=9, color=S.INK, va="top")
    S.brand(fig)
    pdf.savefig(fig)
    plt.close(fig)


def page_upstream_sensitivity_finding(pdf):
    site = "colville"
    green = np.load(GREEN / f"green_functions_{site}.npz")
    g_cub = green["DIC__cub_DIC"]
    g_clb = green["DIC__clb_DIC"]
    # near-boundary grab station (see docs/osse_forte/colville_grab.csv: station at
    # distance_km=3.9 of EL=4.15 km -- cell = round(3.9/0.1) = 39)
    cell_upstream = 39

    fig, ax = plt.subplots(figsize=(11, 7.5))
    fig.subplots_adjust(top=0.78, bottom=0.26)
    fig.suptitle("Finding 2 of 2: the RIVER boundary does not go dark", x=0.05,
                  ha="left", fontsize=14, weight="bold", color="#1e7a4a")
    fig.text(0.05, 0.895,
              f"{S.LABEL[site]}: the SAME sensitivity test, but now for a 10% increase in\n"
              "the RIVER's own upstream dissolved-CO2 level, measured at the grab-sampling\n"
              "station closest to the upstream boundary (3.9 of this river's 4.15 km).",
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
              "FORTE window, because this station sits only 250 m from the boundary itself:\n"
              "there is almost no estuary mixing or chemistry left to erase the signal\n"
              "before it reaches the sensor. This is why grab sampling -- despite learning\n"
              "nothing about the ocean boundary during this season (Finding 1) -- still\n"
              "ends up the single most valuable channel overall (Part 4): it is directly\n"
              "measuring the river boundary, not inferring it indirectly.",
              fontsize=9, color=S.INK, va="top")
    S.brand(fig)
    pdf.savefig(fig)
    plt.close(fig)


REACH_COLOR = {"mouth": "#2a78d6", "mid": "#eb9e34", "head": "#1e7a4a"}
REACH_LABEL = {"mouth": "mouth third", "mid": "middle third", "head": "head third"}


def page_spatial_profile_rate(pdf):
    site = "colville"
    green = np.load(GREEN / f"green_functions_{site}.npz")
    M = green["TOC__rate_aerdeg_mouth"].shape[1]
    dist = np.arange(M) * 0.1
    bounds = _reach_bounds(M)

    fig, ax = plt.subplots(figsize=(11, 7.5))
    fig.subplots_adjust(top=0.80, bottom=0.26)
    fig.suptitle("Along-estuary sensitivity: a local RATE perturbation stays local",
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
              "three thirds apart (see the reach-by-reach ranking two pages on).",
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
    fig.suptitle("Along-estuary sensitivity: a local LOADING perturbation, same test",
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
              "lower in the ranking ahead (page 18) because its pattern is more entangled\n"
              "with discharge (estimated at the same time), not because it is weaker.",
              fontsize=9, color=S.INK, va="top")
    S.brand(fig)
    pdf.savefig(fig)
    plt.close(fig)


def page_hovmoller_sensitivity(pdf):
    site = "colville"
    green = np.load(GREEN / f"green_functions_{site}.npz")
    g = green["TOC__rate_aerdeg_mid"][:274]   # (274 days, M cells)
    M = g.shape[1]
    dist = np.arange(M) * 0.1
    bounds = _reach_bounds(M)

    fig, ax = plt.subplots(figsize=(11, 8.3))
    fig.subplots_adjust(top=0.78, bottom=0.24, right=0.92)
    fig.suptitle("The full picture: sensitivity across BOTH space and time at once",
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
    ax.text(156, dist[-1] * 0.96, "FORTE campaign\nstarts here (Jun 1)",
             fontsize=8, color=S.INK, va="top")
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
    colors = {"grab": "#1e7a4a", "mooring_inshore": S.WARN, "airborne": S.ACCENT,
              "pace": S.MUTED}

    fig, axes = plt.subplots(1, 2, figsize=(12, 6.2))
    fig.subplots_adjust(top=0.78, bottom=0.26, wspace=0.3)
    fig.suptitle("Which channel resolves WHERE: reach-by-reach value", x=0.05,
                  ha="left", fontsize=14, weight="bold")
    fig.text(0.05, 0.89,
              "The same information as the ranking (Part 4), but split by reach instead of\n"
              "totaled -- this is what Part 3's along-estuary figures predict should happen.",
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
# PART 4 -- the ranking, explained from scratch
# =====================================================================================
def page_what_is_dfs(pdf):
    fig, ax = plt.subplots(figsize=(11, 8.5))
    ax.axis("off")
    ax.set_title("How \"which measurement is most valuable\" is actually scored",
                  loc="left", fontsize=14, weight="bold", pad=14)
    fig.text(0.08, 0.80,
              "For each synthetic measurement (Part 2), combined with how sensitive the\n"
              "model is to each boundary input at that place and time (Part 3), and how\n"
              "precise the instrument actually is, this analysis computes how much that\n"
              "single measurement would have narrowed down our uncertainty about each of\n"
              "the 25 inputs tested (page 12 -- ocean, river, and along-estuary). Adding\n"
              "this up correctly across every measurement a channel takes gives a single\n"
              "score per channel, called DEGREES OF FREEDOM FOR SIGNAL (DFS) -- a standard\n"
              "information-theory score used throughout atmospheric and ocean data\n"
              "assimilation.\n\n"
              "A score of 0 means a channel taught us nothing; a score of 25 (the number\n"
              "of inputs tested) would mean it fully pinned down every single one. In\n"
              "practice, scores land somewhere in between, and the split BETWEEN inputs\n"
              "(which ones got constrained, which didn't) matters as much as the total --\n"
              "a channel can score well by fully nailing down a FEW inputs (as grab\n"
              "sampling does, for the upstream and along-estuary ones) rather than\n"
              "partially touching all 25.\n\n"
              "ONE IMPORTANT CORRECTION APPLIED: a sensor recording every 6 minutes for\n"
              "months doesn't really deliver thousands of independent pieces of evidence\n"
              "-- consecutive readings are almost identical. This analysis explicitly\n"
              "corrects for that (an \"effective sample size\" adjustment), so a flood of\n"
              "repetitive readings is not mistaken for a flood of new information.",
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
    fig.suptitle("Final ranking: which measurement actually taught us the most",
                  x=0.05, ha="left", fontsize=13, weight="bold")
    S.tidy(ax, grid="x")
    fig.text(0.05, 0.135,
              "Grab sampling (highlighted) comes out on top -- not from the ocean boundary\n"
              "(Finding 1: cut off all season) but from the RIVER's own upstream boundary,\n"
              "measured almost directly by its near-headwater stations (Finding 2), PLUS a\n"
              "further boost from along-estuary (reach-local) sensitivity (Part 3, pages\n"
              "15-18). Airborne and PACE pick up smaller versions of both effects.",
              fontsize=8, color=S.MUTED, style="italic", va="top")
    S.brand(fig)
    pdf.savefig(fig)
    plt.close(fig)


def page_comparison(pdf):
    # Values from the original, full-year (2007/2008/2009 analog-year), 11-control
    # (marine+discharge only, no upstream) exercise, read directly from
    # docs/osse/observation_value_ranking_combined.csv so this table can never
    # silently drift from that report's own numbers.
    orig_path = ROOT / "docs" / "osse" / "observation_value_ranking_combined.csv"
    orig = {}
    if orig_path.exists():
        with open(orig_path) as f:
            for r in csv.DictReader(f):
                orig[r["channel"]] = float(r["dfs"])

    rows = _ranking_combined()
    forte_full = {r["channel"]: float(r["dfs"]) for r in rows}
    forte_marine_only = {r["channel"]: float(r["dfs_marine_discharge_only_11ctrl"])
                          for r in rows}
    forte_boundary_only = {r["channel"]: float(r["dfs_boundary_only_19ctrl"])
                            for r in rows}

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
        elif r == 1:  # grab sampling row
            cell.set_facecolor("#fdf0ec" if c == 2 else "#eef6ee" if c == 3
                                else "#eaf1fb" if c == 4 else "#ffffff")
    ax.set_title("Three different effects, disentangled: timing, upstream, along-estuary",
                  loc="left", fontsize=13, weight="bold", pad=12)
    fig.text(0.02, 0.21,
              "Reading grab sampling's row left to right: (A -> B) TIMING ALONE crashes its\n"
              "value 5.62 -> 1.00 -- the summer window cannot see the ocean boundary (Finding\n"
              "1). (B -> C) The 8 upstream controls restore it to 8.99/19 -- grab was ALWAYS\n"
              "also measuring the river boundary, column (B) just never asked (Finding 2).\n"
              "(C -> D) The 6 along-estuary controls lift it further to 12.81/25 -- grab's 5\n"
              "stations also resolve WHERE along the channel things happen (pages 15-18), a\n"
              "question boundary-only testing cannot even pose. Airborne shows the same\n"
              "3-step pattern; moorings/drifters/SWOT barely move in C or D, since they are\n"
              "not well-positioned to sense either structure.",
              fontsize=8.0, color=S.MUTED, style="italic", va="top")
    S.brand(fig)
    pdf.savefig(fig)
    plt.close(fig)


def page_conclusions(pdf):
    fig = plt.figure(figsize=(11, 8.5))
    fig.text(0.07, 0.94, "Conclusions and practical implications", fontsize=15,
              weight="bold", color=S.INK)
    fig.text(0.07, 0.875,
              "1. WITHIN the campaign window (Jun 1 - Sep 30), no instrument -- including\n"
              "   grab sampling -- can learn anything about the OCEAN boundary. The spring\n"
              "   flood has already cut off that connection by the time the campaign\n"
              "   starts, for every channel equally. The ONE exception is a mooring placed\n"
              "   directly in the ocean, which reads that boundary directly.\n\n"
              "2. But grab sampling is still the single MOST valuable channel overall --\n"
              "   it also, independently, measures the RIVER's own upstream boundary, and\n"
              "   its headwater-end stations read it almost directly (page 14). Airborne\n"
              "   and PACE pick up a smaller version via optical DOC/chlorophyll proxies.\n\n"
              "3. A THIRD effect, new here: grab and airborne also resolve WHERE along the\n"
              "   channel things happen, not just conditions at the two ends. A local\n"
              "   change in degradation rate or organic-carbon loading leaves a spatially\n"
              "   LOCALIZED signature (pages 15-18) that multi-station/continuous coverage\n"
              "   can tell apart reach by reach -- a single fixed mooring cannot, and is\n"
              "   blind to loading changes entirely (it measures no organic carbon).\n\n"
              "4. Practical reading: grab sampling's real strength is characterizing what\n"
              "   the RIVER is carrying and WHERE along its length -- not the ocean. If the\n"
              "   OCEAN boundary matters too, extend monitoring into the cold season (the\n"
              "   excluded \"Phase 0\") or rely on the offshore mooring.\n\n"
              "5. One limitation of this analysis, not the campaign: 7 of 8 upstream and\n"
              "   all 6 along-estuary controls are tested as simple step changes today,\n"
              "   not time-varying or continuously-distributed inputs -- a natural next\n"
              "   step once real campaign data can supply a finer picture.",
              fontsize=9.0, color=S.INK, va="top")
    fig.text(0.07, 0.145,
              "CAVEATS: these results come from a model, calibrated against real rivers but\n"
              "not identical to them; the analysis assumes one \"typical\" year behaves like\n"
              "every year; and a single, equal 20% uncertainty was assumed for every\n"
              "boundary input, since no better estimate exists for most of them. These\n"
              "affect the exact numbers, not the three main findings, which rest on simple,\n"
              "robust mechanisms (freshwater flushes out salt water; a sensor close to a\n"
              "boundary -- or reach -- reads it almost directly) not sensitive to these choices.",
              fontsize=8.5, color=S.MUTED, style="italic", va="top")
    S.brand(fig)
    pdf.savefig(fig)
    plt.close(fig)


# =====================================================================================
def main():
    out = OUT / "forte_campaign_report.pdf"
    OUT.mkdir(parents=True, exist_ok=True)
    with PdfPages(out) as pdf:
        page_title(pdf)
        page_executive_summary(pdf)

        page_what_is_nsrad(pdf)
        page_campaign_plan(pdf)
        page_why_climatology(pdf)

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
