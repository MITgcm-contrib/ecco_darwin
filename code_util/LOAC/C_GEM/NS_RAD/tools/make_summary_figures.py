"""
Summary figures for the NS-RAD four-river configuration.

Writes docs/ns_rad_model_summary.pdf (multi-page).

IMPORTANT: these describe the model SETUP -- forcings, geometry, mixing and data
provenance. They are NOT model results. No full simulation has been run, so nothing
here shows FCO2 or any state variable.

Usage:  python tools/make_summary_figures.py
"""

import math
import os
import sys
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import nsrad_style as _S  # shared NS-RAD publication/presentation style
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.patches import Rectangle

ROOT = Path(__file__).resolve().parent.parent
CODE = ROOT / "code"
FORC = ROOT / "forcing"
OUT = ROOT / "docs" / "ns_rad_model_summary.pdf"
sys.path.insert(0, str(CODE))

# Validated categorical palette: slots 1,2,3,7 of the reference theme.
# Checked with the skill validator, light mode, --pairs all:
#   worst CVD dE 9.2 (deutan), worst normal-vision dE 16.3, all checks PASS.
# Aqua sits below 3:1 contrast on the light surface, so the relief rule applies --
# every series carries a direct label, never colour alone.
SITES = ["colville", "kuparuk", "sagavanirktok", "canning"]
LABEL = {"colville": "Colville", "kuparuk": "Kuparuk",
         "sagavanirktok": "Sagavanirktok", "canning": "Canning"}
C = {"colville": "#2a78d6", "kuparuk": "#eb6834",
     "sagavanirktok": "#1baf7a", "canning": "#4a3aa7"}

INK, INK2, MUTED = "#0b0b0b", "#52514e", "#8a8880"
GRID, SURFACE = "#e6e5e0", "#fcfcfb"
WARN = "#e34948"

plt.rcParams.update({
    "figure.facecolor": SURFACE, "axes.facecolor": SURFACE,
    "savefig.facecolor": SURFACE,
    "font.size": 8, "axes.titlesize": 9, "axes.labelsize": 8,
    "axes.edgecolor": GRID, "axes.linewidth": 0.8,
    "axes.labelcolor": INK2, "text.color": INK,
    "xtick.color": INK2, "ytick.color": INK2,
    "xtick.labelsize": 7, "ytick.labelsize": 7,
    "legend.frameon": False, "legend.fontsize": 7.5,
    "grid.color": GRID, "grid.linewidth": 0.6,
    "lines.linewidth": 1.6, "lines.solid_capstyle": "round",
})
_S.apply()
_S.install_autoscale(1.2)  # embed fonts + presentation-size bump


def tidy(ax, grid_axis="y"):
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    for s in ("left", "bottom"):
        ax.spines[s].set_color(GRID)
    ax.grid(True, axis=grid_axis, alpha=0.7)
    ax.set_axisbelow(True)


def load(name):
    return np.genfromtxt(open(FORC / name, encoding="utf-8-sig"), delimiter=",")


MONTH0 = [0, 31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334]
MONTHL = list("JFMAMJJASOND")


def month_axis(ax):
    ax.set_xticks(MONTH0)
    ax.set_xticklabels(MONTHL)
    ax.set_xlim(0, 364)


def site_config(site):
    """Import config/init fresh for one site (module-level state, so re-import)."""
    for m in ("config", "variables", "init_module", "fun_module", "file_module"):
        sys.modules.pop(m, None)
    os.environ["CGEM_SITE"] = site
    import warnings
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        import config
        import init_module
    return config, init_module



def main():
    OUT.parent.mkdir(exist_ok=True)
    Q = {s: load(f"{s}_river_discharge_2022_m3sec.csv") for s in SITES}
    t_new = load("river_watertemp_2022_degC.csv")
    air = load("airtemp_2022_degC.csv")

    cfg = {}
    for s in SITES:
        c, im = site_config(s)
        cfg[s] = dict(
            B_lb=c.B_lb, B_ub=c.B_ub, L=c.L_FLARE, D=c.DEPTH_lb, EL=c.EL,
            M=c.M, DELXI=c.DELXI, chezy=c.Chezy_lb, dispmax=c.DISP_MAX,
            width=[im.width_at(i * c.DELXI) for i in range(c.M + 1)],
        )

    with PdfPages(OUT) as pdf:
        page_forcings(pdf, Q, t_new, air)
        page_geometry(pdf, cfg)
        page_dispersion(pdf, cfg, Q)
        page_provenance(pdf)

    print(f"wrote {OUT.relative_to(ROOT)}  ({OUT.stat().st_size/1024:.0f} kB, 4 pages)")


# ---------------------------------------------------------------- page 1
def page_forcings(pdf, Q, t_new, air):
    fig = plt.figure(figsize=(11, 8.5))
    fig.suptitle("NS-RAD — external forcings, 2022",
                 x=0.06, ha="left", fontsize=13, weight="bold")
    fig.text(0.06, 0.945, "Model configuration, not results. Per-river discharge and the "
             "modelled river-temperature forcing.", color=INK2, fontsize=8.5)

    gs = fig.add_gridspec(2, 1, height_ratios=[1.3, 1.0],
                          left=0.07, right=0.83, top=0.90, bottom=0.09, hspace=0.42)

    # (a) discharge
    ax = fig.add_subplot(gs[0])
    for s in SITES:
        ax.plot(Q[s], color=C[s], lw=1.5)
        pk = int(np.argmax(Q[s]))
        ax.plot([pk], [Q[s][pk]], "o", color=C[s], ms=4.5,
                mec=SURFACE, mew=1.2, zorder=5)
    ax.set_yscale("log")
    ax.set_ylabel("discharge  (m³ s⁻¹, log)")
    month_axis(ax)
    tidy(ax)
    ends = sorted(SITES, key=lambda s: -Q[s][-40:].mean())
    lo, hi = ax.get_ylim()
    slots = np.logspace(np.log10(hi * 0.55), np.log10(lo * 3.0), len(ends))
    for slot, s in zip(slots, ends):
        yv = Q[s][-40:].mean()
        ax.annotate(LABEL[s], xy=(364, yv), xytext=(378, slot), color=C[s],
                    fontsize=8, va="center", annotation_clip=False, weight="medium",
                    arrowprops=dict(arrowstyle="-", color=C[s], lw=0.7, alpha=0.45,
                                    shrinkA=0, shrinkB=2))
    ax.set_title("Daily discharge — dots mark the annual peak (the spring freshet)",
                 loc="left", color=INK)

    # (b) temperature
    ax = fig.add_subplot(gs[1])
    ax.axhline(0, color=GRID, lw=0.8)
    ax.plot(air, color=MUTED, lw=1.0, label="air (NDBC PRDA2)")
    ax.plot(t_new, color="#2a78d6", lw=1.8, label="modelled river temperature")
    ax.set_ylabel("temperature  (°C)")
    month_axis(ax)
    tidy(ax)
    ax.legend(loc="upper left", ncol=1, labelcolor=INK2)
    ax.set_title("River-temperature forcing (from air-temperature regression, "
                 "constrained ≥ 0 °C)", loc="left", color=INK)

    pdf.savefig(fig)
    plt.close(fig)


# ---------------------------------------------------------------- page 2
def page_geometry(pdf, cfg):
    fig = plt.figure(figsize=(11, 8.5))
    fig.suptitle("Channel geometry — SWORD v17c widths, USGS survey depths",
                 x=0.06, ha="left", fontsize=13, weight="bold")
    fig.text(0.06, 0.945, "Width law changed from whole-domain exponential to "
             "flare + prismatic (config.WIDTH_MODEL).", color=INK2, fontsize=8.5)

    gs = fig.add_gridspec(2, 2, left=0.07, right=0.96, top=0.90, bottom=0.08,
                          hspace=0.38, wspace=0.22, height_ratios=[1, 1])

    # (a) width profiles
    ax = fig.add_subplot(gs[0, 0])
    for s in SITES:
        c = cfg[s]
        x = np.arange(c["M"] + 1) * c["DELXI"] / 1000
        ax.plot(x, c["width"], color=C[s], lw=1.8)
        ax.axvline(c["L"] / 1000, color=C[s], lw=0.8, ls=":", alpha=0.55)
        ax.annotate(LABEL[s], xy=(x[-1], c["width"][-1]), xytext=(x[-1] + 0.4, c["width"][-1]),
                    color=C[s], fontsize=7.5, va="center", annotation_clip=False)
    ax.set_yscale("log")
    ax.set_xlabel("distance upstream from mouth  (km)")
    ax.set_ylabel("channel width  (m, log)")
    ax.set_xlim(0, 30)
    tidy(ax)
    ax.set_title("Width: converges over the flare, then prismatic\n"
                 "(dotted = flare breakpoint L_FLARE)", loc="left", color=INK)

    # (b) flare vs old exponential, Colville
    ax = fig.add_subplot(gs[0, 1])
    c = cfg["colville"]
    x = np.arange(c["M"] + 1) * c["DELXI"] / 1000
    LC = 1.0 / -(1.0 / c["EL"] * math.log(c["B_ub"] / c["B_lb"]))
    expo = [c["B_lb"] * math.exp(-(i * c["DELXI"]) / LC) for i in range(c["M"] + 1)]
    ax.plot(x, expo, color=MUTED, lw=1.5, ls=(0, (4, 2)))
    ax.plot(x, c["width"], color=C["colville"], lw=2.0)
    ax.annotate("flare + prismatic\n(AIC-selected)", xy=(18, c["width"][90]),
                xytext=(13, c["width"][90] * 1.9), color=C["colville"], fontsize=7.5)
    ax.annotate("old: exponential\nover whole domain", xy=(20, expo[100]),
                xytext=(9.5, 1010), color=MUTED, fontsize=7.5,
                arrowprops=dict(arrowstyle="-", color=MUTED, lw=0.7, alpha=0.6))
    ax.set_xlabel("distance upstream from mouth  (km)")
    ax.set_ylabel("channel width  (m)")
    ax.set_xlim(0, 28)
    tidy(ax)
    ax.set_title("Colville — the two width laws compared", loc="left", color=INK)

    # (c) depth: old vs new
    ax = fig.add_subplot(gs[1, 0])
    xs = np.arange(4)
    ax.bar(xs, [cfg[s]["D"] for s in SITES], width=0.55,
           color=[C[s] for s in SITES], zorder=3)
    ax.axhline(15.0, color=WARN, lw=1.4, ls=(0, (4, 2)), zorder=4)
    ax.text(-0.42, 15.55, "shipped value: 15 m — cited to the Sagavanirktok gauge,"
            " whose own surveys median 0.88 m", color=WARN, fontsize=7,
            va="bottom", ha="left")
    for i, s in enumerate(SITES):
        ax.text(i, cfg[s]["D"] + 0.35, f"{cfg[s]['D']:.2f}", ha="center",
                fontsize=7.5, color=INK)
    ax.set_xticks(xs)
    ax.set_xticklabels([LABEL[s] for s in SITES], fontsize=7.5)
    ax.set_ylabel("channel depth  (m)")
    ax.set_ylim(0, 17)
    tidy(ax)
    ax.set_title("Depth from at-a-station hydraulic geometry D = c·Q^f\n"
                 "evaluated at each river's open-water mean discharge", loc="left", color=INK)

    # (d) geometry table
    ax = fig.add_subplot(gs[1, 1])
    ax.axis("off")
    rows = [["", "B_lb", "B_ub", "L_flare", "depth", "surveys"]]
    nsurv = {"colville": "208", "kuparuk": "323", "sagavanirktok": "189", "canning": "28"}
    for s in SITES:
        c = cfg[s]
        rows.append([LABEL[s], f"{c['B_lb']:.0f} m", f"{c['B_ub']:.0f} m",
                     f"{c['L']/1000:.1f} km", f"{c['D']:.2f} m", nsurv[s]])
    tb = ax.table(cellText=rows[1:], colLabels=rows[0], loc="center", cellLoc="right")
    tb.auto_set_font_size(False)
    tb.set_fontsize(7.5)
    tb.scale(1, 1.55)
    for (r, cix), cell in tb.get_celld().items():
        cell.set_edgecolor(GRID)
        cell.set_linewidth(0.6)
        if r == 0:
            cell.set_text_props(color=INK2, weight="bold")
            cell.set_facecolor(SURFACE)
        else:
            cell.set_facecolor(SURFACE)
            if cix == 0:
                cell.set_text_props(color=C[SITES[r - 1]], weight="medium")
                cell._loc = "left"
    ax.set_title("Configured geometry per site", loc="left", color=INK, y=0.98)
    ax.text(0, 0.06, "Sagavanirktok B_ub and L_flare are BORROWED from Canning —\n"
                     "SWORD shows it widening upstream (its own flare fit R² = −0.65).\n"
                     "Canning depth rests on only 28 surveys.",
            transform=ax.transAxes, fontsize=7, color=INK2)

    pdf.savefig(fig)
    plt.close(fig)


# ---------------------------------------------------------------- page 3
def page_dispersion(pdf, cfg, Q):
    fig = plt.figure(figsize=(11, 8.5))
    fig.suptitle("Longitudinal dispersion — Seo & Cheong (1998) river form",
                 x=0.06, ha="left", fontsize=13, weight="bold")
    fig.text(0.06, 0.945, "Computed per timestep from local width, depth and velocity. "
             "The tide-dominated estuary formulation collapses to zero under these "
             "shallow observation-based geometries.", color=INK2, fontsize=8.5)

    gs = fig.add_gridspec(2, 2, left=0.07, right=0.96, top=0.90, bottom=0.08,
                          hspace=0.38, wspace=0.22)

    # open-water mean discharge (ice gate open days 140-303) for the dispersion estimate
    ow = np.zeros(365, bool)
    ow[140:304] = True
    Qow = {s: float(Q[s][ow].mean()) for s in SITES}

    def seo(W, H, chezy, q):
        U = q / (W * H)
        us = math.sqrt(9.81) * U / chezy
        return 5.915 * (W / H) ** 0.620 * (U / us) ** 1.428 * H * us

    def savenije(c, q):
        EL, B_lb, D_lb, DELXI = c["EL"], c["B_lb"], c["D"], c["DELXI"]
        LC = 1.0 / -(1.0 / EL * math.log(c["B_ub"] / B_lb))
        K = 4.38 * D_lb ** 0.36 * B_lb ** -0.21 * LC ** -0.14
        N = math.pi * q / (D_lb * B_lb)
        D0 = 26 * math.sqrt(N * 9.81) * D_lb ** 1.5
        beta = (K * LC * q) / (D0 * B_lb * D_lb)
        return [max(D0 * (1 - beta * (math.exp((i * DELXI) / LC) - 1)), 0.0)
                for i in range(c["M"] + 1)]

    # (a) new dispersion
    ax = fig.add_subplot(gs[0, 0])
    for s in SITES:
        c = cfg[s]
        x = np.arange(c["M"] + 1) * c["DELXI"] / 1000
        k = [seo(max(w, 1e-3), c["D"], c["chezy"], Qow[s]) for w in c["width"]]
        ax.plot(x, k, color=C[s], lw=1.8)
        ax.annotate(LABEL[s], xy=(x[-1], k[-1]), xytext=(x[-1] + 0.4, k[-1]),
                    color=C[s], fontsize=7.5, va="center", annotation_clip=False)
    ax.axhline(cfg["colville"]["dispmax"], color=WARN, lw=1.2, ls=(0, (4, 2)))
    ax.text(0.4, cfg["colville"]["dispmax"] * 1.06,
            f"DISP_MAX = {cfg['colville']['dispmax']:.0f} m² s⁻¹  (numerical ceiling)",
            color=WARN, fontsize=7)
    ax.set_xlabel("distance upstream from mouth  (km)")
    ax.set_ylabel("dispersion K  (m² s⁻¹, log)")
    ax.set_xlim(0, 30)
    ax.set_yscale("log")
    ax.set_ylim(200, 3000)
    tidy(ax)
    ax.set_title("NEW — Seo & Cheong (1998), from local W, H, U\n"
                 "non-zero across the whole domain, well under the ceiling",
                 loc="left", color=INK)

    # (b) old dispersion collapse
    ax = fig.add_subplot(gs[0, 1])
    rowtxt = []
    for s in SITES:
        c = cfg[s]
        x = np.arange(c["M"] + 1) * c["DELXI"] / 1000
        d = savenije(c, Qow[s])
        ax.plot(x, d, color=C[s], lw=1.8)
        nz = sum(1 for v in d if v > 0)
        rowtxt.append((s, nz * c["DELXI"] / 1000))
    ymax = ax.get_ylim()[1]
    ax.text(9.5, ymax * 0.95, "dispersion reaches zero after:", fontsize=7.5,
            color=INK2, va="top")
    for k, (s, km) in enumerate(sorted(rowtxt, key=lambda r: -r[1])):
        ax.plot([10.1], [ymax * (0.845 - 0.075 * k)], "s", color=C[s], ms=5,
                clip_on=False)
        ax.text(10.9, ymax * (0.845 - 0.075 * k),
                f"{LABEL[s]} — {km:.1f} km", fontsize=7.5, color=INK2, va="center")
    ax.set_xlabel("distance upstream from mouth  (km)")
    ax.set_ylabel("dispersion K  (m² s⁻¹)")
    ax.set_xlim(0, 30)
    tidy(ax)
    ax.set_title("OLD — Van der Burgh / Savenije, same geometry\n"
                 "clamped to zero within a few hundred metres", loc="left", color=INK)

    # (c) formula comparison
    ax = fig.add_subplot(gs[1, 0])
    names = ["Fischer\n(1979)", "Seo & Cheong\n(1998)"]
    fisch, seoc, wh = [], [], []
    for s in SITES:
        c = cfg[s]
        W, H, q = c["B_lb"], c["D"], Qow[s]
        U = q / (W * H)
        us = math.sqrt(9.81) * U / c["chezy"]
        fisch.append(0.011 * U ** 2 * W ** 2 / (H * us))
        seoc.append(seo(W, H, c["chezy"], q))
        wh.append(W / H)
    xs = np.arange(4)
    ax.bar(xs - 0.19, fisch, 0.34, color=MUTED, zorder=3, label=names[0])
    ax.bar(xs + 0.19, seoc, 0.34, color=[C[s] for s in SITES], zorder=3, label=names[1])
    ax.axhline(cfg["colville"]["dispmax"], color=WARN, lw=1.2, ls=(0, (4, 2)))
    ax.text(3.45, cfg["colville"]["dispmax"] * 1.5, "DISP_MAX", color=WARN, fontsize=7,
            ha="right")
    ax.set_yscale("log")
    ax.set_xticks(xs)
    ax.set_xticklabels([f"{LABEL[s]}\nW/H={w:.0f}" for s, w in zip(SITES, wh)], fontsize=7)
    ax.set_ylabel("K at the mouth  (m² s⁻¹, log)")
    tidy(ax)
    ax.legend(loc="upper right", labelcolor=INK2, handlelength=1.2)
    ax.set_title("Why Seo & Cheong, not Fischer\n"
                 "Fischer assumes W/H ≈ 10–100 and blows past the ceiling on wide channels",
                 loc="left", color=INK)

    # (d) notes
    ax = fig.add_subplot(gs[1, 1])
    ax.axis("off")
    ax.set_title("What changed and why", loc="left", color=INK, y=0.97)
    txt = (
        "Shear velocity from the model's own friction closure:\n"
        "    u*  =  √g · U / Chezy\n"
        "so no channel slope is needed — SWORD's slopes at these\n"
        "delta reaches are unusable (zeros, and values as absurd\n"
        "as 2.5 m/m). Since U/u* = Chezy/√g identically, the\n"
        "velocity ratio is set by friction alone.\n\n"
        "Numerical ceiling. schemes_module.disp_sch builds\n"
        "    r[i] = 1 − (c1+c2)·DELTI/(8·DELXI²)\n"
        "and r < 0 makes the Crank–Nicolson right-hand side\n"
        "oscillatory, so K ≤ 4·DELXI²/DELTI = 2133 m² s⁻¹ at the\n"
        "shipped DELTI = 75 s, DELXI = 200 m. The cap is enforced\n"
        "in code and reports once per run. It does not currently\n"
        "bind for any site.\n\n"
        "Root cause. C-GEM's width law and its dispersion closure\n"
        "both derive from Savenije's tide-dominated alluvial-\n"
        "estuary theory. This configuration has AMPL = 0.2 m and\n"
        "distance = 1 — negligible tides, essentially no saline\n"
        "intrusion. These are rivers with short delta flares, not\n"
        "funnel estuaries, so the framework was being applied\n"
        "outside its regime. The width misfit and the collapsing\n"
        "dispersion were one problem surfacing twice."
    )
    ax.text(0, 0.90, txt, transform=ax.transAxes, fontsize=7.4, va="top",
            color=INK2, linespacing=1.55, family="DejaVu Sans")

    pdf.savefig(fig)
    plt.close(fig)


# ---------------------------------------------------------------- page 4
def page_provenance(pdf):
    fig = plt.figure(figsize=(11, 8.5))
    fig.suptitle("Data provenance and confidence", x=0.06, ha="left",
                 fontsize=13, weight="bold")
    fig.text(0.06, 0.945, "How well constrained each input is, per river. "
             "Read this before comparing sites.", color=INK2, fontsize=8.5)

    # 0 observed, 1 observed-with-caveat, 2 modelled/regional, 3 reconstructed/borrowed
    GRADE = {
        "Discharge":        {"colville": 1, "kuparuk": 0, "sagavanirktok": 1, "canning": 3},
        "Water temp.":      {"colville": 2, "kuparuk": 2, "sagavanirktok": 2, "canning": 2},
        "Width (SWORD)":    {"colville": 0, "kuparuk": 0, "sagavanirktok": 3, "canning": 0},
        "Depth (USGS)":     {"colville": 1, "kuparuk": 0, "sagavanirktok": 1, "canning": 1},
        "Boundary chem.":   {"colville": 3, "kuparuk": 3, "sagavanirktok": 3, "canning": 3},
    }
    NOTE = {
        "Discharge": {"colville": "gauged 400 km upstream", "kuparuk": "gauge near tidewater",
                      "sagavanirktok": "gauge ~130 km inland", "canning": "no 2022 record — donor-scaled"},
        "Water temp.": {k: "regional air–water fit" for k in SITES},
        "Width (SWORD)": {"colville": "104 nodes", "kuparuk": "149 nodes",
                          "sagavanirktok": "borrowed shape", "canning": "136 nodes"},
        "Depth (USGS)": {"colville": "208 surveys, upstream", "kuparuk": "323 surveys",
                         "sagavanirktok": "189 surveys, upstream", "canning": "only 28 surveys"},
        "Boundary chem.": {k: "placeholder — unsourced" for k in SITES},
    }
    # sequential-by-severity: one hue light->dark would imply order; these are
    # categories of provenance, so use the reserved status ramp semantics instead.
    FILL = {0: "#1baf7a", 1: "#eda100", 2: "#eb6834", 3: "#e34948"}
    KEY = {0: "observed", 1: "observed, with caveat", 2: "modelled / regional",
           3: "reconstructed or borrowed"}
    # distinct in-cell word per grade -- the colour must never be the only cue
    SHORT = {0: "observed", 1: "observed*", 2: "modelled", 3: "inferred"}

    ax = fig.add_axes([0.20, 0.30, 0.70, 0.54])
    rows = list(GRADE)
    for r, k in enumerate(rows):
        for c_, s in enumerate(SITES):
            g = GRADE[k][s]
            ax.add_patch(Rectangle((c_ + 0.03, r + 0.03), 0.94, 0.94,
                                   color=FILL[g], alpha=0.90, lw=0))
            ax.text(c_ + 0.5, r + 0.62, SHORT[g], ha="center", va="center",
                    fontsize=7.5, color="white", weight="bold")
            ax.text(c_ + 0.5, r + 0.32, NOTE[k][s], ha="center", va="center",
                    fontsize=6.2, color="white", alpha=0.95)
    ax.set_xlim(0, 4)
    ax.set_ylim(0, len(rows))
    ax.set_xticks(np.arange(4) + 0.5)
    ax.set_xticklabels([LABEL[s] for s in SITES], fontsize=9)
    for lbl, s in zip(ax.get_xticklabels(), SITES):
        lbl.set_color(C[s])
    ax.set_yticks(np.arange(len(rows)) + 0.5)
    ax.set_yticklabels(rows, fontsize=8.5, color=INK)
    ax.invert_yaxis()
    ax.tick_params(length=0)
    for sp in ax.spines.values():
        sp.set_visible(False)
    ax.xaxis.tick_top()

    hs = [Rectangle((0, 0), 1, 1, color=FILL[i]) for i in range(4)]
    ax.legend(hs, [KEY[i] for i in range(4)], loc="upper left",
              bbox_to_anchor=(0, -0.06), ncol=4, labelcolor=INK2, handlelength=1.1)

    fig.text(0.06, 0.20,
             "Weakest site — Canning: reconstructed discharge (donor-scaled from the Hulahula, no 2022 record), "
             "borrowed temperature, and a depth\nresting on only 28 surveys. Its geometry is observed, but "
             "everything driving it is inferred.\n\n"
             "Largest remaining gap — boundary chemistry: all four sites carry identical riverine DOC, alkalinity "
             "and nutrients. This is now the main\nun-sourced input, and it is the one most directly tied to the "
             "carbonate system.\n\n"
             "This document covers configuration and provenance only. For model output (fields, seasonal cycles, "
             "the prognostic ice cover) see\nns_rad_diagnostics.pdf; for model-vs-observation see "
             "ns_rad_validation.pdf; for channel geometry see ns_rad_geometry.pdf.",
             fontsize=8, color=INK2, linespacing=1.6, va="top")

    pdf.savefig(fig)
    plt.close(fig)


if __name__ == "__main__":
    main()
