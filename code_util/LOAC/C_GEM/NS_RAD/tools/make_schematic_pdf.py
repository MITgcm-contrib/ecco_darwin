"""
Model schematic: data inputs, the per-timestep core, and outputs, for NS-RAD.

Page 1  land-to-ocean framework: the geographic configuration -- headwaters -> reactive
        estuary -> Beaufort Sea, the boundary/forcing sources, and the four rivers'
        ACTUAL per-site geometry & discharge (read live from sites/<name>.py)
Page 2  data-flow overview: inputs (with provenance) -> model core -> outputs
Page 3  per-timestep call graph + transported-species registry + numerics/config

Pages 2-3 are structural diagrams; page 1's per-river numbers are read live from the
site modules and forcing files (so it reflects the real configuration). Hand-laid in
matplotlib on a 0-100 coordinate grid.

Usage:  python tools/make_schematic_pdf.py
"""
import json
import os
import sys
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import (FancyBboxPatch, FancyArrowPatch, Polygon, Rectangle,
                                Circle, PathPatch)
from matplotlib.path import Path as MplPath
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.text import Text as MplText
from matplotlib.backends.backend_pdf import PdfPages

ROOT = Path(__file__).resolve().parent.parent
CODE = ROOT / "code"
FORC = ROOT / "forcing"


def _perf():
    """Measured runtime metrics (tools/bench_timing.sh), or None if not yet benchmarked."""
    try:
        return json.load(open(ROOT / "docs" / "performance_metrics.json"))
    except Exception:
        return None
OUT = ROOT / "docs" / "ns_rad_model_schematic.pdf"

# --- publication design system -------------------------------------------------
# Refined, cohesive, colourblind-aware palette; muted domain accents with pale tints.
INK, INK2, MUTED = "#1a1a1a", "#5c5c5c", "#9a9a96"
GRID, SURFACE, WARN = "#dcdcd7", "#ffffff", "#c0392b"
# input-group accents + matching pale fills (validated categorical hues, softened)
HYD, GEO, MET, CHEM = "#3a6ea5", "#2e8b6f", "#cf7233", "#6a5aa0"
HYD_F, GEO_F, MET_F, CHEM_F = "#eef4fa", "#eaf5f0", "#fcefe4", "#f0edf8"
CORE, CORE_EC = "#f4f4f2", "#c7c7c1"      # neutral core-process box
ICEB, ICE_F, ICE_E = "#3f7bb0", "#dceaf4", "#9fc0da"
OCEAN, OCEAN_D = "#2f6f9e", "#244f70"
OUTC = "#1a1a1a"
MONO = "Menlo"
LAB = {"colville": "Colville", "kuparuk": "Kuparuk",
       "sagavanirktok": "Sagav.", "canning": "Canning"}

plt.rcParams.update({
    "figure.facecolor": SURFACE, "savefig.facecolor": SURFACE, "figure.dpi": 150,
    # DejaVu Sans: the one available sans with COMPLETE glyph coverage (subscripts,
    # superscripts, arrows, units) — so the whole figure stays in a single typeface
    # with no inconsistent per-glyph fallback. Embedded as TrueType below.
    "font.family": "sans-serif",
    "font.sans-serif": ["DejaVu Sans"],
    "font.size": 8, "text.color": INK,
    "pdf.fonttype": 42, "ps.fonttype": 42,   # embed TrueType (journal requirement)
    "axes.edgecolor": GRID, "patch.linewidth": 0.9,
})


def box(ax, x, y, w, h, title, lines=(), fc=CORE, ec=None, tc=INK, fs=8,
        title_w="bold", lw=0.9, title_c=None):
    ec = ec or CORE_EC
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.3,rounding_size=1.1",
                                fc=fc, ec=ec, lw=lw, mutation_aspect=0.6, zorder=2))
    ax.text(x + w / 2, y + h - 2.1, title, ha="center", va="top", fontsize=fs,
            weight=title_w, color=title_c or tc, zorder=3)
    for i, ln in enumerate(lines):
        ax.text(x + w / 2, y + h - 4.7 - i * 2.5, ln, ha="center", va="top",
                fontsize=fs - 1.7, color=INK2 if tc == INK else tc, zorder=3)


def arrow(ax, x1, y1, x2, y2, color=MUTED, lw=1.1, style="-|>"):
    ax.add_patch(FancyArrowPatch((x1, y1), (x2, y2), arrowstyle=style,
                                 mutation_scale=10, color=color, lw=lw,
                                 shrinkA=1, shrinkB=1,
                                 capstyle="round", joinstyle="round", zorder=4))


def header(ax, title, subtitle=None, accent=INK, fig_label=None):
    """Journal-style page header: title, a short accent rule, and a grey subtitle."""
    ax.text(2, 97.6, title, fontsize=16, weight="bold", color=INK, va="top")
    ax.plot([2.2, 15.5], [92.9, 92.9], color=accent, lw=2.6, solid_capstyle="round",
            zorder=5)
    if subtitle:
        ax.text(2, 91.6, subtitle, fontsize=8.7, color=INK2, va="top")
    if fig_label:
        ax.text(98, 97.4, fig_label, fontsize=8.2, color=MUTED, va="top", ha="right",
                style="italic")


def caption(ax, text):
    ax.plot([2, 98], [3.7, 3.7], color=GRID, lw=0.8)
    ax.text(2, 2.7, text, fontsize=7.0, color=MUTED, va="top")


def plabel(ax, x, y, letter):
    """Bold panel letter, journal convention."""
    ax.text(x, y, letter, fontsize=12, weight="bold", color=INK, va="top", ha="left")


def canvas(fig, rect=(0, 0, 1, 1)):
    ax = fig.add_axes(rect)
    ax.set_xlim(0, 100); ax.set_ylim(0, 100)
    ax.axis("off")
    return ax


# Presentation sizing. Every text size in this figure is set explicitly for a dense,
# journal-style layout; for slides we scale ALL of them up uniformly just before
# saving (the graphics stay put, only the type grows). PAGE_SCALE is per page —
# page 1 (the headline framework slide) can take the most; the denser reference
# pages take a bit less so nothing overflows. Override globally with CGEM_FONT_SCALE.
_ENV_SCALE = float(os.environ.get("CGEM_FONT_SCALE", "1.0"))
PAGE_SCALE = {"framework": 1.45, "dataflow": 1.28, "detail": 1.30}


def scale_fonts(fig, k):
    """Multiply the point size of every text object in the figure by k."""
    k *= _ENV_SCALE
    if k == 1.0:
        return
    for t in fig.findobj(match=MplText):
        t.set_fontsize(t.get_fontsize() * k)


# ------------------------------------------------------------------ page 1
def page_dataflow(pdf):
    fig = plt.figure(figsize=(14, 9))
    ax = canvas(fig)
    header(ax, "Data inputs, model core, and outputs", accent=OCEAN,
           fig_label="NS-RAD · Fig. 2")

    # ---------- INPUTS (top band) ----------
    ax.text(2, 89, "DATA INPUTS  (with provenance)", fontsize=9.5, weight="bold", color=INK2)
    iy, ih = 70, 16
    box(ax, 2, iy, 22, ih, "HYDROLOGY", fc=HYD_F, ec=HYD, tc=HYD, fs=8.5, lines=(
        "Discharge — USGS gauges, per river",
        "cfs→m³/s (0.3048³); Canning",
        "reconstructed from Hulahula ×2.97",
        "Geometry via SWORD v17c widths",
    ))
    box(ax, 26, iy, 22, ih, "CHANNEL GEOMETRY", fc=GEO_F, ec=GEO, tc=GEO, fs=8.5, lines=(
        "Width — SWORD v17c nodes / n_chan",
        "flare + prismatic (not exponential)",
        "Depth — USGS ADCP surveys,",
        "hydraulic geometry D = c·Q^f (1–2 m)",
    ))
    box(ax, 50, iy, 24, ih, "METEOROLOGY", fc=MET_F, ec=MET, tc=MET, fs=8.5, lines=(
        "Wind, solar, air T — NDBC PRDA2 buoy",
        "Sea T (marine bnd) — PRDA2 buoy",
        "River T (riverine bnd) — air-T",
        "regression + USGS 00010 obs (Sag/Kup)",
        "Humidity — Deadhorse Airport ISD",
    ))
    box(ax, 76, iy, 22, ih, "CHEMISTRY", fc=CHEM_F, ec=CHEM, tc=CHEM, fs=8.5, lines=(
        "Atmospheric pCO₂ — Barrow",
        "Kuparuk — Arctic LTER (ref. reach);",
        "Colville, Sag — WQP grabs;",
        "Canning — placeholder",
    ))

    # ---------- CORE (middle) ----------
    ax.text(2, 65, "MODEL CORE", fontsize=9.5, weight="bold", color=INK2)
    cy, ch, cw = 44, 15, 13.2
    step = 15.9
    xs = [2 + i * step for i in range(6)]   # 6 core processes, in loop order
    ICEBLUE = "#dcebf7"
    box(ax, xs[0], cy, cw, ch, "HYDRODYNAMICS", fc=CORE, ec=MUTED, fs=7.6, lines=(
        "hyd(): iterate", "coeff_a→tridag", "→conv→update",
        "~44 its/step (numba)"))
    box(ax, xs[1], cy, cw, ch, "DISPERSION", fc=CORE, ec=MUTED, fs=7.6, lines=(
        "river_dispersion()", "Seo & Cheong 1998", "from local W,H,U",
        "(replaces Savenije)"))
    box(ax, xs[2], cy, cw, ch, "TRANSPORT", fc=CORE, ec=MUTED, fs=7.6, lines=(
        "transport(): per-", "species TVD advec.", "+ dispersion",
        "13 fields incl. T"))
    box(ax, xs[3], cy, cw, ch, "HEAT BUDGET", fc=CORE, ec=MUTED, fs=7.6, lines=(
        "heat_budget():", "Q_sw+lw+sens+lat", "warms T; freeze",
        "energy→ice (skip ice)"))
    box(ax, xs[4], cy, cw, ch, "RIVER ICE", fc=ICEBLUE, ec="#2a6aa8", fs=7.6, lines=(
        "ice_step(): grow/", "melt ice_thickness,", "freshet breakup,",
        "bottom-fast cap"))
    box(ax, xs[5], cy, cw, ch, "BIOGEO + SED", fc=CORE, ec=MUTED, fs=7.6, lines=(
        "NPP, resp, nitrif,", "carbonate→pH, FCO₂,", "+CDOM photo, CH₄/N₂O,",
        "benthic; sed(): erode"))
    for i in range(len(xs) - 1):
        arrow(ax, xs[i] + cw, cy + ch / 2, xs[i + 1], cy + ch / 2, color=INK2, lw=1.3)

    # ice couplings note under the core (ice is now a first-class process box above)
    ax.add_patch(FancyBboxPatch((17, 36.3), 80, 5.4,
                 boxstyle="round,pad=0.2,rounding_size=1", fc="#eef5fb", ec="#2a6aa8", lw=1.0,
                 mutation_aspect=0.5))
    ax.text(57, 39.0, "ICE COUPLINGS (ice_frac) — insulates the heat budget · shuts O₂/CO₂ "
            "gas exchange · dims under-ice PAR · conserves (not zeroes) state year-round · "
            "erosion×(1−ice_frac)",
            ha="center", va="center", fontsize=6.6, color="#2a6aa8")

    # arrows inputs -> core (representative couplings)
    arrow(ax, 13, iy, xs[0] + cw / 2, cy + ch, color=HYD)     # discharge -> hyd
    arrow(ax, 37, iy, xs[2] + cw / 2, cy + ch, color=GEO)     # geometry -> transport
    arrow(ax, 62, iy, xs[3] + cw / 2, cy + ch, color=MET)     # met -> heat budget + ice
    arrow(ax, 87, iy, xs[5] + cw / 2, cy + ch, color=CHEM)    # chem -> biogeo

    # ---------- OUTPUTS (bottom) ----------
    ax.text(2, 31, "OUTPUTS", fontsize=9.5, weight="bold", color=INK2)
    oy, oh = 12, 15
    box(ax, 2, oy, 30, oh, "STATE FIELDS  (NetCDF/.dat)", fc="#f4f4f2", ec=INK2, lines=(
        "time × distance matrices, one per",
        "species: S, T, DIC, ALK, pH, O2,",
        "nutrients, SPM, DIA … (13 fields)",
        "+ ice_thickness, ice_frac (prognostic)"))
    box(ax, 35, oy, 30, oh, "PROCESS RATES  (NetCDF/.dat)", fc="#f4f4f2", ec=INK2, lines=(
        "NPP, aer_deg, denit, nit, O2_ex,",
        "NEM, phy_death …",
        "compact NetCDF default;",
        "CGEM_OUTPUT=dat|both for .dat"))
    box(ax, 68, oy, 30, oh, "AIR–SEA CO₂ FLUX", fc="#0b0b0b", ec="#0b0b0b", tc="#ffffff",
        lines=("FCO₂  (>0 = outgassing; mol/kg",
               "carbonate, v2). RESULT: near",
               "equilibrium, small mixed-sign",
               "flux — sensitive to bnd chem"))
    arrow(ax, 17, cy, 17, oy + oh, color=INK2)
    arrow(ax, 50, cy, 50, oy + oh, color=INK2)
    arrow(ax, 85, cy, 83, oy + oh, color=INK2)

    # ---------- PERFORMANCE (measured wall-clock) ----------
    pm = _perf()
    ax.add_patch(FancyBboxPatch((2, 5.4), 96, 4.6,
                 boxstyle="round,pad=0.2,rounding_size=1", fc="#eef5ef", ec=GEO, lw=1.0,
                 mutation_aspect=0.5))
    if pm:
        order = ["colville", "kuparuk", "sagavanirktok", "canning"]
        ps = pm["per_site_sec"]
        per = "   ".join(f"{LAB.get(s, s)} {ps[s] / 60:.1f}" for s in order if s in ps)
        vals = [ps[s] for s in order if s in ps]
        lo, hi = min(vals) / 60, max(vals) / 60
        head = (f"PERFORMANCE (measured):  per river {lo:.0f}–{hi:.0f} min   ·   "
                f"all four in PARALLEL {pm['parallel_total_sec'] / 60:.0f} min   ·   "
                f"serial {pm['serial_total_sec'] / 60:.0f} min")
        ax.text(4, 8.4, head, fontsize=7.8, color=GEO, weight="bold", va="center")
        ax.text(4, 6.5, f"per-river wall-clock (min):   {per}        "
                f"[{pm['config']}]", fontsize=6.6, color=INK2, va="center")
    else:
        ax.text(4, 7.4, "PERFORMANCE: run tools/bench_timing.sh to populate "
                "docs/performance_metrics.json", fontsize=7.2, color=GEO, va="center")

    caption(ax, "Figure 2 | Data flow for the four rivers (2022 climatology; Δt = 75 s; "
            "136 cells over ~27 km). Boxes are coloured by data domain; grey = model "
            "process; arrows show data / flow direction. See docs/ and CLAUDE.md.")
    scale_fonts(fig, PAGE_SCALE["dataflow"])
    pdf.savefig(fig, bbox_inches="tight", pad_inches=0.15)
    plt.close(fig)


def _staggered_grid(ax, x0, y0, w):
    """Schematic of the staggered grid: odd i → concentration nodes (C), even i →
    velocity nodes (U). An ellipsis marks the elided interior; M is even (velocity)."""
    # (label, parity)  parity 1 = odd/concentration, 0 = even/velocity, None = gap
    slots = [("1", 1), ("2", 0), ("3", 1), ("4", 0), ("5", 1), ("6", 0), ("7", 1),
             ("…", None), ("M–1", 1), ("M", 0)]
    dx = w / len(slots)
    ax.plot([x0 + dx * 0.5, x0 + w - dx * 0.5], [y0, y0], color="#ccccc6", lw=1.0, zorder=2)
    for k, (lab, parity) in enumerate(slots):
        x = x0 + (k + 0.5) * dx
        if parity is None:
            ax.text(x, y0, "…", fontsize=11, color=MUTED, ha="center", va="center")
            continue
        if parity == 1:
            ax.plot(x, y0, "o", ms=10, mfc=HYD_F, mec=HYD, mew=1.2, zorder=3)
            ax.text(x, y0, "C", fontsize=6.2, ha="center", va="center", color=HYD,
                    weight="bold", zorder=4)
        else:
            ax.plot(x, y0, "s", ms=10, mfc=MET_F, mec=MET, mew=1.1, zorder=3)
            ax.text(x, y0, "U", fontsize=6.2, ha="center", va="center", color=MET,
                    weight="bold", zorder=4)
        ax.text(x, y0 - 2.5, lab, fontsize=5.6, ha="center", color=MUTED)


# ------------------------------------------------------------------ page 3
def page_detail(pdf):
    fig = plt.figure(figsize=(14, 9))
    ax = canvas(fig)
    header(ax, "Numerical structure — call graph, species, and grid", accent=OCEAN,
           fig_label="NS-RAD · Fig. 3")

    # ---- call graph (left) as a shaded code block ----
    ax.text(2, 88.5, "PER-TIMESTEP CALL GRAPH", fontsize=9.5, weight="bold", color=INK2)
    ax.add_patch(FancyBboxPatch((2, 57.5), 50, 28.5,
                 boxstyle="round,pad=0.4,rounding_size=1.0", fc="#f6f7f5", ec=GRID,
                 lw=0.9, mutation_aspect=0.5, zorder=1))
    steps = [
        "for t in 0 … MAXT step Δt:",
        "  exfread(...) ×N   →  forcings interpolated to t",
        "  river_dispersion(Qr)   →  K[i] from local W,H,U",
        "  hyd(t, Qr)   →  [coeff_a→tridag→conv→update]*",
        "  transport(t)   →  per sp: openbound → tvd → disp_sch",
        "  heat_budget(t, air, wind, solar, RH)   →  fluxes on T",
        "  ice_step(t, …)   →  grow / melt / freshet break-up",
        "  if t > WARMUP:",
        "      biogeo(t, …, pCO2, I0)   →  rates, pH, FCO₂, CH₄/N₂O",
        "      sed(t, previousdays)   →  erode / deposit SPM",
    ]
    for i, s in enumerate(steps):
        head = s.strip().startswith(("for", "if"))
        ax.text(4, 84.2 - i * 2.62, s, fontsize=7.5, family=MONO,
                color=INK2 if head else INK, zorder=3)

    # ---- species registry (right) ----
    ax.text(56, 88.5, "TRANSPORTED SPECIES", fontsize=9.5, weight="bold", color=INK2)
    ax.text(80, 88.4, "variables.v, env = 1", fontsize=7.4, color=MUTED, va="center")
    species = [
        ("S", "salinity (PSU)", 0), ("T", "water temperature (°C)", 0),
        ("DIC", "dissolved inorganic carbon", 0), ("ALK", "alkalinity", 0),
        ("pH", "diagnosed from S / DIC / ALK (env = 0)", 0), ("O2", "dissolved oxygen", 0),
        ("NO3", "nitrate", 0), ("NH4", "ammonium", 0), ("PO4", "phosphate", 0),
        ("dSi", "dissolved silica", 0), ("DIA", "diatoms", 0), ("TOC", "labile DOC", 0),
        ("RDOC", "refractory + chromophoric DOC", 1),
        ("CH4", "dissolved methane", 1), ("N2O", "dissolved nitrous oxide", 1),
        ("SPM", "suspended particulate matter", 0),
    ]
    for i, (sp, desc, new) in enumerate(species):
        y = 84.4 - i * 1.95
        col = MET if new else HYD
        ax.text(57, y, sp, fontsize=7.6, weight="bold", color=col, family=MONO, va="center")
        ax.text(66, y, desc + ("   ← Arctic ext." if new else ""), fontsize=7.1,
                color=(MET if new else INK2), va="center")

    # ---- staggered-grid inset (fills the middle band) ----
    ax.add_patch(FancyBboxPatch((2, 41.5), 96, 12.5,
                 boxstyle="round,pad=0.3,rounding_size=1.0", fc="#fbfbfa", ec=GRID,
                 lw=0.9, mutation_aspect=0.35, zorder=1))
    ax.text(4, 52.2, "STAGGERED GRID", fontsize=9.5, weight="bold", color=INK2)
    ax.text(24, 52.1, "odd indices → concentrations & cross-sections (C, D) · "
            "even indices → velocities & fluxes (U, fl)", fontsize=7.2, color=INK2,
            va="center")
    _staggered_grid(ax, 6, 46.0, 66)
    # legend
    ax.plot(76, 47.4, "o", ms=9, mfc=HYD_F, mec=HYD, mew=1.2)
    ax.text(78, 47.4, "concentration node", fontsize=6.6, color=INK2, va="center")
    ax.plot(76, 44.6, "s", ms=9, mfc=MET_F, mec=MET, mew=1.1)
    ax.text(78, 44.6, "velocity node", fontsize=6.6, color=INK2, va="center")

    # ---- numerics / status box (bottom) ----
    ax.add_patch(FancyBboxPatch((2, 5.5), 96, 32.5, boxstyle="round,pad=0.4,rounding_size=1.2",
                 fc="#f6f6f4", ec=GRID, lw=0.9, mutation_aspect=0.32, zorder=1))
    ax.text(4, 35.4, "NUMERICS, PERFORMANCE & STATUS", fontsize=9.5, weight="bold",
            color=INK2)
    notes = [
        "Hydrodynamics: implicit; hyd iterates to convergence (while rsum != 2.0, an exact-float test). Shallow depths → ~44 iterations/step.",
        "Transport: TVD advection (parity-aware) + Crank–Nicolson dispersion; dispersion capped at DISP_MAX = 2133 m² s⁻¹ for numerical stability.",
        "Performance: hot loops (hydro kernels, density stack, biogeo, transport, sed) are @njit-compiled → ~14–15 min per site; all four in parallel ~15 min (2-yr run, 14 cores).",
        "State: global mutable NumPy arrays in variables.py, mutated in place — never rebind. State is conserved year-round (prognostic ice model); no winter zeroing.",
        "Beyond the shipped model: per-site config, observed forcings, transported T + surface heat budget, prognostic river ice (freeze / melt / freshet break-up), NetCDF output.",
        "Arctic biogeochemistry extension: refractory / chromophoric DOC + CDOM photomineralisation, CH₄ & N₂O cycling + gas exchange, benthic (SOD) DIC / alkalinity efflux, distributed lateral loading (docs/arctic_biogeochemistry.md).",
        "Not yet done: ice-draft reduction of the cross-section (Tier 3); prognostic sediment-OM pool (benthic is a first-order SOD closure); lateral loading adds solute mass without growing discharge; humidity assumed.",
    ]
    for i, n in enumerate(notes):
        y = 31.6 - i * 4.05
        ax.text(4.5, y, "▪", fontsize=7, color=OCEAN, va="center")
        ax.text(6.4, y, n, fontsize=7.4, color=INK2, va="center")

    caption(ax, "Figure 3 | Numerical structure. The per-timestep call graph, the 16 "
            "transported state variables (incl. the Arctic-extension tracers RDOC / CH₄ / "
            "N₂O), and the staggered grid on which they are solved. See CLAUDE.md, "
            "docs/arctic_biogeochemistry.md and docs/performance.md.")
    scale_fonts(fig, PAGE_SCALE["detail"])
    pdf.savefig(fig, bbox_inches="tight", pad_inches=0.15)
    plt.close(fig)


# ------------------------------------------------------------------ page: framework
RIVC = {"colville": "#2a78d6", "kuparuk": "#eb6834",
        "sagavanirktok": "#1baf7a", "canning": "#4a3aa7"}
# short, honest constraint flag per river (the caveats documented in CLAUDE.md)
CAVEAT = {
    "colville": "Q gauged far upstream (Umiat) → understates mouth",
    "kuparuk": "best-constrained: gauge near tidewater + LTER chem",
    "sagavanirktok": "Q understates; upstream flare borrowed from Canning",
    "canning": "Q reconstructed (Hulahula ×2.97); chem still placeholder",
}


def _site_cfg():
    """Read the ACTUAL per-river configuration live from sites/<name>.py + forcings."""
    if str(CODE) not in sys.path:
        sys.path.insert(0, str(CODE))
    import importlib
    out = []
    for s in ("colville", "kuparuk", "sagavanirktok", "canning"):
        m = importlib.import_module(f"sites.{s}")
        q = np.abs(np.genfromtxt(open(FORC / getattr(m, "DISCHARGE_FILE"),
                                      encoding="utf-8-sig"), delimiter=","))
        B = m.BOUNDARIES
        out.append(dict(
            key=s, label=getattr(m, "LABEL"), EL=getattr(m, "EL"),
            B_lb=getattr(m, "B_lb"), B_ub=getattr(m, "B_ub"),
            L_FLARE=getattr(m, "L_FLARE", 7000), D=getattr(m, "DEPTH_lb"),
            Q=float(q.mean()), gauge=getattr(m, "GAUGE", None),
            chem=getattr(m, "BOUNDARY_CHEM_SOURCE", None),
            ALK=B["ALK"][1], DIC=B["DIC"][1], TOC=B["TOC"][1]))
    return out


def _funnel(x_up, x_mouth, yc, hw_up, hw_mouth, flare_span):
    """Plan-view estuary channel: prismatic width hw_up from x_up, flaring OPEN to
    hw_mouth over the last `flare_span` before the mouth (right). Returns polygon pts."""
    xf = x_mouth - flare_span
    return [(x_up, yc + hw_up), (xf, yc + hw_up), (x_mouth, yc + hw_mouth),
            (x_mouth, yc - hw_mouth), (xf, yc - hw_up), (x_up, yc - hw_up)]


def _grad_fill(ax, pts, c_bot, c_top, zorder=2):
    """Fill an arbitrary polygon with a smooth vertical gradient (c_bot→c_top)."""
    xs = [p[0] for p in pts]
    ys = [p[1] for p in pts]
    cmap = LinearSegmentedColormap.from_list("g", [c_bot, c_top])
    im = ax.imshow(np.linspace(0, 1, 256).reshape(-1, 1),
                   extent=[min(xs), max(xs), min(ys), max(ys)], origin="lower",
                   cmap=cmap, aspect="auto", zorder=zorder, interpolation="bilinear")
    im.set_clip_path(MplPath(pts + [pts[0]]), ax.transData)
    return im


def page_framework(pdf):
    cfg = _site_cfg()
    fig = plt.figure(figsize=(14, 9))
    ax = canvas(fig)
    header(ax, "A land-to-ocean modelling framework for North Slope rivers",
           "Four Arctic rivers resolved as 1-D reactive estuaries from Brooks Range "
           "headwaters to the Beaufort Sea. Per-river configuration read live from the "
           "model source.", accent=OCEAN, fig_label="NS-RAD · Fig. 1")

    # ================= (a) ATMOSPHERE =================
    ax.add_patch(FancyBboxPatch((2, 81.5), 96, 8.0,
                 boxstyle="round,pad=0.2,rounding_size=0.8", fc="#f2f7fc", ec="#cfe0ef",
                 lw=0.9, mutation_aspect=0.5, zorder=1))
    plabel(ax, 3, 89.0, "a")
    ax.add_patch(Circle((9.4, 85.5), 1.7, fc="#f3c34a", ec="#e0a92a", lw=0.9, zorder=2))
    for k in range(8):  # sun rays
        a = k * np.pi / 4
        ax.plot([9.4 + 2.1 * np.cos(a), 9.4 + 2.9 * np.cos(a)],
                [85.5 + 2.1 * np.sin(a), 85.5 + 2.9 * np.sin(a)],
                color="#e0a92a", lw=0.9, zorder=2)
    ax.text(14.5, 87.9, "ATMOSPHERE", fontsize=9.5, weight="bold", color=MET)
    ax.text(37.0, 87.85, "shared regional forcing", fontsize=8.0, color=INK2, va="center")
    ax.text(14.5, 85.0, "Wind · Solar · Air temperature — NDBC PRDA2 (Prudhoe Bay)      "
            "Humidity — Deadhorse ISD      Atmospheric pCO₂ — Barrow", fontsize=7.6,
            color=INK2)
    ax.text(14.5, 82.9, "drives the surface heat budget, gas exchange and ice — applied at "
            "every grid cell, every timestep", fontsize=6.9, color=MUTED, style="italic")

    yc = 60.5
    x_up, x_mouth = 18.5, 82
    hw_up, hw_mouth, flare = 4.4, 11.0, 24
    # forcing rains DOWN (met), CO2 rises UP (offset so they never collide)
    for x in (30, 44, 58, 72):
        arrow(ax, x, 81.3, x, yc + _hw_at(x, x_up, x_mouth, hw_up, hw_mouth, flare) + 1.2,
              color=MET, lw=1.0)
    for x in (37, 51, 65):
        hw = _hw_at(x, x_up, x_mouth, hw_up, hw_mouth, flare)
        arrow(ax, x, yc + hw + 1.2, x, 81.3, color=INK, lw=1.0)
    ax.text(80.5, 78.0, "air–sea CO₂\n(FCO₂ ↑ outgassing)", fontsize=6.7, color=INK,
            ha="center", weight="bold", linespacing=1.3)

    # ================= (b) LAND -> OCEAN cross-section =================
    plabel(ax, 3, 71.5, "b")
    # LAND: layered Brooks Range (back hazy, front darker) with snow caps + tundra
    for (bx, bh, fc, ec) in [(5, 12, "#cdc7ba", "#b7b0a1"), (10, 15.5, "#bdb6a6", "#a49c88"),
                             (14.5, 11, "#c6bfae", "#aca492")]:
        ax.add_patch(Polygon([(bx - 3.4, 52), (bx, 52 + bh), (bx + 3.6, 52)], closed=True,
                     fc=fc, ec=ec, lw=0.7, zorder=1))
        ax.add_patch(Polygon([(bx - 1.0, 52 + bh - 3.2), (bx, 52 + bh),
                     (bx + 1.05, 52 + bh - 3.2)], closed=True, fc="#f4f6f8",
                     ec="none", zorder=2))  # snow cap
    ax.add_patch(Rectangle((2, 49.2), 16.5, 3.0, fc="#e3e8d5", ec="#cbd1b6", lw=0.6, zorder=1))
    ax.text(10, 50.6, "tundra", fontsize=6.4, color="#6f7757", ha="center", style="italic")
    ax.text(10.2, 70.3, "HEADWATERS", fontsize=8.3, weight="bold", color="#6f6650",
            ha="center")
    ax.text(10.2, 68.1, "Brooks Range · tundra", fontsize=6.7, color=MUTED, ha="center")

    # CHANNEL: estuary funnel, gradient water, opening to the sea (mouth on the right)
    ch = _funnel(x_up, x_mouth, yc, hw_up, hw_mouth, flare)
    _grad_fill(ax, ch, "#a9d2ec", "#dceff9", zorder=2)
    ax.add_patch(Polygon(ch, closed=True, fc="none", ec="#5f9ec6", lw=1.1, zorder=5))
    for gx in np.linspace(x_up + 1.5, x_mouth - 1.2, 24):  # faint staggered-grid ticks
        hw = _hw_at(gx, x_up, x_mouth, hw_up, hw_mouth, flare)
        ax.plot([gx, gx], [yc - hw + 0.3, yc + hw - 0.3], color="#7fb2d6", lw=0.25,
                alpha=0.55, zorder=3)
    # landfast/river ice slab over the upstream reach
    ax.add_patch(Polygon([(x_up, yc + hw_up), (x_up + 33, yc + hw_up),
                 (x_up + 33, yc + hw_up + 1.7), (x_up, yc + hw_up + 1.7)], closed=True,
                 fc="#eef5fb", ec=ICE_E, lw=0.6, hatch="///", zorder=6))
    ax.text(x_up + 0.5, yc + hw_up + 3.2, "landfast / river ice", fontsize=6.5,
            color="#3d6f99", ha="left", zorder=7)
    # process labels inside the water (centred in the channel)
    xm = (x_up + x_mouth) / 2
    ax.text(xm, yc + 2.0, "1-D reactive estuary", fontsize=9.4, weight="bold",
            color="#1b567f", ha="center", zorder=7)
    ax.text(xm, yc - 0.4, "TVD advection + dispersion · heat budget · prognostic ice · "
            "carbonate → pH · air–sea CO₂", fontsize=6.8, color="#276a97", ha="center",
            zorder=7)
    ax.text(xm, yc - 2.7, "Δt = 75 s · 136 cells · staggered grid · 2-yr climatology",
            fontsize=6.3, color="#5c86a3", ha="center", zorder=7)

    # OCEAN: Beaufort Sea (right), gradient + waves
    ocean = [(82, 49), (98, 49), (98, 71.5), (82, 71.5)]
    _grad_fill(ax, ocean, OCEAN_D, "#4a86ad", zorder=1)
    ax.add_patch(Polygon(ocean, closed=True, fc="none", ec=OCEAN_D, lw=0.9, zorder=5))
    ax.text(90.2, 68.4, "BEAUFORT SEA", fontsize=8.6, weight="bold", color="white",
            ha="center", zorder=6)
    for wy in (53.5, 56.5, 59.5):
        ax.plot(np.linspace(83.2, 96.8, 60),
                wy + 0.45 * np.sin(np.linspace(0, 6 * np.pi, 60)), color="#bcd8ea",
                lw=0.7, alpha=0.8, zorder=6)
    # discreet distance scale bar under the channel, clear of the boundary arrows
    _sb = 33.0
    ax.plot([_sb, _sb + 12.3], [50.2, 50.2], color=INK2, lw=1.1, solid_capstyle="butt")
    for xx in (_sb, _sb + 12.3):
        ax.plot([xx, xx], [49.8, 50.6], color=INK2, lw=1.1)
    ax.text(_sb + 6.15, 48.6, "5 km", fontsize=6.0, color=INK2, ha="center")

    # RIVERINE boundary (left, cub) — one clean arrow into the upstream face
    ax.text(3, 45.6, "RIVERINE BOUNDARY", fontsize=7.6, weight="bold", color=HYD)
    ax.text(3, 43.7, "upstream (cub)", fontsize=6.6, color=MUTED, style="italic")
    for i, (lab, col) in enumerate([("Discharge — USGS gauges", HYD),
                                    ("River T — air-T regression + USGS 00010", MET),
                                    ("River chemistry — WQP / Arctic LTER", GEO)]):
        yy = 41.1 - i * 2.4
        ax.add_patch(Circle((4.0, yy), 0.42, fc=col, ec="none", zorder=6))
        ax.text(5.2, yy, lab, fontsize=6.7, color=INK2, ha="left", va="center")
    arrow(ax, 21.0, 44.4, x_up + 2.0, yc - hw_up - 0.6, color=HYD, lw=1.4)

    # MARINE boundary (right, clb) — one clean arrow into the mouth face
    ax.text(97, 45.6, "MARINE BOUNDARY", fontsize=7.6, weight="bold", color=OCEAN,
            ha="right")
    ax.text(97, 43.7, "downstream (clb)", fontsize=6.6, color=MUTED, style="italic",
            ha="right")
    for i, lab in enumerate(["Sea T — PRDA2 buoy", "Salinity", "DIC / ALK — seawater",
                             "Tides — NOAA CO-OPS harmonic"]):
        yy = 41.1 - i * 2.4
        ax.add_patch(Circle((96.0, yy), 0.42, fc=OCEAN, ec="none", zorder=6))
        ax.text(94.8, yy, lab, fontsize=6.7, color=INK2, ha="right", va="center")
    arrow(ax, 90.0, 44.4, x_mouth - 2.0, yc - hw_mouth - 0.6, color=OCEAN, lw=1.4)

    # ================= (c) THE FOUR RIVERS =================
    plabel(ax, 3, 32.5, "c")
    ax.text(6.5, 32.3, "THE FOUR RIVERS", fontsize=9.5, weight="bold", color=INK)
    ax.text(31.0, 32.25, "actual configuration, read live from the model source",
            fontsize=7.6, color=INK2, va="center")
    cw, gap, x0 = 23.2, 1.6, 3.2
    Bmax = max(c["B_lb"] for c in cfg)
    for j, c in enumerate(cfg):
        x = x0 + j * (cw + gap)
        col = RIVC[c["key"]]
        ax.add_patch(FancyBboxPatch((x, 3.2), cw, 24.4,
                     boxstyle="round,pad=0.3,rounding_size=0.9", fc="#fcfcfb",
                     ec="#d7d7d2", lw=0.9, mutation_aspect=0.7, zorder=2))
        ax.add_patch(Rectangle((x, 24.7), cw, 2.9, fc=col, ec="none", zorder=3))
        ax.add_patch(FancyBboxPatch((x, 3.2), cw, 24.4,
                     boxstyle="round,pad=0.3,rounding_size=0.9", fc="none",
                     ec="#d7d7d2", lw=0.9, mutation_aspect=0.7, zorder=4))
        ax.text(x + cw / 2, 26.1, c["label"], fontsize=9, weight="bold", color="white",
                ha="center", va="center", zorder=5)
        # scaled funnel glyph (real B_lb : B_ub proportions)
        gy = 20.7
        hw_m = 2.3 * (c["B_lb"] / Bmax) + 0.55
        hw_u = 2.3 * (c["B_ub"] / Bmax) + 0.3
        gpts = _funnel(x + 3.2, x + cw - 3.2, gy, hw_u, hw_m,
                       (cw - 6.4) * min(1.0, c["L_FLARE"] / c["EL"]) * 3.4)
        _grad_fill(ax, gpts, "#cfe4f2", "#eaf3fa", zorder=3)
        ax.add_patch(Polygon(gpts, closed=True, fc="none", ec=col, lw=0.9, zorder=5))
        rows = [
            ("mouth width", f"{c['B_lb']:.0f} m"),
            ("upstream width", f"{c['B_ub']:.0f} m"),
            ("depth · domain", f"{c['D']:.2f} m · 27 km"),
            ("mean discharge", f"{c['Q']:.0f} m³ s⁻¹"),
            ("river ALK / DIC", f"{c['ALK']:.0f} / {c['DIC']:.0f}"),
        ]
        for i, (k, val) in enumerate(rows):
            yy = 15.0 - i * 2.05
            ax.text(x + 1.7, yy, k, fontsize=6.2, color=MUTED, ha="left", va="center")
            ax.text(x + cw - 1.7, yy, val, fontsize=6.4, color=INK, ha="right",
                    va="center", weight="medium")
        ax.plot([x + 1.7, x + cw - 1.7], [5.3, 5.3], color="#e4e4df", lw=0.6)
        ax.text(x + cw / 2, 4.2, CAVEAT[c["key"]], fontsize=5.4, color=col,
                ha="center", va="center", style="italic")

    caption(ax, "Figure 1 | The NS-RAD framework. (a) Shared regional "
            "meteorology forces every cell; (b) each river is a 1-D reactive estuary "
            "between a riverine (cub) and marine (clb) boundary, with prognostic landfast "
            "ice; (c) the four rivers' measured geometry, discharge and boundary "
            "chemistry (funnel glyphs to scale). A fifth synthetic site, idealized, is a "
            "time-varying verification fixture (docs/idealized_verification).")
    scale_fonts(fig, PAGE_SCALE["framework"])
    pdf.savefig(fig, bbox_inches="tight", pad_inches=0.15)
    plt.close(fig)


def _hw_at(x, x_up, x_mouth, hw_up, hw_mouth, flare_span):
    """Channel half-width at plan-view x, matching _funnel()."""
    xf = x_mouth - flare_span
    if x <= xf:
        return hw_up
    return hw_up + (hw_mouth - hw_up) * (x - xf) / flare_span


def main():
    OUT.parent.mkdir(exist_ok=True)
    with PdfPages(OUT) as pdf:
        page_framework(pdf)
        page_dataflow(pdf)
        page_detail(pdf)
    print(f"wrote {OUT.relative_to(ROOT)} ({OUT.stat().st_size/1024:.0f} kB, 3 pages)")


if __name__ == "__main__":
    main()
