#!/usr/bin/env python3
"""
Summary-diagnostics PDF for the IDEALIZED verification experiment (sites/idealized.py).

Reads the full seasonal run (default runs/idealized/output.nc, written by
tools/verify_idealized.py --full) and the analytic forcings, and lays out a focused
report:

  Page 1  the experiment: analytic forcings + the verification-check results table
          (the table IS the output of tools/verify_idealized, so the PDF and the
          harness can never disagree)
  Page 2  time-varying BOUNDARY verification -- the model boundary cell overlaid on
          the prescribed forcing series, for every forced clb/cub (the headline)
  Page 3  seasonal evolution: distance x time Hovmollers of the core fields + ice
  Page 4  the freshet / ice story + carbonate: hydrograph vs ice, salt intrusion,
          pH and air-sea CO2 through the season

This is a MODEL-machinery test on synthetic input, not a river result.

Usage:
    python tools/make_idealized_verification_pdf.py
    python tools/make_idealized_verification_pdf.py --run runs/idealized
"""
import argparse
import sys
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

ROOT = Path(__file__).resolve().parent.parent
FORCINGS = ROOT / "forcing"
OUT = ROOT / "docs" / "idealized_verification.pdf"
DELXI_KM = 0.2
DAY = 86400.0

ACCENT = "#7a3a8c"
INK, INK2, MUTED = "#0b0b0b", "#52514e", "#8a8880"
GRID, SURFACE, OK, BAD = "#e6e5e0", "#fcfcfb", "#1a8f4c", "#cc3333"
plt.rcParams.update({
    "figure.facecolor": SURFACE, "axes.facecolor": SURFACE, "savefig.facecolor": SURFACE,
    "font.size": 8, "axes.titlesize": 9, "axes.labelsize": 8, "axes.edgecolor": GRID,
    "axes.linewidth": 0.8, "axes.labelcolor": INK2, "text.color": INK,
    "xtick.color": INK2, "ytick.color": INK2, "xtick.labelsize": 7, "ytick.labelsize": 7,
    "legend.frameon": False, "legend.fontsize": 7, "grid.color": GRID, "grid.linewidth": 0.6,
})
MONTH0 = [0, 31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334]
MONTHL = list("JFMAMJJASOND")

# forced boundaries: (species, end, file, label, cell-picker)  cell: M for cub, 1 for clb
FORCED = [
    ("TOC", "cub", "idealized_TOC_cub_river.csv", "river DOC (TOC)"),
    ("ALK", "cub", "idealized_ALK_cub_river.csv", "river alkalinity"),
    ("DIC", "cub", "idealized_DIC_cub_river.csv", "river DIC"),
    ("NO3", "cub", "idealized_NO3_cub_river.csv", "river nitrate"),
    ("S",   "clb", "idealized_S_clb_marine.csv", "marine salinity"),
]


def tidy(ax, both=True):
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    ax.grid(True, axis="both" if both else "y", alpha=0.55)
    ax.set_axisbelow(True)


def months(ax, tmax):
    ax.set_xticks([m for m in MONTH0 if m <= tmax])
    ax.set_xticklabels([MONTHL[i] for i, m in enumerate(MONTH0) if m <= tmax])
    ax.set_xlim(0, tmax)


def forcing(fn):
    return np.genfromtxt(open(FORCINGS / fn, encoding="utf-8-sig"), delimiter=",", dtype=float)


def interp_forcing(fn, t_days):
    arr = forcing(fn)
    axis = np.linspace(0, 365, 365)
    tt = np.where(t_days > 365, t_days - 365, t_days)
    return np.interp(tt, axis, arr)


class NC:
    def __init__(self, path):
        from netCDF4 import Dataset
        self.d = Dataset(path)
        self.t = np.asarray(self.d.variables["time"][:]) / DAY

    def has(self, v):
        return v in self.d.variables

    def f(self, v):
        return np.ma.filled(self.d.variables[v][:], np.nan)  # (time, x)


# ---------------------------------------------------------------------------
def page_experiment(pdf, nc):
    fig = plt.figure(figsize=(8.5, 11))
    fig.suptitle("Idealized North Slope river — verification experiment",
                 x=0.5, y=0.975, fontsize=15, weight="bold", color=ACCENT)
    fig.text(0.5, 0.952, "a synthetic fixture that exercises the whole model on a known, "
             "analytic, time-varying input", ha="center", fontsize=8.5, color=INK2)

    intro = (
        "Not a river. Round-number geometry (30 km, 2 m deep, 400→100 m flare), analytic "
        "meteorology and discharge, and — the point of the fixture — TIME-VARYING boundary "
        "chemistry driven from forcing files via the general BOUNDARY_FORCING mechanism. "
        "The spring freshet is large enough to trigger hydraulic ice break-up; the air "
        "temperature crosses 0 °C so the prognostic ice model freezes up and melts out. "
        "Everything is generated from one PARAMS block (tools/build_idealized_forcings.py), "
        "so the experiment is reproducible. Run it with tools/verify_idealized.py."
    )
    fig.text(0.07, 0.905, intro, ha="left", va="top", fontsize=8.2, wrap=True,
             color=INK, linespacing=1.5)

    # --- a curated set of forcings across the top ---
    fpanels = [
        ("idealized_discharge_m3sec.csv", "discharge (m³ s⁻¹)"),
        ("idealized_airtemp_degC.csv", "air temperature (°C)"),
        ("idealized_solar_Wm2.csv", "shortwave (W m⁻²)"),
    ]
    for i, (fn, lab) in enumerate(fpanels):
        ax = fig.add_axes([0.09 + i * 0.31, 0.68, 0.23, 0.14])
        y = forcing(fn)
        ax.plot(np.arange(365), y, color=ACCENT, lw=1.5)
        if "airtemp" in fn:
            ax.axhline(0, color=BAD, lw=0.8, ls="--")
        ax.set_title(lab, fontsize=7.6)
        months(ax, 364)
        tidy(ax)

    # --- verification results table (straight from the harness) ---
    fig.text(0.07, 0.635, "Verification checks", fontsize=10.5, weight="bold", color=INK)
    fig.text(0.07, 0.618, "produced by tools/verify_idealized.py --full on this run; "
             "the harness is the source of truth", fontsize=7.3, color=INK2)

    sys.path.insert(0, str(ROOT / "tools"))
    import verify_idealized as V
    import io, contextlib
    chk = V.Checks()
    with contextlib.redirect_stdout(io.StringIO()):
        V.check_output(chk, str(nc.path), days=int(nc.t.max()) + 1, warmup=30, full=True)
    rows = chk.rows

    # two columns of results
    n = len(rows)
    half = (n + 1) // 2
    cols = [rows[:half], rows[half:]]
    y0, dy = 0.595, 0.0165
    for ci, col in enumerate(cols):
        x = 0.07 + ci * 0.47
        for k, (name, ok, detail) in enumerate(col):
            yy = y0 - k * dy
            fig.text(x, yy, "✓" if ok else "✗", color=OK if ok else BAD,
                     fontsize=8.5, weight="bold")
            txt = name if len(name) <= 46 else name[:44] + "…"
            fig.text(x + 0.018, yy, txt, fontsize=6.6, color=INK, va="baseline")

    npass = sum(1 for _, ok, _ in rows if ok)
    fig.text(0.07, 0.055, f"{npass}/{n} checks passed"
             + ("  — ALL PASSED" if npass == n else "  — SOME FAILED"),
             fontsize=9.5, weight="bold", color=OK if npass == n else BAD)
    fig.text(0.07, 0.035, "Model machinery on synthetic input — not a river result.",
             fontsize=7.2, color=MUTED, style="italic")
    pdf.savefig(fig)
    plt.close(fig)


def page_boundaries(pdf, nc):
    M = nc.d.dimensions["x"].size - 1
    fig, axes = plt.subplots(3, 2, figsize=(8.5, 11))
    fig.suptitle("Time-varying boundary verification", x=0.5, y=0.97,
                 fontsize=14, weight="bold", color=ACCENT)
    fig.text(0.5, 0.945, "model boundary cell (points) overlaid on the prescribed "
             "forcing series (line) — the mechanism under test", ha="center",
             fontsize=8.3, color=INK2)
    t = nc.t
    axes = axes.ravel()
    for ax, (sp, end, fn, lab) in zip(axes, FORCED):
        cell = M if end == "cub" else 1
        model = nc.f(sp)[:, cell]
        target = interp_forcing(fn, t)
        ax.plot(t, target, color=INK, lw=2.0, label="prescribed forcing", zorder=2)
        ax.plot(t, model, color=ACCENT, lw=0.0, marker=".", ms=1.6, alpha=0.6,
                label=f"model cell {'M (upstream)' if end=='cub' else '1 (mouth)'}", zorder=3)
        if end == "cub":
            dev = np.nanmax(np.abs(model - target))
            note = f"max dev {dev:.1f} / range {np.ptp(target):.0f}"
        else:
            cc = np.corrcoef(model, target)[0, 1]
            note = f"daily-mean r = {cc:.2f}  (freshet flushes the mouth fresh)"
        ax.set_title(f"{sp}.{end} — {lab}", fontsize=8.6, color=INK)
        ax.text(0.03, 0.06, note, transform=ax.transAxes, fontsize=6.8, color=INK2)
        ax.set_ylabel("mmol m⁻³" if sp != "S" else "PSU", fontsize=7)
        months(ax, t.max())
        tidy(ax)
        ax.legend(loc="upper right", fontsize=6.3)
    axes[-1].axis("off")
    fig.text(0.55, 0.20,
             "Upstream (cub) cells are set almost directly from the forcing, so they\n"
             "track it near-exactly. The marine (clb) salinity cell is not bounded by\n"
             "the sea value: the ~400 m³/s freshet flushes the salt intrusion out to\n"
             "sea (cell → 0 PSU), so it is verified by correlation, not tight match.",
             fontsize=7.6, color=INK, va="top", linespacing=1.6)
    fig.tight_layout(rect=[0, 0, 1, 0.93])
    pdf.savefig(fig)
    plt.close(fig)


def page_hovmoller(pdf, nc):
    picks = [("S", "salinity (PSU)", "viridis"),
             ("T", "temperature (°C)", "inferno"),
             ("TOC", "DOC (mmol m⁻³)", "YlGn"),
             ("DIC", "DIC (mmol m⁻³)", "plasma"),
             ("ice_thickness", "ice thickness (m)", "Blues"),
             ("FCO2", "air–sea CO₂ (>0 outgas)", "magma_r")]
    picks = [p for p in picks if nc.has(p[0])]
    fig, axes = plt.subplots(3, 2, figsize=(8.5, 11))
    fig.suptitle("Seasonal evolution — distance × time", x=0.5, y=0.97,
                 fontsize=14, weight="bold", color=ACCENT)
    t = nc.t
    for ax, (v, lab, cm) in zip(axes.ravel(), picks):
        F = nc.f(v)[:, 1:]
        x = np.arange(F.shape[1]) * DELXI_KM
        step = max(1, len(t) // 600)
        im = ax.pcolormesh(t[::step], x, F[::step].T, cmap=cm, shading="auto")
        cb = fig.colorbar(im, ax=ax, pad=0.02, fraction=0.046)
        cb.ax.tick_params(labelsize=6)
        cb.set_label(lab, fontsize=6.8)
        ax.set_ylabel("distance from mouth (km)", fontsize=7)
        months(ax, t.max())
    for ax in axes.ravel()[len(picks):]:
        ax.axis("off")
    fig.tight_layout(rect=[0, 0, 1, 0.94])
    pdf.savefig(fig)
    plt.close(fig)


def page_arctic(pdf, nc):
    """The Arctic biogeochemistry extension: the new tracers and their fluxes."""
    t = nc.t
    M = nc.d.dimensions["x"].size - 1
    fig, axes = plt.subplots(3, 2, figsize=(8.5, 11))
    fig.suptitle("Arctic biogeochemistry extension — tracers & fluxes", x=0.5, y=0.975,
                 fontsize=14, weight="bold", color=ACCENT)
    fig.text(0.5, 0.952, "refractory/chromophoric DOC + photomineralisation · CH₄ & N₂O · "
             "benthic efflux · distributed lateral loading", ha="center", fontsize=8,
             color=INK2)

    # Hovmöllers of the new concentration tracers
    hov = [("RDOC", "refractory DOC (mmol m⁻³)", "YlOrBr"),
           ("CH4", "dissolved CH₄ (mmol m⁻³)", "viridis"),
           ("photo", "CDOM photomineralisation (mmol m⁻³ s⁻¹)", "magma_r"),
           ("N2O", "dissolved N₂O (mmol m⁻³)", "cividis")]
    for k, (ax, (v, lab, cm)) in enumerate(zip(axes.ravel()[:4], hov)):
        if not nc.has(v):
            ax.axis("off"); continue
        F = nc.f(v)[:, 1:]
        x = np.arange(F.shape[1]) * DELXI_KM
        step = max(1, len(t) // 600)
        im = ax.pcolormesh(t[::step], x, F[::step].T, cmap=cm, shading="auto")
        cb = fig.colorbar(im, ax=ax, pad=0.02, fraction=0.046)
        cb.ax.tick_params(labelsize=6); cb.set_label(lab, fontsize=6.6)
        ax.set_ylabel("distance from mouth (km)", fontsize=7)
        ax.set_title(f"({'abcd'[k]}) {v}", fontsize=8.4, loc="left")
        months(ax, t.max())

    day = 86400.0
    # (e) channel-mean Arctic carbon rates vs time -> the seasonal turn-on
    ax = axes[2, 0]
    for v, lab, col in [("photo", "CDOM photomineralisation", ACCENT),
                        ("rdoc_ox", "refractory DOC oxidation", "#c98a2b"),
                        ("ch4_ox", "methane oxidation", "#1b7a4b"),
                        ("sod", "benthic respiration (SOD)", INK)]:
        if nc.has(v):
            r = np.nanmean(nc.f(v)[:, 1:], axis=1) * day     # mmol m⁻³ d⁻¹
            ax.plot(t, r, lw=1.3, color=col, label=lab)
    ax.set_ylabel("rate (mmol m⁻³ d⁻¹)", fontsize=7.5)
    ax.set_title("(e) Arctic C rates — channel mean", fontsize=8.6)
    months(ax, t.max()); tidy(ax); ax.legend(fontsize=6.2, loc="upper left")

    # (f) greenhouse-gas air–water fluxes vs time (>0 = outgassing)
    ax = axes[2, 1]
    for v, lab, col in [("ch4_ex", "CH₄ air–water flux", "#1b7a4b"),
                        ("n2o_ex", "N₂O air–water flux (×20)", "#8a4fbf")]:
        if nc.has(v):
            f = np.nanmean(nc.f(v)[:, 1:], axis=1) * day     # mmol m⁻³ d⁻¹
            if "N₂O" in lab:
                f = f * 20.0
            ax.plot(t, f, lw=1.3, color=col, label=lab)
    ax.axhline(0, color=MUTED, lw=0.7, ls="--")
    ax.set_ylabel("flux (mmol m⁻³ d⁻¹, >0 outgas)", fontsize=7.5)
    ax.set_title("(f) CH₄ / N₂O outgassing — channel mean", fontsize=8.6)
    months(ax, t.max()); tidy(ax); ax.legend(fontsize=6.4, loc="upper left")

    fig.tight_layout(rect=[0, 0.05, 1, 0.94])
    fig.text(0.5, 0.025, "Arctic extension. (a–d) the new transported tracers and the "
             "photomineralisation rate (distance × time); (e) channel-mean carbon rates — "
             "photochemistry and methane oxidation switch on with open water, benthic "
             "respiration runs year-round; (f) CH₄ / N₂O outgassing, zero under ice, "
             "venting at break-up. Magnitudes are literature-plausible "
             "(docs/arctic_biogeochemistry.md); the idealized river runs warm (to 20 °C), "
             "inflating the Q₁₀ rates vs a real Arctic river.", ha="center", va="bottom",
             fontsize=6.6, color=MUTED, wrap=True)
    pdf.savefig(fig)
    plt.close(fig)


def page_story(pdf, nc):
    M = nc.d.dimensions["x"].size - 1
    t = nc.t
    fig, axes = plt.subplots(2, 2, figsize=(8.5, 11))
    fig.suptitle("Freshet, ice and carbonate through the season", x=0.5, y=0.97,
                 fontsize=14, weight="bold", color=ACCENT)

    # (a) discharge hydrograph + ice thickness at an upstream cell
    ax = axes[0, 0]
    q = interp_forcing("idealized_discharge_m3sec.csv", t)
    ax.plot(t, q, color=ACCENT, lw=1.6, label="discharge")
    ax.set_ylabel("discharge (m³ s⁻¹)", color=ACCENT, fontsize=7.5)
    ax.set_title("(a) hydrograph vs ice: freshet clears the ice", fontsize=8.6)
    if nc.has("ice_thickness"):
        ice_up = nc.f("ice_thickness")[:, M - 6]
        ax2 = ax.twinx()
        ax2.plot(t, ice_up, color="#2a6fb0", lw=1.6, label="ice (upstream)")
        ax2.set_ylabel("ice thickness (m)", color="#2a6fb0", fontsize=7.5)
        ax2.spines["top"].set_visible(False)
    months(ax, t.max())
    tidy(ax, both=False)

    # (b) salt intrusion: S at mouth + a few cells inland
    ax = axes[0, 1]
    for cell, lab, col in [(1, "mouth", INK), (10, "2 km", ACCENT), (25, "5 km", "#c98a2b")]:
        ax.plot(t, nc.f("S")[:, cell], lw=1.3, color=col, label=lab)
    ax.set_ylabel("salinity (PSU)", fontsize=7.5)
    ax.set_title("(b) salt intrusion: freshet flush + summer freshening", fontsize=8.6)
    months(ax, t.max())
    tidy(ax)
    ax.legend(fontsize=6.5, loc="upper right")

    # (c) pH through the season at mouth + mid-channel
    ax = axes[1, 0]
    for cell, lab, col in [(1, "mouth", INK), (M // 2, "mid-channel", ACCENT), (M, "upstream", "#c98a2b")]:
        ax.plot(t, nc.f("pH")[:, cell], lw=1.3, color=col, label=lab)
    ax.set_ylabel("pH", fontsize=7.5)
    ax.set_title("(c) pH — carbonate solve stays physical", fontsize=8.6)
    ax.axhspan(2, 5, color=BAD, alpha=0.05)
    months(ax, t.max())
    tidy(ax)
    ax.legend(fontsize=6.5, loc="lower right")

    # (d) air-sea CO2 flux, spatial mean, open water
    ax = axes[1, 1]
    if nc.has("FCO2"):
        fco2 = nc.f("FCO2")[:, 1:]
        # all-NaN rows (fully ice-covered / pre-warmup) -> plot as 0, no warning
        finite = np.isfinite(fco2).any(axis=1)
        mean = np.zeros(len(t))
        mean[finite] = np.nanmean(fco2[finite], axis=1)
        ax.plot(t, mean, lw=1.4, color=ACCENT)
        ax.axhline(0, color=INK2, lw=0.7, ls="--")
        ax.set_ylabel("air–sea CO₂ flux (>0 outgas)", fontsize=7.5)
    ax.set_title("(d) air–sea CO₂ (channel mean)", fontsize=8.6)
    months(ax, t.max())
    tidy(ax)

    fig.tight_layout(rect=[0, 0, 1, 0.94])
    pdf.savefig(fig)
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--run", default=str(ROOT / "runs" / "idealized"))
    ap.add_argument("--out", default=str(OUT))
    args = ap.parse_args()
    ncpath = Path(args.run) / "output.nc"
    if not ncpath.exists():
        sys.exit(f"no output.nc at {ncpath}. Run: python tools/verify_idealized.py --full")
    nc = NC(ncpath)
    nc.path = ncpath
    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    with PdfPages(args.out) as pdf:
        page_experiment(pdf, nc)
        page_boundaries(pdf, nc)
        page_hovmoller(pdf, nc)
        page_arctic(pdf, nc)
        page_story(pdf, nc)
    print(f"wrote {args.out}  ({Path(args.out).stat().st_size/1e3:.0f} KB, 5 pages)")


if __name__ == "__main__":
    main()
