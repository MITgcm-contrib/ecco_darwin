"""
Channel-geometry document: the idealized C-GEM geometry vs the observed data, per river.

Page 1  along-channel WIDTH profiles: the flare model actually used, the per-channel
        SWORD width it derives from, and the delta-mouth distributary-sum widening.
Page 2  geometry + hydraulics table (widths, depth, mouth velocity, W/H, dispersion)
        and the data provenance / SWORD-pending status.

The raw SWORD node-width profiles (per-channel and braided sum) will be overlaid here
once the SWORD v17 extraction lands (docs/sword_widths.json); until then the seaward
mouth uses an INTERIM morphological distributary count (N_CHAN_MOUTH in sites/*.py).

Usage:  python tools/make_geometry_pdf.py
"""
import json
import math
import importlib.util
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import nsrad_style as _S  # shared NS-RAD publication/presentation style
from matplotlib.backends.backend_pdf import PdfPages

ROOT = Path(__file__).resolve().parent.parent
SITEDIR = ROOT / "code" / "sites"
OUT = ROOT / "docs" / "ns_rad_geometry.pdf"
SWORD_JSON = ROOT / "docs" / "sword_widths.json"   # filled by the SWORD extraction

SITES = ["colville", "kuparuk", "sagavanirktok", "canning"]
LABEL = {"colville": "Colville", "kuparuk": "Kuparuk",
         "sagavanirktok": "Sagavanirktok", "canning": "Canning"}
C = {"colville": "#2a78d6", "kuparuk": "#eb6834",
     "sagavanirktok": "#1baf7a", "canning": "#4a3aa7"}
INK, INK2, MUTED = "#0b0b0b", "#52514e", "#8a8880"
GRID, SURFACE, WARN = "#e6e5e0", "#fcfcfb", "#e34948"
plt.rcParams.update({
    "figure.facecolor": SURFACE, "axes.facecolor": SURFACE, "savefig.facecolor": SURFACE,
    "font.size": 8, "axes.titlesize": 9, "axes.labelsize": 8, "axes.edgecolor": GRID,
    "axes.linewidth": 0.8, "axes.labelcolor": INK2, "text.color": INK,
    "xtick.color": INK2, "ytick.color": INK2, "xtick.labelsize": 7, "ytick.labelsize": 7,
    "legend.frameon": False, "grid.color": GRID, "grid.linewidth": 0.6,
})
_S.apply()
_S.install_autoscale(1.2)  # embed fonts + presentation-size bump

# 2022 open-water mean discharge per river [m^3/s] (for the mouth velocity)
Q_MEAN = {"colville": 239.0, "kuparuk": 63.8, "sagavanirktok": 47.6, "canning": 45.4}


def load_site(name):
    """Import sites/<name>.py (as a package, for its relative imports) and pull geometry."""
    import sys
    p = str(ROOT / "code")
    if p not in sys.path:
        sys.path.insert(0, p)
    m = importlib.import_module(f"sites.{name}")
    per = getattr(m, "B_lb_perchan", m.B_lb)
    nch = getattr(m, "N_CHAN_MOUTH", 1)
    return dict(B_lb=m.B_lb, B_ub=m.B_ub, L_FLARE=m.L_FLARE, EL=m.EL,
                DEPTH=m.DEPTH_lb, B_lb_perchan=per, N_CHAN=nch)


def flare(s, B_lb, B_ub, L_FLARE):
    """C-GEM flare width law: exponential B_lb->B_ub over 0..L_FLARE, then prismatic."""
    if s <= 0:
        return float(B_lb)
    if s >= L_FLARE:
        return float(B_ub)
    LC = -float(L_FLARE) / math.log(float(B_ub) / float(B_lb))
    return B_lb * math.exp(-s / LC)


def sword():
    try:
        return json.load(open(SWORD_JSON))
    except Exception:
        return {}


def page_width(pdf, geo, sw):
    fig = plt.figure(figsize=(11, 8.5))
    fig.suptitle("Channel width — idealized C-GEM geometry vs SWORD", x=0.055,
                 ha="left", fontsize=13, weight="bold")
    fig.text(0.055, 0.935, "Line = model flare (mouth width → prismatic upstream). "
             "Black = SWORD per-channel width (median + IQR). Triangle = SWORD delta-mouth "
             "distributary sum, the conveyance width the flare starts from (B_lb).",
             wrap=True,
             color=INK2, fontsize=8.0)
    gs = fig.add_gridspec(2, 2, left=0.07, right=0.97, top=0.90, bottom=0.09,
                          hspace=0.34, wspace=0.2)
    for k, s in enumerate(SITES):
        ax = fig.add_subplot(gs[k // 2, k % 2])
        g = geo[s]
        x = np.linspace(0, 27.0, 200)
        xm = x * 1000.0
        w_model = [flare(v, g["B_lb"], g["B_ub"], g["L_FLARE"]) for v in xm]
        ax.plot(x, w_model, color=C[s], lw=1.8, zorder=4,
                label=f"model flare (mouth {g['B_lb']:.0f} m → {g['B_ub']:.0f} m)")
        d = sw.get(s, {})
        # SWORD per-channel width: every node (faint scatter) + robust median line & IQR
        # band (main stem, non-trib). Braided raw totals are divided by n_chan_mod, so
        # this is the single-channel width the prismatic B_ub should match.
        if d.get("scatter_dist_km"):
            ax.scatter(d["scatter_dist_km"], d["scatter_width"], s=5, color=MUTED,
                       alpha=0.35, edgecolor="none", zorder=1,
                       label=f"SWORD per-channel nodes ({len(d['scatter_width'])})")
        if d.get("dist_km") and d.get("width_perchan_med"):
            dk = np.array(d["dist_km"]); wm = np.array(d["width_perchan_med"])
            ax.fill_between(dk, d["width_perchan_p25"], d["width_perchan_p75"],
                            color=MUTED, alpha=0.18, lw=0, zorder=2)
            ax.plot(dk, wm, ".-", ms=3, lw=1.0, color=INK, alpha=0.8, zorder=3,
                    label="SWORD per-channel median")
        # delta-mouth conveyance = distributary SUM, marked at x=0 (this is B_lb)
        if d.get("delta_sum_m"):
            ax.plot([0], [d["delta_sum_m"]], "v", ms=8, color=C[s], mec="white",
                    mew=0.8, zorder=5, label=f"SWORD delta sum ({d['delta_sum_m']} m)")
        ax.set_title(LABEL[s], loc="left", color=C[s], weight="medium")
        ax.set_xlabel("distance from mouth (km)", fontsize=7.5)
        ax.set_ylabel("width (m)", fontsize=7.5)
        ax.set_xlim(0, 27)
        for sp in ("top", "right"):
            ax.spines[sp].set_visible(False)
        ax.grid(True, alpha=0.5); ax.set_axisbelow(True)
        ax.legend(loc="upper right", fontsize=6.2, labelcolor=INK2)
    pdf.savefig(fig)
    plt.close(fig)


def page_table(pdf, geo, sw):
    fig = plt.figure(figsize=(11, 8.5))
    fig.suptitle("Channel geometry & hydraulics — parameters and provenance", x=0.055,
                 ha="left", fontsize=13, weight="bold")
    ax = fig.add_axes([0.05, 0.52, 0.9, 0.38]); ax.axis("off")
    cols = ["River", "Mouth width (model)", "Per-channel", "N_chan", "Upstream B_ub",
            "Flare", "Depth", "Mouth U", "W/H"]
    xc = [0.0, 0.13, 0.28, 0.38, 0.46, 0.58, 0.66, 0.74, 0.85]
    for x, c in zip(xc, cols):
        ax.text(x, 0.98, c, fontsize=7.2, weight="bold", color=INK2, va="top")
    for i, s in enumerate(SITES):
        g = geo[s]
        U = Q_MEAN[s] / (g["B_lb"] * g["DEPTH"])
        wh = g["B_lb"] / g["DEPTH"]
        widened = g["N_CHAN"] > 1
        vals = [LABEL[s], f"{g['B_lb']:.0f} m", f"{g['B_lb_perchan']:.0f} m",
                f"×{g['N_CHAN']}" if widened else "1 (single-ch.)",
                f"{g['B_ub']:.0f} m", f"{g['L_FLARE']/1000:.0f} km",
                f"{g['DEPTH']:.2f} m", f"{U:.2f} m/s", f"{wh:.0f}"]
        y = 0.86 - i * 0.11
        for x, v in zip(xc, vals):
            ax.text(x, y, v, fontsize=7.2, va="top", color=C[s])

    ax2 = fig.add_axes([0.05, 0.06, 0.9, 0.42]); ax2.axis("off")
    swstatus = ("LOADED — raw node widths overlaid on page 1." if sw
                else "PENDING — the SWORD v17 bundle (1.83 GB) is downloading; raw node "
                     "widths and the exact braided sum will be overlaid and the mouth "
                     "widths refined when it lands.")
    notes = [
        "WIDTH is width-averaged (1-D model, no lateral/vertical structure). Derived from SWORD v17c per-channel "
        "node widths (~200 m spacing) over the model domain, then a flare law: exponential convergence over "
        "L_FLARE down to a prismatic upstream width B_ub. Depth from at-a-station hydraulic geometry D = c·Q^f on "
        "USGS ADCP surveys, at each river's 2022 open-water mean discharge.",
        "DELTA-MOUTH WIDTH. The seaward boundary of a delta spans all active distributary channels, which diverge "
        "to separate sea outlets and convey the total discharge in PARALLEL — so the conveyance/salt-exchange "
        "width is the distributary SUM, not the per-channel width. Using per-channel width with total Q "
        "over-estimated the mouth velocity ~n_chan-fold and flushed the salt intrusion to ~0. Kuparuk and "
        "Sagavanirktok (braided/deltaic mouths) are therefore widened to the distributary sum; Colville and "
        "Canning are single-channel at the mouth (SWORD n_chan≈1) and unchanged.",
        "INTERIM vs SWORD. The widened mouths currently use per-channel × N_CHAN_MOUTH, a morphological "
        "distributary count (Kuparuk ×4 for Gwydyr Bay, Sag ×4). These are placeholders: the exact raw SWORD "
        "node-width sum replaces them once the extraction completes (sites/*.py SWORD_MOUTH_SUM).",
        f"SWORD data status: {swstatus}",
        "EFFECT. The widening lowers the Kuparuk/Sag mouth velocity from ~0.47–0.50 to ~0.10–0.12 m/s (comparable "
        "to Colville) and raises the Seo–Cheong dispersion (higher W/H), lengthening the salt intrusion toward the "
        "observed brackish grabs. It is a first-order change to velocity and dispersion, so it also shifts the "
        "carbon flux — see the validation PDF for the salinity and FCO₂ impact.",
    ]
    for i, n in enumerate(notes):
        ax2.text(0.0, 0.96 - i * 0.20, "•", transform=ax2.transAxes, color=MUTED, fontsize=9)
        ax2.text(0.02, 0.96 - i * 0.20, n, transform=ax2.transAxes, fontsize=7.2,
                 color=INK2, va="top")
    pdf.savefig(fig)
    plt.close(fig)


def main():
    OUT.parent.mkdir(exist_ok=True)
    geo = {s: load_site(s) for s in SITES}
    sw = sword()
    with PdfPages(OUT) as pdf:
        page_width(pdf, geo, sw)
        page_table(pdf, geo, sw)
    print(f"wrote {OUT.relative_to(ROOT)} ({OUT.stat().st_size/1024:.0f} kB, 2 pages)")


if __name__ == "__main__":
    main()
