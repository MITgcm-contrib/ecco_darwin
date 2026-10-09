"""Nature-journal figure style for matplotlib.

    import sys, os
    sys.path.insert(0, "<dir containing nature_style.py>")
    from nature_style import *          # applies the style on import

    fig, axs = plt.subplots(1, 2, figsize=size(DOUBLE, 70))
    ...
    panel(axs[0], "a"); panel(axs[1], "b")
    cbar(fig, axs[1], im, "Chl (mg m$^{-3}$)")
    save(fig, "outdir/fig1_chl")        # writes fig1_chl.pdf + fig1_chl.png

Nature artwork guide: 89 mm single / 183 mm double column, max 170 mm tall;
sans-serif 5-7 pt at final size; bold lowercase panel labels 8 pt; lines
0.25-1 pt; RGB; vector PDF with embedded fonts.
"""
import os

import matplotlib
from matplotlib import font_manager

if os.environ.get("MPLBACKEND") is None and not os.environ.get("DISPLAY"):
    matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

MM = 1 / 25.4
SINGLE, ONEHALF, DOUBLE = 89 * MM, 120 * MM, 183 * MM
MAX_H = 170 * MM

# Okabe-Ito: color-blind-safe categorical colors, in fixed order
OKABE_ITO = ["#0072B2", "#D55E00", "#009E73", "#CC79A7", "#E69F00", "#56B4E9", "#F0E442", "#000000"]
OBS, MOD, GRAY = "#0072B2", "#D55E00", "#7F7F7F"   # observations, model, reference
INK, LAND = "#222222", "#D9D4CC"


def _font():
    # Pleiades has no Arial; copies live in /nobackup/dcarrol2/fonts (2026-10-07).
    for d in (os.environ.get("NATURE_FONT_DIR", ""), "/nobackup/dcarrol2/fonts"):
        if d and os.path.isdir(d):
            for f in sorted(os.listdir(d)):
                if f.lower().endswith((".ttf", ".otf")):
                    try:
                        font_manager.fontManager.addfont(os.path.join(d, f))
                    except Exception:
                        pass
    have = {f.name for f in font_manager.fontManager.ttflist}
    for name in ("Arial", "Helvetica", "Liberation Sans", "DejaVu Sans"):
        if name in have:
            return name
    return "sans-serif"


STYLE = {
    "font.family": _font(), "font.size": 6.5, "axes.titlesize": 7, "axes.labelsize": 6.5,
    "xtick.labelsize": 6, "ytick.labelsize": 6, "legend.fontsize": 6,
    "axes.linewidth": 0.5, "xtick.major.width": 0.5, "ytick.major.width": 0.5,
    "xtick.minor.width": 0.4, "ytick.minor.width": 0.4,
    "xtick.major.size": 2.5, "ytick.major.size": 2.5, "xtick.minor.size": 1.5, "ytick.minor.size": 1.5,
    "xtick.direction": "out", "ytick.direction": "out", "lines.linewidth": 1.0, "lines.markersize": 3,
    "axes.edgecolor": INK, "axes.labelcolor": INK, "xtick.color": INK, "ytick.color": INK,
    "axes.spines.top": False, "axes.spines.right": False, "axes.grid": False,
    "axes.prop_cycle": matplotlib.cycler(color=OKABE_ITO),
    "image.cmap": "viridis", "legend.frameon": False,
    "pdf.fonttype": 42, "ps.fonttype": 42, "svg.fonttype": "none",
    "savefig.dpi": 600, "figure.dpi": 150, "mathtext.default": "regular",
}
# mathtext in the same face, so $\\mathit{Genus}$ and $^{-1}$ match the body text
# (default mathtext italic is DejaVu Sans)
STYLE.update({"mathtext.fontset": "custom", "mathtext.rm": STYLE["font.family"],
              "mathtext.it": STYLE["font.family"] + ":italic",
              "mathtext.bf": STYLE["font.family"] + ":bold"})
plt.rcParams.update(STYLE)


def size(width, height_mm):
    """figsize tuple from a width in inches (SINGLE/ONEHALF/DOUBLE) and a height in mm."""
    h = height_mm * MM
    if h > MAX_H + 1e-9:
        raise ValueError("Nature figures are at most 170 mm tall (got %.0f mm)" % height_mm)
    return (width, h)


def panel(ax, s, x=-0.02, y=1.02):
    """Bold lowercase panel label (a, b, c) at 8 pt, top left outside the axes."""
    ax.text(x, y, s, transform=ax.transAxes, fontsize=8, fontweight="bold",
            va="bottom", ha="right", color=INK)


def boxed(ax):
    """Restore top/right spines, e.g. for maps and images."""
    ax.spines["top"].set_visible(True)
    ax.spines["right"].set_visible(True)


def cbar(fig, ax_or_axes, mappable, label, **kw):
    """Thin colorbar with a units label, styled to match the axes."""
    cb = fig.colorbar(mappable, ax=ax_or_axes, shrink=kw.pop("shrink", 0.9),
                      pad=kw.pop("pad", 0.02), aspect=kw.pop("aspect", 22), **kw)
    cb.set_label(label)
    cb.outline.set_linewidth(0.5)
    cb.ax.tick_params(width=0.5, length=2)
    return cb


def save(fig, path, png=True):
    """Write <path>.pdf (vector, fonts embedded) and <path>.png (600 dpi); close the figure."""
    path = os.path.splitext(path)[0]
    d = os.path.dirname(path)
    if d:
        os.makedirs(d, exist_ok=True)
    fig.savefig(path + ".pdf", bbox_inches="tight", pad_inches=0.02)
    if png:
        fig.savefig(path + ".png", bbox_inches="tight", pad_inches=0.02)
    plt.close(fig)
    return path


__all__ = ["plt", "MM", "SINGLE", "ONEHALF", "DOUBLE", "MAX_H", "OKABE_ITO", "OBS", "MOD", "GRAY",
           "INK", "LAND", "STYLE", "size", "panel", "boxed", "cbar", "save"]
