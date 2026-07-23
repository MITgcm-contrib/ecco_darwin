"""
Shared publication / presentation style for all NS-RAD figures.

One place for the look of every PDF: embedded TrueType fonts (journal requirement),
a single complete-coverage typeface (DejaVu Sans — has subscripts/superscripts/arrows),
presentation-readable sizes, a cohesive palette, and the per-river colours used across
the diagnostics / validation / schematic figures.

Usage in a figure tool:

    import nsrad_style as S
    S.apply()                       # set rcParams (call once, before plotting)
    ...
    S.brand(fig)                    # NS-RAD footer tag (optional, per page)

`apply(scale=...)` multiplies every font size, so a tool can dial presentation size up.
`CGEM_FONT_SCALE` in the environment scales globally on top of that.
"""
import os
import matplotlib.pyplot as plt

# --- palette ---------------------------------------------------------------
INK, INK2, MUTED = "#1a1a1a", "#5c5c5c", "#9a9a96"
GRID, SURFACE, WARN = "#dcdcd7", "#ffffff", "#c0392b"
ACCENT = "#2f6f9e"                     # NS-RAD ocean-blue accent
# per-river colours (consistent across diagnostics / validation / movies / schematic)
RIVC = {"colville": "#2a78d6", "kuparuk": "#eb6834",
        "sagavanirktok": "#1baf7a", "canning": "#4a3aa7"}
LABEL = {"colville": "Colville", "kuparuk": "Kuparuk",
         "sagavanirktok": "Sagavanirktok", "canning": "Canning"}

_ENV_SCALE = float(os.environ.get("CGEM_FONT_SCALE", "1.0"))


def apply():
    """Set the shared look: embedded TrueType fonts, DejaVu Sans, NS-RAD colours. Does
    NOT change font SIZES — each tool keeps its tuned relative sizing; the uniform
    presentation bump is applied at save time by install_autoscale()."""
    plt.rcParams.update({
        "figure.facecolor": SURFACE, "axes.facecolor": SURFACE, "savefig.facecolor": SURFACE,
        "figure.dpi": 150,
        # DejaVu Sans: the one available sans with COMPLETE glyph coverage, embedded as
        # TrueType so the PDFs are journal-submission ready and render identically anywhere.
        "font.family": "sans-serif", "font.sans-serif": ["DejaVu Sans"],
        "pdf.fonttype": 42, "ps.fonttype": 42,
        "text.color": INK, "axes.labelcolor": INK2, "axes.edgecolor": GRID,
        "axes.linewidth": 0.9, "xtick.color": INK2, "ytick.color": INK2,
        "legend.frameon": False, "grid.color": GRID, "grid.linewidth": 0.6,
        "lines.solid_capstyle": "round",
    })


_installed = [False]


def install_autoscale(scale=1.28):
    """Monkeypatch PdfPages.savefig so EVERY text object (explicit-size or rcParams)
    is enlarged by `scale` just before the page is written — the uniform, presentation-
    readable bump, applied once per figure. Idempotent per figure and per process."""
    global _installed
    k = scale * _ENV_SCALE
    if _installed[0] or k == 1.0:
        return
    _installed[0] = True
    from matplotlib.backends.backend_pdf import PdfPages
    from matplotlib.text import Text
    _orig = PdfPages.savefig

    def savefig(self, figure=None, **kw):
        fig = figure if (figure is not None and hasattr(figure, "findobj")) else plt.gcf()
        for t in fig.findobj(match=Text):
            if not getattr(t, "_nsrad_scaled", False):
                t.set_fontsize(t.get_fontsize() * k)
                t._nsrad_scaled = True
        return _orig(self, figure, **kw)

    PdfPages.savefig = savefig


def tidy(ax, grid="y"):
    """Despine and add a light grid — the shared axis look."""
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)
    if grid:
        ax.grid(True, axis=grid, alpha=0.55)
    ax.set_axisbelow(True)


def brand(fig, extra=""):
    """Small NS-RAD tag at the figure foot (bottom-right)."""
    txt = "NS-RAD — North Slope River-Aquatic-Delta Model" + (f"  ·  {extra}" if extra else "")
    fig.text(0.99, 0.006, txt, ha="right", va="bottom", fontsize=6.8, color=MUTED,
             style="italic")
