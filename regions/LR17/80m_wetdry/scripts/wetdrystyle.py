"""Shared figure style for the pkg/wetdry figures and movies (Nature format;
matches ~/.claude/skills/nature-figures: no in-figure titles, the caption
carries the message; exponents and subscripts in mathtext, e.g. m s$^{-1}$,
kg m$^{-3}$, NO$_3$: Arial here lacks the Unicode ⁻ ₂ ₃ glyphs; ², ³, µ, ° are fine).

Nature figure guide: single column 89 mm, double column 183 mm, max height
170 mm; sans-serif text 5-7 pt; panel letters 8 pt bold lower case; lines
>= 0.25 pt; vector output preferred, raster at >= 300 dpi. Colours are
colourblind-safe (Okabe-Ito categorical, single-hue sequential, two-hue
diverging with a neutral midpoint).

    import wetdrystyle as ws
    ws.use()
    fig, ax = plt.subplots(figsize=ws.size(89, 60))
    ...
    ws.letter(ax, 'a')
    ws.save(fig, 'name')          # figures/<topic>/name.pdf + .png (600 dpi)
"""
import os

import matplotlib
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
FIGDIR = os.path.join(ROOT, 'figures')     # subfolders wetdry/, sediment/, colville/
MOVDIR = os.path.join(ROOT, 'movies')      # subfolders wetdry/, sediment/, colville/
MM = 1 / 25.4
SINGLE, ONEHALF, DOUBLE = 89, 120, 183          # mm

# Okabe-Ito, fixed order (never cycled)
OKABE = ['#0072B2', '#D55E00', '#009E73', '#CC79A7', '#E69F00',
         '#56B4E9', '#F0E442', '#000000']
BLUE, VERMILION, GREEN, PURPLE, ORANGE, SKY = OKABE[:6]
INK, INK2, MUTED, GRID = '#000000', '#404040', '#8a8a8a', '#d9d9d9'
DRY, BED, LAND, ICE = '#e3d9c6', '#b9ad94', '#7a7466', '#ffffff'
WATER = '#cfe3f3'

# Colormaps: cmocean (perceptually uniform, color-blind-safe), registered
# under the old names so the scripts need no change. TOPO stays a two-hue
# diverging map (blue below MSL, brown above) centered on MSL.
import cmocean  # noqa: E402


def _cmo(name, base):
    c = getattr(cmocean.cm, base).copy()
    c.name = name
    return c


SEQ = _cmo('wetdry_seq', 'deep')            # depth, speed: light -> dark blue
DIV = _cmo('wetdry_div', 'balance')         # signed differences
TOPO = LinearSegmentedColormap.from_list(         # bed: blue below MSL,
    'wetdry_topo', ['#08306b', '#2171b5', '#6baed6', '#c6dbef', '#f5f0e6',
                 '#dfc27d', '#bf812d', '#8c510a'])  # brown above
SEQ_SED = _cmo('wetdry_sed', 'turbid')      # sediment
SEQ_ICE = _cmo('wetdry_ice', 'ice')         # ice: dark (water) -> white
SEQ_TR = _cmo('wetdry_tr', 'dense')         # tracers
SEQ_S = _cmo('wetdry_s', 'haline')          # salinity

def use(movie=False):
    """Set rcParams. movie=True scales text for a 1920x1080 frame."""
    s = 1.0 if not movie else 1.9
    matplotlib.rcParams.update({
        'font.family': 'sans-serif',
        'font.sans-serif': ['Arial', 'Helvetica', 'DejaVu Sans'],
        'font.size': 7 * s, 'axes.titlesize': 7 * s,
        'axes.labelsize': 7 * s, 'xtick.labelsize': 6 * s,
        'ytick.labelsize': 6 * s, 'legend.fontsize': 6 * s,
        'axes.linewidth': 0.5 * s, 'lines.linewidth': 1.0 * s,
        'xtick.major.width': 0.5 * s, 'ytick.major.width': 0.5 * s,
        'xtick.minor.width': 0.4 * s, 'ytick.minor.width': 0.4 * s,
        'xtick.major.size': 2.5 * s, 'ytick.major.size': 2.5 * s,
        'xtick.minor.size': 1.5 * s, 'ytick.minor.size': 1.5 * s,
        'xtick.direction': 'out', 'ytick.direction': 'out',
        'axes.spines.top': False, 'axes.spines.right': False,
        'axes.edgecolor': INK, 'axes.labelcolor': INK,
        'xtick.color': INK, 'ytick.color': INK, 'text.color': INK,
        'axes.titleweight': 'normal', 'axes.titlepad': 3 * s,
        'axes.labelpad': 2 * s,
        'axes.prop_cycle': matplotlib.cycler(color=OKABE),
        'legend.frameon': False, 'legend.handlelength': 1.6,
        'figure.dpi': 150, 'savefig.dpi': 600,
        'savefig.bbox': 'tight', 'savefig.pad_inches': 0.02,
        'pdf.fonttype': 42, 'ps.fonttype': 42, 'svg.fonttype': 'none',
        'mathtext.default': 'regular',
        'image.cmap': 'wetdry_seq',
        'hatch.linewidth': 0.3 * s, 'hatch.color': INK2,
    })
    if movie:                     # frames are 19.2 x 10.8 in: 100 dpi = 1920 x 1080
        matplotlib.rcParams['savefig.dpi'] = 100
        matplotlib.rcParams['savefig.bbox'] = 'standard'
    for c in (SEQ, DIV, TOPO, SEQ_SED, SEQ_ICE, SEQ_TR, SEQ_S):
        try:
            matplotlib.colormaps.register(c)
        except ValueError:
            pass


def size(w_mm, h_mm):
    return (w_mm * MM, h_mm * MM)


def letter(ax, s, x=-0.02, y=1.02, movie=False):
    """Nature panel letter: 8 pt bold lower case, top left, outside."""
    ax.text(x, y, s, transform=ax.transAxes, fontsize=8 * (1.9 if movie else 1),
            fontweight='bold', va='bottom', ha='right')


panel = letter                     # nature-figures name


def cbar(fig, mappable, ax, label, **kw):
    kw.setdefault('fraction', 0.05)
    kw.setdefault('pad', 0.02)
    cb = fig.colorbar(mappable, ax=ax, **kw)
    cb.set_label(label)
    cb.outline.set_linewidth(0.4)
    cb.ax.tick_params(width=0.4, length=2)
    return cb


def topic(name):
    """Topic folder of a figure, still or movie called name:
    colville_* -> colville; LR17_* -> LR17; *mangrove* -> mangrove;
    fig5_*, fig6_*, *sediment* and *wave* -> sediment;
    everything else (wetting and drying) -> wetdry."""
    if name.startswith('colville'):
        return 'colville'
    if name.startswith('LR17'):
        return 'LR17'
    if 'mangrove' in name:
        return 'mangrove'
    if name.startswith(('fig5_', 'fig6_')) or 'sediment' in name or 'wave' in name:
        return 'sediment'
    return 'wetdry'


def figdir(name):
    """figures/<topic>/ for a figure or still called name (created)."""
    d = os.path.join(FIGDIR, topic(name))
    os.makedirs(d, exist_ok=True)
    return d


def movdir(name):
    """movies/<topic>/ for a movie called name (created)."""
    d = os.path.join(MOVDIR, topic(name))
    os.makedirs(d, exist_ok=True)
    return d


def save(fig, name, png_dpi=600):
    """Write figures/<topic>/<name>.pdf (vector) and .png (topic: figdir)."""
    d = figdir(name)
    for ext in ('pdf', 'png'):
        fig.savefig(os.path.join(d, f'{name}.{ext}'),
                    dpi=png_dpi if ext == 'png' else None)
    plt.close(fig)
    print('wrote', os.path.join(d, name + '.{pdf,png}'))
