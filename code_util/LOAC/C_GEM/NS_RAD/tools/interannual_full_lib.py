"""
Loader library for the FULL-FORCING interannual variant (runs/interannual_full/<site>/,
the *_interannual_full site variants -- see sites/colville_interannual_full.py and
CLAUDE.md -> "Interannual forcing" -> "Full-forcing interannual (2005-2023)").

This is a DELIBERATELY LEAN, parallel counterpart to make_interannual_pdf.py -- not a
parametrized version of it, and not a full PDF-report builder. It exists purely to
give tools/make_publication_figures_full.py the same loader surface
(`load`, `freshet_metrics`, `area_normalized_annual_series`,
`annual_openwater_mean_series`, `_linregress_np`, `_trend_annotation`, `_site_delxi`,
`SITES`) that tools/make_publication_figures.py gets from `make_interannual_pdf as I`,
pointed at the 2005-2023/19-year record instead of 1980-2023/44-year. The original
make_interannual_pdf.py is NOT touched -- it still builds docs/ns_rad_interannual.pdf
from runs/interannual/ exactly as before; this is a separate, independent module so
that work carries zero risk of regressing it.

WHY A SEPARATE MODULE, NOT A PARAMETRIZED ORIGINAL: make_interannual_pdf.py is ~1000
lines, mostly PDF-page-rendering functions (page_forcing, page_discharge_validation,
page_freshet, etc.) that make_publication_figures.py never calls at all. Threading a
variant flag through the whole file risked subtle breakage in the PDF-report path for
no benefit; most of the functions actually reused by make_publication_figures.py
(area_normalized_annual_series, annual_openwater_mean_series) already derive their
year count from the data itself (`t.max()`), not a hardcoded constant, so they needed
NO changes at all -- copied verbatim. Only `_forcing_daily`, `freshet_metrics`, and the
module-level path/year constants actually needed adapting.
"""
import math
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent.parent
RUN = ROOT / "runs" / "interannual_full"
FORC = ROOT / "forcing"

SITES = ["colville", "kuparuk", "sagavanirktok"]  # canning excluded, same reason as always

YEAR_START = 2005
YEAR_END = 2023
N_YEARS = YEAR_END - YEAR_START + 1  # 19
WARMUP_DAYS = 365
BREAKUP_Q_FACTOR = 3.0  # matches config.py's hydraulic-breakup multiplier
MOLAR_MASS = {"C": 12.011, "N": 14.007, "N2O": 44.013, "O2": 31.998}  # copied verbatim
# from make_interannual_pdf.py -- used by make_publication_figures_full.py's Fig 4
# unit conversions (mmol -> g), fully generic, no year-count dependency.

_CACHE = {}


def _forcing_daily(site, kind):
    """One river's full-forcing-variant discharge/TOC CSV (2005-2023, sliced from the
    same 1980-2023 PWBM series make_interannual_pdf.py uses -- see
    tools/build_interannual_met.py's docstring for the window's full derivation) as
    (calendar_year_float, values)."""
    if kind == "discharge":
        p = FORC / f"{site}_river_discharge_interannual_{YEAR_START}-{YEAR_END}_m3sec.csv"
    else:
        p = FORC / f"{site}_toc_interannual_{YEAR_START}-{YEAR_END}_mmolC_m3.csv"
    if not p.exists():
        return None, None
    vals = np.genfromtxt(open(p, encoding="utf-8-sig"), delimiter=",")
    years = YEAR_START + np.arange(vals.size) / 365.0
    return years, vals


def load(site, var):
    """Load a field (time [days], values[cells 1..M]) from runs/interannual_full/<site>/
    output.nc, memoised."""
    key = (site, var)
    if key in _CACHE:
        return _CACHE[key]
    ncp = RUN / site / "output.nc"
    if not ncp.exists():
        _CACHE[key] = (None, None)
        return _CACHE[key]
    from netCDF4 import Dataset
    d = Dataset(ncp)
    if var not in d.variables:
        _CACHE[key] = (None, None)
        return _CACHE[key]
    t = np.asarray(d.variables["time"][:]) / 86400.0
    F = np.asarray(d.variables[var][:])
    _CACHE[key] = (t, F[:, 1:])
    return _CACHE[key]


def _site_delxi(site):
    """This site's DELXI [m] -- geometry is inherited unchanged by the
    _interannual_full variant from the regular site, same as the original
    _interannual variant, so importing the regular site's config is equivalent."""
    import importlib, sys, os, warnings
    p = str(ROOT / "code")
    if p not in sys.path:
        sys.path.insert(0, p)
    for m in ("config", "sites"):
        for mod in list(sys.modules):
            if mod == m or mod.startswith(m + "."):
                del sys.modules[mod]
    os.environ["CGEM_SITE"] = site
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        import config as _cfg
    return _cfg.DELXI


def area_normalized_annual_series(site, var):
    """Area-normalized ANNUAL total of a volumetric rate for every full post-warmup
    calendar year -- copied verbatim from make_interannual_pdf.py: `n_years` is
    already derived from `t.max()`, not a hardcoded year count, so this needed no
    change beyond the module-level YEAR_START it closes over."""
    t, F = load(site, var)
    _, D = load(site, "depth")
    _, W = load(site, "width")
    if F is None or D is None or W is None:
        return np.array([]), np.array([])
    w0 = W[0]
    n_years = int((t.max() - WARMUP_DAYS) // 365)
    years, vals = [], []
    for i in range(n_years):
        lo, hi = WARMUP_DAYS + i * 365, WARMUP_DAYS + (i + 1) * 365
        mask = (t >= lo) & (t < hi)
        if not mask.any():
            continue
        areal = F[mask] * D[mask]
        spatial_mean = np.sum(areal * w0, axis=1) / np.sum(w0)
        tt = t[mask] * 86400.0
        finite = np.isfinite(spatial_mean)
        spatial_mean, tt = spatial_mean[finite], tt[finite]
        if len(tt) < 2:
            continue
        years.append(YEAR_START + 1 + i)
        vals.append(float(np.trapezoid(spatial_mean, tt)))
    return np.array(years), np.array(vals)


def annual_openwater_mean_series(site, var):
    """Per-calendar-year spatial+temporal mean of `var` over open-water cells, for
    every full post-warmup calendar year -- copied verbatim, same reason as above."""
    t, F = load(site, var)
    _, ICE = load(site, "ice_frac")
    if F is None:
        return np.array([]), np.array([])
    n_years = int((t.max() - WARMUP_DAYS) // 365)
    years, vals = [], []
    for i in range(n_years):
        lo, hi = WARMUP_DAYS + i * 365, WARMUP_DAYS + (i + 1) * 365
        mask = (t >= lo) & (t < hi)
        if not mask.any():
            continue
        Fy = F[mask]
        if ICE is not None and ICE.shape == F.shape:
            m = float(np.nanmean(np.where(ICE[mask] < 0.5, Fy, np.nan)))
        else:
            m = float(np.nanmean(Fy))
        if np.isfinite(m):
            years.append(YEAR_START + 1 + i)
            vals.append(m)
    return np.array(years), np.array(vals)


def freshet_metrics(site):
    """Per-calendar-year spring-freshet characteristics straight from the daily
    discharge/TOC forcing -- same method as make_interannual_pdf.py, bound to this
    variant's 2005-2023/19-year window instead of 1980-2023/44."""
    q_yrs, q = _forcing_daily(site, "discharge")
    _, toc = _forcing_daily(site, "toc")
    if q_yrs is None:
        return {}
    year, peak_doy, peak_q, flashiness, duration, doc_load = [], [], [], [], [], []
    for y in range(YEAR_START, YEAR_START + N_YEARS):
        m = (q_yrs >= y) & (q_yrs < y + 1)
        if not m.any():
            continue
        qy = q[m]
        mean_q = qy.mean()
        i_peak = int(np.argmax(qy))
        year.append(y)
        peak_doy.append(i_peak)
        peak_q.append(float(qy[i_peak]))
        flashiness.append(float(qy[i_peak] / mean_q) if mean_q > 0 else np.nan)
        duration.append(float(np.sum(qy > BREAKUP_Q_FACTOR * mean_q)))
        if toc is not None:
            tocy = toc[m]
            mass_mmolC = np.trapezoid(qy * tocy, dx=86400.0)
            doc_load.append(mass_mmolC * 12.011 / 1e9)
        else:
            doc_load.append(np.nan)
    return dict(year=np.array(year), peak_doy=np.array(peak_doy, float),
                peak_q=np.array(peak_q), flashiness=np.array(flashiness),
                duration=np.array(duration), doc_load_tC=np.array(doc_load))


def _linregress_np(x, y):
    """OLS slope/intercept/r/two-sided-p for the slope -- copied verbatim from
    make_interannual_pdf.py (fully generic, no year-count dependency)."""
    x = np.asarray(x, float)
    y = np.asarray(y, float)
    ok = np.isfinite(x) & np.isfinite(y)
    x, y = x[ok], y[ok]
    n = x.size
    if n < 3:
        return np.nan, np.nan, np.nan, np.nan, n
    xm, ym = x.mean(), y.mean()
    sxx = np.sum((x - xm) ** 2)
    syy = np.sum((y - ym) ** 2)
    sxy = np.sum((x - xm) * (y - ym))
    if sxx == 0 or syy == 0:
        return np.nan, np.nan, np.nan, np.nan, n
    slope = sxy / sxx
    intercept = ym - slope * xm
    dof = n - 2
    resid = y - (slope * x + intercept)
    s2 = np.sum(resid ** 2) / dof if dof > 0 else np.nan
    se_slope = np.sqrt(s2 / sxx) if s2 == s2 else np.nan
    tstat = slope / se_slope if se_slope and se_slope > 0 else np.nan
    p = 2.0 * (1.0 - 0.5 * (1.0 + math.erf(abs(tstat) / math.sqrt(2.0)))) if tstat == tstat else np.nan
    r = sxy / np.sqrt(sxx * syy)
    return slope, intercept, r, p, n


def _trend_annotation(ax, x, y, unit_per_decade, color, row, fmt="{:+.2f}"):
    """Fit `_linregress_np`, draw a colored trend line, annotate -- copied verbatim."""
    slope, intercept, r, p, n = _linregress_np(x, y)
    if n >= 3 and np.isfinite(slope):
        xs = np.array([min(x), max(x)])
        ax.plot(xs, slope * xs + intercept, "--", color=color, lw=1.0, alpha=0.8)
        sig = "*" if p < 0.05 else ""
        ax.text(0.02, 0.03 + 0.075 * row,
                f"{fmt.format(slope * 10)} {unit_per_decade}/decade{sig}  "
                f"(r={r:.2f}, p={p:.2f}, n={n})",
                transform=ax.transAxes, fontsize=6.0, color=color, va="bottom")
    return slope, intercept, r, p, n
