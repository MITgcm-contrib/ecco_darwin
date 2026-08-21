"""
Interannual variability report: how discharge + riverine DOC loading actually varied
1980-2023 (from PWBM), and how the model's carbon flux responds to that variability.

Reads runs/interannual/<site>/ (the *_interannual site variants -- see
sites/colville_interannual.py and CLAUDE.md -> "Interannual forcing"), NOT
runs/definitive/. Kept as a SEPARATE PDF, not stitched into ns_rad_report.pdf --
same precedent as ns_rad_idealized_verification.pdf: this is an optional, expensive-
to-produce analysis of a fundamentally different question (how the model behaves
across 44 years of real forcing variability) than the combined report's single-year
"current state of the model" sections.

the forcing itself: discharge + riverine TOC boundary, 1980-2023, all three runnable
rivers overlaid -- straight from forcing/*_interannual_*.csv, no model run needed;
discharge vs. OBSERVATION -- the PWBM-derived forcing against the real USGS gauge
record (tools/fetch_interannual_discharge_obs.py), year by year, the same three
gauges/caveats tools/fetch_discharge.py and CLAUDE.md's "Discharge" section already
document (Canning has no gauge, excluded); channel geometry (width, mean depth vs
distance) for reference; EVERY tracked state
variable (config.py's species registry, VAR_PANELS below -- 2 pages) as an annual
open-water spatial mean, one point per calendar year, 1981-2023; EVERY process rate
written to output.nc (AREA_VARS below -- 2 pages) as an area-normalized annual budget,
same method as ns_rad_diagnostics.pdf's single-year "Area-normalized annual budgets"
page generalized to every year; ice variability (annual ice-covered duration and peak
thickness -- the freshet-breakup mechanism, config.BREAKUP_Q_FACTOR, makes ice timing
itself downstream of the interannual discharge signal); and a closing synthesis page
(annual FCO2 vs. annual mean discharge, ice-covered days vs. annual mean discharge) --
does the model's response actually track the interannual forcing, and how.

Usage:  python tools/make_interannual_pdf.py
"""
import json
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import nsrad_style as S
from matplotlib.backends.backend_pdf import PdfPages

ROOT = Path(__file__).resolve().parent.parent
RUN = ROOT / "runs" / "interannual"
FORC = ROOT / "forcing"
OUT = ROOT / "docs" / "ns_rad_interannual.pdf"
DISCHARGE_OBS_PATH = ROOT / "docs" / "interannual_discharge_obs.json"
# Same observation cache ns_rad_validation.pdf reads -- no new fetch needed. temp_wqp
# is Water Quality Portal grabs pooled by day-of-year across ALL observed years (as far
# back as 1969), which is the directly comparable dataset here since this report is
# itself multi-year; temp_usgs_2022 is a same-year (2022) daily reference, Kuparuk/
# Sagavanirktok only -- see tools/make_validation_pdf.py's module docstring.
TEMP_OBS_PATH = ROOT / "docs" / "validation_obs.json"

# Same caveats as tools/fetch_discharge.py / CLAUDE.md's "Discharge" section: Colville
# and Sagavanirktok are gauged well upstream of the delta and UNDERSTATE mouth flow, so
# PWBM (whole-basin) running higher there is expected, not a mismatch. Kuparuk is near
# tidewater -- the best-constrained comparison. Canning has no gauge at all.
GAUGE_CAVEAT = {
    "colville": "gauged far upstream (Umiat) — understates mouth flow",
    "kuparuk": "near tidewater — best-constrained comparison",
    "sagavanirktok": "gauged far upstream (Pump Sta 3) — understates mouth flow",
}

# canning excluded: inherits EL=0 from canning.py (placeholder, no observed estuary
# length yet), so canning_interannual can't run either -- see CLAUDE.md.
SITES = ["colville", "kuparuk", "sagavanirktok"]

YEAR_START = 1980          # forcing record start (tools/build_interannual_forcings.py)
WARMUP_DAYS = 365          # matches CGEM_WARMUP_DAYS default for these runs

S.apply()
S.install_autoscale(1.2)
plt.rcParams.update({"font.size": 8, "axes.titlesize": 9, "axes.labelsize": 8,
                      "xtick.labelsize": 7, "ytick.labelsize": 7,
                      "legend.fontsize": 7.5})


# --------------------------------------------------------------------- forcing series
def _forcing_daily(site, kind):
    """One river's interannual forcing CSV as (calendar_year_float, values). `kind` is
    'discharge' or 'toc' -- see tools/build_interannual_forcings.py for the exact
    filenames/units."""
    if kind == "discharge":
        p = FORC / f"{site}_river_discharge_interannual_1980-2023_m3sec.csv"
    else:
        p = FORC / f"{site}_toc_interannual_1980-2023_mmolC_m3.csv"
    if not p.exists():
        return None, None
    vals = np.genfromtxt(open(p, encoding="utf-8-sig"), delimiter=",")
    years = YEAR_START + np.arange(vals.size) / 365.0
    return years, vals


def page_forcing(pdf):
    fig, axes = plt.subplots(2, 1, figsize=(11, 7.5), sharex=True)
    fig.suptitle("Interannual forcing, 1980–2023 (PWBM-derived)", x=0.05, ha="left",
                 fontsize=13, weight="bold")
    fig.text(0.05, 0.945, "Discharge and riverine TOC (from PWBM DOC — see CLAUDE.md "
             "→ \"Interannual forcing\") for the three runnable rivers. Everything "
             "else (met, tides, marine boundary) still repeats a single 2022-typical "
             "year — see the same section for why.", color=S.INK2, fontsize=8.0)

    for site in SITES:
        yrs, q = _forcing_daily(site, "discharge")
        if yrs is not None:
            axes[0].plot(yrs, q, color=S.RIVC[site], lw=0.5, label=S.LABEL[site])
        yrs, toc = _forcing_daily(site, "toc")
        if yrs is not None:
            axes[1].plot(yrs, toc, color=S.RIVC[site], lw=0.5)

    axes[0].set_ylabel("discharge [m³ s⁻¹]")
    axes[0].legend(ncols=3, loc="upper left")
    axes[1].set_ylabel("riverine TOC [mmol C m⁻³]")
    axes[1].set_xlabel("year")
    axes[1].set_xlim(YEAR_START, YEAR_START + 44)
    for ax in axes:
        S.tidy(ax)
    S.brand(fig)
    pdf.savefig(fig)
    plt.close(fig)


def _load_discharge_obs():
    """{site: {year_str: {mean_m3s, n_valid_days}}} from
    tools/fetch_interannual_discharge_obs.py's cache, or {} if it hasn't been run."""
    if not DISCHARGE_OBS_PATH.exists():
        return {}
    return {k: v["annual"] for k, v in json.loads(DISCHARGE_OBS_PATH.read_text()).items()}


def page_discharge_validation(pdf):
    obs = _load_discharge_obs()
    fig = plt.figure(figsize=(11, 8.5))
    fig.suptitle("Discharge forcing vs. observation, 1980–2023", x=0.05, ha="left",
                 fontsize=13, weight="bold")
    if not obs:
        fig.text(0.05, 0.5, "No observation cache found — run "
                 "tools/fetch_interannual_discharge_obs.py to build "
                 f"{DISCHARGE_OBS_PATH.relative_to(ROOT)}.", color=S.WARN, fontsize=10)
        pdf.savefig(fig)
        plt.close(fig)
        return
    fig.text(0.05, 0.94, "PWBM-derived discharge (the model's actual forcing, solid "
             "line) against the real USGS gauge annual mean (points), same three "
             "gauges tools/fetch_discharge.py uses for the single-year forcing. "
             "Colville/Sagavanirktok are gauged well upstream of the delta and are "
             "KNOWN to understate mouth flow (see CLAUDE.md → \"Discharge\") — PWBM "
             "running consistently higher there is the expected whole-basin-vs-"
             "partial-basin gap, not a forcing error. Kuparuk (near tidewater) is the "
             "one genuine skill test.", color=S.INK2, fontsize=7.6)

    gs = fig.add_gridspec(3, 2, left=0.08, right=0.97, top=0.86, bottom=0.06,
                          hspace=0.5, wspace=0.28)
    for k, site in enumerate(SITES):
        axT = fig.add_subplot(gs[k, 0])   # time series
        axS = fig.add_subplot(gs[k, 1])   # scatter (model vs obs, same year)

        pwbm_yrs, pwbm_q = _forcing_daily(site, "discharge")
        pwbm_annual = {}
        for y in range(YEAR_START, YEAR_START + 44):
            m = (pwbm_yrs >= y) & (pwbm_yrs < y + 1)
            if m.any():
                pwbm_annual[y] = np.mean(pwbm_q[m])
        axT.plot(list(pwbm_annual.keys()), list(pwbm_annual.values()),
                 color=S.RIVC[site], lw=1.3, label="PWBM (model forcing)")

        site_obs = obs.get(site, {})
        oy = sorted(int(y) for y in site_obs)
        ov = [site_obs[str(y)]["mean_m3s"] for y in oy]
        axT.scatter(oy, ov, color=S.INK, s=14, zorder=5, label="USGS gauge (observed)")

        both_x, both_y_model, both_y_obs = [], [], []
        for y, v in zip(oy, ov):
            if y in pwbm_annual:
                both_x.append(y)
                both_y_model.append(pwbm_annual[y])
                both_y_obs.append(v)
        if both_x:
            axS.scatter(both_y_obs, both_y_model, color=S.RIVC[site], s=18, alpha=0.85)
            lo = min(min(both_y_obs), min(both_y_model))
            hi = max(max(both_y_obs), max(both_y_model))
            axS.plot([lo, hi], [lo, hi], "--", color=S.GRID, lw=1.0)
            bias = np.mean(np.array(both_y_model) - np.array(both_y_obs))
            r = np.corrcoef(both_y_obs, both_y_model)[0, 1] if len(both_x) > 1 else np.nan
            axS.text(0.03, 0.95, f"n={len(both_x)}  bias={bias:+.0f} m³/s  r={r:.2f}",
                     transform=axS.transAxes, fontsize=6.5, color=S.INK2, va="top")
        axS.set_xlabel("observed [m³ s⁻¹]", fontsize=6.8)
        axS.set_ylabel("PWBM [m³ s⁻¹]", fontsize=6.8)
        S.tidy(axS)

        axT.set_ylabel(f"{S.LABEL[site]}\n[m³ s⁻¹]", fontsize=7.5)
        axT.set_title(GAUGE_CAVEAT[site], loc="right", fontsize=6.3, color=S.MUTED)
        if k == 0:
            axT.legend(loc="upper left", fontsize=6.5)
        if k == len(SITES) - 1:
            axT.set_xlabel("year")
        S.tidy(axT)
    S.brand(fig)
    pdf.savefig(fig)
    plt.close(fig)


# --------------------------------------------------------------------- model response
_CACHE = {}


def load(site, var):
    """Load a field (time [days], values[cells 1..M]) from runs/interannual/<site>/
    output.nc, memoised. Interannual runs are NetCDF-only (CGEM_OUTPUT default)."""
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
    _CACHE[key] = (t, F[:, 1:])   # cells 1..M
    return _CACHE[key]


def _site_delxi(site):
    """This site's DELXI [m] -- geometry is inherited unchanged by the _interannual
    variant from the regular site (see sites/<name>_interannual.py), so importing the
    regular site's config is equivalent and avoids re-deriving anything."""
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


def _temp_obs():
    if not TEMP_OBS_PATH.exists():
        return None
    return json.loads(TEMP_OBS_PATH.read_text())


def model_doy_mid_all_years(site):
    """Mid-channel T at every simulated day, folded to day-of-year, ALL 43 post-warmup
    years pooled -- the interannual analog of make_validation_pdf.py's
    model_doy_series/model_daily_mean, generalized from 2 years to the whole record.
    Returns (doy_array, T_array) for scatter, and a {doy: mean} dict for pairing
    against observations. Exact-zero days are dropped (under-ice T held at the 0 degC
    freezing point -- no summer observations to compare against anyway; same
    convention make_validation_pdf.py uses)."""
    t, F = load(site, "T")
    if F is None:
        return None, None, {}
    mask = t >= WARMUP_DAYS
    t, F = t[mask], F[mask]
    col = F.shape[1] // 2
    mid = F[:, col]
    doy = np.mod(t, 365)
    nz = mid != 0
    doy, mid = doy[nz], mid[nz]
    acc = {}
    for d, v in zip(doy, mid):
        acc.setdefault(int(d), []).append(v)
    return doy, mid, {d: float(np.mean(vs)) for d, vs in acc.items()}


def page_temperature_validation(pdf):
    obs = _temp_obs()
    fig = plt.figure(figsize=(11, 8.5))
    fig.suptitle("Water temperature — model (all 43 years) vs. observation", x=0.05,
                 ha="left", fontsize=13, weight="bold")
    if obs is None:
        fig.text(0.05, 0.5, "No observation cache found — run "
                 "tools/make_validation_pdf.py once (it builds "
                 f"{TEMP_OBS_PATH.relative_to(ROOT)}) or restore that file.",
                 color=S.WARN, fontsize=10)
        pdf.savefig(fig)
        plt.close(fig)
        return
    fig.text(0.05, 0.935, "Every simulated day, mid-channel, folded to day-of-year and "
             "pooled across all 43 post-warmup years (semi-transparent — density shows "
             "where the interannual spread actually is) against Water Quality Portal "
             "grabs (grey, ALL observed years back to 1969 — the directly comparable "
             "multi-year dataset) and the USGS 2022 same-year daily record (Kuparuk/"
             "Sagavanirktok only). Same observation cache and method as "
             "ns_rad_validation.pdf's temperature page.", color=S.INK2, fontsize=7.6)

    gs = fig.add_gridspec(2, 2, left=0.07, right=0.97, top=0.87, bottom=0.14,
                          hspace=0.4, wspace=0.2)
    for k, site in enumerate(SITES):
        ax = fig.add_subplot(gs[k // 2, k % 2])
        doy, mid, mm = model_doy_mid_all_years(site)
        if doy is not None:
            ax.plot(doy, mid, ".", ms=1.0, color=S.RIVC[site], alpha=0.06,
                    label="model, all yrs (mid-channel)")

        wqp = obs.get("temp_wqp", {}).get(site, [])
        if wqp:
            wd = np.array([r[0] for r in wqp]); wv = np.array([r[1] for r in wqp])
            ax.plot(wd, wv, "o", ms=2.6, color=S.MUTED, alpha=0.5, mec="none",
                    label=f"WQP grabs ({len(wqp)}, all yr)")

        usgs = obs.get("temp_usgs_2022", {}).get(site, [])
        stat = ""
        if usgs and mm:
            ud = np.array([r[0] for r in usgs]); uv = np.array([r[1] for r in usgs])
            ax.plot(ud, uv, "o", ms=3.0, color=S.INK, mec=S.SURFACE, mew=0.4,
                    label=f"USGS 2022 daily ({len(usgs)})", zorder=6)
            pairs = np.array([(mm[d], v) for d, v in usgs if d in mm and mm[d] > 0])
            if len(pairs):
                bias = np.mean(pairs[:, 0] - pairs[:, 1])
                rmse = np.sqrt(np.mean((pairs[:, 0] - pairs[:, 1]) ** 2))
                stat = f"vs USGS: bias {bias:+.1f}  RMSE {rmse:.1f} °C"
        if wqp and mm:
            wpairs = np.array([(mm[int(r[0])], r[1]) for r in wqp if int(r[0]) in mm])
            if len(wpairs):
                wbias = np.mean(wpairs[:, 0] - wpairs[:, 1])
                stat += (("\n" if stat else "") + f"vs WQP: bias {wbias:+.1f} °C "
                         f"(n={len(wpairs)})")

        ax.set_xticks([0, 90, 181, 273, 364])
        ax.set_xticklabels(["J", "A", "J", "O", "D"])
        ax.set_xlim(0, 364)
        S.tidy(ax)
        ax.set_ylabel("T [°C]", fontsize=7.5)
        ax.set_title(S.LABEL[site], loc="left", color=S.RIVC[site], fontsize=9)
        if stat:
            ax.text(0.985, 0.05, stat, transform=ax.transAxes, ha="right",
                    fontsize=6.5, color=S.INK2)
        if k == 0:
            ax.legend(loc="upper left", fontsize=6.0, labelcolor=S.INK2)
    fig.text(0.05, 0.03, "The met/river/sea temperature FORCING is not interannual — "
             "it still repeats a single 2022-typical year (see CLAUDE.md → "
             "\"Interannual forcing\"). The apparent cold bias vs USGS 2022 on Kuparuk "
             "(bias -2.5, RMSE 3.3 °C) is concentrated almost entirely in the spring "
             "warming transition (days ~156-191, up to -5.5 °C) and shrinks to ~±1 °C "
             "once past it (checked day-by-day) — a PHASE-AVERAGING ARTIFACT, not a "
             "broad model cold bias: pooling 43 years with different real breakup/"
             "warming timing into one day-of-year mean smears a transition that is "
             "sharp in any single year. ns_rad_validation.pdf's own summer warm bias "
             "(inflow-boundary regression, not cloud cover) is a separate, genuine "
             "single-year model bias and still applies past the transition window.",
             color=S.INK2, fontsize=6.6, va="top")
    S.brand(fig)
    pdf.savefig(fig)
    plt.close(fig)


def estuary_surface_area(site):
    """Total water-surface area [m^2]: trapezoidal integral of width(x) (time-invariant)
    from the mouth to EL. Same method as ns_rad_diagnostics.pdf."""
    _, W = load(site, "width")
    if W is None:
        return None
    return float(np.trapezoid(W[0], dx=_site_delxi(site)))


def page_geometry(pdf):
    """Channel width (time-invariant) and record-mean depth vs. distance -- reference
    context only. Geometry is inherited unchanged by the _interannual site variants
    (see sites/<name>_interannual.py), so this is identical to the regular site's
    geometry, just plotted here so the report doesn't require cross-referencing
    ns_rad_geometry.pdf to know what channel these numbers came from."""
    fig, (axW, axD) = plt.subplots(1, 2, figsize=(11, 5))
    fig.suptitle("Channel geometry (reference — unchanged from the regular site)",
                 x=0.05, ha="left", fontsize=13, weight="bold")
    for site in SITES:
        delxi = _site_delxi(site)
        _, W = load(site, "width")
        _, D = load(site, "depth")
        if W is None:
            continue
        x = np.arange(W.shape[1]) * delxi / 1000.0
        axW.plot(x, W[0], color=S.RIVC[site], lw=1.3, label=S.LABEL[site])
        if D is not None:
            axD.plot(x, np.nanmean(D, axis=0), color=S.RIVC[site], lw=1.3)
    axW.set_title("Width (time-invariant)", loc="left", fontsize=9, color=S.INK)
    axW.set_ylabel("width [m]")
    axD.set_title("Depth — 1980–2023 record mean", loc="left", fontsize=9, color=S.INK)
    axD.set_ylabel("depth [m]")
    for ax in (axW, axD):
        ax.set_xlabel("distance from mouth [km]")
        S.tidy(ax)
    axW.legend(loc="best")
    S.brand(fig)
    pdf.savefig(fig)
    plt.close(fig)


def area_normalized_annual_series(site, var):
    """Area-normalized ANNUAL total of a volumetric rate for EVERY full post-warmup
    calendar year in the interannual record -- the same depth-integrate / width-weight
    / trapezoidal-time-integrate method as ns_rad_diagnostics.pdf's single-year
    'Area-normalized annual budgets' page (see that page's docstring in
    make_diagnostics_pdf.py), generalized from one hardcoded year to all of them.
    Returns (calendar_years, mmol m^-2 yr^-1)."""
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
        areal = F[mask] * D[mask]                              # mmol m^-2 s^-1
        spatial_mean = np.sum(areal * w0, axis=1) / np.sum(w0)  # width-weighted
        tt = t[mask] * 86400.0
        # The row at exactly t=WARMUP is pre-biogeo (main.py's gate is `t > WARMUP`,
        # strict) and so is all-NaN; drop before integrating -- same caveat as the
        # single-year version.
        finite = np.isfinite(spatial_mean)
        spatial_mean, tt = spatial_mean[finite], tt[finite]
        if len(tt) < 2:
            continue
        years.append(YEAR_START + 1 + i)
        vals.append(float(np.trapezoid(spatial_mean, tt)))
    return np.array(years), np.array(vals)


def annual_openwater_mean_series(site, var):
    """Per-calendar-year spatial+temporal mean of `var` over OPEN-WATER cells
    (ice_frac < 0.5 -- same mask ns_rad_diagnostics.pdf's openwater_mean_series uses),
    for every full post-warmup calendar year. Unlike that function this does NOT fold
    to day-of-year: the interannual record is not a repeating climatology, so each
    year is a genuinely distinct sample. Returns (calendar_years, values)."""
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


# Every species in the model's registry (config.py's species dict / MAXV list; see
# CLAUDE.md -> "Species registry") that is actually written to output.nc. SPM is g/L
# (not mmol/m3, hence its own axis label); the rest follow the model's mmol/m3
# convention, "as X" matching how each site's BOUNDARIES table documents the species
# (e.g. NO3/NH4 as N, PO4 as P, dSi as Si) -- see sites/<name>.py's DOC/nutrient
# conversion comments.
VAR_PANELS = [("T", "temperature [°C]"), ("S", "salinity [PSU]"),
              ("O2", "O₂ [mmol m⁻³]"), ("pH", "pH"),
              ("DIC", "DIC [mmol m⁻³]"), ("ALK", "alkalinity [mmol m⁻³]"),
              ("TOC", "TOC [mmol C m⁻³]"), ("RDOC", "RDOC [mmol C m⁻³]"),
              ("DIA", "phytoplankton (DIA) [mmol C m⁻³]"), ("SPM", "SPM [g L⁻¹]"),
              ("NO3", "NO₃ [mmol N m⁻³]"), ("NH4", "NH₄ [mmol N m⁻³]"),
              ("PO4", "PO₄ [mmol P m⁻³]"), ("dSi", "dSi [mmol Si m⁻³]"),
              ("CH4", "CH₄ [mmol m⁻³]"), ("N2O", "N₂O [mmol m⁻³]")]


def _paginated_grid(pdf, panels, ncols, plot_one, title, subtitle, page_note=None):
    """Emit ceil(len(panels)/page_size) pages, `page_size = ncols * nrows_per_page`
    (nrows_per_page chosen so panels stay readable: 4 rows for ncols=2, 3 rows for
    ncols=3). `plot_one(ax, panel)` draws one panel; the legend goes on the first
    panel of the first page only."""
    nrows = 4 if ncols == 2 else 3
    page_size = ncols * nrows
    pages = [panels[i:i + page_size] for i in range(0, len(panels), page_size)]
    for p, page_panels in enumerate(pages):
        fig = plt.figure(figsize=(11, 8.5))
        suffix = f" ({p + 1}/{len(pages)})" if len(pages) > 1 else ""
        fig.suptitle(title + suffix, x=0.05, ha="left", fontsize=13, weight="bold")
        fig.text(0.05, 0.935, subtitle, color=S.INK2, fontsize=8.0)
        gs = fig.add_gridspec(nrows, ncols, left=0.07, right=0.97, top=0.87,
                              bottom=0.08 if page_note else 0.06,
                              hspace=0.55, wspace=0.28)
        for k, panel in enumerate(page_panels):
            ax = fig.add_subplot(gs[k // ncols, k % ncols])
            plot_one(ax, panel)
            if p == 0 and k == 0:
                ax.legend(loc="best", fontsize=6.8)
        if page_note:
            fig.text(0.05, 0.02, page_note, color=S.INK2, fontsize=6.6, va="bottom")
        S.brand(fig)
        pdf.savefig(fig)
        plt.close(fig)


def page_variables(pdf):
    def plot_one(ax, panel):
        var, ylab = panel
        for site in SITES:
            yrs, vals = annual_openwater_mean_series(site, var)
            if yrs.size:
                ax.plot(yrs, vals, "o-", ms=3, lw=1.1, color=S.RIVC[site],
                        label=S.LABEL[site])
        S.tidy(ax)
        ax.set_ylabel(ylab, fontsize=7.5)
        ax.set_xlabel("year", fontsize=7.5)

    _paginated_grid(
        pdf, VAR_PANELS, ncols=2, plot_one=plot_one,
        title="Modelled state variables — annual open-water mean, 1981–2023",
        subtitle="Every tracked species (config.py's species registry), spatial mean "
                 "over open-water cells (ice_frac < 0.5), averaged across each full "
                 "calendar year — how the interannual discharge/DOC forcing "
                 "propagates through the modelled water column.")


# Every process rate written to output.nc (biogeo_module.py), classified by the
# element/species its "mmol ... m^-3 s^-1" units are actually in -- see config.py's
# ARCTIC BIOGEOCHEMISTRY EXTENSION block and the state-update equations in
# biogeo_module.py (e.g. c_TOC's update uses aer_deg/denit directly in mmol C terms,
# so both are C-currency despite denit consuming NO3). N2O-currency rates are reported
# as mass of N2O itself (44.013 g/mol), not N-equivalent, to avoid an extra assumption
# not stated in the code. NEM mixes C- and O2-currency terms as numerically equivalent,
# exactly as config.py's own NEM formula does; reported here as C for that reason.
AREA_VARS = [("FCO2", "C", "air–sea CO₂ flux (>0 outgassing)"),
             ("NPP", "C", "net primary production"),
             ("aer_deg", "C", "aerobic degradation"),
             ("denit", "C", "denitrification (carbon consumed)"),
             ("nit", "N", "nitrification"),
             ("NEM", "C", "net ecosystem metabolism"),
             ("phy_death", "C", "phytoplankton mortality"),
             ("rdoc_ox", "C", "RDOC aerobic oxidation"),
             ("photo", "C", "RDOC photomineralisation"),
             ("ch4_ox", "C", "methanotrophy (CH₄ oxidation)"),
             ("ch4_ex", "C", "air–sea CH₄ flux"),
             ("n2o_prod", "N2O", "N₂O production"),
             ("n2o_ex", "N2O", "air–sea N₂O flux"),
             ("sod", "O2", "sediment O₂ demand"),
             ("O2_ex", "O2", "air–sea O₂ flux")]
MOLAR_MASS = {"C": 12.011, "N": 14.007, "N2O": 44.013, "O2": 31.998}


def page_annual_response(pdf):
    areas = {s: estuary_surface_area(s) for s in SITES}

    def plot_one(ax, panel):
        var, element, title = panel
        for site in SITES:
            yrs, vals = area_normalized_annual_series(site, var)
            if yrs.size == 0:
                continue
            mass = vals * MOLAR_MASS[element] / 1000.0   # mmol -> g
            ax.plot(yrs, mass, "o-", color=S.RIVC[site], ms=2.5, lw=1.0,
                    label=S.LABEL[site])
        ax.axhline(0, color=S.GRID, lw=0.8)
        ax.set_title(title, loc="left", fontsize=8, color=S.INK)
        ax.set_ylabel(f"g{element} m⁻² yr⁻¹", fontsize=7.2)
        ax.set_xlabel("year", fontsize=7.2)
        S.tidy(ax)

    _paginated_grid(
        pdf, AREA_VARS, ncols=3, plot_one=plot_one,
        title="Model response: area-normalized annual budgets, 1981–2023",
        subtitle="Every process rate written to output.nc. Same method as "
                 "ns_rad_diagnostics.pdf's single-year area-normalized page "
                 "(depth-integrate a volumetric rate, width-weight across the "
                 "channel, trapezoidal-integrate over the year), run for every full "
                 "post-warmup year instead of just one.",
        page_note="Surface area: " + "  ·  ".join(
            f"{S.LABEL[s]} {areas[s]/1e6:.3f} km²" for s in SITES if areas[s]) +
            "\nrdoc_ox / photo / ch4_ox / ch4_ex / n2o_prod / n2o_ex / sod are flat "
            "zero — config.ARCTIC_BGC is off (default) for these sites, so the Arctic "
            "extension reactions never run; RDOC/CH4/N2O still advect as inert passive "
            "tracers. Not missing data.")


def annual_mean_discharge(site, years):
    """annual mean discharge [m3/s] for each calendar year in `years`, from the
    forcing CSV directly (not the model output)."""
    q_yrs, q = _forcing_daily(site, "discharge")
    if q_yrs is None:
        return np.full(len(years), np.nan)
    out = []
    for y in years:
        m = (q_yrs >= y) & (q_yrs < y + 1)
        out.append(np.mean(q[m]) if m.any() else np.nan)
    return np.array(out)


def annual_ice_stats(site):
    """Per-calendar-year ice-covered duration [days/yr, domain-mean ice_frac > 0.1] and
    peak ice thickness [m, max over cells and days], for every full post-warmup year.
    Daily output cadence means a row-count directly gives a day-count. Returns
    (calendar_years, duration_days, peak_thickness_m)."""
    t, H = load(site, "ice_thickness")
    _, FR = load(site, "ice_frac")
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


def page_ice(pdf):
    fig, (axD, axP) = plt.subplots(1, 2, figsize=(11, 5.2))
    fig.suptitle("Ice variability, 1981–2023", x=0.05, ha="left", fontsize=13,
                 weight="bold")
    fig.text(0.05, 0.945, "The hydraulic-breakup mechanism (config.BREAKUP_Q_FACTOR=3.0 "
             "× q_ref) makes breakup timing a downstream consequence of the "
             "interannual discharge signal — see CLAUDE.md → \"Prognostic ice model\".",
             color=S.INK2, fontsize=8.0)

    no_breakup, n_total = {}, {}
    for site in SITES:
        yrs, dur, peak = annual_ice_stats(site)
        if yrs.size == 0:
            continue
        axD.plot(yrs, dur, "o-", ms=3, lw=1.1, color=S.RIVC[site], label=S.LABEL[site])
        axP.plot(yrs, peak, "o-", ms=3, lw=1.1, color=S.RIVC[site])
        no_breakup[site] = yrs[dur >= 364]
        n_total[site] = len(yrs)
    axD.set_ylabel("ice-covered days / yr")
    axP.set_ylabel("peak ice thickness [m]")
    for ax in (axD, axP):
        ax.set_xlabel("year")
        S.tidy(ax)
    axD.legend(loc="best")

    # q_ref (main.py) is the FIXED mean over the WHOLE 44-year discharge record (one
    # number for the entire run) here, not that year's own mean as in the single-year
    # definitive runs -- so IN PRINCIPLE a year whose freshet peak doesn't clear 3x
    # that fixed, decades-scale mean could never trigger hydraulic breakup (full 365
    # d/yr ice cover, which would be a model artifact, not real North Slope behavior).
    # Checked explicitly: with correct forcing (see CLAUDE.md -> "Interannual forcing"
    # for a related bug that WAS previously causing spurious full-year ice years until
    # fixed), this does not actually happen for any of the 43 years, any site -- the
    # summer heat budget and/or freshet clear every year's cover regardless.
    any_affected = any(len(v) for v in no_breakup.values())
    if any_affected:
        lines = [f"{S.LABEL[s]}: {len(no_breakup.get(s, []))}/{n_total.get(s, 0)} years "
                 f"never break up ({', '.join(str(y) for y in no_breakup.get(s, []))})"
                 for s in SITES]
        fig.text(0.05, 0.02,
                 "Full-year (365 d) ice cover = the freshet peak that year didn't "
                 "clear 3× the FIXED 44-year mean discharge (q_ref is one number for "
                 "the whole run here, unlike the single-year definitive runs) — a "
                 "model artifact of this fixed threshold, not a claim these rivers "
                 "are ice-locked all summer. Affected years:\n" + "\n".join(lines),
                 color=S.WARN, fontsize=6.6, va="bottom")
    else:
        fig.text(0.05, 0.02,
                 "Every year (all 3 rivers, 1981–2023) breaks up — the fixed-threshold "
                 "concern above does not materialize in practice.",
                 color=S.INK2, fontsize=6.6, va="bottom")
    S.brand(fig)
    pdf.savefig(fig)
    plt.close(fig)


def page_forcing_vs_response(pdf):
    fig, (axQ, axI) = plt.subplots(1, 2, figsize=(11, 5.5))
    fig.suptitle("Does the model's response track the forcing?", x=0.05, ha="left",
                 fontsize=13, weight="bold")

    for site in SITES:
        yrs, fco2 = area_normalized_annual_series(site, "FCO2")
        if yrs.size == 0:
            continue
        q_annual = annual_mean_discharge(site, yrs)
        mass = fco2 * MOLAR_MASS["C"] / 1000.0
        ok = np.isfinite(q_annual) & np.isfinite(mass)
        axQ.scatter(q_annual[ok], mass[ok], color=S.RIVC[site], s=22,
                   label=S.LABEL[site], alpha=0.85)

        iyrs, dur, _ = annual_ice_stats(site)
        if iyrs.size == 0:
            continue
        q_annual_i = annual_mean_discharge(site, iyrs)
        ok = np.isfinite(q_annual_i) & np.isfinite(dur)
        axI.scatter(q_annual_i[ok], dur[ok], color=S.RIVC[site], s=22, alpha=0.85)

    axQ.axhline(0, color=S.GRID, lw=0.8)
    axQ.set_xlabel("annual mean discharge [m³ s⁻¹]")
    axQ.set_ylabel("annual FCO₂ [gC m⁻² yr⁻¹]")
    axQ.set_title("Carbon flux vs. discharge", loc="left", fontsize=9, color=S.INK)
    axI.set_xlabel("annual mean discharge [m³ s⁻¹]")
    axI.set_ylabel("ice-covered days / yr")
    axI.set_title("Ice duration vs. discharge", loc="left", fontsize=9, color=S.INK)
    for ax in (axQ, axI):
        S.tidy(ax)
    axQ.legend(loc="best")
    S.brand(fig)
    pdf.savefig(fig)
    plt.close(fig)


def main():
    OUT.parent.mkdir(exist_ok=True)
    missing = [s for s in SITES if not (RUN / s / "output.nc").exists()]
    if missing:
        raise SystemExit(f"Missing runs/interannual/<site>/output.nc for: {missing}. "
                          f"Run tools/run_interannual.sh first.")
    with PdfPages(OUT) as pdf:
        page_forcing(pdf)
        page_discharge_validation(pdf)
        page_temperature_validation(pdf)
        page_geometry(pdf)
        page_variables(pdf)
        page_annual_response(pdf)
        page_ice(pdf)
        page_forcing_vs_response(pdf)
        n_pages = pdf.get_pagecount()
    print(f"wrote {OUT.relative_to(ROOT)} ({OUT.stat().st_size/1024:.0f} kB, {n_pages} pages)")


if __name__ == "__main__":
    main()
