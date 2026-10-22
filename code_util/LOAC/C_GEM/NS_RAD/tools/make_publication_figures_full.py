"""
FULL-FORCING counterpart to tools/make_publication_figures.py -- same four figure
sets, same method, but sourced from runs/interannual_full/<site>/output.nc (2005-2023,
19 post-warmup years, EVERY forcing category genuinely time-varying -- see
sites/colville_interannual_full.py and CLAUDE.md -> "Interannual forcing" ->
"Full-forcing interannual (2005-2023)") instead of runs/interannual/ (1980-2023,
discharge/DOC the only genuinely multi-year forcings). Written to
docs/publication_figures_full/, a separate directory -- the original script and its
docs/publication_figures/ output are completely untouched.

Fig 3's mouth tidal-aliasing calibration needs a 6-min-cadence reference from the SAME
2005-2023 forcing window (runs/definitive/ is a repeated-2022-climatology run, not
drawn from this period at all, so it cannot serve that role here). tools/
run_definitive_full.sh built runs/definitive_full/<site>/ for exactly this purpose: the
*_interannual_full site variants run at 6-min cadence for the first 730 days (real
2005-2006) of their already-built forcing -- no new forcing pipeline, just a different
CGEM_MAXT_DAYS/CGEM_TS. This script's `_dload()` is a small inline loader pointed at
that run (the full `make_diagnostics_pdf.py` dependency isn't needed -- only 4 call
sites, same pattern as interannual_full_lib.py's relationship to make_interannual_pdf.py).

Reuses tools/interannual_full_lib.py (as I) for the interannual loader/trend-analysis
surface, exactly as the original reuses make_interannual_pdf.py.

WHY THIS CANNOT BE A MECHANICAL YEAR-SUBSTITUTION OF THE ORIGINAL: in the 44-yr
record, every solute boundary except T and TOC is held CONSTANT, and water
temperature/wind have no interannual equivalent at all (shared single-year
climatology). In the full-forcing variant this is no longer true: the MARINE (mouth,
`clb`) boundary is now genuinely time-varying for every solute (S, DIC, ALK, NO3, NH4,
PO4, dSi, O2, TOC -- via ECCO-Darwin v5, 2005-2023), and water temperature (river AND,
per-river, sea) and wind are now genuinely multi-year and hourly-resolved. Only the
UPSTREAM (river, `cub`) boundary remains constant for every solute except T and TOC --
same limitation as before, just narrower in scope (marine side only used to be the one
exception; now it's the rule). Fig 1b and the Hovmoller/filtering caveat text are
rewritten below to reflect this accurately, not just relabeled with new years.

Four figure sets (file-for-file parallel to the original, same numbering):
  1. Forcing overview       -- fig1a (decadal discharge+TOC, 2005-2023), fig1b
                               (seasonal day-of-year climatology, ALL SIX panels now
                               genuinely multi-year: discharge, river water temp,
                               sea water temp [now per-river], wind [now hourly],
                               DOC -- see CONST_CAVEAT_FULL below)
  2. Hovmoller maps         -- fig2_hovmoller_<var>, 2006-2023, daily
  3. Estuarine filtering    -- fig3a/b + CSV, 2006-2023 mean, calibrated against
     mass balance             runs/definitive_full/ (real 2005-2006, 6-min)
  4. Freshet's role in      -- fig4a/b, 2006-2023 (18 post-warmup years)
     biogeochemistry

Usage:  python tools/make_publication_figures_full.py
"""
import csv
import sys
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import nsrad_style as S

sys.path.insert(0, str(Path(__file__).resolve().parent))
import interannual_full_lib as I   # noqa: E402  (full-forcing interannual loader)

ROOT = Path(__file__).resolve().parent.parent
OUT = ROOT / "docs" / "publication_figures_full"
SITES = ["colville", "kuparuk", "sagavanirktok"]   # canning excluded: EL=0, non-runnable

S.apply()   # idempotent


# =====================================================================================
# Small inline loader for runs/definitive_full/ -- the 6-min, real-2005-2006 run used
# ONLY to calibrate Fig 3's mouth tidal-aliasing correction (4 call sites below, same
# role make_diagnostics_pdf.D plays in the original script). Not a module import: the
# full diagnostics-PDF machinery isn't needed for 4 call sites.
# =====================================================================================
DRUN = ROOT / "runs" / "definitive_full"
_DCACHE = {}


def _dload(site, var):
    key = (site, var)
    if key in _DCACHE:
        return _DCACHE[key]
    ncp = DRUN / site / "output.nc"
    if not ncp.exists():
        _DCACHE[key] = (None, None)
        return _DCACHE[key]
    from netCDF4 import Dataset
    d = Dataset(ncp)
    if var not in d.variables:
        _DCACHE[key] = (None, None)
        return _DCACHE[key]
    t = np.asarray(d.variables["time"][:]) / 86400.0
    F = np.asarray(d.variables[var][:])
    _DCACHE[key] = (t, F[:, 1:])
    return _DCACHE[key]


def _save(fig, name, caption=None):
    OUT.mkdir(parents=True, exist_ok=True)
    if caption:
        fig.subplots_adjust(bottom=max(fig.subplotpars.bottom, 0.18))
        fig.text(0.05, 0.032, caption, ha="left", va="bottom", fontsize=6.0,
                  color=S.MUTED, style="italic", wrap=True)
    S.brand(fig)
    path = OUT / f"{name}.pdf"
    fig.savefig(path)
    plt.close(fig)
    print(f"  wrote {path.relative_to(ROOT)} ({path.stat().st_size / 1024:.0f} kB)")


CONST_CAVEAT_FULL = ("Upstream (river) boundaries are CONSTANT for every solute "
                      "except T and TOC -- same limitation as the original 44-yr "
                      "record. Marine (mouth) boundaries, by contrast, are "
                      "genuinely time-varying here for EVERY solute (S, DIC, ALK, "
                      "NO3, NH4, PO4, dSi, O2, TOC -- ECCO-Darwin v5, 2005-2023), "
                      "which is new relative to the original figures (only T was "
                      "marine-time-varying there). Structure seen here reflects "
                      "BOTH the model's own internal dynamics acting on a constant "
                      "upstream source AND a genuinely time-varying ocean signal "
                      "propagating in from the mouth. See CLAUDE.md -> \"Interannual "
                      "forcing\" -> \"Full-forcing interannual (2005-2023)\".")
SALT_CAVEAT = (" This is an ADVECTIVE-flux-only estimate. For S (salinity) "
               "specifically, dispersive/tidal-pumping transport -- the dominant "
               "salt-intrusion mechanism on this coast (see CLAUDE.md -> \"Width and "
               "dispersion\") -- is NOT captured, so S's entry likely understates "
               "true salt exchange; the biogeochemical solutes are advection-"
               "dominated and less affected by this limitation.")


# =====================================================================================
# Fig 1 -- forcing overview
# =====================================================================================
def fig1a_forcing_decadal():
    fig, axes = plt.subplots(2, 1, figsize=(9, 6), sharex=True)
    fig.suptitle(f"River forcing, {I.YEAR_START}–{I.YEAR_END}", x=0.05, ha="left",
                 fontsize=13, weight="bold")
    fig.text(0.05, 0.94, "Discharge and upstream DOC (as TOC), the full 19-yr "
             "full-forcing window (see 1b for every other forcing, also genuinely "
             "multi-year here).", color=S.INK2, fontsize=8.0)
    for site in SITES:
        yrs, q = I._forcing_daily(site, "discharge")
        if yrs is not None:
            axes[0].plot(yrs, q, color=S.RIVC[site], lw=0.5, label=S.LABEL[site])
        yrs, toc = I._forcing_daily(site, "toc")
        if yrs is not None:
            axes[1].plot(yrs, toc, color=S.RIVC[site], lw=0.5)
    axes[0].set_ylabel("discharge [m³ s⁻¹]")
    axes[0].legend(ncols=3, loc="upper left", fontsize=8)
    axes[1].set_ylabel("upstream DOC, as TOC [mmol C m⁻³]")
    axes[1].set_xlabel("year")
    axes[1].set_xlim(I.YEAR_START, I.YEAR_END + 1)
    for ax in axes:
        S.tidy(ax)
    _save(fig, "fig1a_forcing_decadal",
          f"Source: tools/build_interannual_forcings.py, sliced to the full-forcing "
          f"window ({I.YEAR_START}-{I.YEAR_END}) from Michael's daily PWBM product. "
          "No new model run needed for this figure.")


def _forc_full(fname):
    p = ROOT / "forcing" / fname
    return np.genfromtxt(open(p, encoding="utf-8-sig"), delimiter=",")


def _climatology_values(values, n_per_day=1):
    """Day-of-year climatology (doy[0..364], mean, p10, p90) across the full
    2005-2023 full-forcing record. `values` may already be daily (n_per_day=1) or
    hourly (n_per_day=24 -- reshape to a daily mean first: 19*8760 hourly rows =
    19*365 days x 24 hours exactly, no leap days, matching
    tools/build_interannual_met_hourly.py). The record is exactly N_YEARS*365 daily
    rows with no leap days, so day-of-year = index % 365 recovers the calendar day
    cleanly, same as the original _climatology()."""
    if n_per_day > 1:
        values = values.reshape(-1, n_per_day).mean(axis=1)
    doy_of_row = np.arange(values.size) % 365
    mean = np.empty(365); p10 = np.empty(365); p90 = np.empty(365)
    for d in range(365):
        vals = values[doy_of_row == d]
        mean[d] = np.mean(vals)
        p10[d] = np.percentile(vals, 10)
        p90[d] = np.percentile(vals, 90)
    return np.arange(365), mean, p10, p90


def _climatology(site, kind):
    _, v = I._forcing_daily(site, kind)
    return _climatology_values(v)


def fig1b_forcing_seasonal():
    fig, axes = plt.subplots(5, 1, figsize=(9, 11.5), sharex=True)
    fig.subplots_adjust(top=0.90)
    fig.suptitle(f"River forcing, seasonal climatology ({I.YEAR_START}–"
                 f"{I.YEAR_END})", x=0.05, ha="left", fontsize=13, weight="bold")
    fig.text(0.05, 0.965, "Day-of-year mean (line) + 10-90th percentile (shaded), "
             f"all {I.N_YEARS} years of the full-forcing record -- EVERY panel here "
             "is genuinely multi-year, unlike the original 44-yr figures (where "
             "water temp/wind were a shared single-year climatology). River water "
             "temp and wind are a shared regional record across all 3 rivers (same "
             "source station); sea temperature is per-river (ECCO-Darwin v5, each "
             "river's own nearest-wet-cell point).",
             color=S.INK2, fontsize=7.6, va="top")
    axQ, axTr, axTs, axW, axC = axes

    DISCHARGE_FLOOR = 1e-2   # m3/s; log-scale display floor, see original fig1b
    for site in SITES:
        doy, mean, p10, p90 = _climatology(site, "discharge")
        axQ.fill_between(doy, np.maximum(p10, DISCHARGE_FLOOR),
                          np.maximum(p90, DISCHARGE_FLOOR), color=S.RIVC[site],
                          alpha=0.2, linewidth=0)
        axQ.plot(doy, np.maximum(mean, DISCHARGE_FLOOR), color=S.RIVC[site], lw=1.3,
                 label=S.LABEL[site])
    axQ.set_yscale("log")
    axQ.set_ylabel("discharge [m³ s⁻¹]")
    axQ.legend(loc="upper left", ncols=3, fontsize=7.5)

    river_t = _forc_full("river_temp_interannual_2005-2023_degC.csv")
    doy, mean, p10, p90 = _climatology_values(river_t)
    axTr.fill_between(doy, p10, p90, color="#2a78d6", alpha=0.2, linewidth=0)
    axTr.plot(doy, mean, color="#2a78d6", lw=1.3)
    axTr.axhline(0, color=S.GRID, lw=0.8)
    axTr.set_ylabel("river water temp\n(upstream bnd) [°C]")

    for site in SITES:
        sea_t = _forc_full(f"{site}_t_marine_interannual_2005-2023.csv")
        doy, mean, p10, p90 = _climatology_values(sea_t)
        axTs.fill_between(doy, p10, p90, color=S.RIVC[site], alpha=0.2, linewidth=0)
        axTs.plot(doy, mean, color=S.RIVC[site], lw=1.3, label=S.LABEL[site])
    axTs.axhline(0, color=S.GRID, lw=0.8)
    axTs.set_ylabel("sea water temp\n(marine bnd) [°C]")
    axTs.legend(loc="upper left", fontsize=7.5, ncols=3)

    wind = _forc_full("windspeed_hourly_interannual_2005-2023_msec.csv")
    doy, mean, p10, p90 = _climatology_values(wind, n_per_day=24)
    axW.fill_between(doy, p10, p90, color=S.MUTED, alpha=0.2, linewidth=0)
    axW.plot(doy, mean, color=S.MUTED, lw=1.3)
    axW.set_ylabel("wind speed [m s⁻¹]\n(shared, all rivers)")

    for site in SITES:
        doy, mean, p10, p90 = _climatology(site, "toc")
        axC.fill_between(doy, p10, p90, color=S.RIVC[site], alpha=0.2, linewidth=0)
        axC.plot(doy, mean, color=S.RIVC[site], lw=1.3, label=S.LABEL[site])
    axC.set_ylabel("upstream DOC, as TOC\n[mmol C m⁻³]")
    axC.set_xlabel("day of year")
    axC.legend(loc="upper right", fontsize=7.5, ncols=3)

    for ax in axes:
        S.tidy(ax)
        ax.set_xlim(0, 364)
        ax.set_xticks([0, 90, 181, 273, 364])
        ax.set_xticklabels(["J", "A", "J", "O", "D"])
    _save(fig, "fig1b_forcing_seasonal",
          f"Day-of-year mean + 10-90th percentile across all {I.N_YEARS} years of "
          "forcing/*_interannual_2005-2023_*.csv (tools/build_interannual_*.py). "
          "Wind is hourly, aggregated to a daily mean before the day-of-year fold "
          "(24 hours/day exactly, no leap days). No new model run needed.")


# =====================================================================================
# Fig 2 -- Hovmoller maps
# =====================================================================================
HOVMOLLER_VARS = [
    ("FCO2",  True,  "RdBu_r",  "diverging", "air–sea CO₂ flux [mmol C m⁻³ s⁻¹] (>0 outgassing)"),
    ("TOC",   False, "YlOrBr",  "sequential", "DOC, as TOC [mmol C m⁻³]"),
    ("DIC",   False, "PuBuGn",  "sequential", "DIC [mmol C m⁻³]"),
    ("REMIN", True,  "magma",   "sequential", "remineralization = aer_deg+denit [mmol C m⁻³ s⁻¹]"),
    ("NPP",   True,  "Greens",  "sequential", "primary production [mmol C m⁻³ s⁻¹]"),
]


def _load_var(loader, site, name):
    if name == "REMIN":
        t1, F1 = loader(site, "aer_deg")
        t2, F2 = loader(site, "denit")
        if F1 is None or F2 is None:
            return None, None
        return t1, F1 + F2
    return loader(site, name)


def _hovmoller_prep(name, site):
    t, F = _load_var(I.load, site, name)
    if F is None:
        return None
    mask = t >= I.WARMUP_DAYS   # post-warmup, 2006-2023
    tt, Fm = t[mask], F[mask]
    x = I.YEAR_START + tt / 365.0
    step = max(1, len(x) // 3000)   # ~daily-scale render cost cap over 19 years
    x, Fm = x[::step], Fm[::step]
    delxi_km = I._site_delxi(site) / 1000.0
    y = np.arange(F.shape[1]) * delxi_km
    return x, y, Fm


def _hovmoller_figure(name, is_rate, cmap, scale_kind, label):
    panels = {}
    for site in SITES:
        res = _hovmoller_prep(name, site)
        if res is not None:
            panels[site] = res
    if not panels:
        return

    pooled = np.concatenate([Fm.ravel() for _, _, Fm in panels.values()])
    pooled = pooled[np.isfinite(pooled)]
    if scale_kind == "diverging":
        vext = np.nanpercentile(np.abs(pooled), 99) if pooled.size else 1.0
        vmin, vmax = -vext, vext
    else:
        vmin = 0.0
        vmax = np.nanpercentile(pooled, 99) if pooled.size else 1.0
        if vmax <= vmin:
            vmax = vmin + 1.0

    fig, axes = plt.subplots(len(SITES), 1, figsize=(9, 8.5), sharex=True)
    fig.subplots_adjust(bottom=0.18)
    fig.suptitle(f"{label} — distance × time ({I.YEAR_START + 1}–"
                 f"{I.YEAR_END}, daily)", x=0.05, ha="left", fontsize=12.5,
                 weight="bold")
    im = None
    for ax, site in zip(axes, SITES):
        if site not in panels:
            ax.axis("off")
            continue
        x, y, Fm = panels[site]
        im = ax.pcolormesh(x, y, Fm.T, cmap=cmap, vmin=vmin, vmax=vmax,
                            shading="auto", rasterized=True)
        ax.set_ylabel(f"{S.LABEL[site]}\ndistance from\nmouth [km]", fontsize=7.5,
                       color=S.RIVC[site])
        S.tidy(ax, grid=None)
    axes[-1].set_xlabel("year")
    if im is not None:
        cb = fig.colorbar(im, ax=list(axes), pad=0.02, fraction=0.04, extend=(
            "both" if scale_kind == "diverging" else "max"))
        cb.ax.tick_params(labelsize=7)
    fig.text(0.05, 0.032, CONST_CAVEAT_FULL, ha="left", va="bottom", fontsize=5.8,
              color=S.MUTED, style="italic", wrap=True)
    S.brand(fig)
    path = OUT / f"fig2_hovmoller_{name.lower()}.pdf"
    OUT.mkdir(parents=True, exist_ok=True)
    fig.savefig(path)
    plt.close(fig)
    print(f"  wrote {path.relative_to(ROOT)} ({path.stat().st_size / 1024:.0f} kB)")


def make_hovmoller_figures():
    for name, is_rate, cmap, scale_kind, label in HOVMOLLER_VARS:
        _hovmoller_figure(name, is_rate, cmap, scale_kind, label)


# =====================================================================================
# Fig 3 -- estuarine filtering mass balance
# =====================================================================================
TRANSPORTED_SPECIES = ["DIA", "dSi", "NO3", "NH4", "PO4", "O2", "TOC", "RDOC", "CH4",
                        "N2O", "S", "SPM", "DIC", "ALK", "T"]

SPECIES_SPEC = {
    "DIA": (12.011, "mmol"), "dSi": (28.085, "mmol"), "NO3": (14.007, "mmol"),
    "NH4": (14.007, "mmol"), "PO4": (30.974, "mmol"), "O2": (31.998, "mmol"),
    "TOC": (12.011, "mmol"), "RDOC": (12.011, "mmol"), "CH4": (12.011, "mmol"),
    "N2O": (44.013, "mmol"), "DIC": (12.011, "mmol"), "ALK": (None, "mmol_eq"),
    "SPM": (None, "gL"), "S": (None, "psu"), "T": (None, "degC"),
}


def _boundary_series_definitive(site, name, end, k=2):
    """(t_seconds, Q_seaward[m3/s], C[native units]) from runs/definitive_full/, the
    first post-warmup year (day 365-730, i.e. real 2006), native 6-min cadence. Same
    reconstruction/sign-flip method as the original script's
    _boundary_series_definitive -- see that docstring for the full derivation."""
    t, Vel = _dload(site, "velocity")
    _, Dep = _dload(site, "depth")
    _, W = _dload(site, "width")
    _, C = _dload(site, name)
    if any(x is None for x in (Vel, Dep, W, C)):
        return None
    M = Vel.shape[1]
    col = k if end == "mouth" else M - 1 - k
    mask = (t >= I.WARMUP_DAYS) & (t < I.WARMUP_DAYS + 365.0)
    tt = t[mask] * 86400.0
    Q = -(Vel[mask, col] * Dep[mask, col] * W[mask, col])
    Cc = C[mask, col]
    finite = np.isfinite(Q) & np.isfinite(Cc)
    tt, Q, Cc = tt[finite], Q[finite], Cc[finite]
    if len(tt) < 2:
        return None
    return tt, Q, Cc


def _boundary_series_interannual(site, name, end, k=2):
    """(t_days, Q_seaward[m3/s], C[native units]) from runs/interannual_full/, ALL
    post-warmup years (2006-2023) at once, native daily cadence."""
    t, Vel = I.load(site, "velocity")
    _, Dep = I.load(site, "depth")
    _, W = I.load(site, "width")
    _, C = I.load(site, name)
    if any(x is None for x in (Vel, Dep, W, C)):
        return None
    M = Vel.shape[1]
    col = k if end == "mouth" else M - 1 - k
    mask = t >= I.WARMUP_DAYS
    tt = t[mask]
    Q = -(Vel[mask, col] * Dep[mask, col] * W[mask, col])
    Cc = C[mask, col]
    finite = np.isfinite(Q) & np.isfinite(Cc)
    tt, Q, Cc = tt[finite], Q[finite], Cc[finite]
    if len(tt) < 2:
        return None
    return tt, Q, Cc


MOUTH_ALIAS_STRIDE = 240   # 86400s / 360s -- see original script's derivation
_CORR_CACHE = {}


def mouth_bias_correction(site, name, k=2):
    key = (site, name, k)
    if key in _CORR_CACHE:
        return _CORR_CACHE[key]
    res = _boundary_series_definitive(site, name, "mouth", k)
    if res is None:
        _CORR_CACHE[key] = None
        return None
    tt, Q, C = res
    _, kind = SPECIES_SPEC[name]
    flux = Q * C * 1000.0 if kind == "gL" else Q * C
    true_mean = float(np.trapezoid(flux, tt) / (tt[-1] - tt[0]))
    daily_mean = float(np.mean(flux[::MOUTH_ALIAS_STRIDE]))
    corr = true_mean / daily_mean if daily_mean and np.isfinite(daily_mean) else None
    _CORR_CACHE[key] = corr
    return corr


def annual_flux_loads(site, name, end, k=2, apply_correction=False):
    res = _boundary_series_interannual(site, name, end, k)
    if res is None:
        return np.array([]), np.array([])
    tt, Q, C = res
    _, kind = SPECIES_SPEC[name]
    flux = Q * C * 1000.0 if kind == "gL" else Q * C
    if apply_correction:
        corr = mouth_bias_correction(site, name, k)
        if corr is not None:
            flux = flux * corr
    n_years = int((tt.max() - I.WARMUP_DAYS) // 365)
    years, vals = [], []
    for i in range(n_years):
        lo, hi = I.WARMUP_DAYS + i * 365, I.WARMUP_DAYS + (i + 1) * 365
        m = (tt >= lo) & (tt < hi)
        if m.sum() < 2:
            continue
        years.append(I.YEAR_START + 1 + i)
        vals.append(float(np.trapezoid(flux[m], tt[m] * 86400.0)))
    return np.array(years), np.array(vals)


def _report(name, raw):
    if raw is None:
        return None, ""
    mm, kind = SPECIES_SPEC[name]
    if kind == "gL":
        return raw / 1e6, "t/yr"
    if kind == "psu":
        return raw / 1e3, "t/yr"
    if kind == "degC":
        return raw / 1e6, "°C·hm³/yr"
    if kind == "mmol_eq":
        return raw / 1e6, "kmol/yr"
    return raw * mm / 1e9, "t/yr"


def validate_discharge_reconstruction():
    print(f"Discharge reconstruction sanity check (reconstructed vs. forcing, "
          f"{I.YEAR_START + 1}-{I.YEAR_END} mean):")
    for site in SITES:
        res = _boundary_series_interannual(site, "S", "upstream", k=2)
        if res is None:
            print(f"  {site:14s} -- no data"); continue
        _, Q, _ = res
        recon = float(np.mean(Q))
        yrs, q = I._forcing_daily(site, "discharge")
        forc = float(np.mean(q[yrs >= I.YEAR_START + 1]))
        pct = 100.0 * (recon - forc) / forc if forc else float("nan")
        print(f"  {site:14s} reconstructed={recon:9.2f}  forcing={forc:9.2f}  "
              f"diff={pct:+.1f}%")

    print("Mouth tidal-aliasing correction (definitive_full-run calibration, "
          "true 6-min mean / daily-sampled mean):")
    for site in SITES:
        line = [f"{site:14s}"]
        for name in ("TOC", "DIC", "O2", "NO3"):
            corr = mouth_bias_correction(site, name, k=2)
            line.append(f"{name}={corr:.3f}" if corr is not None else f"{name}=n/a")
        print("  " + "  ".join(line))


def filtering_mass_balance():
    rows = []
    for site in SITES:
        for name in TRANSPORTED_SPECIES:
            for k in (0, 2):
                corr = mouth_bias_correction(site, name, k)
                yrs_in, raw_in_series = annual_flux_loads(site, name, "upstream", k)
                yrs_out, raw_out_series = annual_flux_loads(site, name, "mouth", k,
                                                              apply_correction=True)
                in_by_year = dict(zip(yrs_in, raw_in_series))
                out_by_year = dict(zip(yrs_out, raw_out_series))
                for y in sorted(set(yrs_in) & set(yrs_out)):
                    val_in, unit = _report(name, in_by_year[y])
                    val_out, _ = _report(name, out_by_year[y])
                    net = val_out - val_in
                    pct = net / abs(val_in) * 100.0 if val_in != 0.0 else None
                    rows.append(dict(site=site, species=name, k=k, year=int(y),
                                      input_load=val_in, output_load=val_out,
                                      net_effect=net, pct_change=pct, unit=unit,
                                      mouth_correction=corr))
    return rows


def summarize_by_site_species(rows, k=2):
    out = []
    for site in SITES:
        for name in TRANSPORTED_SPECIES:
            srows = [r for r in rows if r["site"] == site and r["species"] == name
                     and r["k"] == k]
            if not srows:
                continue
            val_in = float(np.mean([r["input_load"] for r in srows]))
            val_out = float(np.mean([r["output_load"] for r in srows]))
            net = val_out - val_in
            pct = net / abs(val_in) * 100.0 if val_in != 0.0 else None
            out.append(dict(site=site, species=name, input_load=val_in,
                             output_load=val_out, net_effect=net, pct_change=pct,
                             unit=srows[0]["unit"], n_years=len(srows)))
    return out


def write_filtering_csv(rows):
    OUT.mkdir(parents=True, exist_ok=True)
    path = OUT / "filtering_mass_balance.csv"
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["site", "species", "k", "year", "input_load",
                                            "output_load", "net_effect", "pct_change",
                                            "unit", "mouth_correction"])
        w.writeheader()
        for r in rows:
            w.writerow(r)
    print(f"  wrote {path.relative_to(ROOT)} ({len(rows)} rows: "
          f"{len(SITES)} sites x {len(TRANSPORTED_SPECIES)} species x 2 k-values x "
          f"~{I.N_YEARS - 1} years)")


def _filtering_bar_figures(summary):
    yr_range = f"{I.YEAR_START + 1}–{I.YEAR_END}"
    fig, axes = plt.subplots(1, len(SITES), figsize=(11, 6.5), sharey=True)
    fig.suptitle(f"Estuarine filtering — annual load in vs. out ({yr_range} mean)",
                 x=0.05, ha="left", fontsize=13, weight="bold")
    for ax, site in zip(axes, SITES):
        srows = [r for r in summary if r["site"] == site]
        y = np.arange(len(srows))
        vin = [r["input_load"] for r in srows]
        vout = [r["output_load"] for r in srows]
        ax.barh(y + 0.2, vin, height=0.35, color=S.MUTED, label="input (upstream)")
        ax.barh(y - 0.2, vout, height=0.35, color=S.RIVC[site], label="output (mouth)")
        ax.set_yticks(y)
        ax.set_yticklabels([f"{r['species']} [{r['unit']}]" for r in srows], fontsize=6.8)
        ax.set_title(S.LABEL[site], loc="left", color=S.RIVC[site], fontsize=10)
        ax.axvline(0, color=S.GRID, lw=0.8)
        ax.set_xscale("symlog", linthresh=1.0)
        S.tidy(ax, grid="x")
        ax.legend(loc="lower right", fontsize=6.5)
    _save(fig, "fig3a_filtering_loads",
          f"Input = upstream boundary, output = mouth, both 2 cells in from the "
          f"boundary-forced cell, mean over all post-warmup years ({yr_range}) of "
          "the 19-yr full-forcing record. Mouth loads are corrected for the ~12h-"
          "tide daily-snapshot aliasing bias (mouth_bias_correction -- calibrated "
          "once against runs/definitive_full/'s true 6-min output, real 2005-2006). "
          "Symlog x-axis (linear near 0). See filtering_mass_balance.csv for every "
          "year and the k=0 vs k=2 check. " + CONST_CAVEAT_FULL + SALT_CAVEAT)

    fig, axes = plt.subplots(1, len(SITES), figsize=(11, 6.5), sharey=True)
    fig.suptitle(f"Estuarine filtering — % retained (−) / added (+), "
                 f"{yr_range} mean", x=0.05, ha="left", fontsize=13, weight="bold")
    for ax, site in zip(axes, SITES):
        srows = [r for r in summary if r["site"] == site]
        y = np.arange(len(srows))
        vin = [r["input_load"] for r in srows]
        vout = [r["output_load"] for r in srows]
        pct = [r["pct_change"] if r["pct_change"] is not None else np.nan for r in srows]
        degenerate = [abs(vi) < 1e-3 * abs(vo) for vi, vo in zip(vin, vout)]
        pct_clip = np.clip([0.0 if d else p for d, p in zip(degenerate, pct)], -300, 300)
        colors = [S.WARN if (p == p and p > 0) else S.ACCENT for p in pct]
        ax.barh(y, pct_clip, color=colors, height=0.6)
        for yi, (p, pc, d) in enumerate(zip(pct, pct_clip, degenerate)):
            if d:
                ax.text(0, yi, "  in≈0", fontsize=5.5, va="center", ha="left", color=S.INK2)
            elif p == p and abs(p) > 300:
                ax.text(pc, yi, f"  {p:+.0f}%" if pc > 0 else f"{p:+.0f}%  ", fontsize=5.5,
                        va="center", ha="left" if pc > 0 else "right", color=S.INK2)
        ax.set_yticks(y)
        ax.set_yticklabels([r["species"] for r in srows], fontsize=6.8)
        ax.set_title(S.LABEL[site], loc="left", color=S.RIVC[site], fontsize=10)
        ax.axvline(0, color=S.GRID, lw=0.8)
        ax.set_xlim(-320, 320)
        S.tidy(ax, grid="x")
    _save(fig, "fig3b_filtering_percent",
          f"Percent change from input (upstream) to output (mouth) load, mean over "
          f"{yr_range}, clipped at ±300% for display; 'in≈0' marks a near-"
          "zero input base (pct_change is not meaningful there -- see the CSV for "
          "the raw value). Orange = net source (added in the estuary), blue = net "
          "sink (removed). " + CONST_CAVEAT_FULL + SALT_CAVEAT)


# =====================================================================================
# Fig 4 -- role of the freshet in biogeochemistry
# =====================================================================================
FRESHET_RESPONSE_VARS = [
    ("NPP", "rate", "C", "net primary production"),
    ("REMIN", "rate", "C", "remineralization (aer_deg+denit)"),
    ("FCO2", "rate", "C", "air–sea CO₂ flux"),
    ("DIC", "conc", None, "DIC"),
    ("TOC", "conc", None, "DOC (as TOC)"),
    ("NO3", "conc", None, "NO₃"),
]


def _annual_series(site, name, kind):
    if name == "REMIN":
        ya, va = I.area_normalized_annual_series(site, "aer_deg")
        yd, vd = I.area_normalized_annual_series(site, "denit")
        if ya.size == 0 or yd.size == 0 or not np.array_equal(ya, yd):
            return np.array([]), np.array([])
        return ya, va + vd
    if kind == "rate":
        return I.area_normalized_annual_series(site, name)
    return I.annual_openwater_mean_series(site, name)


def fig4a_freshet_trends():
    yr_range = f"{I.YEAR_START + 1}–{I.YEAR_END}"
    fig = plt.figure(figsize=(11, 8.5))
    fig.suptitle("Does freshet duration/intensity shape estuarine biogeochemistry?",
                 x=0.05, ha="left", fontsize=13, weight="bold")
    fig.text(0.05, 0.945, "Freshet duration (days > BREAKUP_Q_FACTOR× that year's "
             f"mean Q) vs. annual budgets/means of each response variable, {yr_range}, "
             "one point per site-year. Rate variables (NPP/REMIN/FCO2) are area-"
             "normalized annual budgets (gC m⁻² yr⁻¹); "
             "concentrations are annual open-water spatial means.",
             color=S.INK2, fontsize=7.6)
    gs = fig.add_gridspec(2, 3, left=0.07, right=0.97, top=0.87, bottom=0.15,
                          hspace=0.5, wspace=0.32)
    for k, (name, kind, element, title) in enumerate(FRESHET_RESPONSE_VARS):
        ax = fig.add_subplot(gs[k // 3, k % 3])
        for site in SITES:
            fm = I.freshet_metrics(site)
            if not fm:
                continue
            fyear_to_i = {y: i for i, y in enumerate(fm["year"])}
            yrs, vals = _annual_series(site, name, kind)
            if kind == "rate" and element:
                vals = vals * I.MOLAR_MASS[element] / 1000.0
            x, y2 = [], []
            for y, v in zip(yrs, vals):
                if y in fyear_to_i:
                    x.append(fm["duration"][fyear_to_i[y]])
                    y2.append(v)
            if not x:
                continue
            ax.scatter(x, y2, color=S.RIVC[site], s=16, alpha=0.8)
            _, _, r, p, n = I._linregress_np(x, y2)
            if n >= 3:
                ax.text(0.02, 0.03 + 0.09 * SITES.index(site),
                        f"{S.LABEL[site][:4]}: r={r:.2f} (p={p:.2f})",
                        transform=ax.transAxes, fontsize=5.8, color=S.RIVC[site])
        ax.axhline(0, color=S.GRID, lw=0.8) if kind == "rate" else None
        ax.set_title(title, loc="left", fontsize=8.5, color=S.INK)
        ylab = f"g{element} m⁻² yr⁻¹" if kind == "rate" else "mmol m⁻³"
        ax.set_ylabel(ylab, fontsize=7.2)
        ax.set_xlabel(f"freshet duration [days > {I.BREAKUP_Q_FACTOR:.0f}× mean Q]",
                       fontsize=7)
        S.tidy(ax)
    _save(fig, "fig4a_freshet_biogeochemistry_trends",
          f"FORTE Hypothesis 1: spring-freshet timing/intensity modulate carbon/"
          f"nutrient fluxes and air-sea CO2 flux. Uses runs/interannual_full/ "
          f"({I.N_YEARS - 1} post-warmup years). " + CONST_CAVEAT_FULL)


def freshet_window_means(site, name, kind):
    if name == "REMIN":
        t, Fa = I.load(site, "aer_deg")
        _, Fd = I.load(site, "denit")
        F = None if (Fa is None or Fd is None) else Fa + Fd
    else:
        t, F = I.load(site, name)
    if F is None:
        return np.array([]), np.array([]), np.array([])

    if kind == "rate":
        _, Dep = I.load(site, "depth")
        _, W = I.load(site, "width")
        w0 = W[0]
        areal = F * Dep
        daily = np.sum(areal * w0, axis=1) / np.sum(w0)
    else:
        _, ICE = I.load(site, "ice_frac")
        if ICE is not None and ICE.shape == F.shape:
            daily = np.nanmean(np.where(ICE < 0.5, F, np.nan), axis=1)
        else:
            daily = np.nanmean(F, axis=1)

    q_yrs, q = I._forcing_daily(site, "discharge")
    n_years = int((t.max() - I.WARMUP_DAYS) // 365)
    years, fmean, bmean = [], [], []
    for i in range(n_years):
        lo, hi = I.WARMUP_DAYS + i * 365, I.WARMUP_DAYS + (i + 1) * 365
        mask = (t >= lo) & (t < hi)
        if not mask.any():
            continue
        d_year = daily[mask]
        yr = I.YEAR_START + 1 + i
        qmask = (q_yrs >= yr) & (q_yrs < yr + 1)
        qy = q[qmask]
        n = min(len(qy), len(d_year))
        if n < 2:
            continue
        qy, d_year = qy[:n], d_year[:n]
        mean_q = qy.mean()
        if mean_q <= 0:
            continue
        fmask = qy > I.BREAKUP_Q_FACTOR * mean_q
        finite = np.isfinite(d_year)
        if not (fmask & finite).any() or not ((~fmask) & finite).any():
            continue
        years.append(yr)
        fmean.append(float(np.nanmean(d_year[fmask & finite])))
        bmean.append(float(np.nanmean(d_year[(~fmask) & finite])))
    return np.array(years), np.array(fmean), np.array(bmean)


def fig4b_freshet_window_comparison():
    yr_range = f"{I.YEAR_START + 1}–{I.YEAR_END}"
    fig = plt.figure(figsize=(11, 8.5))
    fig.suptitle(f"Freshet-window vs. baseflow-window biogeochemistry, {yr_range}",
                 x=0.05, ha="left", fontsize=13, weight="bold")
    fig.text(0.05, 0.945, "Per calendar year: mean value on freshet days (Q > "
             f"{I.BREAKUP_Q_FACTOR:.0f}× that year's mean) minus the mean on "
             "all other days. Positive = elevated during the freshet.",
             color=S.INK2, fontsize=8.0)
    gs = fig.add_gridspec(2, 3, left=0.07, right=0.97, top=0.87, bottom=0.15,
                          hspace=0.5, wspace=0.32)
    for k, (name, kind, element, title) in enumerate(FRESHET_RESPONSE_VARS):
        ax = fig.add_subplot(gs[k // 3, k % 3])
        for site in SITES:
            yrs, fm, bm = freshet_window_means(site, name, kind)
            if yrs.size == 0:
                continue
            diff = fm - bm
            if kind == "rate" and element:
                diff = diff * I.MOLAR_MASS[element] / 1000.0 * 86400.0
            ax.plot(yrs, diff, "o-", ms=2.5, lw=0.9, color=S.RIVC[site],
                    label=S.LABEL[site])
        ax.axhline(0, color=S.GRID, lw=0.8)
        ax.set_title(title, loc="left", fontsize=8.5, color=S.INK)
        ylab = f"g{element} m⁻² d⁻¹" if kind == "rate" else "mmol m⁻³"
        ax.set_ylabel(f"freshet − baseflow\n[{ylab}]", fontsize=7.0)
        ax.set_xlabel("year", fontsize=7.2)
        S.tidy(ax)
        if k == 0:
            ax.legend(loc="best", fontsize=6.5)
    _save(fig, "fig4b_freshet_window_comparison",
          "New calculation (not a trapezoidal annual budget -- freshet days are "
          "non-contiguous, so a per-day-subset mean is used instead). Uses runs/"
          "interannual_full/. " + CONST_CAVEAT_FULL)


# =====================================================================================
def main():
    OUT.mkdir(parents=True, exist_ok=True)
    print("Fig 1 -- forcing overview")
    fig1a_forcing_decadal()
    fig1b_forcing_seasonal()

    print("Fig 2 -- Hovmoller maps")
    make_hovmoller_figures()

    print("Fig 3 -- estuarine filtering mass balance")
    validate_discharge_reconstruction()
    rows = filtering_mass_balance()
    write_filtering_csv(rows)
    _filtering_bar_figures(summarize_by_site_species(rows, k=2))

    print("Fig 4 -- role of the freshet")
    fig4a_freshet_trends()
    fig4b_freshet_window_comparison()

    n = len(list(OUT.glob("*.pdf"))) + len(list(OUT.glob("*.csv")))
    print(f"done: {n} files in {OUT.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
