"""
Standalone, publication-ready NS-RAD figures -- one self-contained PDF per figure (not
a bundled multi-page report like ns_rad_diagnostics.pdf/ns_rad_interannual.pdf), meant
to drop directly into a manuscript or supplement. Written to
docs/publication_figures/.

Reads ONLY existing output -- runs/interannual/<site>/output.nc (1980-2023 daily-
forced record, 43 post-warmup years) and forcing/*.csv are the basis for EVERY
figure. runs/definitive/<site>/output.nc (2022 climatology, native 6-min cadence) is
used in exactly one place: as a CALIBRATION reference for Fig 3's mouth tidal-
aliasing correction (see below) -- it is not a data source for any figure directly.
No new model runs.

Reuses tools/make_diagnostics_pdf.py (as D: the definitive-run loader D.load, used
only for the Fig-3 calibration) and tools/make_interannual_pdf.py (as I: the
interannual loader I.load, I._forcing_daily, I.freshet_metrics,
I.area_normalized_annual_series, I.annual_openwater_mean_series, I._linregress_np,
I._site_delxi) rather than reimplementing any of that.

Four figure sets:
  1. Forcing overview       -- fig1a (decadal discharge+TOC, the full 44-yr record),
                               fig1b (discharge+DOC sliced to one representative
                               year from that SAME interannual record; water temp/
                               wind are the shared single-year climatology, which
                               has no interannual equivalent in this model)
  2. Hovmoller maps         -- fig2_hovmoller_<var>, var in {FCO2, TOC, DIC,
                               REMIN=aer_deg+denit, NPP}, distance x time,
                               1981-2023, daily (the interannual record's native
                               cadence -- no sub-daily/tidal-resolved version)
  3. Estuarine filtering    -- fig3a/b (loads in/out, % retained-added) + a CSV, for
     mass balance             every env=1 transported species (15; pH is env=0/
                               diagnosed, excluded), MEANED over 1981-2023 from
                               runs/interannual/. The upstream (input) side has no
                               tidal signal and needs no correction (verified in
                               validate_discharge_reconstruction). The mouth
                               (output) side is a daily SNAPSHOT of a quantity that
                               genuinely oscillates on the ~12h tidal cycle, which
                               would otherwise bias the estimate by a real,
                               systematic (NOT random/cancelling) 1-13% depending on
                               the river -- measured, not assumed (see
                               mouth_bias_correction). The correction factor is
                               calibrated ONCE per (site, species) against
                               runs/definitive's true 6-min output (both run
                               families share the identical deterministic tidal
                               forcing and t=0 clock, so this is an exact match of
                               the aliasing mechanism, not an approximation) and
                               applied to every interannual year -- the one
                               remaining assumption is that the correction ratio
                               itself is insensitive to that year's actual discharge
                               magnitude, which is calibrated against only the
                               2022-forcing year and not independently re-verified
                               per year.
  4. Freshet's role in      -- fig4a (freshet duration/intensity vs. annual
     biogeochemistry           NPP/remin/FCO2/DIC/TOC/NO3 budgets, trend-annotated,
                               1981-2023), fig4b (freshet-window vs. baseflow-window
                               mean, per year)

IMPORTANT INTERPRETIVE CAVEAT, repeated on the relevant figures: every solute
boundary except T and TOC (both genuinely time-varying in the interannual record) is
held CONSTANT. Any seasonal or decadal structure seen elsewhere in a Hovmoller/mass-
balance figure for the other solutes (DIC, NO3, NH4, PO4, dSi, O2, SPM, RDOC, CH4,
N2O, ALK) is produced by the MODEL'S OWN internal dynamics (ice gating of reactions/
gas-exchange, temperature-dependent rates, discharge-driven residence time, tidal
mixing) acting on a constant river source, not by a real observed time-varying river
concentration.

Usage:  python tools/make_publication_figures.py
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
import make_diagnostics_pdf as D   # noqa: E402  (definitive-run loader, GEOM, forcing)
import make_interannual_pdf as I   # noqa: E402  (interannual loader, freshet, trends)

ROOT = Path(__file__).resolve().parent.parent
OUT = ROOT / "docs" / "publication_figures"
SITES = ["colville", "kuparuk", "sagavanirktok"]   # canning excluded: EL=0, non-runnable

S.apply()   # idempotent; D/I already called this + install_autoscale(1.2) at import


def _save(fig, name, caption=None):
    OUT.mkdir(parents=True, exist_ok=True)
    if caption:
        # Figures built with bare plt.subplots() (no explicit gridspec bottom=) leave
        # only a thin default margin, which the caption would otherwise collide with
        # (axis tick labels) or overflow past (figure edge). subplots_adjust is a
        # no-op for gridspec-positioned axes (they bypass subplotpars), so this only
        # affects the plain-subplots figures that actually need it.
        fig.subplots_adjust(bottom=max(fig.subplotpars.bottom, 0.18))
        # y=0.032, not 0.006: S.brand()'s footer sits at y=0.006, va="bottom" -- a
        # multi-line caption anchored at the SAME y would have its last line collide
        # with the brand text (both anchor their bottom edge to the same baseline).
        fig.text(0.05, 0.032, caption, ha="left", va="bottom", fontsize=6.0,
                  color=S.MUTED, style="italic", wrap=True)
    S.brand(fig)
    path = OUT / f"{name}.pdf"
    fig.savefig(path)
    plt.close(fig)
    print(f"  wrote {path.relative_to(ROOT)} ({path.stat().st_size / 1024:.0f} kB)")


CONST_CAVEAT_SHORT = ("Upstream boundaries (except T, and TOC in the 44-yr record) "
                       "are CONSTANT -- structure shown here is the model's own "
                       "internal dynamics acting on a constant source, not a real "
                       "time-varying river signal. See CLAUDE.md -> \"Publication "
                       "figure suite\".")
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
    fig.suptitle("River forcing, 1980–2023", x=0.05, ha="left", fontsize=13,
                 weight="bold")
    fig.text(0.05, 0.94, "Discharge and upstream DOC (as TOC) -- the two forcings "
             "that are genuinely multi-year in this model (see 1b for the rest).",
             color=S.INK2, fontsize=8.0)
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
    axes[1].set_xlim(1980, 2024)
    for ax in axes:
        S.tidy(ax)
    _save(fig, "fig1a_forcing_decadal",
          "Source: tools/build_interannual_forcings.py, from Michael's daily PWBM "
          "product (1980-2023). No model run needed.")


def _forc_shared(fname):
    """A SHARED (not per-river) single-year forcing CSV -- the same file both run
    families read, since met forcing/tides are not interannual (see CLAUDE.md ->
    'Interannual forcing')."""
    p = ROOT / "forcing" / fname
    return np.genfromtxt(open(p, encoding="utf-8-sig"), delimiter=",")


def _climatology(site, kind):
    """Day-of-year climatology across the FULL 1980-2023 interannual record (all 44
    years, not just one) for 'discharge' or 'toc'. Returns (doy[0..364], mean, p10,
    p90). The record is exactly 44*365 daily values with no leap days (see
    tools/build_interannual_forcings.py), so `day-of-year = index % 365` recovers the
    calendar day cleanly with no reindexing needed."""
    _, v = I._forcing_daily(site, kind)
    doy_of_row = np.arange(v.size) % 365
    mean = np.empty(365); p10 = np.empty(365); p90 = np.empty(365)
    for d in range(365):
        vals = v[doy_of_row == d]
        mean[d] = np.mean(vals)
        p10[d] = np.percentile(vals, 10)
        p90[d] = np.percentile(vals, 90)
    return np.arange(365), mean, p10, p90


def fig1b_forcing_seasonal():
    fig, axes = plt.subplots(4, 1, figsize=(9, 9.5), sharex=True)
    fig.subplots_adjust(top=0.87)   # clear the 2-line subtitle above axQ's legend
    fig.suptitle("River forcing, seasonal climatology (1980–2023)",
                 x=0.05, ha="left", fontsize=13, weight="bold")
    # Explicit \n breaks, not wrap=True: S.install_autoscale() rescales font size at
    # savefig time (after this call), so matplotlib's own wrap-width calculation (done
    # at call time, against the PRE-scale font) mismatches the actually-rendered text
    # and can shove wrapped lines up into the title -- seen directly when this was
    # tried.
    fig.text(0.05, 0.955, "Discharge/DOC: day-of-year mean (line) + 10-90th percentile "
             "(shaded), all 44 years of the interannual record.\n"
             "Water temp/wind: the shared single-year climatology both run families "
             "use -- met/tide forcing has no interannual equivalent.",
             color=S.INK2, fontsize=7.6, va="top")
    axQ, axT, axW, axC = axes

    DISCHARGE_FLOOR = 1e-2   # m3/s; purely a log-scale display floor -- some winter
                              # calendar days are frozen solid (Q=0.0) in EVERY one
                              # of the 44 years, and log(0) cannot be plotted; the
                              # underlying data/CSV is untouched
    for site in SITES:
        doy, mean, p10, p90 = _climatology(site, "discharge")
        axQ.fill_between(doy, np.maximum(p10, DISCHARGE_FLOOR),
                          np.maximum(p90, DISCHARGE_FLOOR), color=S.RIVC[site], alpha=0.2,
                          linewidth=0)
        axQ.plot(doy, np.maximum(mean, DISCHARGE_FLOOR), color=S.RIVC[site], lw=1.3,
                 label=S.LABEL[site])
    axQ.set_yscale("log")
    axQ.set_ylabel("discharge [m³ s⁻¹]")
    axQ.legend(loc="upper left", ncols=3, fontsize=7.5)

    axT.plot(_forc_shared("river_watertemp_2022_degC.csv"), color="#2a78d6", lw=1.3,
              label="river (upstream bnd)")
    axT.plot(_forc_shared("watertemp.csv"), color=S.WARN, lw=1.1, label="sea (marine bnd)")
    axT.axhline(0, color=S.GRID, lw=0.8)
    axT.set_ylabel("water temperature [°C]")
    axT.legend(loc="upper left", fontsize=7.5)

    axW.plot(_forc_shared("windspeed.csv"), color=S.MUTED, lw=1.0)
    axW.set_ylabel("wind speed [m s⁻¹]  (shared, all rivers)")

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
          "Discharge/DOC: day-of-year mean + 10-90th percentile across all 44 years "
          "of forcing/*_interannual_1980-2023_*.csv (tools/build_interannual_"
          "forcings.py). Water temp/wind: forcing/river_watertemp_2022_degC.csv, "
          "watertemp.csv, windspeed.csv (shared regional record, not per-river, not "
          "interannual). No model run needed.")


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
    """`name` is an output.nc field, or 'REMIN' = aer_deg + denit (both carbon-
    currency in this model's own c_TOC update -- see CLAUDE.md)."""
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
    mask = t >= 365.0   # post-warmup, 1981-2023 -- consistent with every other
                          # interannual-record page in make_interannual_pdf.py
    tt, Fm = t[mask], F[mask]
    x = 1980.0 + tt / 365.0
    step = max(1, len(x) // 3000)   # ~daily-scale render cost cap over 43 years
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
    fig.subplots_adjust(bottom=0.18)   # room for the caption below the x-axis labels
    fig.suptitle(f"{label} — distance × time (1981–2023, daily)", x=0.05, ha="left",
                 fontsize=12.5, weight="bold")
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
    fig.text(0.05, 0.032, CONST_CAVEAT_SHORT, ha="left", va="bottom", fontsize=5.8,
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
# env=1 transported species (see code/variables.py `names` + code/init_module.py's
# initialize_substance, which sets env=1 for everything except pH). pH is env=0
# (diagnosed from S/DIC/ALK), so it is excluded -- there is no boundary flux to speak
# of for a diagnosed, non-transported field.
TRANSPORTED_SPECIES = ["DIA", "dSi", "NO3", "NH4", "PO4", "O2", "TOC", "RDOC", "CH4",
                        "N2O", "S", "SPM", "DIC", "ALK", "T"]

# (molar_mass_g_per_mol_or_None, unit_kind). "mmol" species have a defensible molar
# mass and are converted to tonnes/yr; ALK has no defensible mass conversion (it is a
# charge-equivalent, not literally carbon) and is reported in kmol/yr instead; SPM is
# already g/L (not mmol/m3); S (salinity, PSU) is converted via a salt-mass proxy
# (rho=1000 kg/m3 approx, fine at these near-fresh salinities); T is NOT a solute --
# it is reported as a net advective heat-TRANSPORT proxy (degC*hm3/yr), not a mass,
# and is heavily modulated by the model's own atmospheric heat exchange along the
# channel (heat_module.py), unlike the biogeochemical species. CH4 and N2O follow the
# same currency convention make_interannual_pdf.py's AREA_VARS already uses (CH4 as
# carbon, N2O as the mass of N2O itself, not N-equivalent).
SPECIES_SPEC = {
    "DIA": (12.011, "mmol"), "dSi": (28.085, "mmol"), "NO3": (14.007, "mmol"),
    "NH4": (14.007, "mmol"), "PO4": (30.974, "mmol"), "O2": (31.998, "mmol"),
    "TOC": (12.011, "mmol"), "RDOC": (12.011, "mmol"), "CH4": (12.011, "mmol"),
    "N2O": (44.013, "mmol"), "DIC": (12.011, "mmol"), "ALK": (None, "mmol_eq"),
    "SPM": (None, "gL"), "S": (None, "psu"), "T": (None, "degC"),
}


def _boundary_series_definitive(site, name, end, k=2):
    """(t_seconds, Q_seaward[m3/s], C[native units]) from runs/definitive/, year 2
    (day 365-730, post-warmup), native 6-min cadence, `k` grid cells in from the
    `end` boundary ('mouth' or 'upstream') -- avoids sampling exactly at the
    openbound-forced cell. Q is RECONSTRUCTED as velocity*depth*width
    (DEPTH[i]=D[i]/B[i] exactly, per code/init_module.py, so this equals the model's
    own cross-section-integrated flow) and SIGN-FLIPPED so positive = seaward:
    main.py negates the discharge forcing before using it to set the upstream
    velocity BC (`Qr = -Qr*Q_FRACTION`; `U[M] = Qr/D[M]`), so the model's raw
    velocity sign is landward-positive; flipping it here makes both boundaries read
    as "positive = net transport toward the ocean", matching the discharge forcing's
    own sign convention (see validate_discharge_reconstruction). Used ONLY to
    calibrate the mouth tidal-aliasing correction below -- the reported figures are
    sourced from runs/interannual/ (_boundary_series_interannual)."""
    t, Vel = D.load(site, "velocity")
    _, Dep = D.load(site, "depth")
    _, W = D.load(site, "width")
    _, C = D.load(site, name)
    if any(x is None for x in (Vel, Dep, W, C)):
        return None
    M = Vel.shape[1]
    col = k if end == "mouth" else M - 1 - k
    mask = (t >= 365.0) & (t < 730.0)
    tt = t[mask] * 86400.0
    Q = -(Vel[mask, col] * Dep[mask, col] * W[mask, col])
    Cc = C[mask, col]
    finite = np.isfinite(Q) & np.isfinite(Cc)
    tt, Q, Cc = tt[finite], Q[finite], Cc[finite]
    if len(tt) < 2:
        return None
    return tt, Q, Cc


def _boundary_series_interannual(site, name, end, k=2):
    """(t_days, Q_seaward[m3/s], C[native units]) from runs/interannual/, ALL post-
    warmup years (1981-2023) at once, native daily cadence. Same reconstruction and
    sign-flip as the definitive-run version above."""
    t, Vel = I.load(site, "velocity")
    _, Dep = I.load(site, "depth")
    _, W = I.load(site, "width")
    _, C = I.load(site, name)
    if any(x is None for x in (Vel, Dep, W, C)):
        return None
    M = Vel.shape[1]
    col = k if end == "mouth" else M - 1 - k
    mask = t >= 365.0
    tt = t[mask]
    Q = -(Vel[mask, col] * Dep[mask, col] * W[mask, col])
    Cc = C[mask, col]
    finite = np.isfinite(Q) & np.isfinite(Cc)
    tt, Q, Cc = tt[finite], Q[finite], Cc[finite]
    if len(tt) < 2:
        return None
    return tt, Q, Cc


# 86400s (interannual daily save) / 360s (definitive native cadence, CGEM_TS=12 x
# DELTI=30s) = 240 exactly. Both run families share the identical deterministic tidal
# forcing and t=0 clock, so subsampling the definitive run's TRUE 6-min output at this
# exact stride reproduces precisely the tidal phase the interannual run's daily save
# actually lands on -- letting the aliasing bias be MEASURED, not just caveated about.
MOUTH_ALIAS_STRIDE = 240
_CORR_CACHE = {}


def mouth_bias_correction(site, name, k=2):
    """Multiplicative correction for the tidal-aliasing bias a daily snapshot
    introduces at the MOUTH only (the upstream boundary has no tidal signal --
    verified in validate_discharge_reconstruction, so it needs no correction).
    correction = true_6min_annual_mean / daily_sampled_annual_mean, both computed
    from the SAME definitive-run year-2 series. Multiply this onto the interannual
    run's (aliased) daily-sampled mouth flux. Memoised (same value reused across all
    43 interannual years). Returns None if the definitive-run comparison is
    unavailable, or the daily-sampled mean is ~0 (correction undefined)."""
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
    """(calendar_years, per-year time-integral of Q_seaward*C in the species' native
    per-second basis), from runs/interannual/, one value per full post-warmup year
    1981-2023. `apply_correction`: multiply the mouth-side flux by
    mouth_bias_correction before integrating (only meaningful/used for end='mouth')."""
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
    n_years = int((tt.max() - 365.0) // 365)
    years, vals = [], []
    for i in range(n_years):
        lo, hi = 365 + i * 365, 365 + (i + 1) * 365
        m = (tt >= lo) & (tt < hi)
        if m.sum() < 2:
            continue
        years.append(1980 + 1 + i)
        vals.append(float(np.trapezoid(flux[m], tt[m] * 86400.0)))
    return np.array(years), np.array(vals)


def _report(name, raw):
    """Convert a raw flux-integral to a human-scale unit. Returns (value, unit)."""
    if raw is None:
        return None, ""
    mm, kind = SPECIES_SPEC[name]
    if kind == "gL":
        return raw / 1e6, "t/yr"            # g -> t
    if kind == "psu":
        return raw / 1e3, "t/yr"            # kg -> t (salt-mass proxy)
    if kind == "degC":
        return raw / 1e6, "°C·hm³/yr"   # m3 -> hm3; NOT a mass
    if kind == "mmol_eq":
        return raw / 1e6, "kmol/yr"         # mmol -> kmol
    return raw * mm / 1e9, "t/yr"           # mmol*g/mol = mg; mg -> t is /1e9


def validate_discharge_reconstruction():
    """Sanity check the flux-reconstruction method on the interannual record:
    reconstructed mean Q at the upstream boundary (k=2), averaged over all
    post-warmup years, should closely match the interannual discharge forcing's own
    all-years mean (the file the model's boundary condition literally reads) --
    confirms there is no tidal contamination at this end (unlike the mouth, which is
    why mouth_bias_correction below exists)."""
    print("Discharge reconstruction sanity check (reconstructed vs. forcing, "
          "1981-2023 mean):")
    for site in SITES:
        res = _boundary_series_interannual(site, "S", "upstream", k=2)
        if res is None:
            print(f"  {site:14s} -- no data"); continue
        _, Q, _ = res
        recon = float(np.mean(Q))
        yrs, q = I._forcing_daily(site, "discharge")
        forc = float(np.mean(q[yrs >= 1981]))
        pct = 100.0 * (recon - forc) / forc if forc else float("nan")
        print(f"  {site:14s} reconstructed={recon:9.2f}  forcing={forc:9.2f}  "
              f"diff={pct:+.1f}%")

    print("Mouth tidal-aliasing correction (definitive-run calibration, "
          "true 6-min mean / daily-sampled mean):")
    for site in SITES:
        line = [f"{site:14s}"]
        for name in ("TOC", "DIC", "O2", "NO3"):
            corr = mouth_bias_correction(site, name, k=2)
            line.append(f"{name}={corr:.3f}" if corr is not None else f"{name}=n/a")
        print("  " + "  ".join(line))


def filtering_mass_balance():
    """One row per (site, species, k) -- input/output loads, net effect, %change."""
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
    """Collapse the per-year rows (k=2 only) to one row per (site, species): mean
    across 1981-2023 of input/output load, and net/pct derived from THOSE means (not
    a mean-of-ratios, which would be sensitive to whatever year happens to have a
    tiny input load)."""
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
          "~43 years)")


def _filtering_bar_figures(summary):
    fig, axes = plt.subplots(1, len(SITES), figsize=(11, 6.5), sharey=True)
    fig.suptitle("Estuarine filtering — annual load in vs. out (1981–2023 mean)",
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
          "Input = upstream boundary, output = mouth, both 2 cells in from the "
          "boundary-forced cell, mean over all post-warmup years (1981-2023) of the "
          "44-yr interannual record. Mouth loads are corrected for the ~12h-tide "
          "daily-snapshot aliasing bias (mouth_bias_correction -- calibrated once "
          "against the definitive run's true 6-min output; see CLAUDE.md -> "
          "\"Publication figure suite\"). Symlog x-axis (linear near 0). See "
          "filtering_mass_balance.csv for every year and the k=0 vs k=2 check. "
          + CONST_CAVEAT_SHORT + SALT_CAVEAT)

    fig, axes = plt.subplots(1, len(SITES), figsize=(11, 6.5), sharey=True)
    fig.suptitle("Estuarine filtering — % retained (−) / added (+), 1981–2023 mean",
                 x=0.05, ha="left", fontsize=13, weight="bold")
    for ax, site in zip(axes, SITES):
        srows = [r for r in summary if r["site"] == site]
        y = np.arange(len(srows))
        vin = [r["input_load"] for r in srows]
        vout = [r["output_load"] for r in srows]
        pct = [r["pct_change"] if r["pct_change"] is not None else np.nan for r in srows]
        # A near-zero input load (e.g. S's upstream boundary is exactly 0 -- pure
        # freshwater) makes pct_change numerically explode (division by floating-
        # point noise, e.g. ~1e-6). The absolute noise floor differs by river/grid
        # (Sagavanirktok's is ~1e-5, Colville's ~1e-7), so flag degenerate RELATIVE
        # to that species' own output magnitude, not a fixed absolute threshold --
        # input is negligible if it's not even 0.1% of what comes out the mouth.
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
          "Percent change from input (upstream) to output (mouth) load, mean over "
          "1981-2023, clipped at ±300% for display; 'in≈0' marks a near-zero input "
          "base (pct_change is not meaningful there -- see the CSV for the raw "
          "value). Orange = net source (added in the estuary), blue = net sink "
          "(removed). " + CONST_CAVEAT_SHORT + SALT_CAVEAT)


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
    fig = plt.figure(figsize=(11, 8.5))
    fig.suptitle("Does freshet duration/intensity shape estuarine biogeochemistry?",
                 x=0.05, ha="left", fontsize=13, weight="bold")
    fig.text(0.05, 0.945, "Freshet duration (days > BREAKUP_Q_FACTOR× that year's "
             "mean Q -- I.freshet_metrics' own definition) vs. annual budgets/means "
             "of each response variable, 1981–2023, one point per site-year. Rate "
             "variables (NPP/REMIN/FCO2) are area-normalized annual budgets "
             "(gC m⁻² yr⁻¹); concentrations are annual open-water spatial means.",
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
          "FORTE Hypothesis 1: spring-freshet timing/intensity modulate carbon/"
          "nutrient fluxes and air-sea CO2 flux. Uses runs/interannual/ (43 "
          "post-warmup years). " + CONST_CAVEAT_SHORT)


def freshet_window_means(site, name, kind):
    """Per calendar year: mean of `name` during freshet days (Q > BREAKUP_Q_FACTOR x
    that year's own mean daily discharge -- I.freshet_metrics' definition) vs. the
    rest of the year. Rate variables are depth-integrated + width-weighted per day
    (matching I.area_normalized_annual_series' spatial method) but NOT trapezoidal-
    time-integrated -- freshet days are not contiguous, so a per-day-subset MEAN
    avoids any integration-over-gaps pitfall. Concentration variables use the same
    ice-aware open-water spatial mean as I.annual_openwater_mean_series. Returns
    (years, freshet_mean, baseflow_mean)."""
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
    n_years = int((t.max() - 365.0) // 365)
    years, fmean, bmean = [], [], []
    for i in range(n_years):
        lo, hi = 365 + i * 365, 365 + (i + 1) * 365
        mask = (t >= lo) & (t < hi)
        if not mask.any():
            continue
        d_year = daily[mask]
        yr = 1980 + 1 + i
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
    fig = plt.figure(figsize=(11, 8.5))
    fig.suptitle("Freshet-window vs. baseflow-window biogeochemistry, 1981–2023",
                 x=0.05, ha="left", fontsize=13, weight="bold")
    fig.text(0.05, 0.945, "Per calendar year: mean value on freshet days (Q > "
             f"{I.BREAKUP_Q_FACTOR:.0f}× that year's mean) minus the mean on all "
             "other days. Positive = elevated during the freshet.",
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
                diff = diff * I.MOLAR_MASS[element] / 1000.0 * 86400.0  # per s -> per day
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
          "interannual/. " + CONST_CAVEAT_SHORT)


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
