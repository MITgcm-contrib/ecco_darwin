"""
Diagnostic plots of modelled fields for all four North Slope rivers.

Reads the definitive run (runs/definitive/<site>/): transported T + heat budget + observed-blended T
boundary. Every panel overlays the four rivers so they can be compared directly.

Page 1  seasonal cycle (spatial-mean, open water) of the core state variables
Page 2  longitudinal profiles — a mid-summer snapshot down each channel
Page 3  distance x time Hovmoller of temperature and salinity, per river
Page 4  process-rate seasonal cycles + a data-quality panel (known defects)

These are MODEL fields, not validation. For model-vs-observation see
ns_rad_validation.pdf.

Usage:  python tools/make_diagnostics_pdf.py
"""
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import nsrad_style as _S  # shared NS-RAD publication/presentation style
from matplotlib.backends.backend_pdf import PdfPages

ROOT = Path(__file__).resolve().parent.parent
RUN = ROOT / "runs" / "definitive"
OUT = ROOT / "docs" / "ns_rad_diagnostics.pdf"

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
    "legend.frameon": False, "legend.fontsize": 7.5, "grid.color": GRID,
    "grid.linewidth": 0.6, "lines.solid_capstyle": "round",
})
_S.apply()
_S.install_autoscale(1.2)  # embed fonts + presentation-size bump
MONTH0 = [0, 31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334]
MONTHL = list("JFMAMJJASOND")
DELXI_KM = 0.2


def tidy(ax, both=False):
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    ax.grid(True, axis="both" if both else "y", alpha=0.6)
    ax.set_axisbelow(True)


def month_axis(ax):
    ax.set_xticks(MONTH0)
    ax.set_xticklabels(MONTHL)
    ax.set_xlim(0, 364)


_CACHE = {}


def load(site, var):
    """Load a field (time, values[cells 1..M]), memoised. Prefers NetCDF (output.nc);
    falls back to the legacy .dat with pandas' fast C parser."""
    key = (site, var)
    if key in _CACHE:
        return _CACHE[key]
    res = _load_nc(RUN / site, var)
    if res is None:
        res = _load_dat(RUN / site, var)
    _CACHE[key] = res
    return res


def _load_nc(sitedir, var):
    ncp = sitedir / "output.nc"
    if not ncp.exists():
        return None
    from netCDF4 import Dataset
    d = Dataset(ncp)
    if var not in d.variables:
        return (None, None)
    t = np.asarray(d.variables["time"][:]) / 86400.0
    F = np.asarray(d.variables[var][:])   # (time, x=0..M)
    return (t, F[:, 1:])                  # cells 1..M, matching the .dat convention


def _load_dat(sitedir, var):
    p = sitedir / f"{var}.dat"
    if not p.exists():
        return (None, None)
    try:
        import pandas as pd
        d = pd.read_csv(p, sep="\t", header=None).to_numpy()
    except Exception:
        d = np.genfromtxt(p, delimiter="\t")
    if np.isnan(d[:, -1]).all():
        d = d[:, :-1]
    return (d[:, 0] / 86400.0, d[:, 1:])


def openwater_mean_series(site, var):
    """Spatial-mean over OPEN-WATER cells vs day-of-year (both model years).

    With the prognostic ice model, winter state is CONSERVED (non-zero), so the old
    'non-zero == active' heuristic no longer isolates open water. Mask by ice_frac
    where it exists (open = ice_frac < 0.5), falling back to the non-zero heuristic
    for legacy runs that zeroed winter state."""
    t, F = load(site, var)
    if F is None:
        return None, None
    _, ICE = load(site, "ice_frac")
    if ICE is not None and ICE.shape == F.shape:
        open_mask = ICE < 0.5
        Fm = np.where(open_mask, F, np.nan)
        series = np.nanmean(Fm, axis=1)                 # NaN where fully ice-covered
    else:
        active = np.abs(F).sum(axis=1) > 0
        series = np.where(active, np.nanmean(np.where(F != 0, F, np.nan), axis=1), np.nan)
    doy = np.mod(t, 365)
    o = np.argsort(doy)
    return doy[o], series[o]


FORC = ROOT / "forcing"
DISCH = {"colville": "colville_river_discharge_2022_m3sec.csv",
         "kuparuk": "kuparuk_river_discharge_2022_m3sec.csv",
         "sagavanirktok": "sagavanirktok_river_discharge_2022_m3sec.csv",
         "canning": "canning_river_discharge_2022_m3sec.csv"}


def _forc(fname):
    return np.genfromtxt(open(FORC / fname, encoding="utf-8-sig"), delimiter=",")


def _boundaries(site):
    import importlib, sys
    p = str(ROOT / "code")
    if p not in sys.path:
        sys.path.insert(0, p)
    m = importlib.import_module(f"sites.{site}")
    return getattr(m, "BOUNDARIES", None) or importlib.import_module("sites._baseline").BOUNDARIES


# --------------------------------------------------- page: boundary conditions & forcing
def page_boundaries(pdf):
    fig = plt.figure(figsize=(11, 8.5))
    fig.suptitle("Boundary conditions & forcing", x=0.05, ha="left",
                 fontsize=13, weight="bold")
    fig.text(0.05, 0.935, "External drivers of the model. Discharge and river temperature "
             "are per-river; meteorology and marine pCO₂ are shared regional records. "
             "Boundary concentrations at the table foot.", color=INK2, fontsize=8.3)
    gs = fig.add_gridspec(2, 3, left=0.07, right=0.97, top=0.90, bottom=0.42,
                          hspace=0.5, wspace=0.32)

    # (a) discharge per river
    ax = fig.add_subplot(gs[0, 0])
    for s in SITES:
        ax.plot(_forc(DISCH[s]), color=C[s], lw=1.1, label=LABEL[s])
    ax.set_yscale("log"); month_axis(ax); tidy(ax)
    ax.set_ylabel("discharge (m³ s⁻¹)", fontsize=7.3)
    ax.set_title("Discharge — USGS gauges (per river)", loc="left", fontsize=8)
    ax.legend(loc="lower center", ncol=2, fontsize=5.6, labelcolor=INK2)

    # (b) temperature: river forcing + sea boundary
    ax = fig.add_subplot(gs[0, 1])
    ax.plot(_forc("river_watertemp_2022_degC.csv"), color="#2a78d6", lw=1.3,
            label="river (upstream bnd)")
    ax.plot(_forc("watertemp.csv"), color=WARN, lw=1.1, label="sea (marine bnd)")
    ax.axhline(0, color=GRID, lw=0.8); month_axis(ax); tidy(ax)
    ax.set_ylabel("temperature (°C)", fontsize=7.3)
    ax.set_title("Water-temperature boundaries", loc="left", fontsize=8)
    ax.legend(loc="upper left", fontsize=5.8, labelcolor=INK2)

    # (c) wind + solar (shared)
    ax = fig.add_subplot(gs[0, 2])
    ax.plot(_forc("windspeed.csv"), color=MUTED, lw=1.0)
    ax.set_ylabel("wind (m s⁻¹)", fontsize=7.3, color=MUTED)
    month_axis(ax); tidy(ax)
    ax2 = ax.twinx()
    ax2.plot(_forc("solarradiation.csv"), color="#eb9b34", lw=1.0)
    ax2.set_ylabel("solar (W m⁻²)", fontsize=7.3, color="#eb9b34")
    ax2.spines["top"].set_visible(False)
    ax.set_title("Wind & solar — NDBC PRDA2 (shared)", loc="left", fontsize=8)

    # (d) pCO2 (shared marine)
    ax = fig.add_subplot(gs[1, 0])
    ax.plot(_forc("pCO2_Barrow_2022.csv"), color="#4a3aa7", lw=1.2)
    month_axis(ax); tidy(ax)
    ax.set_ylabel("pCO₂ (µatm)", fontsize=7.3)
    ax.set_title("Atmospheric pCO₂ — Barrow (shared)", loc="left", fontsize=8)

    # (e) humidity (shared)
    ax = fig.add_subplot(gs[1, 1])
    ax.plot(_forc("relhum_2022_frac.csv") * 100, color="#1baf7a", lw=1.1)
    month_axis(ax); tidy(ax)
    ax.set_ylabel("rel. humidity (%)", fontsize=7.3)
    ax.set_title("Humidity — Deadhorse Airport (shared)", loc="left", fontsize=8)

    # (f) air temperature (shared, drives heat budget + ice)
    ax = fig.add_subplot(gs[1, 2])
    ax.plot(_forc("airtemp_2022_degC.csv"), color=INK2, lw=1.0)
    ax.axhline(0, color=GRID, lw=0.8); month_axis(ax); tidy(ax)
    ax.set_ylabel("air temp (°C)", fontsize=7.3)
    ax.set_title("Air temperature — PRDA2 (shared)", loc="left", fontsize=8)

    # boundary-condition table (downstream marine / upstream river) per site
    axT = fig.add_axes([0.05, 0.05, 0.92, 0.32]); axT.axis("off")
    axT.set_title("Boundary concentrations — downstream (marine) | upstream (river), per site",
                  loc="left", fontsize=9, color=INK)
    species = ["S", "T", "DIC", "ALK", "pH", "NO3", "NH4", "PO4", "TOC", "O2"]
    units = {"S": "PSU", "T": "°C", "pH": "", "O2": "mmol/m³"}
    xc = [0.0] + list(np.linspace(0.14, 0.98, len(species)))
    axT.text(xc[0], 0.95, "Site / bnd", fontsize=6.8, weight="bold", color=INK2, va="top")
    for x, sp in zip(xc[1:], species):
        axT.text(x, 0.95, sp, fontsize=6.8, weight="bold", color=INK2, va="top", ha="center")
    row = 0.88
    for s in SITES:
        B = _boundaries(s)
        for j, (bnd, idx) in enumerate([("marine", 0), ("river", 1)]):
            y = row - (0 if j == 0 else 0.045)
            axT.text(xc[0], y, f"{LABEL[s][:9]} {bnd}", fontsize=6.2,
                     color=C[s], va="top")
            for x, sp in zip(xc[1:], species):
                val = B[sp][idx] if sp in B else float("nan")
                axT.text(x, y, f"{val:.4g}", fontsize=6.2, color=C[s], va="top", ha="center")
        row -= 0.11
    axT.text(0.0, row + 0.02, "Kuparuk river chemistry = Arctic LTER headwater "
             "(natural reach); Colville/Sag/Canning river chemistry still placeholder. "
             "S/T marine boundary shared; river S = 0.", fontsize=6.2, color=INK2, va="top")
    pdf.savefig(fig)
    plt.close(fig)


# ------------------------------------------------------------- page 1: seasonal cycle
def page_seasonal(pdf):
    fig = plt.figure(figsize=(11, 8.5))
    fig.suptitle("Modelled state variables — seasonal cycle (spatial mean, open water)",
                 x=0.05, ha="left", fontsize=13, weight="bold")
    fig.text(0.05, 0.935, "Definitive run: transported T + heat budget + prognostic ice. "
             "Series masked to OPEN-WATER cells (ice_frac < 0.5). All four rivers overlaid.",
             color=INK2, fontsize=8.3)
    panels = [("T", "temperature (°C)"), ("S", "salinity (PSU)"),
              ("DIC", "DIC (mmol m⁻³)"), ("pH", "pH"),
              ("O2", "O₂ (mmol m⁻³)"), ("NO3", "nitrate (mmol m⁻³)"),
              ("TOC", "TOC (mmol m⁻³)"), ("SPM", "SPM (mg L⁻¹)")]
    gs = fig.add_gridspec(4, 2, left=0.07, right=0.97, top=0.90, bottom=0.07,
                          hspace=0.55, wspace=0.18)
    for k, (var, ylab) in enumerate(panels):
        ax = fig.add_subplot(gs[k // 2, k % 2])
        for s in SITES:
            doy, ser = openwater_mean_series(s, var)
            if ser is not None:
                ax.plot(doy, ser, color=C[s], lw=1.2, label=LABEL[s])
        month_axis(ax)
        tidy(ax)
        ax.set_ylabel(ylab, fontsize=7.5)
        if k == 0:
            ax.legend(loc="upper left", ncol=2, fontsize=6.3, labelcolor=INK2)
        if var == "pH":
            ax.text(0.5, 0.06, "pH now physical (5–9.5); the inherited negative-H⁺ solver "
                    "defect is fixed (see p.6)", transform=ax.transAxes, ha="center",
                    fontsize=5.8, color=INK2)
    pdf.savefig(fig)
    plt.close(fig)


# ------------------------------------------------------------- page 2: profiles
def page_profiles(pdf):
    fig = plt.figure(figsize=(11, 8.5))
    fig.suptitle("Longitudinal profiles — mid-summer snapshot (≈ day 210)",
                 x=0.05, ha="left", fontsize=13, weight="bold")
    fig.text(0.05, 0.935, "Value vs distance from the mouth, down each channel. "
             "Mouth = 0 km (marine boundary), head ≈ 27 km (riverine boundary).",
             color=INK2, fontsize=8.3)
    panels = [("T", "temperature (°C)"), ("S", "salinity (PSU)"),
              ("DIC", "DIC (mmol m⁻³)"), ("O2", "O₂ (mmol m⁻³)"),
              ("NO3", "nitrate (mmol m⁻³)"), ("SPM", "SPM (mg L⁻¹)")]
    gs = fig.add_gridspec(3, 2, left=0.07, right=0.97, top=0.90, bottom=0.07,
                          hspace=0.42, wspace=0.18)
    for k, (var, ylab) in enumerate(panels):
        ax = fig.add_subplot(gs[k // 2, k % 2])
        for s in SITES:
            t, F = load(s, var)
            if F is None:
                continue
            doy = np.mod(t, 365)
            act = np.abs(F).sum(axis=1) > 0
            cand = np.where(act & (np.abs(doy - 210) < 8))[0]
            if len(cand) == 0:
                continue
            prof = np.nanmean(F[cand], axis=0)
            x = np.arange(F.shape[1]) * DELXI_KM
            ax.plot(x, prof, color=C[s], lw=1.3, label=LABEL[s])
        tidy(ax, both=True)
        ax.set_ylabel(ylab, fontsize=7.5)
        ax.set_xlabel("distance from mouth (km)", fontsize=7.5)
        ax.set_xlim(0, 27)
        if k == 0:
            ax.legend(loc="upper right", fontsize=6.3, labelcolor=INK2)
    pdf.savefig(fig)
    plt.close(fig)


# ------------------------------------------------------------- page 3: Hovmoller
def page_hovmoller(pdf):
    fig = plt.figure(figsize=(11, 8.5))
    fig.suptitle("Temperature & salinity — distance × time (Hovmöller), per river",
                 x=0.05, ha="left", fontsize=13, weight="bold")
    fig.text(0.05, 0.935, "Colour = field value; x = day of year, y = distance from "
             "mouth. State is conserved year-round; the winter T≈0 °C band is the ice "
             "season (water held at freezing), not missing data.",
             color=INK2, fontsize=8.3)
    gs = fig.add_gridspec(4, 2, left=0.08, right=0.9, top=0.90, bottom=0.06,
                          hspace=0.5, wspace=0.22)
    for row, s in enumerate(SITES):
        for col, (var, cmap, vlab) in enumerate(
                [("T", "inferno", "T (°C)"), ("S", "viridis", "S (PSU)")]):
            ax = fig.add_subplot(gs[row, col])
            t, F = load(s, var)
            if F is None:
                continue
            doy = np.mod(t, 365)
            yr1 = t <= 365
            d1 = doy[yr1]
            # No zero-masking: state is conserved year-round, so the winter field is a
            # REAL value (T held at 0 °C under ice), not missing. Masking sum==0 used to
            # white out that ice season and read as a data gap.
            Fm = F[yr1]
            o = np.argsort(d1)
            d1, Fm = d1[o], Fm[o]
            # downsample time to ~daily -- 35k pcolormesh columns render slowly and
            # look identical at this figure size
            step = max(1, len(d1) // 730)
            d1, Fm = d1[::step], Fm[::step]
            x = np.arange(F.shape[1]) * DELXI_KM
            im = ax.pcolormesh(d1, x, Fm.T, cmap=cmap, shading="auto")
            month_axis(ax)
            ax.set_ylabel(f"{LABEL[s]}\nkm", fontsize=7, color=C[s])
            if row == 0:
                ax.set_title(vlab, loc="left", fontsize=8.5, color=INK)
            cb = fig.colorbar(im, ax=ax, pad=0.02, fraction=0.05)
            cb.ax.tick_params(labelsize=6)
    pdf.savefig(fig)
    plt.close(fig)


# ------------------------------------------------------------- page: prognostic ice
def page_ice(pdf):
    fig = plt.figure(figsize=(11, 8.5))
    fig.suptitle("Prognostic river ice — seasonal cycle & distance × time",
                 x=0.05, ha="left", fontsize=13, weight="bold")
    fig.text(0.05, 0.935, "New ice_module: freeze-up, ~1–1.3 m bottom-fast winter cover, "
             "abrupt hydraulic breakup at the freshet, summer open water, refreeze.",
             color=INK2, fontsize=8.3)

    # top: channel-mean ice thickness + areal fraction, all rivers overlaid
    axT = fig.add_axes([0.07, 0.60, 0.42, 0.30])
    axF = fig.add_axes([0.56, 0.60, 0.42, 0.30])
    any_ice = False
    for s in SITES:
        t, H = load(s, "ice_thickness")
        _, FR = load(s, "ice_frac")
        if H is None:
            continue
        any_ice = True
        doy = np.mod(t, 365)
        o = np.argsort(doy)
        axT.plot(doy[o], np.nanmean(H, axis=1)[o], color=C[s], lw=1.2, label=LABEL[s])
        if FR is not None:
            axF.plot(doy[o], np.nanmean(FR, axis=1)[o], color=C[s], lw=1.2, label=LABEL[s])
    for ax, ylab in ((axT, "mean ice thickness (m)"), (axF, "mean ice fraction (0–1)")):
        month_axis(ax); tidy(ax); ax.set_ylabel(ylab, fontsize=7.5)
    # legend in the summer open-water gap (ice = 0, ~Jun–Sep) to clear the curves
    axT.legend(loc="center", bbox_to_anchor=(0.56, 0.5), ncol=2, fontsize=6.3,
               labelcolor=INK2, columnspacing=1.0, handletextpad=0.4)
    if not any_ice:
        axT.text(0.5, 0.5, "no ice fields in run (ICE_MODEL off?)", transform=axT.transAxes,
                 ha="center", color=WARN)

    # bottom: per-river Hovmoller of ice thickness
    gs = fig.add_gridspec(2, 2, left=0.08, right=0.9, top=0.50, bottom=0.06,
                          hspace=0.45, wspace=0.22)
    for k, s in enumerate(SITES):
        ax = fig.add_subplot(gs[k // 2, k % 2])
        t, H = load(s, "ice_thickness")
        if H is None:
            ax.axis("off"); continue
        doy = np.mod(t, 365)
        yr = t <= 365
        d1, H1 = doy[yr], H[yr]
        o = np.argsort(d1)
        step = max(1, len(d1) // 730)
        d1, H1 = d1[o][::step], H1[o][::step]
        x = np.arange(H.shape[1]) * DELXI_KM
        im = ax.pcolormesh(d1, x, H1.T, cmap="Blues", shading="auto", vmin=0)
        month_axis(ax)
        ax.set_ylabel(f"{LABEL[s]}\nkm", fontsize=7, color=C[s])
        cb = fig.colorbar(im, ax=ax, pad=0.02, fraction=0.05)
        cb.ax.tick_params(labelsize=6)
        cb.set_label("ice (m)", fontsize=6)
    pdf.savefig(fig)
    plt.close(fig)


# ------------------------------------------------------------- page 4: rates + QC
def page_rates_qc(pdf):
    fig = plt.figure(figsize=(11, 8.5))
    fig.suptitle("Process rates & data-quality diagnostics", x=0.05, ha="left",
                 fontsize=13, weight="bold")
    rates = [("NPP", "net primary prod. (mmol C m⁻³ s⁻¹)"),
             ("aer_deg", "aerobic degradation"), ("nit", "nitrification"),
             ("FCO2", "air–sea CO₂ flux (>0 outgassing)")]
    gs = fig.add_gridspec(3, 2, left=0.07, right=0.97, top=0.90, bottom=0.30,
                          hspace=0.5, wspace=0.18)
    for k, (var, ylab) in enumerate(rates):
        ax = fig.add_subplot(gs[k // 2, k % 2])
        for s in SITES:
            doy, ser = openwater_mean_series(s, var)
            if ser is not None:
                ax.plot(doy, ser, color=C[s], lw=1.1, label=LABEL[s])
        month_axis(ax)
        tidy(ax)
        ax.set_ylabel(ylab, fontsize=7)
        if k == 0:
            ax.legend(loc="upper left", ncol=2, fontsize=6.3, labelcolor=INK2)

    # QC / known-defect panel
    ax = fig.add_axes([0.07, 0.06, 0.9, 0.2])
    ax.axis("off")
    ax.set_title("Data-quality notes — read alongside these fields", loc="left",
                 color=INK, fontsize=9)
    notes = [
        "pH is now physical (5–9.5, all cells inside [2,12]). The inherited Follows-2006 solver used to return a "
        "huge/negative H⁺ at ice-out (when the gate leaves DIC/ALK tiny-positive), giving pH < 0; this is FIXED "
        "with a DIC/ALK floor and a physical-range bound. (The borate-alkalinity unit mismatch remains, but is "
        "negligible for these near-fresh rivers.)",
        "Salinity is ≈0 through the domain in open water — the saline zone is outside the configured reach "
        "(distance=1, weak tides). Do not read the S panels as estuarine mixing.",
        "Prognostic ice model (ice_module): under ice, state is CONSERVED and keeps reacting, so DIC builds up "
        "under the winter cover and vents at the freshet breakup (~day 150). Gas exchange (FCO2, O2) is scaled to "
        "zero under ice, and under-ice PAR is attenuated through the slab. Toggle off with CGEM_ICE=off.",
        "Kuparuk river chemistry uses Arctic LTER headwater values; the other three carry a placeholder — so "
        "absolute nutrient/DIC/TOC levels are not comparable across rivers.",
        "These are the raw model fields. Only temperature has same-year in-situ validation "
        "(see ns_rad_validation.pdf); the rest are unconstrained by data.",
    ]
    for i, n in enumerate(notes):
        ax.text(0.0, 0.95 - i * 0.2, "•", transform=ax.transAxes, color=MUTED, fontsize=9)
        ax.text(0.02, 0.95 - i * 0.2, n, transform=ax.transAxes, fontsize=7.2,
                color=INK2, va="top")
    pdf.savefig(fig)
    plt.close(fig)


def main():
    OUT.parent.mkdir(exist_ok=True)
    with PdfPages(OUT) as pdf:
        page_boundaries(pdf)
        page_seasonal(pdf)
        page_profiles(pdf)
        page_hovmoller(pdf)
        page_ice(pdf)
        page_rates_qc(pdf)
    print(f"wrote {OUT.relative_to(ROOT)} ({OUT.stat().st_size/1024:.0f} kB, 6 pages)")


if __name__ == "__main__":
    main()
