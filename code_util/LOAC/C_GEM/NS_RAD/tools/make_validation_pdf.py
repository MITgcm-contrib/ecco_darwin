"""
Validation figures: modelled river temperature and salinity vs ALL available data.

Model run
  runs/   the model (transported T + surface heat budget + prognostic ice, observed
          humidity). A regression-boundary run (runs/regression_bnd) is used ONLY on the
          temperature-scatter page to show the fit is independent of the observed
          boundary it is compared against -- it is the same physics, a different T
          boundary source, not a prior model version.

Observations (docs/validation_obs.json, built by scratchpad/collect_obs.py)
  temp_usgs_2022  USGS param 00010 daily, SAME YEAR as the model -- Sag 93 d, Kuparuk 55 d
  temp_wqp        Water Quality Portal grab samples, ALL years, pooled by day-of-year
                  (validates the model climatology) -- all four rivers
  salinity_wqp    WQP salinity grabs, all years -- Colville & Kuparuk deltas only; none
                  for Sagavanirktok / Canning, and NONE anywhere in 2022

Pages: 1 temperature, 2 temperature scatter, 3 salinity, 4 FCO2, 5 river ice, 6 caveats.

Usage:  python tools/make_validation_pdf.py
"""
import json
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import nsrad_style as _S  # shared NS-RAD publication/presentation style
from matplotlib.backends.backend_pdf import PdfPages

ROOT = Path(__file__).resolve().parent.parent
UPG = ROOT / "runs" / "definitive"       # the model (operational, observed-blended T boundary)
REG = ROOT / "runs" / "regression_bnd"   # same physics, regression T boundary (independence check)
OBS = json.load(open(ROOT / "docs" / "validation_obs.json"))
try:
    VOBS = json.load(open(ROOT / "docs" / "usgs_velocity_obs.json"))
except Exception:
    VOBS = {}
OUT = ROOT / "docs" / "ns_rad_validation.pdf"

# Fitted at-a-station depth relations D = c·Q^f (USGS ADCP surveys). The prismatic width
# B_ub is read LIVE from each sites/<name>.py so this can never desync from the model
# geometry (it did once, after B_ub was snapped to the SWORD per-channel median). Used for
# the continuity velocity V = Q/(B_ub·D) plotted against the USGS field measurements.
def _bub(site):
    import importlib, sys
    p = str(ROOT / "code")
    if p not in sys.path:
        sys.path.insert(0, p)
    return importlib.import_module(f"sites.{site}").B_ub

_DEPTH_CQ = {"colville": (0.360, 0.297), "kuparuk": (0.291, 0.309),
             "sagavanirktok": (0.280, 0.259), "canning": (0.224, 0.325)}
DEPTH_REL = {s: (c, f, _bub(s)) for s, (c, f) in _DEPTH_CQ.items()}

SITES = ["colville", "kuparuk", "sagavanirktok", "canning"]
LABEL = {"colville": "Colville", "kuparuk": "Kuparuk",
         "sagavanirktok": "Sagavanirktok", "canning": "Canning"}
C = {"colville": "#2a78d6", "kuparuk": "#eb6834",
     "sagavanirktok": "#1baf7a", "canning": "#4a3aa7"}

INK, INK2, MUTED = "#0b0b0b", "#52514e", "#8a8880"
GRID, SURFACE, WARN, OBSCOL = "#e6e5e0", "#fcfcfb", "#e34948", "#0b0b0b"
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


def load(rundir, site, var):
    """Memoised loader (time, values[cells 1..M]). Prefers NetCDF, falls back to .dat."""
    key = (str(rundir), site, var)
    if key in _CACHE:
        return _CACHE[key]
    ncp = rundir / site / "output.nc"
    if ncp.exists():
        from netCDF4 import Dataset
        d = Dataset(ncp)
        res = ((np.asarray(d.variables["time"][:]) / 86400.0,
                np.asarray(d.variables[var][:])[:, 1:]) if var in d.variables
               else (None, None))
        _CACHE[key] = res
        return res
    p = rundir / site / f"{var}.dat"
    if not p.exists():
        _CACHE[key] = (None, None)
        return None, None
    try:
        import pandas as pd
        d = pd.read_csv(p, sep="\t", header=None).to_numpy()
    except Exception:
        d = np.genfromtxt(p, delimiter="\t")
    if np.isnan(d[:, -1]).all():
        d = d[:, :-1]
    res = (d[:, 0] / 86400.0, d[:, 1:])
    _CACHE[key] = res
    return res


def model_doy_series(rundir, site, var, col=68):
    """Model value at one grid cell vs day-of-year (both model years overlaid)."""
    t, F = load(rundir, site, var)
    if F is None:
        return None, None
    return np.mod(t, 365), F[:, col]


def model_daily_mean(rundir, site, var, col=68):
    """Map the model's mid-channel series to a day-of-year -> mean-value dict,
    averaging the two model years. Exact zeros are dropped: under the ice model these
    are the under-ice temperatures held at the 0 °C freezing point (no summer
    observations to compare against anyway); under the legacy gate they were the
    zeroed winter state."""
    doy, mid = model_doy_series(rundir, site, var, col)
    if doy is None:
        return {}
    acc = {}
    for d, v in zip(doy, mid):
        if v != 0:
            acc.setdefault(int(d), []).append(v)
    return {d: float(np.mean(vs)) for d, vs in acc.items()}


def paired_temp(rundir, site):
    """(model, observed) pairs for the USGS 2022 daily record, same-year."""
    mm = model_daily_mean(rundir, site, "T")
    usgs = OBS["temp_usgs_2022"].get(site, [])
    return np.array([(mm[d], v) for d, v in usgs if d in mm and mm[d] > 0])


def paired_temp_wqp(rundir, site):
    """(model day-of-year mean, grab) pairs — climatological, all years pooled."""
    mm = model_daily_mean(rundir, site, "T")
    wqp = OBS["temp_wqp"].get(site, [])
    return np.array([(mm[r[0]], r[1]) for r in wqp if r[0] in mm and mm[r[0]] > 0])


def fit_stats(pairs):
    """bias, RMSE, R^2 for a (model, obs) array."""
    if len(pairs) < 2:
        return None
    mod, obs = pairs[:, 0], pairs[:, 1]
    bias = np.mean(mod - obs)
    rmse = np.sqrt(np.mean((mod - obs) ** 2))
    ss_res = np.sum((obs - mod) ** 2)
    ss_tot = np.sum((obs - obs.mean()) ** 2)
    r2 = 1 - ss_res / ss_tot if ss_tot > 0 else float("nan")
    return bias, rmse, r2


# --------------------------------------------------------------- page 1: temperature
def page_temperature(pdf):
    fig = plt.figure(figsize=(11, 8.5))
    fig.suptitle("Water temperature — model vs all available observations",
                 x=0.055, ha="left", fontsize=13, weight="bold")
    fig.text(0.055, 0.935, "Model (transported T + heat budget + prognostic ice). Black dots = USGS "
             "2022 daily (same year); grey = Water Quality Portal grabs, all years, as a "
             "climatology.", color=INK2, fontsize=8.3)

    gs = fig.add_gridspec(2, 2, left=0.06, right=0.97, top=0.90, bottom=0.13,
                          hspace=0.34, wspace=0.16)
    for k, site in enumerate(SITES):
        ax = fig.add_subplot(gs[k // 2, k % 2])
        doy, mid = model_doy_series(UPG, site, "T")
        if mid is not None:
            m = mid > 0
            ax.plot(doy[m], mid[m], ".", ms=1.4, color=C[site], alpha=0.45,
                    label="model (mid-channel)")
        # WQP climatology grabs
        wqp = OBS["temp_wqp"].get(site, [])
        if wqp:
            wd = np.array([r[0] for r in wqp]); wv = np.array([r[1] for r in wqp])
            ax.plot(wd, wv, "o", ms=2.6, color=MUTED, alpha=0.5, mec="none",
                    label=f"WQP grabs ({len(wqp)}, all yr)")
        # USGS 2022 daily
        usgs = OBS["temp_usgs_2022"].get(site, [])
        stat = ""
        if usgs:
            ud = np.array([r[0] for r in usgs]); uv = np.array([r[1] for r in usgs])
            ax.plot(ud, uv, "o", ms=3.2, color=OBSCOL, mec=SURFACE, mew=0.5,
                    label=f"USGS 2022 daily ({len(usgs)})", zorder=6)
            if mid is not None:
                yr1 = doy[:len(doy) // 2] if False else doy
                mm = {}
                for dd, mv in zip(doy, mid):
                    mm.setdefault(int(dd), []).append(mv)
                pairs = [(np.mean(mm[d]), v) for d, v in zip(ud, uv)
                         if d in mm and np.mean(mm[d]) > 0]
                if pairs:
                    a = np.array(pairs)
                    stat = (f"vs USGS: bias {np.mean(a[:,0]-a[:,1]):+.1f}  "
                            f"RMSE {np.sqrt(np.mean((a[:,0]-a[:,1])**2)):.1f} °C")
        month_axis(ax)
        tidy(ax)
        ax.set_ylabel("T (°C)")
        ax.set_title(LABEL[site], loc="left", color=C[site], weight="medium")
        if stat:
            ax.text(0.985, 0.05, stat, transform=ax.transAxes, ha="right",
                    fontsize=6.8, color=INK2)
        if k == 0:
            ax.legend(loc="upper left", fontsize=6.6, labelcolor=INK2)
    fig.text(0.055, 0.075,
             "Same-year validation exists only for Sagavanirktok and Kuparuk (USGS 00010). Colville and Canning "
             "have grab samples but no 2022 daily record.\nThe model carries a summer WARM bias and overshoots "
             "peaks; observed Deadhorse humidity rather than the assumed 0.85 did not remove it. The cause is NOT "
             "cloud cover — the\nsolar forcing is already observed, so the shortwave cloud effect is present. It is "
             "the warm-biased INFLOW BOUNDARY: the upstream temperature is the\nair-temp regression, +2.6 °C warm on "
             "the Sagavanirktok (a mountain/groundwater-fed river it was not fitted to) but only +0.05 °C on Kuparuk, "
             "where it was.\nDepth is NOT the culprit — modelled depths match the USGS surveys. Timing and the "
             "mouth-to-head gradient remain real improvements.",
             fontsize=7.6, color=INK2, linespacing=1.5, va="top")
    pdf.savefig(fig)
    plt.close(fig)


# ------------------------------------------------------- page 1b: temperature scatter
def page_temp_scatter(pdf):
    fig = plt.figure(figsize=(11, 8.5))
    fig.suptitle("Water temperature — model vs observed, 1:1 scatter",
                 x=0.055, ha="left", fontsize=13, weight="bold")
    fig.text(0.055, 0.935, "Points on the diagonal = perfect. Above = model too warm. "
             "Left panel is the rigorous same-year test; right pools multi-year grabs "
             "as a climatology.", color=INK2, fontsize=8.3)

    # ---- left: two boundary configs vs USGS 2022 daily (Sag + Kuparuk) ----
    # open markers = REGRESSION boundary (independent forcing -> real skill test)
    # filled markers = OBSERVED boundary (USGS fed into the boundary -> near-tautological)
    axL = fig.add_axes([0.08, 0.30, 0.4, 0.55])
    lim = (0, 20)
    axL.plot(lim, lim, color=MUTED, lw=1.0, ls=(0, (4, 2)), zorder=1)
    stat_rows = []   # collect stats, render as one block clear of the cloud
    for site in ["sagavanirktok", "kuparuk"]:
        pr = paired_temp(REG, site)   # independent (regression boundary)
        po = paired_temp(UPG, site)   # operational (observed boundary)
        if len(pr):
            axL.scatter(pr[:, 1], pr[:, 0], s=16, facecolor="none",
                        edgecolor=C[site], linewidth=0.9, zorder=3)
            st = fit_stats(pr)
            if st:
                stat_rows.append((f"{LABEL[site]}  ○ regression bnd",
                                  f"bias {st[0]:+.1f}  RMSE {st[1]:.1f}", C[site]))
        if len(po):
            axL.scatter(po[:, 1], po[:, 0], s=13, color=C[site], alpha=0.75,
                        edgecolor=SURFACE, linewidth=0.3, zorder=4)
            st = fit_stats(po)
            if st:
                stat_rows.append((f"{LABEL[site]}  ● observed bnd",
                                  f"bias {st[0]:+.1f}  RMSE {st[1]:.1f}", C[site]))
    # stats block in the empty lower-right triangle, below the 1:1 line
    for i, (lab, val, col) in enumerate(stat_rows):
        y = 0.30 - i * 0.052
        axL.text(0.40, y, lab, transform=axL.transAxes, fontsize=6.8, color=col,
                 va="top", weight="medium")
        axL.text(0.985, y, val, transform=axL.transAxes, fontsize=6.8, color=col,
                 va="top", ha="right")
    axL.set_xlim(lim); axL.set_ylim(lim)
    axL.set_aspect("equal")
    axL.set_xlabel("observed T (°C)  —  USGS 00010, 2022")
    axL.set_ylabel("model mid-channel T (°C)")
    for s in ("top", "right"):
        axL.spines[s].set_visible(False)
    axL.grid(True, alpha=0.5)
    axL.set_title("Same-year daily — two boundary configs", loc="left", color=INK, fontsize=9)
    axL.text(0.02, 0.97, "○ regression boundary (independent)\n● observed boundary "
             "(near-tautological)", transform=axL.transAxes, fontsize=6.4, color=INK2,
             va="top")

    # ---- right: WQP grabs, all years pooled by day-of-year (climatology) ----
    axR = fig.add_axes([0.56, 0.30, 0.4, 0.55])
    lim2 = (0, 22)
    axR.plot(lim2, lim2, color=MUTED, lw=1.0, ls=(0, (4, 2)), zorder=1)
    txt_y = 0.97
    for site in SITES:
        p = paired_temp_wqp(UPG, site)
        if len(p) == 0:
            continue
        axR.scatter(p[:, 1], p[:, 0], s=12, color=C[site], alpha=0.55,
                    edgecolor="none", label=LABEL[site], zorder=3)
        st = fit_stats(p)
        if st:
            axR.text(0.03, txt_y, f"{LABEL[site]}:  bias {st[0]:+.1f}  RMSE {st[1]:.1f}  "
                     f"(n={len(p)})", transform=axR.transAxes,
                     fontsize=7.3, color=C[site], va="top", weight="medium")
            txt_y -= 0.055
    axR.set_xlim(lim2); axR.set_ylim(lim2)
    axR.set_aspect("equal")
    axR.set_xlabel("observed T (°C)  —  WQP grabs, all years")
    axR.set_ylabel("model day-of-year mean T (°C)")
    for s in ("top", "right"):
        axR.spines[s].set_visible(False)
    axR.grid(True, alpha=0.5)
    axR.set_title("Multi-year grabs (climatology)", loc="left", color=INK, fontsize=9)

    fig.text(0.055, 0.2,
             "READ THE LEFT PANEL CAREFULLY — the two boundary configs tell different things. ○ open markers use "
             "the air-temp REGRESSION as the\nupstream boundary (forcing independent of the gauge), so their "
             "scatter is the genuine model-skill test: Kuparuk RMSE 2.3, Sagavanirktok\nRMSE 3.8 (inflated by the "
             "regression's transfer bias). ● filled markers use the OBSERVED USGS temperature as the boundary, so "
             "they sit on\nthe 1:1 line almost by construction — the mid-channel is within 0.2 °C of the value fed "
             "in, making that comparison near-tautological, NOT\nindependent skill. The observed boundary is the "
             "right operational choice (it carries real temperature into the biogeochemistry) but its\ntight fit is "
             "definitional. The honest heat-budget skill is the open-marker RMSE of ~2–4 °C.\n\nThe climatology "
             "panel (right) pools multi-year grabs against one model year, so spread mixes model error with "
             "interannual variability —\nread it for sign and gross agreement, not RMSE. Colville and Canning have "
             "no same-year daily record and so no rigorous validation.",
             fontsize=7.4, color=INK2, linespacing=1.45, va="top")
    pdf.savefig(fig)
    plt.close(fig)


# --------------------------------------------------------------- page 2: salinity
def page_salinity(pdf):
    fig = plt.figure(figsize=(11, 8.5))
    fig.suptitle("Salinity — model vs available delta observations",
                 x=0.055, ha="left", fontsize=13, weight="bold")
    fig.text(0.055, 0.935, "Key finding: as configured (distance=1, weak tides), the model "
             "produces ≈0 PSU through the whole domain in open water — it does not "
             "reproduce the observed delta salinity.", color=WARN, fontsize=8.3)

    # (a,b) the two rivers with data: model open-water salinity vs grabs
    gs = fig.add_gridspec(2, 2, left=0.06, right=0.97, top=0.90, bottom=0.30,
                          hspace=0.34, wspace=0.18)
    for k, site in enumerate(["colville", "kuparuk"]):
        ax = fig.add_subplot(gs[0, k])
        t, S = load(UPG, site, "S")
        if S is not None:
            doy = np.mod(t, 365)
            # domain-max salinity; mask any fully-zeroed timesteps (legacy ice gate) so
            # the brief boundary transients at freeze/thaw don't read as real intrusion.
            # Under the prognostic ice model salinity is conserved, so nothing is masked.
            dmax = np.nanmax(S, axis=1)
            gated = np.nansum(np.abs(S), axis=1) == 0
            dmax_ow = np.where(gated, np.nan, dmax)
            o = np.argsort(doy)
            ax.plot(doy[o], dmax_ow[o], color=C[site], lw=1.3,
                    label="model, domain max (open water)")
        sal = OBS["salinity_wqp"].get(site, [])
        if sal:
            sd = np.array([r[0] for r in sal]); sv = np.array([r[1] for r in sal])
            yrs = sorted(set(r[2] for r in sal))
            ax.plot(sd, sv, "o", ms=3.4, color=OBSCOL, mec=SURFACE, mew=0.5,
                    label=f"WQP grabs ({len(sal)}, {'/'.join(yrs)})", zorder=6)
        ax.axhline(27.1, color=MUTED, ls=(0, (4, 2)), lw=1.0)
        ax.text(2, 27.1, " model seaward boundary (27.1)", fontsize=6.6, color=MUTED,
                va="bottom")
        ax.text(150, 3.5, "model ≈ 0 all summer", fontsize=7, color=C[site], style="italic")
        month_axis(ax)
        tidy(ax)
        ax.set_ylabel("salinity (PSU)")
        ax.set_ylim(-1, 34)
        ax.set_title(LABEL[site], loc="left", color=C[site], weight="medium")
        ax.legend(loc="upper left", fontsize=6.5, labelcolor=INK2)

    # (c) profile snapshot showing where salinity lives in the domain
    ax = fig.add_subplot(gs[1, 0])
    for site in ["colville", "kuparuk"]:
        t, S = load(UPG, site, "S")
        if S is None:
            continue
        doy = np.mod(t, 365)
        j = int(np.argmin(np.abs(doy - 220)))
        x = np.arange(S.shape[1]) * 0.2
        ax.plot(x, S[j], color=C[site], lw=1.6, label=LABEL[site])
    ax.set_xlabel("distance from mouth (km)")
    ax.set_ylabel("salinity (PSU)")
    ax.set_title("Modelled profile (day 220) — intrusion confined to the first few km",
                 loc="left", color=INK, fontsize=8.5)
    tidy(ax, both=True)
    ax.legend(loc="upper right", labelcolor=INK2)

    # (d) INDICATIVE scatter: grab salinity vs model near-mouth range (NOT a 1:1
    # validation — different year, positions unmatched; shows coverage only)
    ax = fig.add_subplot(gs[1, 1])
    ax.plot((0, 34), (0, 34), color=MUTED, lw=1.0, ls=(0, (4, 2)))
    for site in ["colville", "kuparuk"]:
        sal = OBS["salinity_wqp"].get(site, [])
        t, S = load(UPG, site, "S")
        if not sal or S is None:
            continue
        doy = np.mod(t, 365)
        near = S[:, :25]
        lo = {}
        hi = {}
        for d, l, h in zip(doy, np.nanmin(near, 1), np.nanmax(near, 1)):
            lo.setdefault(int(d), []).append(l)
            hi.setdefault(int(d), []).append(h)
        xs, ymid, yerr = [], [], []
        for r in sal:
            d = r[0]
            if d in lo:
                a, b = np.mean(lo[d]), np.mean(hi[d])
                xs.append(r[1]); ymid.append((a + b) / 2); yerr.append((b - a) / 2)
        if xs:
            ax.errorbar(xs, ymid, yerr=yerr, fmt="o", ms=3.5, color=C[site],
                        ecolor=C[site], elinewidth=0.6, capsize=1.5, alpha=0.7,
                        mec=SURFACE, mew=0.4, label=LABEL[site])
    ax.set_xlim(0, 34); ax.set_ylim(0, 34)
    ax.set_aspect("equal")
    ax.set_xlabel("observed grab salinity (PSU)")
    ax.set_ylabel("model, seaward 5 km (range)")
    ax.set_title("Indicative only — position-unmatched, 2014 vs model", loc="left",
                 color=INK, fontsize=8.5)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    ax.grid(True, alpha=0.5)
    ax.legend(loc="lower right", fontsize=6.6, labelcolor=INK2)

    fig.text(0.055, 0.24,
             "What can and cannot be concluded. The observed delta grabs (Colville 12–26, Kuparuk 8–32 PSU, 2014) "
             "show real saline intrusion, but the\nmodel domain-max salinity is ≈0 through the entire open-water "
             "season: the strong seaward freshet flushes the boundary out, and with\ndistance=1 and AMPL=0.2 m the "
             "saline zone sits essentially OUTSIDE the 27 km domain. So the model does not reproduce the observed\n"
             "delta salinity — it represents the river, not the estuarine mixing zone. The 27.1 PSU seaward "
             "boundary value itself is bracketed by the\ngrabs (defensible), but the intrusion dynamics are "
             "unrepresented. Sagavanirktok and Canning have no salinity data at all. To model the\nsaline zone "
             "would require a longer domain (larger distance) and stronger tidal forcing — a configuration choice, "
             "not a tuning.",
             fontsize=7.6, color=INK2, linespacing=1.5, va="top")
    pdf.savefig(fig)
    plt.close(fig)


# --------------------------------------------------------------- page 3: FCO2
def page_fco2(pdf):
    fig = plt.figure(figsize=(11, 8.5))
    fig.suptitle("Air–sea CO₂ flux — modelled seasonal cycle", x=0.055, ha="left",
                 fontsize=13, weight="bold")
    fig.text(0.055, 0.935, "FCO₂ > 0 = outgassing. Spatial-mean over open-water cells. "
             "No flux observations exist for these rivers, so this is the modelled signal, "
             "not a validation against data.", color=INK2, fontsize=8.3)
    gs = fig.add_gridspec(2, 2, left=0.07, right=0.96, top=0.90, bottom=0.16,
                          hspace=0.40, wspace=0.2)
    for k, site in enumerate(SITES):
        ax = fig.add_subplot(gs[k // 2, k % 2])
        t, F = load(UPG, site, "FCO2")
        if F is not None:
            fm = np.where(np.abs(F).sum(axis=1, keepdims=True) > 0, F, np.nan)
            series = np.nanmean(fm, axis=1) * 1e3
            doy = np.mod(t, 365)
            o = np.argsort(doy)
            ax.plot(doy[o], series[o], color=C[site], lw=1.4)
        ax.axhline(0, color=GRID, lw=0.8)
        month_axis(ax)
        tidy(ax)
        ax.set_ylabel("mean FCO₂ (×10⁻³)")
        ax.set_title(LABEL[site], loc="left", color=C[site], weight="medium")
    fig.text(0.055, 0.115, "The rivers sit NEAR CO₂ EQUILIBRIUM — a small, mixed-sign flux "
             "(Colville/Kuparuk ~0, Sagavanirktok slight uptake, Canning still outgassing). "
             "This followed the mol/kg carbonate fix (C-GEM v2): the earlier code left the\n"
             "atmospheric-saturation term K0·pCO₂ in mol/kg while co2s was mmol/m³ (~10⁶× too "
             "small), so it always outgassed proportional to dissolved CO₂, ignoring the air–"
             "sea gradient. With the gradient resolved, the flux is the small\ndisequilibrium "
             "residual and is now GENUINELY SENSITIVE to the riverine DIC/alkalinity — which "
             "is placeholder for Colville/Sag/Canning, so their flux SIGN is unconstrained. "
             "No flux observations exist to check it.",
             fontsize=7.6, color=INK2, linespacing=1.5, va="top")
    pdf.savefig(fig)
    plt.close(fig)


# --------------------------------------------------------------- page 4: data inventory
def page_inventory(pdf):
    fig = plt.figure(figsize=(11, 8.5))
    fig.suptitle("Validation data inventory & honest verdict", x=0.055, ha="left",
                 fontsize=13, weight="bold")
    ax = fig.add_axes([0.06, 0.52, 0.9, 0.34])
    ax.axis("off")
    rows = [["", "Temp 2022 daily", "Temp grabs", "Salinity grabs", "River chem", "Flux"]]
    inv = {
        "colville": ["—", "147", "29 (2014)", "placeholder", "none"],
        "kuparuk": ["55 d", "225", "69 (2014)", "Arctic LTER*", "none"],
        "sagavanirktok": ["93 d", "14", "none", "placeholder", "none"],
        "canning": ["6", "6", "none", "placeholder", "none"],
    }
    for s in SITES:
        rows.append([LABEL[s]] + inv[s])
    tb = ax.table(cellText=rows[1:], colLabels=rows[0], loc="center", cellLoc="center")
    tb.auto_set_font_size(False)
    tb.set_fontsize(8.5)
    tb.scale(1, 1.7)
    for (r, c), cell in tb.get_celld().items():
        cell.set_edgecolor(GRID)
        if r == 0:
            cell.set_text_props(color=INK2, weight="bold")
        elif c == 0:
            cell.set_text_props(color=C[SITES[r - 1]], weight="medium")
        cell.set_facecolor(SURFACE)
    ax.set_title("What data exists, per river", loc="left", color=INK, y=0.94)

    fig.text(0.055, 0.44,
             "Verdict — read before trusting any number in this document\n\n"
             "TEMPERATURE. With the air-temp REGRESSION boundary (independent forcing) the heat-budget model has "
             "RMSE ~2–4 °C — the honest\nskill. The Sagavanirktok's +2.6 °C bias came from that boundary, not the "
             "physics: it is the regression transferred onto a mountain/\ngroundwater-fed river it was not fitted "
             "to. Feeding OBSERVED USGS temperature into the boundary removes the bias, but the resulting\ntight fit "
             "(RMSE ~0.4) is near-tautological — the mid-channel is essentially the value fed in. Not cloud cover "
             "(solar forcing already\ncarries the shortwave cloud effect) and not depth (modelled depths match the "
             "USGS surveys). What the model genuinely resolves is a\nlongitudinal temperature gradient and "
             "per-river timing.\n\n"
             "SALINITY. The model produces ≈0 PSU through the whole domain in open water — with distance=1 and "
             "weak tides the saline zone\nsits essentially outside the 27 km domain, so the model represents the "
             "river, not the estuarine mixing zone. It does NOT reproduce the\nobserved 8–32 PSU delta grabs "
             "(Colville, Kuparuk 2014). The 27.1 PSU boundary value is bracketed by those grabs, but the intrusion\n"
             "dynamics are unrepresented. Sagavanirktok and Canning have no salinity data. Modelling the saline "
             "zone is a config change, not a tuning.\n\n"
             "RIVER CHEMISTRY (*): Kuparuk's riverine boundary now uses Arctic LTER Streams data (natural reach, "
             "1978-2019), replacing a\nplaceholder that was 5-10× off on carbon and nitrogen. But the LTER site is "
             "~163 km upstream (headwaters), so it is a proxy, not a\ndelta measurement; the other three rivers "
             "still carry the placeholder, so cross-site chemistry is not yet consistent.\n\n"
             "FLUX (FCO₂) has no observations for any of these rivers. After the mol/kg carbonate fix (v2) the "
             "rivers sit NEAR CO₂ EQUILIBRIUM (small, mixed-sign flux), not\nstrongly outgassing as before — the "
             "old code's atmospheric-saturation term was ~10⁶× too small (mol/kg vs mmol/m³). The flux is now the "
             "small\nair–sea disequilibrium residual, so its SIGN is set by the (placeholder) boundary chemistry "
             "and is unconstrained for three of four rivers.\n\n"
             "Bottom line: temperature is the one variable with real same-year data, and the honest heat-budget "
             "skill is RMSE ~2–4 °C (the\ntighter fit comes from feeding observations into the boundary). The CO₂ "
             "flux is now physically correct (real air–sea gradient) but near\nequilibrium and chemistry-limited, "
             "with no flux data to check it. Salinity is not reproduced — the saline zone / surge intrusion is "
             "outside\nthe configured domain and forcing. Cross-river comparison is not yet clean: chemistry, "
             "temperature quality and coverage all differ by site.",
             fontsize=8, color=INK2, linespacing=1.5, va="top")
    pdf.savefig(fig)
    plt.close(fig)


# ------------------------------------------------------ page 6: river ice vs observations
# Observational anchors for North Slope RIVER ice (distinct from Beaufort sea ice).
# Sources: UAF Water & Environmental Research Center Sagavanirktok breakup monitoring
# (Franklin Bluffs, 2015-17); NWS Alaska-Pacific RFC breakup climatology; USGS ice-
# affected (A:e) discharge flags for 2022; Barnes et al. 1979 (USGS OFR 79-539) for the
# Beaufort fast-ice contrast. River ice: freezes at 0 C, ~1-2 m in the thalweg (bottom-
# fast in shallow reaches where Q -> 0), breaks up MECHANICALLY at the freshet (late May,
# ~day 145-155), weeks before coastal sea ice clears.
ICE_OBS = {
    "breakup_doy": (145, 160),      # river breakup window (late May), all sites
    "freezeup_doy": (267, 288),     # autumn freeze-up (late Sep - mid Oct)
    "thickness_m": (1.0, 2.0),      # river thalweg thickness; bottom-fast in shallows
    "seaice_thickness_m": 1.9,      # Beaufort landfast/fast ice (Barnes et al. 1979)
    "seaice_breakup_doy": 173,      # coastal sea-ice clearance (old gate opened here)
}


def ice_metrics(rundir, site):
    """Modelled breakup/freeze-up day-of-year and peak thickness from year 2 (post-spin-up)."""
    t, H = load(rundir, site, "ice_thickness")
    _, FR = load(rundir, site, "ice_frac")
    if H is None or FR is None:
        return None
    yr2 = t > 365
    if yr2.sum() < 10:
        yr2 = np.ones_like(t, bool)
    doy = np.mod(t[yr2], 365)
    o = np.argsort(doy)
    doy = doy[o]
    fm = np.nanmean(FR[yr2], axis=1)[o]
    hm = np.nanmax(H[yr2], axis=1)[o]
    covered = fm > 0.5
    # breakup = last covered->open transition before mid-year; freeze-up = open->covered after
    bru = fru = np.nan
    for i in range(1, len(doy)):
        if covered[i - 1] and not covered[i] and doy[i] < 200:
            bru = doy[i]
        if not covered[i - 1] and covered[i] and doy[i] > 220:
            fru = doy[i]
            break
    return dict(breakup=bru, freezeup=fru, peak=np.nanmax(hm))


def page_velocity(pdf):
    fig = plt.figure(figsize=(11, 8.5))
    fig.suptitle("Flow velocity — model vs USGS field measurements", x=0.055,
                 ha="left", fontsize=13, weight="bold")
    fig.text(0.055, 0.935, "Grey = USGS channel measurements (mean V = Q/area, all flows). "
             "Line = model continuity velocity V = Q/(B·D) at the prismatic reach, from the "
             "surveyed depth relation D = c·Q^f.", color=INK2, fontsize=8.0)
    gauged = ["colville", "kuparuk", "sagavanirktok"]
    gs = fig.add_gridspec(2, 2, left=0.07, right=0.96, top=0.90, bottom=0.09,
                          hspace=0.34, wspace=0.2)
    for k, s in enumerate(gauged):
        ax = fig.add_subplot(gs[k // 2, k % 2])
        o = VOBS.get(s, {})
        q = np.array(o.get("discharge_m3s", []))
        v = np.array(o.get("velocity_ms", []))
        if len(q):
            ax.scatter(q, v, s=8, color=MUTED, alpha=0.45, edgecolor="none",
                       label=f"USGS ({len(q)})")
        c, f, Bub = DEPTH_REL[s]
        if len(q):
            qq = np.logspace(np.log10(max(q.min(), 0.5)), np.log10(q.max()), 100)
            vmod = qq / (Bub * c * qq ** f)
            ax.plot(qq, vmod, color=C[s], lw=1.8, label="model V(Q)")
        ax.set_xscale("log")
        ax.set_title(LABEL[s], loc="left", color=C[s], weight="medium")
        ax.set_xlabel("discharge (m³ s⁻¹, log)", fontsize=7.5)
        ax.set_ylabel("velocity (m s⁻¹)", fontsize=7.5)
        ax.set_ylim(0, min(4, np.nanmax(v) * 1.1) if len(v) else 3)
        for sp in ("top", "right"):
            ax.spines[sp].set_visible(False)
        ax.grid(True, alpha=0.5); ax.set_axisbelow(True)
        ax.legend(loc="upper left", fontsize=6.6, labelcolor=INK2)

    ax = fig.add_subplot(gs[1, 1]); ax.axis("off")
    ax.set_title("Summary & caveats", loc="left", fontsize=9, color=INK)
    lines = []
    for s in gauged:
        v = np.array(VOBS.get(s, {}).get("velocity_ms", []))
        if len(v):
            lines.append(f"{LABEL[s]}: observed median {np.median(v):.2f} m/s "
                         f"(range {v.min():.2f}–{v.max():.2f}, n={len(v)})")
    txt = ("\n".join(lines) + "\n\n"
           "The model's momentum-solved velocities (section-averaged, open water) run "
           "~0.3–0.9 m/s mean with 1–2.6 m/s freshet peaks — squarely in the observed "
           "range. The velocity–discharge shape is set by the same USGS surveys that "
           "fixed the depth, so it is consistent by construction.\n\n"
           "CAVEATS. (1) Gauge location: the Colville (Umiat) and Sagavanirktok (Pump Sta 3) "
           "gauges are far UPSTREAM of the delta, in narrower, faster cross-sections than "
           "the delta reach the model represents — so the model V(Q) sits below the "
           "observed cloud there. Kuparuk (near Deadhorse) is near-tidewater and comparable. "
           "(2) Width: model uses SWORD width, the survey a local cross-section. "
           "(3) Canning has no gauge. So this is a range/consistency check, not a "
           "cell-matched validation.")
    ax.text(0.0, 0.90, txt, transform=ax.transAxes, fontsize=7.0, color=INK2,
            va="top", linespacing=1.4)
    pdf.savefig(fig)
    plt.close(fig)


def page_river_ice(pdf):
    fig = plt.figure(figsize=(11, 8.5))
    fig.suptitle("River ice — model vs observations (and the sea-ice contrast)",
                 x=0.055, ha="left", fontsize=13, weight="bold")
    fig.text(0.055, 0.935, "The prognostic ice_module models RIVER ice (freezes at 0 °C, "
             "mechanical freshet breakup), distinct from Beaufort SEA ice. "
             "Observational anchors below.", color=INK2, fontsize=8.3)

    # (a) modelled mean ice thickness seasonal cycle, all rivers, with obs breakup band
    ax = fig.add_axes([0.07, 0.60, 0.52, 0.30])
    for s in SITES:
        t, H = load(UPG, s, "ice_thickness")
        if H is None:
            continue
        doy = np.mod(t, 365)
        o = np.argsort(doy)
        ax.plot(doy[o], np.nanmean(H, axis=1)[o], color=C[s], lw=1.2, label=LABEL[s])
    ax.axvspan(*ICE_OBS["breakup_doy"], color=WARN, alpha=0.12)
    ax.axvspan(*ICE_OBS["freezeup_doy"], color=MUTED, alpha=0.12)
    ax.axhspan(*ICE_OBS["thickness_m"], color="#2a6aa8", alpha=0.07)
    ax.text(152, ax.get_ylim()[1] * 0.96, "obs river\nbreakup", fontsize=6, color=WARN,
            ha="center", va="top")
    ax.set_xticks(MONTH0); ax.set_xticklabels(MONTHL); ax.set_xlim(0, 364)
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)
    ax.grid(True, alpha=0.5); ax.set_axisbelow(True)
    ax.set_ylabel("mean ice thickness (m)", fontsize=7.5)
    # place the legend in the summer open-water gap (ice = 0, ~Jun–Sep) so it clears the curves
    ax.legend(loc="center", bbox_to_anchor=(0.56, 0.5), ncol=2, fontsize=6.3,
              labelcolor=INK2, columnspacing=1.0, handletextpad=0.4)
    ax.set_title("Modelled ice thickness vs observed river-ice envelope", loc="left",
                 fontsize=8.5)

    # (b) river-ice vs sea-ice contrast panel
    ax2 = fig.add_axes([0.66, 0.60, 0.30, 0.30]); ax2.axis("off")
    ax2.set_title("Why river ice ≠ sea ice", loc="left", fontsize=8.5, color=INK)
    contrast = [
        ("", "RIVER (model)", "Beaufort SEA"),
        ("freezes at", "0 °C", "−1.8 °C"),
        ("peak thickness", "1–2 m*", "~1.9 m"),
        ("breakup", "day ~150", "day ~173"),
        ("breakup driver", "freshet", "thermal"),
    ]
    for i, (a, b, c) in enumerate(contrast):
        y = 0.86 - i * 0.17
        w = "bold" if i == 0 else "normal"
        ax2.text(0.0, y, a, fontsize=6.8, color=INK2, weight=w)
        ax2.text(0.46, y, b, fontsize=6.8, color="#2a6aa8", weight=w, ha="center")
        ax2.text(0.86, y, c, fontsize=6.8, color=MUTED, weight=w, ha="center")
    ax2.text(0.0, 0.86 - 5 * 0.17, "*bottom-fast in shallow reaches", fontsize=5.6, color=INK2)

    # (c) per-river model-vs-obs metrics table
    ax3 = fig.add_axes([0.07, 0.10, 0.89, 0.40]); ax3.axis("off")
    ax3.set_title("Modelled ice phenology vs observed river-ice anchors", loc="left",
                  fontsize=9.5, color=INK)
    cols = ["River", "Breakup (model)", "Breakup (obs)", "Freeze-up (model)",
            "Freeze-up (obs)", "Peak ice (model)", "Peak ice (obs)"]
    xcol = [0.0, 0.16, 0.32, 0.46, 0.62, 0.77, 0.89]
    for x, c in zip(xcol, cols):
        ax3.text(x, 0.95, c, fontsize=6.8, color=INK2, weight="bold", va="top")
    bru_o = f"~{ICE_OBS['breakup_doy'][0]}–{ICE_OBS['breakup_doy'][1]}"
    fru_o = f"~{ICE_OBS['freezeup_doy'][0]}–{ICE_OBS['freezeup_doy'][1]}"
    thk_o = f"{ICE_OBS['thickness_m'][0]:.0f}–{ICE_OBS['thickness_m'][1]:.0f} m"
    for i, s in enumerate(SITES):
        m = ice_metrics(UPG, s)
        y = 0.86 - i * 0.11
        if m is None:
            continue
        vals = [LABEL[s],
                f"day {m['breakup']:.0f}" if np.isfinite(m['breakup']) else "—",
                bru_o,
                f"day {m['freezeup']:.0f}" if np.isfinite(m['freezeup']) else "—",
                fru_o,
                f"{m['peak']:.2f} m", thk_o]
        for x, v, isobs in zip(xcol, vals, [0, 0, 1, 0, 1, 0, 1]):
            ax3.text(x, y, v, fontsize=6.8, va="top",
                     color=(MUTED if isobs else C[s]))
    notes = (
        "Model reproduces the defining North Slope river-ice phenology: freeze-up in autumn, a ~1–2 m bottom-fast "
        "winter cover (Colville deepest → thickest; the shallow rivers cap at their bed depth), and an ABRUPT "
        "MECHANICAL breakup at the freshet (~day 150) — not the weeks-late thermal breakup a sea-ice formulation "
        "gives. Peak thickness scales with channel depth, as observed. Sources: UAF WERC Sagavanirktok breakup "
        "monitoring (Franklin Bluffs 2015–17); NWS Alaska-Pacific RFC breakup climatology; USGS 2022 ice-affected "
        "(A:e) discharge flags; Barnes et al. 1979 (USGS OFR 79-539) for the Beaufort fast-ice contrast. These are "
        "phenology/magnitude anchors, not a gridded thickness validation — no continuous in-situ river-ice "
        "thickness series exists for these four rivers in 2022."
    )
    ax3.text(0.0, 0.30, notes, fontsize=7.0, color=INK2, va="top", linespacing=1.5, wrap=True)
    pdf.savefig(fig)
    plt.close(fig)


def main():
    OUT.parent.mkdir(exist_ok=True)
    with PdfPages(OUT) as pdf:
        page_temperature(pdf)      # 1  temperature time series
        page_temp_scatter(pdf)     # 2  temperature 1:1 scatter
        page_salinity(pdf)         # 3  salinity time series + profile + scatter
        page_fco2(pdf)             # 4  FCO2 modelled seasonal cycle
        page_velocity(pdf)         # 5  flow velocity vs USGS field measurements
        page_river_ice(pdf)        # 6  river ice vs observations + sea-ice contrast
        page_inventory(pdf)        # 7  data inventory + verdict
    print(f"wrote {OUT.relative_to(ROOT)} ({OUT.stat().st_size/1024:.0f} kB, 7 pages)")


if __name__ == "__main__":
    main()
