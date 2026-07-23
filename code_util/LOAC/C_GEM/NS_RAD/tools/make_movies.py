"""
Animations of the time-space evolution of the major diagnostics, one movie per river.

For a 1-D channel model the natural "movie" is the along-channel PROFILE animating
through time: distance from the mouth on x, field value on y, one frame per output
save. Each river gets a multi-panel GIF (docs/movies/<site>_evolution.gif) with the
core state variables + the prognostic ice cover animating together, plus a companion
Hovmoller GIF with a sweeping time cursor.

Reads the same runs/<site>/output.nc (or legacy .dat) as the PDF tools. Writes mp4
(H.264, yuv420p) by default -- the macOS-native format that QuickTime plays. macOS
QuickTime cannot open the AVI *container* at all (any AVI codec shows "not compatible"),
so AVI is offered only via --avi for tools that specifically need it. --gif for a
portable animated GIF.

Usage:
    python tools/make_movies.py                 # all four rivers, mp4 (QuickTime-native)
    python tools/make_movies.py kuparuk         # one river
    python tools/make_movies.py --avi kuparuk   # AVI (mpeg4) instead of mp4
    python tools/make_movies.py --gif kuparuk   # animated GIF
    python tools/make_movies.py --run /tmp/x kuparuk   # a specific run dir (single site)
"""
import sys
import json
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, PillowWriter

ROOT = Path(__file__).resolve().parent.parent
RUN = ROOT / "runs" / "definitive"
OUTDIR = ROOT / "docs" / "movies"

SITES = ["colville", "kuparuk", "sagavanirktok", "canning"]
LABEL = {"colville": "Colville", "kuparuk": "Kuparuk",
         "sagavanirktok": "Sagavanirktok", "canning": "Canning",
         "idealized": "Idealized (verification fixture)"}
C = {"colville": "#2a78d6", "kuparuk": "#eb6834",
     "sagavanirktok": "#1baf7a", "canning": "#4a3aa7",
     "idealized": "#7a3a8c"}
# `idealized` is the synthetic verification site (sites/idealized.py). It is not in
# SITES, so the default all-sites run ignores it; movie it explicitly against its run
# dir:  python tools/make_movies.py --run runs/idealized idealized
INK, INK2, MUTED = "#0b0b0b", "#52514e", "#8a8880"
GRID, SURFACE = "#e6e5e0", "#fcfcfb"
ICEBLUE = "#7fb2d9"
MEANCOL = "#e34948"      # observation-mean marker
plt.rcParams.update({
    "figure.facecolor": SURFACE, "axes.facecolor": SURFACE, "savefig.facecolor": SURFACE,
    "font.size": 8, "axes.titlesize": 9, "axes.labelsize": 8, "axes.edgecolor": GRID,
    "axes.linewidth": 0.8, "axes.labelcolor": INK2, "text.color": INK,
    "xtick.color": INK2, "ytick.color": INK2, "xtick.labelsize": 7, "ytick.labelsize": 7,
    "legend.frameon": False, "grid.color": GRID, "grid.linewidth": 0.6,
})
DELXI_KM = 0.2

# (variable, label, colormap-for-hovmoller). ice_* are skipped gracefully if absent.
PANELS = [
    ("T",   "temperature (°C)",        "inferno"),
    ("S",   "salinity (PSU)",          "viridis"),
    ("DIC", "DIC (mmol m⁻³)",          "plasma"),
    ("O2",  "O₂ (mmol m⁻³)",           "cividis"),
    ("pH",  "pH",                      "coolwarm"),
    # FCO2 is single-signed (all outgassing, >= 0 under the positive-out convention), so
    # a sequential map is correct -- a divergent RdBu would imply a sign change at its
    # white midpoint that never occurs. magma_r: near-zero flux light -> strong outgassing
    # dark; perceptually uniform, CVD-safe.
    ("FCO2", "air–sea CO₂ (>0 outgas)", "magma_r"),
    ("ice_thickness", "ice thickness (m)", "Blues"),
    ("ice_frac",      "ice fraction (0–1)", "Blues"),
]

_CACHE = {}

# Observations for the scatter overlay (same file the validation PDF uses).
#   temp_usgs_2022  same-year USGS river temperature (Sag, Kuparuk)
#   temp_wqp        Water Quality Portal grabs, all years pooled by day-of-year
#   salinity_wqp    WQP salinity grabs (brackish, near the mouth)
try:
    OBS = json.load(open(ROOT / "docs" / "validation_obs.json"))
except Exception:
    OBS = {}
OBS_WINDOW = 6           # +/- days-of-year a point is shown around the current frame
MIDX_KM = 68 * DELXI_KM  # mid-channel x used by the validation comparison (~13.6 km)


def site_obs(site):
    """Per-variable observations to overlay: {var: (pts[[doy,val],...], x_km, label)}.
    Temperature obs (river gauges + WQP) are placed at mid-channel to match the
    validation methodology; salinity grabs are brackish and placed near the mouth."""
    out = {}
    T = [(d, v) for d, v in OBS.get("temp_usgs_2022", {}).get(site, [])]
    T += [(r[0], r[1]) for r in OBS.get("temp_wqp", {}).get(site, [])]
    if T:
        out["T"] = (np.array(T, float), MIDX_KM, "USGS/WQP obs")
    S = [(r[0], r[1]) for r in OBS.get("salinity_wqp", {}).get(site, [])]
    if S:
        out["S"] = (np.array(S, float), 1.0, "WQP grab")
    return out


def load(sitedir, var):
    key = (str(sitedir), var)
    if key in _CACHE:
        return _CACHE[key]
    res = _load_nc(sitedir, var)
    if res is None or res[1] is None:
        res = _load_dat(sitedir, var)
    _CACHE[key] = res
    return res


def _load_nc(sitedir, var):
    ncp = Path(sitedir) / "output.nc"
    if not ncp.exists():
        return None
    from netCDF4 import Dataset
    d = Dataset(ncp)
    if var not in d.variables:
        return (None, None)
    t = np.asarray(d.variables["time"][:]) / 86400.0
    F = np.asarray(d.variables[var][:])
    return (t, F[:, 1:])


def _load_dat(sitedir, var):
    p = Path(sitedir) / f"{var}.dat"
    if not p.exists():
        return (None, None)
    try:
        import pandas as pd
        d = pd.read_csv(p, sep="\t", header=None).to_numpy()
    except Exception:
        d = np.genfromtxt(p, delimiter="\t")
    if d.ndim == 2 and np.isnan(d[:, -1]).all():
        d = d[:, :-1]
    return (d[:, 0] / 86400.0, d[:, 1:])


def _ylim(F):
    """Robust fixed y-limits from finite values (1–99th pct, padded)."""
    v = F[np.isfinite(F)]
    if v.size == 0:
        return (-1.0, 1.0)
    lo, hi = np.nanpercentile(v, 1), np.nanpercentile(v, 99)
    if hi <= lo:
        hi = lo + 1.0
    pad = 0.08 * (hi - lo)
    return (lo - pad, hi + pad)


def _frames(nt, target=240):
    """Even frame indices, capped so GIFs stay a sane size."""
    if nt <= target:
        return np.arange(nt)
    return np.linspace(0, nt - 1, target).astype(int)


def make_profile_movie(site, sitedir, writer, ext):
    # gather available panels
    data, ylims = {}, {}
    t_ref = None
    for var, lab, _ in PANELS:
        t, F = load(sitedir, var)
        if F is None:
            continue
        data[var] = (t, F)
        ylims[var] = _ylim(F)
        t_ref = t if t_ref is None else t_ref
    if not data:
        print(f"  {site}: no fields found, skipping")
        return
    have = [(v, l) for (v, l, _) in PANELS if v in data]
    ncol = 2
    nrow = int(np.ceil(len(have) / ncol))
    x = np.arange(data[have[0][0]][1].shape[1]) * DELXI_KM

    fig = plt.figure(figsize=(9.5, 2.0 * nrow + 0.8))
    fig.suptitle(f"{LABEL.get(site, site)} — along-channel evolution",
                 x=0.06, ha="left", fontsize=12, weight="bold", color=C.get(site, INK))
    daytxt = fig.text(0.94, 0.965, "", ha="right", fontsize=11, weight="bold", color=INK)
    gs = fig.add_gridspec(nrow, ncol, left=0.08, right=0.97, top=0.90, bottom=0.08,
                          hspace=0.55, wspace=0.24)
    obs = site_obs(site)
    lines, scats, mean_scats = {}, {}, {}
    for k, (var, lab) in enumerate(have):
        ax = fig.add_subplot(gs[k // ncol, k % ncol])
        col = ICEBLUE if var.startswith("ice") else C.get(site, INK)
        (ln,) = ax.plot([], [], color=col, lw=1.6, zorder=2)
        ax.set_xlim(0, x.max())
        # widen y-limits so overlaid observations are never clipped
        ylo, yhi = ylims[var]
        if var in obs and len(obs[var][0]):
            ov = obs[var][0][:, 1]
            ylo, yhi = min(ylo, np.nanmin(ov)), max(yhi, np.nanmax(ov))
            pad = 0.06 * (yhi - ylo + 1e-9)
            ylo, yhi = ylo - pad, yhi + pad
        ax.set_ylim(ylo, yhi)
        ax.set_ylabel(lab, fontsize=7.3)
        if k // ncol == nrow - 1:
            ax.set_xlabel("distance from mouth (km)", fontsize=7.3)
        for sp in ("top", "right"):
            ax.spines[sp].set_visible(False)
        ax.grid(True, alpha=0.5)
        ax.set_axisbelow(True)
        lines[var] = ln
        # observation overlay: individual grabs (open circles) that light up when the
        # movie day matches, plus a filled marker at their MEAN when 2+ coincide.
        if var in obs:
            sc = ax.scatter([], [], s=26, facecolor="none", edgecolor=INK,
                            linewidths=1.1, zorder=3, label=obs[var][2])
            mn = ax.scatter([], [], s=70, marker="D", facecolor=MEANCOL,
                            edgecolor="white", linewidths=0.8, zorder=4, label="obs mean")
            scats[var] = sc
            mean_scats[var] = mn
            ax.legend(loc="upper right", fontsize=5.6, labelcolor=INK2,
                      handletextpad=0.3, borderpad=0.2)

    # common frame axis from the first field's time vector
    t0 = data[have[0][0]][0]
    frames = _frames(len(t0))

    def update(fi):
        arts = []
        for var, _ in have:
            t, F = data[var]
            j = min(fi, F.shape[0] - 1)
            lines[var].set_data(x, F[j])
            arts.append(lines[var])
        day = np.mod(t0[min(fi, len(t0) - 1)], 365)
        # observations within +/- OBS_WINDOW days-of-year of the current frame
        for var, sc in scats.items():
            pts, xk, _ = obs[var]
            near = pts[np.abs(pts[:, 0] - day) <= OBS_WINDOW]
            sc.set_offsets(np.column_stack([np.full(len(near), xk), near[:, 1]])
                           if len(near) else np.empty((0, 2)))
            arts.append(sc)
            # mean marker when 2+ observations coincide at this point/time
            mn = mean_scats[var]
            mn.set_offsets([[xk, float(np.nanmean(near[:, 1]))]] if len(near) >= 2
                           else np.empty((0, 2)))
            arts.append(mn)
        daytxt.set_text(f"day {day:5.1f}")
        arts.append(daytxt)
        return arts

    anim = FuncAnimation(fig, update, frames=frames, blit=False)
    out = OUTDIR / f"{site}_evolution.{ext}"
    anim.save(out, writer=writer, dpi=110)
    plt.close(fig)
    print(f"  wrote {out.relative_to(ROOT)} ({out.stat().st_size/1e6:.1f} MB, {len(frames)} frames)")


def make_hovmoller_movie(site, sitedir, writer, ext):
    """Hovmoller (distance x time) of core fields with a sweeping day cursor."""
    picks = [(v, l, cm) for (v, l, cm) in PANELS
             if v in ("T", "S", "DIC", "FCO2", "ice_thickness") and load(sitedir, v)[1] is not None]
    if not picks:
        return
    nrow = len(picks)
    fig = plt.figure(figsize=(9.0, 1.5 * nrow + 0.8))
    fig.suptitle(f"{LABEL.get(site, site)} — distance × time, with day cursor",
                 x=0.06, ha="left", fontsize=12, weight="bold", color=C.get(site, INK))
    gs = fig.add_gridspec(nrow, 1, left=0.10, right=0.90, top=0.90, bottom=0.09, hspace=0.45)
    cursors = []
    for k, (var, lab, cm) in enumerate(picks):
        ax = fig.add_subplot(gs[k])
        t, F = load(sitedir, var)
        x = np.arange(F.shape[1]) * DELXI_KM
        step = max(1, len(t) // 730)
        im = ax.pcolormesh(t[::step], x, F[::step].T, cmap=cm, shading="auto")
        # y-axis is DISTANCE; the field value is the colour, so it labels the colorbar
        ax.set_ylabel("distance (km)", fontsize=7.0)
        cb = fig.colorbar(im, ax=ax, pad=0.02, fraction=0.045)
        cb.ax.tick_params(labelsize=6)
        cb.set_label(lab, fontsize=6.5)
        cur = ax.axvline(t[0], color=INK, lw=1.2)
        cursors.append((cur, t))
        if k == nrow - 1:
            ax.set_xlabel("day", fontsize=7.3)
    t0 = load(sitedir, picks[0][0])[0]
    frames = _frames(len(t0))

    def update(fi):
        day = t0[min(fi, len(t0) - 1)]
        for cur, _ in cursors:
            cur.set_xdata([day, day])
        return [c for c, _ in cursors]

    anim = FuncAnimation(fig, update, frames=frames, blit=False)
    out = OUTDIR / f"{site}_hovmoller.{ext}"
    anim.save(out, writer=writer, dpi=110)
    plt.close(fig)
    print(f"  wrote {out.relative_to(ROOT)} ({out.stat().st_size/1e6:.1f} MB)")


def main():
    args = [a for a in sys.argv[1:]]
    # Default: mp4 (H.264, yuv420p) -- QuickTime-native on macOS. --avi and --gif override.
    if "--gif" in args:
        writer, ext = PillowWriter(fps=12), "gif"
        args.remove("--gif")
    elif "--avi" in args:
        from matplotlib.animation import FFMpegWriter
        # MPEG-4 Part 2 in an AVI container. NOTE: macOS QuickTime cannot open AVI at
        # all -- use the default mp4 there; AVI is only for tools that require it.
        writer, ext = FFMpegWriter(
            fps=15, codec="mpeg4", extra_args=["-pix_fmt", "yuv420p", "-q:v", "4"]), "avi"
        args.remove("--avi")
    else:
        from matplotlib.animation import FFMpegWriter
        # H.264 + yuv420p + faststart: the broadly-compatible mp4 that QuickTime,
        # Preview, browsers and PowerPoint all play natively.
        writer, ext = FFMpegWriter(
            fps=15, codec="h264",
            extra_args=["-pix_fmt", "yuv420p", "-crf", "23", "-movflags", "+faststart"]), "mp4"
    rundir = None
    if "--run" in args:
        i = args.index("--run")
        rundir = Path(args[i + 1])
        del args[i:i + 2]
    sites = args if args else SITES

    OUTDIR.mkdir(parents=True, exist_ok=True)
    for s in sites:
        sitedir = rundir if rundir is not None else (RUN / s)
        print(f"{LABEL.get(s, s)}: {sitedir}")
        _CACHE.clear()
        make_profile_movie(s, sitedir, writer, ext)
        make_hovmoller_movie(s, sitedir, writer, ext)


if __name__ == "__main__":
    main()
