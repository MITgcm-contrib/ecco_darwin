"""Diagnostic figure: area-weighted mean vertical profile for every one of
the 36 PTRACERS initial conditions, on the real ECCOv4r6 50-level LLC90
grid -- a check that the regrid_to_llc90.py vertical interpolation (23
Steph levels -> 50 ECCOv4r6 levels, see that module's _interp_vertical)
produced sensible profiles, not a check of horizontal structure (that's
plot_ic_surface_fields.py's job).

Reuses plot_ic_surface_fields.TRACERS for the tracer metadata (source
file/kind per tracer) rather than duplicating it -- see that module's own
header for why each tracer is "file" or "derived_unknown".

Weighting: RAC (cell area) * hFacC (wet fraction at each level), NOT
volume (no DRF factor) -- profiles are in the tracer's own per-volume
units, so an area-weighted mean at each level is the right average; a
volume weight would double-count the level thickness already implicit in
delR's variation with depth.
"""
import os
import numpy as np
import MITgcmutils as mu
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import regrid_to_llc90 as r
from plot_ic_surface_fields import TRACERS, IC_DIR

OUT_PNG = "figures/IC_vertical_profiles.png"
REAL_COLOR = "#1b6a8c"   # single hue for real (file-based) profiles
ZERO_COLOR = "#9a9a9a"   # muted gray for zero / derived-zero placeholders


def area_weights(grid_dir=r.GRID_DIR):
    """RAC (cell area, tile-indexed) broadcast against hFacC's per-level
    wet fraction -> (50,13,90,90) area weight, zero on dry cells."""
    import ecco_v4_py as ecco
    rac = ecco.llc_compact_to_tiles(mu.rdmds(grid_dir + "/RAC"), less_output=True)  # (13,90,90)
    _, _, hfacc = r.target_grid(grid_dir)  # (50,13,90,90)
    return rac[None, :, :, :] * hfacc


def load_profile(filename, weights):
    """Area-weighted mean at each of the 50 levels for one IC file."""
    import ecco_v4_py as ecco
    raw = np.fromfile(os.path.join(IC_DIR, filename), dtype=">f4").reshape(50, 1170, 90)
    field = ecco.llc_compact_to_tiles(raw, less_output=True)  # (50,13,90,90)
    wsum = weights.sum(axis=(1, 2, 3))
    fsum = (field * weights).sum(axis=(1, 2, 3))
    return np.where(wsum > 0, fsum / np.where(wsum > 0, wsum, 1.0), np.nan)


def main():
    os.makedirs(os.path.dirname(OUT_PNG), exist_ok=True)
    weights = area_weights()
    depth = r.target_rc()  # (50,) cell-center depths, positive down
    neg_depth = -depth     # plot on a plain (non-inverted) y-axis: surface
                            # (0) at top, deep (negative) at bottom. Avoids
                            # ax.invert_yaxis() entirely -- calling it once
                            # per subplot on a sharey=True axis group toggles
                            # the shared inversion state each time, so an
                            # even number of calls (36, here) nets out to no
                            # inversion at all, which is what actually
                            # produced the upside-down first version.

    fig, axes = plt.subplots(6, 6, figsize=(20, 22), sharey=True)
    fig.suptitle(
        "PTRACERS initial conditions -- area-weighted mean vertical profile "
        "(all 36, real ECCOv4r6 50-level LLC90 grid)",
        fontsize=15, y=0.995,
    )

    cache = {}  # filename -> profile, since c07-c10 share one file
    for ax, (idx, name, units, kind, filename, source, note) in zip(axes.flat, TRACERS):
        if kind == "file":
            if filename not in cache:
                cache[filename] = load_profile(filename, weights)
            profile = cache[filename]
            color, ls = REAL_COLOR, "-"
        else:
            profile = np.zeros_like(depth)
            color, ls = ZERO_COLOR, "--"

        ax.plot(profile, neg_depth, color=color, linestyle=ls, linewidth=2)
        ax.set_title(f"{idx}: {name}", fontsize=9)
        ax.set_xlabel(units, fontsize=7)
        ax.yaxis.set_major_formatter(lambda y, _: f"{abs(y):.0f}")
        ax.tick_params(labelsize=7)
        ax.axvline(0, color="0.85", linewidth=0.8, zorder=0)
        finite = profile[np.isfinite(profile)]
        print(
            f"{idx:2d} {name:16s} [{kind:13s}] "
            f"surface={profile[0]:.4g} bottom={profile[-1]:.4g} "
            f"min={finite.min():.4g} max={finite.max():.4g}"
        )

    axes[0, 0].set_ylabel("Depth [m]", fontsize=9)
    for row in range(6):
        axes[row, 0].set_ylabel("Depth [m]", fontsize=9)

    fig.text(
        0.5, 0.005,
        f"Solid ({REAL_COLOR}): real regridded field, now including c01-c10 (all 10\n"
        "plankton tracers share the same real biomass restart, reverted).\n"
        "Dashed gray: Chl01-06, a placeholder zero -- the real cold-start value (genuine\n"
        "acclimation from real biomass) isn't computable here. See\n"
        "plot_ic_surface_fields.py and data.ptracers for source detail.",
        ha="center", fontsize=8,
    )
    fig.tight_layout(rect=(0, 0.015, 1, 0.98))
    fig.savefig(OUT_PNG, dpi=150)
    print(f"\nDone: {OUT_PNG}")


if __name__ == "__main__":
    main()
