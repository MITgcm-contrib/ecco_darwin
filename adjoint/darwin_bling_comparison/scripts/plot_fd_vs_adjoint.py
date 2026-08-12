#!/usr/bin/env python3
"""Darwin forward-FD vs BLING adjoint sensitivity, dJ/d(tracer IC) at
grid point (j=40,i=97), after matching Darwin's ICs to BLING's own
(eliminating the ~7% ALK / ~2% DIC IC-provenance mismatch found earlier).
Mirrors the style of notes/from_kay/conc_fd_vs_adj.png (GCHP CO2 adjoint
validation): scatter of adjoint sensitivity vs FD sensitivity with a 1:1
reference line.
"""
import os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

OUT_DIR = os.environ.get("FIGURES_OUT_DIR", "../figures")

# dJ/dtracer, both sides converted to consistent per-mol/m3 (or per-eq/m3)
# units at grid point (j=40,i=97). Darwin: full-field cost-function FD,
# post-IC-match rerun. BLING: adxx_ptr* * sqrt(ctrl weight), the raw
# dJ/dtracer implied by MITgcm's genarr3d control-preconditioning.
DARWIN = {
    'DIC': 2.198565e+04,
    'ALK': -1.910626e+04,
    'NO3': -7.601870e+03,
    'O2':  -5.141986e+01,
    'PO4': -5.420134e+05,
    'FeT':  0.0,
}
BLING = {
    'DIC': 2.562473e+04,
    'ALK': -2.152067e+04,
    'NO3': -2.462582e+03,
    'O2':   1.382956e-01,
    'PO4':  2.361490e+04,
    'FeT': -3.050708e+06,
}

DIRECT = ['DIC', 'ALK']       # direct carbonate-chemistry pathway
INDIRECT = ['NO3', 'O2', 'PO4', 'FeT']  # biology-mediated pathway

# Okabe-Ito colorblind-safe categorical palette (Okabe & Ito 2008), fixed
# hue order, one fully-distinct color per tracer -- same assignment used
# in plot_multipoint_comparison.py for consistency across figures.
COLORS = {'DIC': '#0072B2', 'ALK': '#D55E00', 'O2': '#009E73',
          'NO3': '#CC79A7', 'PO4': '#E69F00', 'FeT': '#56B4E9'}

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 5))

# --- Panel 1: DIC/ALK, linear scale, the headline result ---
xs = [DARWIN[t] for t in DIRECT]
ys = [BLING[t] for t in DIRECT]
for t in DIRECT:
    ax1.scatter(DARWIN[t], BLING[t], s=90, color=COLORS[t], zorder=3, label=t)
    ax1.annotate(t, (DARWIN[t], BLING[t]), textcoords="offset points",
                 xytext=(8, 6), fontsize=11, fontweight='bold')

lims = [min(xs + ys) * 1.25, max(xs + ys) * 1.25]
ax1.plot(lims, lims, 'k--', linewidth=1, label='1:1', zorder=1)
ax1.axhline(0, color='0.8', linewidth=0.8, zorder=0)
ax1.axvline(0, color='0.8', linewidth=0.8, zorder=0)
ax1.set_xlim(lims)
ax1.set_ylim(lims)
ax1.set_xlabel(r'Darwin forward-FD  $dJ/d(\mathrm{tracer})$  (per mol m$^{-3}$)')
ax1.set_ylabel(r'BLING adjoint  $dJ/d(\mathrm{tracer})$  (per mol m$^{-3}$)')
ax1.set_title('Direct carbon-chemistry tracers\n(matched initial conditions)')
ax1.legend(loc='lower right', frameon=False)
ax1.set_aspect('equal', adjustable='box')

# --- Panel 2: all 6 tracers, symlog scale, full context ---
all_t = DIRECT + INDIRECT
xs2 = np.array([DARWIN[t] for t in all_t])
ys2 = np.array([BLING[t] for t in all_t])
colors = [COLORS[t] for t in all_t]
ax2.scatter(xs2, ys2, s=90, c=colors, zorder=3)
for t, x, y in zip(all_t, xs2, ys2):
    ax2.annotate(t, (x, y), textcoords="offset points", xytext=(8, 6), fontsize=10)

lin_thresh = 1.0
lims2 = 5e6
ax2.set_xscale('symlog', linthresh=lin_thresh)
ax2.set_yscale('symlog', linthresh=lin_thresh)
ax2.plot([-lims2, lims2], [-lims2, lims2], 'k--', linewidth=1, label='1:1', zorder=1)
ax2.axhline(0, color='0.8', linewidth=0.8, zorder=0)
ax2.axvline(0, color='0.8', linewidth=0.8, zorder=0)
ax2.set_xlim(-lims2, lims2)
ax2.set_ylim(-lims2, lims2)
ax2.set_xlabel(r'Darwin forward-FD  $dJ/d(\mathrm{tracer})$  (symlog)')
ax2.set_ylabel(r'BLING adjoint  $dJ/d(\mathrm{tracer})$  (symlog)')
ax2.set_title('All 6 overlapping tracers')
handles = [plt.Line2D([], [], marker='o', linestyle='', color=COLORS[t], label=t) for t in all_t]
ax2.legend(handles=handles, loc='lower right', frameon=False, fontsize=8, ncol=2)

fig.suptitle('Darwin vs BLING: forward-FD vs adjoint sensitivity of SOCAT cost function\n'
             'grid point (j=40, i=97), single-point IC perturbation, day-15 snapshot',
             fontsize=11)
fig.tight_layout(rect=[0, 0, 1, 0.93])

os.makedirs(OUT_DIR, exist_ok=True)
OUT = f'{OUT_DIR}/darwin_vs_bling_fd_vs_adjoint.png'
fig.savefig(OUT, dpi=150)
print(f"saved to {OUT}")
