#!/usr/bin/env python3
"""38-point x 6-tracer Darwin FD vs BLING adjoint regression plot, Kay-
style (slope/correlation reported per tracer, 1:1 reference line), with a
map of the validation point locations as the first panel.

Requires MITgcmutils (ships with MITgcm under utils/python/MITgcmutils --
add that directory to PYTHONPATH) and cartopy.

Configure via environment variables:
  DARWIN_RUN_DIR  -- Darwin experiment dir with run/, run_multi_*/,
                      multipoint_comparison_results.csv (from
                      compute_multipoint_comparison.py)
  FIGURES_OUT_DIR -- where to write the output PNG
"""
import csv
import json
import os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import numpy as np
from MITgcmutils import rdmds

DARWIN = os.environ.get("DARWIN_RUN_DIR", ".")
OUT_DIR = os.environ.get("FIGURES_OUT_DIR", "../figures")
rows = list(csv.DictReader(open(f'{DARWIN}/multipoint_comparison_results.csv')))

TRACERS = ['DIC', 'ALK', 'O2', 'NO3', 'PO4', 'FeT']
# Okabe-Ito colorblind-safe categorical palette (Okabe & Ito 2008), fixed
# hue order -- each tracer gets a fully distinct color instead of shades
# within a "blue family"/"red family" that are hard to tell apart.
COLORS = {'DIC': '#0072B2', 'ALK': '#D55E00', 'O2': '#009E73',
          'NO3': '#CC79A7', 'PO4': '#E69F00', 'FeT': '#56B4E9'}

fig = plt.figure(figsize=(17, 5.5))
gs = gridspec.GridSpec(1, 3, figure=fig, width_ratios=[1.1, 1, 1])
ax0 = fig.add_subplot(gs[0], projection=ccrs.Robinson(central_longitude=0))
ax1 = fig.add_subplot(gs[1])
ax2 = fig.add_subplot(gs[2])

# --- Panel 0: map of the 38 validation point locations ---
XC = rdmds(f'{DARWIN}/run/XC'); YC = rdmds(f'{DARWIN}/run/YC')
with open(f'{DARWIN}/run_multi_DIC_plus/perturb_log.json') as fh:
    plog = json.load(fh)
lons = np.array([XC[p['j'], p['i']] for p in plog['points']])
lats = np.array([YC[p['j'], p['i']] for p in plog['points']])
lons_180 = np.where(lons > 180, lons - 360, lons)

# "low sensitivity" sites: bottom quartile of |dJ_darwin/d(tracer)| -- these
# are the points where both models' actual carbon-chemistry sensitivity is
# small, so ratio-based agreement metrics get noisy (see readme.txt's
# outlier discussion) even though the absolute disagreement is trivial.
dic_mag = {(r['j'], r['i']): abs(float(r['dJ_darwin'])) for r in rows if r['tracer'] == 'DIC'}
alk_mag = {(r['j'], r['i']): abs(float(r['dJ_darwin'])) for r in rows if r['tracer'] == 'ALK'}
dic_thresh = np.percentile(list(dic_mag.values()), 25)
alk_thresh = np.percentile(list(alk_mag.values()), 25)
keys = [(str(p['j']), str(p['i'])) for p in plog['points']]
low_dic = np.array([dic_mag[k] < dic_thresh for k in keys])
low_alk = np.array([alk_mag[k] < alk_thresh for k in keys])
low_both = low_dic & low_alk
low_dic_only = low_dic & ~low_alk
low_alk_only = low_alk & ~low_dic
normal = ~low_dic & ~low_alk

ax0.set_global()
ax0.add_feature(cfeature.LAND, facecolor='0.85', zorder=1)
ax0.add_feature(cfeature.OCEAN, facecolor='#eaf3fb', zorder=0)
ax0.coastlines(linewidth=0.6, zorder=2)
ax0.gridlines(draw_labels=False, linewidth=0.3, color='0.7', alpha=0.6)
ax0.scatter(lons_180[normal], lats[normal], s=45, c='#b2182b', edgecolor='white',
            linewidth=0.6, zorder=3, transform=ccrs.PlateCarree(), label='normal sensitivity')
ax0.scatter(lons_180[low_dic_only], lats[low_dic_only], s=70, c='#ffee00', edgecolor='black',
            linewidth=0.9, marker='D', zorder=4, transform=ccrs.PlateCarree(),
            label='low DIC sensitivity only')
ax0.scatter(lons_180[low_alk_only], lats[low_alk_only], s=70, c='#66c2ff', edgecolor='black',
            linewidth=0.9, marker='s', zorder=4, transform=ccrs.PlateCarree(),
            label='low ALK sensitivity only')
ax0.scatter(lons_180[low_both], lats[low_both], s=90, c='#9b59b6', edgecolor='black',
            linewidth=0.9, marker='*', zorder=5, transform=ccrs.PlateCarree(),
            label='low DIC + ALK sensitivity')
ax0.legend(loc='lower left', fontsize=7, frameon=True, framealpha=0.9)
ax0.set_title(f'{len(plog["points"])} validation grid points\n(bottom quartile |dJ/d(tracer)| highlighted)', fontsize=11)

# --- Panel 1: DIC + ALK, linear scale, the headline result ---
stats_text = []
for t in ['DIC', 'ALK']:
    sub = [r for r in rows if r['tracer'] == t]
    x = np.array([float(r['dJ_darwin']) for r in sub])
    y = np.array([float(r['dJ_bling']) for r in sub])
    ax1.scatter(x, y, s=36, color=COLORS[t], alpha=0.8, label=f'{t} (n={len(x)})')
    slope, intercept = np.polyfit(x, y, 1)
    corr = np.corrcoef(x, y)[0, 1]
    stats_text.append(f'{t}: slope={slope:.2f}, r={corr:.3f}')

all_x = np.concatenate([[float(r['dJ_darwin']) for r in rows if r['tracer'] == t] for t in ['DIC', 'ALK']])
all_y = np.concatenate([[float(r['dJ_bling']) for r in rows if r['tracer'] == t] for t in ['DIC', 'ALK']])
lims = [min(all_x.min(), all_y.min()) * 1.15, max(all_x.max(), all_y.max()) * 1.15]
ax1.plot(lims, lims, 'k--', linewidth=1, zorder=0, label='1:1')
ax1.axhline(0, color='0.85', linewidth=0.8, zorder=0)
ax1.axvline(0, color='0.85', linewidth=0.8, zorder=0)
ax1.set_xlim(lims); ax1.set_ylim(lims)
ax1.set_aspect('equal', adjustable='box')
ax1.set_xlabel(r'Darwin forward-FD  $dJ/d(\mathrm{tracer})$  (per mol m$^{-3}$)')
ax1.set_ylabel(r'BLING adjoint  $dJ/d(\mathrm{tracer})$  (per mol m$^{-3}$)')
ax1.set_title('DIC + ALK: 38 global points\n(direct carbon-chemistry pathway)')
ax1.legend(loc='upper left', fontsize=9, frameon=False)
ax1.text(0.97, 0.03, '\n'.join(stats_text), transform=ax1.transAxes,
         ha='right', va='bottom', fontsize=9,
         bbox=dict(boxstyle='round', facecolor='white', alpha=0.8, edgecolor='0.7'))

# --- Panel 2: all 6 tracers, symlog scale, full context ---
for t in TRACERS:
    sub = [r for r in rows if r['tracer'] == t]
    x = np.array([float(r['dJ_darwin']) for r in sub])
    y = np.array([float(r['dJ_bling']) for r in sub])
    ax2.scatter(x, y, s=26, color=COLORS[t], alpha=0.7, label=t)

lin_thresh = 1.0
lims2 = 2e6
ax2.set_xscale('symlog', linthresh=lin_thresh)
ax2.set_yscale('symlog', linthresh=lin_thresh)
ax2.plot([-lims2, lims2], [-lims2, lims2], 'k--', linewidth=1, zorder=0, label='1:1')
ax2.axhline(0, color='0.85', linewidth=0.8, zorder=0)
ax2.axvline(0, color='0.85', linewidth=0.8, zorder=0)
ax2.set_xlim(-lims2, lims2); ax2.set_ylim(-lims2, lims2)
ax2.set_xlabel(r'Darwin forward-FD  $dJ/d(\mathrm{tracer})$  (symlog)')
ax2.set_ylabel(r'BLING adjoint  $dJ/d(\mathrm{tracer})$  (symlog)')
ax2.set_title('All 6 tracers, 38 points each\n(228 total point-tracer pairs)')
ax2.legend(loc='lower right', fontsize=8, frameon=False, ncol=2)

fig.suptitle('Darwin vs BLING: forward-FD vs adjoint sensitivity of SOCAT cost function\n'
             '38 globally-distributed grid points, matched initial conditions, +/-5% relative perturbation',
             fontsize=11)
fig.tight_layout(rect=[0, 0, 1, 0.90])

os.makedirs(OUT_DIR, exist_ok=True)
OUT = f'{OUT_DIR}/darwin_vs_bling_multipoint.png'
fig.savefig(OUT, dpi=150)
print(f"saved to {OUT}")
