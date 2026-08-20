#!/usr/bin/env python3
"""Darwin FD vs BLING adjoint sensitivity, nutrients + iron only (O2, NO3,
PO4, FeT) -- companion to plot_multipoint_comparison.py's DIC+ALK headline
panel, pulled out separately since these four tracers span ~9 orders of
magnitude and would otherwise wash out on the same axes.

Map panel and marker-shape-by-category convention are identical to
plot_multipoint_comparison.py (circle=normal, diamond=low DIC sensitivity
only, square=low ALK sensitivity only, star=low DIC+ALK sensitivity) --
see that script for the category derivation.

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

TRACERS = ['O2', 'NO3', 'PO4', 'FeT']
COLORS = {'O2': '#009E73', 'NO3': '#CC79A7', 'PO4': '#E69F00', 'FeT': '#56B4E9'}

fig = plt.figure(figsize=(12, 5.5))
gs = gridspec.GridSpec(1, 2, figure=fig, width_ratios=[1.1, 1])
ax0 = fig.add_subplot(gs[0], projection=ccrs.Robinson(central_longitude=0))
ax1 = fig.add_subplot(gs[1])

# --- Panel 0: map of the 38 validation point locations (same as the main figure) ---
XC = rdmds(f'{DARWIN}/run/XC'); YC = rdmds(f'{DARWIN}/run/YC')
with open(f'{DARWIN}/run_multi_DIC_plus/perturb_log.json') as fh:
    plog = json.load(fh)
lons = np.array([XC[p['j'], p['i']] for p in plog['points']])
lats = np.array([YC[p['j'], p['i']] for p in plog['points']])
lons_180 = np.where(lons > 180, lons - 360, lons)

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

CAT_STYLE = {
    'normal':        dict(marker='o', s=36, edgecolors='none', linewidths=0),
    'low_dic_only':  dict(marker='D', s=55, edgecolors='black', linewidths=0.7),
    'low_alk_only':  dict(marker='s', s=55, edgecolors='black', linewidths=0.7),
    'low_both':      dict(marker='*', s=90, edgecolors='black', linewidths=0.7),
}
cat_lookup = {}
for idx, k in enumerate(keys):
    if low_both[idx]:
        cat_lookup[k] = 'low_both'
    elif low_dic_only[idx]:
        cat_lookup[k] = 'low_dic_only'
    elif low_alk_only[idx]:
        cat_lookup[k] = 'low_alk_only'
    else:
        cat_lookup[k] = 'normal'


def scatter_by_category(ax, sub, color, base_label, zorder_base=2):
    x = np.array([float(r['dJ_darwin']) for r in sub])
    y = np.array([float(r['dJ_bling']) for r in sub])
    cats = np.array([cat_lookup[(r['j'], r['i'])] for r in sub])
    labeled = False
    for cat, style in CAT_STYLE.items():
        mask = cats == cat
        if not mask.any():
            continue
        label = base_label if not labeled else None
        ax.scatter(x[mask], y[mask], color=color, alpha=0.8, label=label,
                   zorder=zorder_base, **style)
        labeled = True
    return x, y


# --- Panel 1: O2, NO3, PO4, FeT, symlog scale ---
# O2's dJ values are ~9 orders of magnitude smaller than FeT's, so its points
# cluster at the origin under everything else -- draw it last / on top
# (highest zorder) so it stays visible instead of being buried.
ZORDER = {'O2': 10, 'NO3': 3, 'PO4': 4, 'FeT': 2}
stats_text = []
for t in TRACERS:
    sub = [r for r in rows if r['tracer'] == t]
    x, y = scatter_by_category(ax1, sub, COLORS[t], f'{t} (n={len(sub)})', zorder_base=ZORDER[t])
    if np.ptp(x) < 1e-6 * max(abs(x).max(), 1):
        stats_text.append(f'{t}: sensitivity ~0 (both models), slope/r undefined')
    else:
        slope, intercept = np.polyfit(x, y, 1)
        corr = np.corrcoef(x, y)[0, 1]
        stats_text.append(f'{t}: slope={slope:.2f}, r={corr:.3f}')

lin_thresh = 1.0
lims = 6e7
ax1.set_xscale('symlog', linthresh=lin_thresh)
ax1.set_yscale('symlog', linthresh=lin_thresh)
ax1.plot([-lims, lims], [-lims, lims], 'k--', linewidth=1, zorder=0, label='1:1')
ax1.axhline(0, color='0.85', linewidth=0.8, zorder=0)
ax1.axvline(0, color='0.85', linewidth=0.8, zorder=0)
ax1.set_xlim(-lims, lims); ax1.set_ylim(-lims, lims)
ax1.set_xlabel(r'Darwin forward-FD  $dJ/d(\mathrm{tracer})$  (symlog)')
ax1.set_ylabel(r'BLING adjoint  $dJ/d(\mathrm{tracer})$  (symlog)')
ax1.set_title('Nutrients + iron: 38 global points\n(O2, NO3, PO4, FeT)')
ax1.legend(loc='upper left', fontsize=8, frameon=False)
ax1.text(0.97, 0.03, '\n'.join(stats_text), transform=ax1.transAxes,
         ha='right', va='bottom', fontsize=8,
         bbox=dict(boxstyle='round', facecolor='white', alpha=0.8, edgecolor='0.7'))

fig.suptitle('Darwin vs BLING: forward-FD vs adjoint sensitivity of SOCAT cost function\n'
             'nutrients + iron only, 38 globally-distributed grid points, +/-5% relative perturbation',
             fontsize=11)
fig.text(0.5, 0.01, 'marker shape matches the sensitivity category shown on the map (circle/diamond/square/star)',
          ha='center', fontsize=8, style='italic', color='0.3')
fig.tight_layout(rect=[0, 0.02, 1, 0.90])

os.makedirs(OUT_DIR, exist_ok=True)
OUT = f'{OUT_DIR}/darwin_vs_bling_nutrients.png'
fig.savefig(OUT, dpi=150)
print(f"saved to {OUT}")
