#!/usr/bin/env python3
"""Darwin FD vs pkg/dic adjoint sensitivity, PO4 and DOP only -- pulled
out separately from plot_multipoint_comparison_dic.py's all-5-tracers
panel since their axes are hard to read squeezed in with DIC/ALK's much
larger range, mirroring darwin_bling_comparison's
plot_multipoint_nutrients.py split-out.

Both tracers' clean-point ranges are modest (a few thousand at most), so
linear axes are used here (unlike the BLING nutrients plot, which needed
symlog for FeT's much larger range).

Corrupted points (see README.md's "A correction" section -- a genuine,
field-specific pkg/dic adjoint instability, confirmed to persist even
when each tracer is the sole registered control, i.e. NOT the
slot-position bug documented in hybrid_darwin_dic/run_dic/data.ctrl) are
shown as grey X's and excluded from the regression, not silently dropped.
"""
import csv
import json
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import numpy as np
import sys

sys.path.insert(0, '/Users/carrolld/Documents/research/ECCO/offline/scripts')
from MITgcmutils import rdmds

DARWIN = '/Users/carrolld/Documents/research/adjoint/MITgcm/verification/global_oce_biogeo_darwin'
COMPARISON = '/Users/carrolld/Documents/research/adjoint/MITgcm/verification/darwin_dic_comparison'
rows = list(csv.DictReader(open(f'{COMPARISON}/results/multipoint_comparison_results_dic.csv')))

TRACERS = ['PO4', 'DOP']
COLORS = {'PO4': '#E69F00', 'DOP': '#009E73'}

fig = plt.figure(figsize=(15, 5.5))
gs = gridspec.GridSpec(1, 3, figure=fig, width_ratios=[1.1, 1, 1])
ax0 = fig.add_subplot(gs[0], projection=ccrs.Robinson(central_longitude=0))
ax1 = fig.add_subplot(gs[1])
ax2 = fig.add_subplot(gs[2])

# --- Panel 0: map of the 38 validation point locations (same categorization as the main figure) ---
XC = rdmds(f'{DARWIN}/run/XC'); YC = rdmds(f'{DARWIN}/run/YC')
with open(f'{DARWIN}/run_multi_DIC_plus/perturb_log.json') as fh:
    plog = json.load(fh)
lons = np.array([XC[p['j'], p['i']] for p in plog['points']])
lats = np.array([YC[p['j'], p['i']] for p in plog['points']])
lons_180 = np.where(lons > 180, lons - 360, lons)

bling_rows = list(csv.DictReader(open(f'{DARWIN}/multipoint_comparison_results.csv')))
dic_mag = {(r['j'], r['i']): abs(float(r['dJ_darwin'])) for r in bling_rows if r['tracer'] == 'DIC'}
alk_mag = {(r['j'], r['i']): abs(float(r['dJ_darwin'])) for r in bling_rows if r['tracer'] == 'ALK'}
dic_thresh = np.percentile(list(dic_mag.values()), 25)
alk_thresh = np.percentile(list(alk_mag.values()), 25)
keys = [(str(p['j']), str(p['i'])) for p in plog['points']]
low_dic = np.array([dic_mag[k] < dic_thresh for k in keys])
low_alk = np.array([alk_mag[k] < alk_thresh for k in keys])
low_both = low_dic & low_alk
low_dic_only = low_dic & ~low_alk
low_alk_only = low_alk & ~low_dic
normal = ~low_dic & ~low_alk

# corrupted-point flags, keyed by (j,i), shared across both tracers for the map overlay
# The 'corrupted' flag was withdrawn (see README.md); no points are excluded.
is_corrupted = np.zeros(len(plog['points']), dtype=bool)

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
ax0.scatter(lons_180[is_corrupted], lats[is_corrupted], s=110, facecolors='none',
            edgecolors='0.2', linewidth=1.5, marker='o', zorder=6,
            transform=ccrs.PlateCarree(), label=None)
ax0.legend(loc='lower left', fontsize=6.5, frameon=True, framealpha=0.9)
ax0.set_title(f'{len(plog["points"])} validation grid points', fontsize=11)

# --- Panels 1-2: PO4 and DOP, linear scale, one tracer per panel ---
for ax, t in zip([ax1, ax2], TRACERS):
    sub_all = [r for r in rows if r['tracer'] == t]
    sub = sub_all
    corrupted = []

    x = np.array([float(r['dJ_darwin']) for r in sub])
    y = np.array([float(r['dJ_dic']) for r in sub])
    ax.scatter(x, y, s=42, color=COLORS[t], alpha=0.85, zorder=3)

    if corrupted:
        cx = [float(r['dJ_darwin']) for r in corrupted]
        # corrupted dJ_dic values are enormous (~1e7-1e8); mark where they sit on
        # the x-axis (Darwin's own FD value is unaffected by the pkg/dic bug)
        # rather than forcing the axes out to their off-scale y-value.
        for xi in cx:
            ax.axvline(xi, color='0.75', linewidth=0.9, linestyle=':', zorder=1)
        ax.text(0.03, 0.97, f'{len(corrupted)} corrupted point(s) excluded\n(off-scale, dJ/d{t} ~1e7-1e8)',
                transform=ax.transAxes, ha='left', va='top', fontsize=7.5, color='0.4',
                bbox=dict(boxstyle='round', facecolor='white', alpha=0.75, edgecolor='0.7'))

    slope, intercept = np.polyfit(x, y, 1)
    corr = np.corrcoef(x, y)[0, 1]

    # x and y are on very different physical scales for these biology-mediated
    # tracers (esp. DOP), so -- unlike the DIC+ALK panel -- axes are scaled
    # independently rather than forced square/shared, or the smaller-range
    # axis collapses to a sliver.
    xpad = 0.15 * (x.max() - x.min())
    ypad = 0.15 * (y.max() - y.min())
    ax.set_xlim(x.min() - xpad, x.max() + xpad)
    ax.set_ylim(y.min() - ypad, y.max() + ypad)
    ax.axhline(0, color='0.85', linewidth=0.8, zorder=0)
    ax.axvline(0, color='0.85', linewidth=0.8, zorder=0)
    ax.set_xlabel(r'Darwin forward-FD  $dJ/d(\mathrm{tracer})$  (per mol m$^{-3}$)')
    ax.set_ylabel(r'pkg/dic adjoint  $dJ/d(\mathrm{tracer})$  (per mol m$^{-3}$)')
    ax.set_title(f'{t}: {len(sub)} clean global points\n(biology-mediated pathway)')
    ax.text(0.97, 0.03, f'slope={slope:.2f}, r={corr:.3f}', transform=ax.transAxes,
            ha='right', va='bottom', fontsize=9,
            bbox=dict(boxstyle='round', facecolor='white', alpha=0.8, edgecolor='0.7'))

fig.suptitle('Darwin vs pkg/dic: forward-FD vs adjoint sensitivity of SOCAT cost function\n'
             'PO4 and DOP only, all 38 points (see README.md)',
             fontsize=11)
fig.tight_layout(rect=[0, 0, 1, 0.90])

OUT = f'{COMPARISON}/figures/darwin_vs_dic_nutrients.png'
fig.savefig(OUT, dpi=150)
print(f"saved to {OUT}")
