#!/usr/bin/env python3
"""38-point x 5-tracer Darwin forward-FD vs pkg/dic adjoint regression
plot, mirroring darwin_bling_comparison/scripts/plot_multipoint_comparison.py:
map of validation points (colored/shaped by DIC+ALK sensitivity category,
reusing the exact same categorization since it's a property of the
physical grid points, not of which adjoint model is being compared),
DIC+ALK headline panel (linear), all-5-tracers panel (symlog).
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

TRACERS = ['DIC', 'ALK', 'PO4', 'DOP', 'O2']
COLORS = {'DIC': '#0072B2', 'ALK': '#D55E00', 'PO4': '#E69F00', 'DOP': '#009E73', 'O2': '#56B4E9'}
ZORDER = {'DIC': 3, 'ALK': 3, 'PO4': 3, 'DOP': 3, 'O2': 10}  # O2 is all-zero -> draw last/on top

fig = plt.figure(figsize=(17, 5.5))
gs = gridspec.GridSpec(1, 3, figure=fig, width_ratios=[1.1, 1, 1])
ax0 = fig.add_subplot(gs[0], projection=ccrs.Robinson(central_longitude=0))
ax1 = fig.add_subplot(gs[1])
ax2 = fig.add_subplot(gs[2])

# --- Panel 0: map of the 38 validation point locations (same categorization as BLING comparison) ---
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
    y = np.array([float(r['dJ_dic']) for r in sub])
    ax1.scatter(x, y, s=36, color=COLORS[t], alpha=0.8, label=f'{t} (n={len(x)})', zorder=ZORDER[t])
    slope, intercept = np.polyfit(x, y, 1)
    corr = np.corrcoef(x, y)[0, 1]
    stats_text.append(f'{t}: slope={slope:.2f}, r={corr:.3f}')

all_x = np.concatenate([[float(r['dJ_darwin']) for r in rows if r['tracer'] == t] for t in ['DIC', 'ALK']])
all_y = np.concatenate([[float(r['dJ_dic']) for r in rows if r['tracer'] == t] for t in ['DIC', 'ALK']])
lims = [min(all_x.min(), all_y.min()) * 1.15, max(all_x.max(), all_y.max()) * 1.15]
ax1.plot(lims, lims, 'k--', linewidth=1, zorder=0, label='1:1')
ax1.axhline(0, color='0.85', linewidth=0.8, zorder=0)
ax1.axvline(0, color='0.85', linewidth=0.8, zorder=0)
ax1.set_xlim(lims); ax1.set_ylim(lims)
ax1.set_aspect('equal', adjustable='box')
ax1.set_xlabel(r'Darwin forward-FD  $dJ/d(\mathrm{tracer})$  (per mol m$^{-3}$)')
ax1.set_ylabel(r'pkg/dic adjoint  $dJ/d(\mathrm{tracer})$  (per mol m$^{-3}$)')
ax1.set_title('DIC + ALK: 38 global points\n(direct carbon-chemistry pathway)')
ax1.legend(loc='upper left', fontsize=9, frameon=False)
ax1.text(0.97, 0.03, '\n'.join(stats_text), transform=ax1.transAxes,
         ha='right', va='bottom', fontsize=9,
         bbox=dict(boxstyle='round', facecolor='white', alpha=0.8, edgecolor='0.7'))

# --- Panel 2: all 5 tracers, symlog scale, full context ---
# All 38 points per tracer are plotted and fitted. An earlier version excluded
# points flagged `abs(adxx) > 1e6` as "corrupted"; that flag was withdrawn --
# it applied one absolute cut to raw control-space fields whose per-tracer
# scale spans ~400x, so it removed 15 of 38 perfectly ordinary DOP points and
# inverted that regression. See README.md. Because OLS on these samples can be
# carried by a single high-leverage point, the leverage-immune Theil-Sen slope
# and Spearman rank r are annotated alongside it.
def _theil_sen(x, y):
    sl = [(y[b] - y[a]) / (x[b] - x[a])
          for a in range(len(x)) for b in range(a + 1, len(x)) if x[b] != x[a]]
    return float(np.median(sl)) if sl else float('nan')


def _spearman(x, y):
    rx = np.argsort(np.argsort(x)).astype(float)
    ry = np.argsort(np.argsort(y)).astype(float)
    return float(np.corrcoef(rx, ry)[0, 1])


stats_text2 = []
for t in TRACERS:
    sub = [r for r in rows if r['tracer'] == t]
    x = np.array([float(r['dJ_darwin']) for r in sub])
    y = np.array([float(r['dJ_dic']) for r in sub])
    ax2.scatter(x, y, s=26, color=COLORS[t], alpha=0.7, label=t, zorder=ZORDER[t])
    if np.ptp(x) < 1e-6 * max(abs(x).max(), 1) or np.ptp(y) == 0:
        stats_text2.append(f'{t}: sensitivity ~0 (both models)')
    else:
        slope = np.polyfit(x, y, 1)[0]
        corr = np.corrcoef(x, y)[0, 1]
        stats_text2.append(
            f'{t}: OLS {slope:.2f}/r={corr:.3f} | TS {_theil_sen(x, y):.2f}/'
            f'rank={_spearman(x, y):.2f} (n={len(sub)})')

lin_thresh = 1.0
lims2 = 1e8
ax2.set_xscale('symlog', linthresh=lin_thresh)
ax2.set_yscale('symlog', linthresh=lin_thresh)
ax2.plot([-lims2, lims2], [-lims2, lims2], 'k--', linewidth=1, zorder=0, label='1:1')
ax2.axhline(0, color='0.85', linewidth=0.8, zorder=0)
ax2.axvline(0, color='0.85', linewidth=0.8, zorder=0)
ax2.set_xlim(-lims2, lims2); ax2.set_ylim(-lims2, lims2)
ax2.set_xlabel(r'Darwin forward-FD  $dJ/d(\mathrm{tracer})$  (symlog)')
ax2.set_ylabel(r'pkg/dic adjoint  $dJ/d(\mathrm{tracer})$  (symlog)')
ax2.set_title('All 5 tracers, 38 points each\n(190 total point-tracer pairs)')
ax2.legend(loc='upper left', fontsize=8, frameon=False, ncol=2)
ax2.text(0.97, 0.03, '\n'.join(stats_text2), transform=ax2.transAxes,
         ha='right', va='bottom', fontsize=7.5,
         bbox=dict(boxstyle='round', facecolor='white', alpha=0.8, edgecolor='0.7'))

fig.suptitle('Darwin vs pkg/dic: forward-FD vs adjoint sensitivity of SOCAT cost function\n'
             '38 globally-distributed grid points, matched initial conditions, +/-5% relative perturbation',
             fontsize=11)
fig.tight_layout(rect=[0, 0, 1, 0.90])

OUT = f'{COMPARISON}/figures/darwin_vs_dic_multipoint.png'
fig.savefig(OUT, dpi=150)
print(f"saved to {OUT}")
