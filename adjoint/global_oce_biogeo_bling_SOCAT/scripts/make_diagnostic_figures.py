# NOTE: RUN/BATHY/OUTDIR below point at the paths used to generate the
# reference figures/*.png shipped with this experiment. Point RUN at your
# own run/grdchk output directory and re-run to regenerate from a fresh
# gradient check. The hardcoded `gradcheck` dict and monthly `fc_contrib`
# values below reflect this specific reference campaign -- for a fresh
# run, pull the adj/fd pairs from your own grdchk log and the monthly
# cost breakdown from output_optim_000.txt instead of hardcoding them.
import os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import cartopy.crs as ccrs

RUN = os.environ.get("BLING_RUN_OPT_DIR", "../run_opt")
BATHY = os.environ.get("BATHY_BIN", "../../tutorial_global_oce_biogeo/input/bathy.bin")
OUTDIR = os.environ.get("FIGURES_OUT_DIR", "../figures")

nx, ny, nz = 128, 64, 15
sNx, sNy = 64, 32
delX, delY = 2.8125, 2.8125
mlon = (np.arange(nx) + 0.5) * delX
mlat = -90.0 + (np.arange(ny) + 0.5) * delY

bathy = np.fromfile(BATHY, dtype='>f4').reshape(ny, nx)
ocean_mask = bathy < 0  # True = ocean

# ---------- helpers ----------

def load_tiled_surface(varname, iterstr="0000000000"):
    """Assemble the 4 (sNx x sNy x nz) tiles into a global (ny,nx) surface (k=0) field."""
    full = np.full((ny, nx), np.nan, dtype='f8')
    tile_map = {(1, 1): (0, 0), (1, 2): (0, 1), (2, 1): (1, 0), (2, 2): (1, 1)}
    for (xt, yt), (bi, bj) in tile_map.items():
        fname = f"{RUN}/{varname}.{iterstr}.{xt:03d}.{yt:03d}.data"
        a = np.fromfile(fname, dtype='>f4').astype('f8').reshape(nz, sNy, sNx)
        surf = a[0]  # k=1 (surface)
        y0, y1 = bj * sNy, (bj + 1) * sNy
        x0, x1 = bi * sNx, (bi + 1) * sNx
        full[y0:y1, x0:x1] = surf
    full[~ocean_mask] = np.nan
    return full

# ============================================================
# Figure 1: Gradient check validation (DIC + ALK)
# ============================================================
gradcheck = {
    'DIC': {
        'points': ['(34,10,1)', '(35,10,1)'],
        'adj': [-459.09582519531, -332.98162841797],
        'fd':  [-459.09637192381, -332.98208691122],
        'relerr': [1.19088099604e-6, 1.37693256197e-6],
        'rms': 1.2872724705626e-6,
    },
    'ALK': {
        'points': ['(34,10,1)', '(35,10,1)'],
        'adj': [407.41873168945, 295.65786743164],
        'fd':  [407.41926095507, 295.65834338428],
        'relerr': [1.2990703963656e-6, 1.6098087924910e-6],
        'rms': 1.4627146411885e-6,
    },
}

fig, axes = plt.subplots(1, 2, figsize=(11, 4.8))

colors = {'DIC': '#3b6fa0', 'ALK': '#c2703d'}

ax = axes[0]
all_vals = []
for var, d in gradcheck.items():
    all_vals += d['adj'] + d['fd']
lims = [min(all_vals) * 1.05, max(all_vals) * 1.05]
ax.plot(lims, lims, color='0.6', lw=1, ls='--', zorder=1, label='1:1')
for var, d in gradcheck.items():
    ax.scatter(d['adj'], d['fd'], s=70, color=colors[var], label=var, zorder=3,
               edgecolor='white', linewidth=0.6)
ax.set_xlabel('Adjoint gradient  (∂fc / ∂xx, TAF)')
ax.set_ylabel('Finite-difference gradient  (centered, ε=1e-5)')
ax.set_title('Adjoint vs. finite-difference gradient')
ax.legend(frameon=False)
ax.spines[['top', 'right']].set_visible(False)

ax = axes[1]
bar_labels = []
bar_vals = []
bar_colors = []
for var, d in gradcheck.items():
    for i, pt in enumerate(d['points']):
        bar_labels.append(f"{var}\n{pt}")
        bar_vals.append(abs(d['relerr'][i]))
        bar_colors.append(colors[var])
xpos = np.arange(len(bar_labels))
ax.bar(xpos, bar_vals, color=bar_colors, width=0.6)
for var, d in gradcheck.items():
    pass
ax.axhline(gradcheck['DIC']['rms'], color=colors['DIC'], lw=1, ls=':', alpha=0.7)
ax.axhline(gradcheck['ALK']['rms'], color=colors['ALK'], lw=1, ls=':', alpha=0.7)
ax.set_yscale('log')
ax.set_xticks(xpos)
ax.set_xticklabels(bar_labels, fontsize=9)
ax.set_ylabel('|1 − fd_grad / adj_grad|  (relative error)')
ax.set_title('Gradient check relative error\n(RMS: DIC=1.29e-6, ALK=1.46e-6)')
ax.spines[['top', 'right']].set_visible(False)
ax.set_ylim(1e-7, 1e-5)

fig.suptitle('Adjoint gradient check — pkg/grdchk, 1-month test config', fontsize=12, y=1.02)
fig.tight_layout()
fig.savefig(f"{OUTDIR}/gradient_check_validation.png", dpi=160, bbox_inches='tight')
plt.close(fig)
print("wrote gradient_check_validation.png")

# ============================================================
# Figure 2: Monthly cost function breakdown (1-year run)
# ============================================================
months = ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec']
fc_contrib = [35041.3311664302, 23205.7783327819, 19837.5719235390, 14098.2944822456,
              14842.3173166944, 16322.0065689015, 16780.8977505925, 15370.0496288447,
              13844.7958051219, 15582.3170741114, 12582.1511987328, 14865.1914731840]
nobs = [2759, 2733, 2695, 2378, 2101, 1931, 2050, 1978, 2012, 2334, 2344, 2564]
fc_total = sum(fc_contrib)

fig, ax = plt.subplots(figsize=(9.5, 4.8))
xpos = np.arange(12)
bar_color = '#3b6fa0'
bars = ax.bar(xpos, fc_contrib, color=bar_color, width=0.65)
ax.set_xticks(xpos)
ax.set_xticklabels(months)
ax.set_ylabel('prof_PCO cost contribution')
ax.set_title(f'Monthly SOCAT pCO2 cost function contribution — 1-year run\ntotal fc = {fc_total:,.1f}  '
             f'(sum of 12 monthly terms, optimcycle=0 / first-guess control vector)')
ax.spines[['top', 'right']].set_visible(False)

ax2 = ax.twinx()
ax2.plot(xpos, nobs, color='#c2703d', marker='o', markersize=5, lw=1.5, label='N observations')
ax2.set_ylabel('N observations (per month)', color='#c2703d')
ax2.tick_params(axis='y', colors='#c2703d')
ax2.spines['top'].set_visible(False)
ax2.spines['right'].set_color('#c2703d')

fig.tight_layout()
fig.savefig(f"{OUTDIR}/monthly_cost_breakdown_1yr.png", dpi=160, bbox_inches='tight')
plt.close(fig)
print("wrote monthly_cost_breakdown_1yr.png")

# ============================================================
# Figure 3: Spatial maps of surface adjoint sensitivity
# ============================================================
fields = [
    ('adxx_theta', 'd(fc) / d(theta_IC)  [surface, K$^{-1}$]'),
    ('adxx_ptr1',  'd(fc) / d(DIC_IC)  [surface]'),
    ('adxx_ptr6',  'd(fc) / d(Fe_IC)  [surface]'),
]

fig = plt.figure(figsize=(16, 4.0))
gs = fig.add_gridspec(1, 3, wspace=0.08, top=0.95, bottom=0.02)
axes = [fig.add_subplot(gs[0, i], projection=ccrs.Robinson()) for i in range(3)]

for ax, (varname, title) in zip(axes, fields):
    data = load_tiled_surface(varname)
    # robust color range: 95th percentile keeps a handful of extreme coastal
    # grid cells from saturating the scale and hiding the broad-scale pattern
    vmax = np.nanpercentile(np.abs(data), 95)
    if not np.isfinite(vmax) or vmax == 0:
        vmax = np.nanmax(np.abs(data)) or 1.0
    mesh = ax.pcolormesh(mlon, mlat, data, transform=ccrs.PlateCarree(),
                          cmap='RdBu_r', vmin=-vmax, vmax=vmax, shading='auto')
    ax.coastlines(linewidth=0.5, color='0.3')
    ax.set_title(title, fontsize=10)
    cb = fig.colorbar(mesh, ax=ax, orientation='horizontal', pad=0.05, shrink=0.85, extend='both')
    cb.ax.tick_params(labelsize=8)

fig.suptitle('Surface adjoint sensitivity of the 1-year pCO2 cost function to initial conditions\n'
             '(color range clipped to the 95th percentile per panel to reveal broad-scale structure)',
             fontsize=12, y=1.14)
fig.savefig(f"{OUTDIR}/adjoint_sensitivity_maps.png", dpi=160, bbox_inches='tight')
plt.close(fig)
print("wrote adjoint_sensitivity_maps.png")
