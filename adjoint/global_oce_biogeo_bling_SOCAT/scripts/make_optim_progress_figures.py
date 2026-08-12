# NOTE: RUN/BATHY/OUTDIR below point at the paths used to generate the
# reference figures/*.png shipped with this experiment. Point RUN at your
# own run directory and re-run to regenerate from a fresh optimization.
# The hardcoded `fc` array in the convergence figure below reflects this
# specific reference campaign (cycles 0-10); for a fresh run, parse fc
# values out of your own run's stdout ("global fc = ...") instead of
# hardcoding them -- see readme.txt, "Reproducing the figures."
import os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
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
ocean_mask = bathy < 0

def load_tiled_surface(varname, iterstr):
    full = np.full((ny, nx), np.nan, dtype='f8')
    tile_map = {(1, 1): (0, 0), (1, 2): (0, 1), (2, 1): (1, 0), (2, 2): (1, 1)}
    for (xt, yt), (bi, bj) in tile_map.items():
        fname = f"{RUN}/{varname}.{iterstr}.{xt:03d}.{yt:03d}.data"
        a = np.fromfile(fname, dtype='>f4').astype('f8').reshape(nz, sNy, sNx)
        surf = a[0]
        y0, y1 = bj * sNy, (bj + 1) * sNy
        x0, x1 = bi * sNx, (bi + 1) * sNx
        full[y0:y1, x0:x1] = surf
    full[~ocean_mask] = np.nan
    return full

# ============================================================
# Figure 1: Cost function convergence (leg 1, cycles 0-10)
# ============================================================
fc = [212372.702721180, 205231.624625958, 204090.983013349, 203085.762665770,
      202190.714912161, 201357.618549722, 200503.021795958, 199551.351134046,
      198944.783752850, 198619.084252134, 198346.427310668]
cycles = np.arange(len(fc))
fc = np.array(fc)
pct_cum = 100 * (fc - fc[0]) / fc[0]
steps = -np.diff(fc)

fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(8.5, 7), sharex=True,
                                gridspec_kw={'height_ratios': [2, 1], 'hspace': 0.08})

ax1.plot(cycles, fc, color='#3b6fa0', marker='o', markersize=6, lw=2, zorder=3)
ax1.set_ylabel('Cost function  fc')
ax1.set_title('M1QN3 offline optimization — cost function trajectory\n'
              f'{fc[0]:,.0f} → {fc[-1]:,.0f}  ({pct_cum[-1]:.2f}% cumulative decrease)')
ax1.spines[['top', 'right']].set_visible(False)
for x, y, p in zip(cycles, fc, pct_cum):
    if x in (0, len(cycles) - 1):
        ax1.annotate(f'{p:+.2f}%' if x > 0 else 'baseline',
                     (x, y), textcoords='offset points', xytext=(0, 10),
                     ha='center', fontsize=9, color='#3b6fa0')
ax1.axvline(10, color='0.6', lw=1, ls='--', zorder=1)
ax1.annotate('m1qn3 hit its iteration cap\n(omode=4) at cycle 10 —\n"leg 2" restarts with a cold Hessian',
             xy=(10, fc[-1]), xycoords='data',
             xytext=(0.62, 0.55), textcoords='axes fraction',
             fontsize=8.5, color='0.35', ha='left', va='center',
             arrowprops=dict(arrowstyle='-', color='0.6', lw=0.8))

ax2.bar(cycles[1:], steps, color='#3b6fa0', width=0.6)
ax2.set_yscale('log')
ax2.set_xlabel('Optimization cycle (model evaluation)')
ax2.set_ylabel('|Δfc| from\nprevious cycle')
ax2.set_xticks(cycles)
ax2.spines[['top', 'right']].set_visible(False)
ax2.axvline(10, color='0.6', lw=1, ls='--', zorder=1)

fig.savefig(f"{OUTDIR}/optimization_convergence.png", dpi=160, bbox_inches='tight')
plt.close(fig)
print("wrote optimization_convergence.png")

# ============================================================
# Figure 2: Spatial pattern of the optimized control adjustments
#           (surface, raw xx_* fields at the final leg-1 iterate)
# ============================================================
fields = [
    ('xx_theta', 'Optimized Δ(theta$_{IC}$)  [surface, K]', '0000000010'),
    ('xx_ptr1',  'Optimized Δ(DIC$_{IC}$)  [surface]',       '0000000010'),
    ('xx_ptr2',  'Optimized Δ(ALK$_{IC}$)  [surface]',       '0000000010'),
]

fig = plt.figure(figsize=(16, 3.6))
gs = fig.add_gridspec(1, 3, wspace=0.08, top=0.98, bottom=0.02)
axes = [fig.add_subplot(gs[0, i], projection=ccrs.Robinson()) for i in range(3)]

for ax, (varname, title, iterstr) in zip(axes, fields):
    data = load_tiled_surface(varname, iterstr)
    vmax = np.nanpercentile(np.abs(data), 99)
    if not np.isfinite(vmax) or vmax == 0:
        vmax = np.nanmax(np.abs(data)) or 1.0
    mesh = ax.pcolormesh(mlon, mlat, data, transform=ccrs.PlateCarree(),
                          cmap='RdBu_r', vmin=-vmax, vmax=vmax, shading='auto')
    ax.coastlines(linewidth=0.5, color='0.3')
    ax.set_title(title, fontsize=10)
    cb = fig.colorbar(mesh, ax=ax, orientation='horizontal', pad=0.05, shrink=0.85, extend='both')
    cb.ax.tick_params(labelsize=8)
    cb.formatter.set_powerlimits((-2, 2))
    cb.update_ticks()

fig.suptitle('Optimized control-vector adjustments after 10 M1QN3 iterations (cycle 10 iterate)\n'
             'these perturbations, added to the first-guess initial conditions, reduced fc by 6.60%',
             fontsize=12, y=1.06)
fig.savefig(f"{OUTDIR}/control_adjustment_maps.png", dpi=160, bbox_inches='tight')
plt.close(fig)
print("wrote control_adjustment_maps.png")
