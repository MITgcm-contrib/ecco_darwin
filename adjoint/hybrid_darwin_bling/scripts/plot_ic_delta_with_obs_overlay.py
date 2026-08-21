#!/usr/bin/env python3
"""Overlay SOCAT month-01 observation locations on the DIC IC delta map,
to show that the hybrid pilot's IC adjustments are concentrated near
actual SOCAT ship-track locations (~34% of ocean grid cells have any
month-01 observation at all) rather than spread uniformly -- the cost
function only has direct gradient information where there's an
observation; adjustments elsewhere come only from weak indirect
advective spread over the 30-day window.

Configure via environment variables:
  HYBRID_DIR          -- hybrid_darwin_bling/ root (default: parent of scripts/)
  DARWIN_INPUT_AD_DIR -- global_oce_biogeo_darwin/input_ad/ (for socat_pco2_clim_month01.nc)
  BLING_PRIOR_DIR     -- global_oce_biogeo_bling_SOCAT/input_ad/ (prior DIC IC)
  FIGURES_OUT_DIR     -- where to write the output PNG

Requires MITgcmutils (ships with MITgcm under utils/python/MITgcmutils --
add that directory to PYTHONPATH before running this script).
"""
import os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import netCDF4 as nc
from MITgcmutils import rdmds

HYBRID = os.environ.get("HYBRID_DIR", "..")
DARWIN_INPUT_AD = os.environ.get("DARWIN_INPUT_AD_DIR", "../../global_oce_biogeo_darwin/input_ad")
BLING_PRIOR = os.environ.get("BLING_PRIOR_DIR", "../../global_oce_biogeo_bling_SOCAT/input_ad")
OUT_DIR = os.environ.get("FIGURES_OUT_DIR", "../figures")

NZ, NY, NX = 15, 64, 128

XC = rdmds(f"{HYBRID}/run_darwin/XC")
YC = rdmds(f"{HYBRID}/run_darwin/YC")
lon_axis = XC[0, :]
lat_axis = YC[:, 0]

ds = nc.Dataset(f"{DARWIN_INPUT_AD}/socat_pco2_clim_month01.nc")
obs_lon = ds.variables["prof_lon"][:]
obs_lat = ds.variables["prof_lat"][:]
obs_lon360 = np.where(obs_lon < 0, obs_lon + 360, obs_lon)

prior = np.fromfile(f"{BLING_PRIOR}/dic_init.bin", dtype=">f4").reshape(NZ, NY, NX)[0]
final = np.fromfile(f"{HYBRID}/run_darwin/dic_darwin_init.bin", dtype=">f4").reshape(NZ, NY, NX)[0] / 1000.0
delta = np.ma.masked_where(prior == 0, final - prior)
dmax = np.abs(delta).max()

fig, ax = plt.subplots(figsize=(11, 5.5))
im = ax.imshow(delta, origin="lower", extent=[0, 360, -90, 90], vmin=-dmax, vmax=dmax,
                cmap="RdBu_r", aspect="auto")
ax.scatter(obs_lon360, obs_lat, s=4, c="black", alpha=0.5, marker=".",
           label=f"SOCAT month-01 obs (n={len(obs_lon)})")
ax.set_title("DIC IC adjustment (optimized - prior) with SOCAT month-01 observation locations overlaid")
ax.set_xlabel("longitude")
ax.set_ylabel("latitude")
ax.legend(loc="lower left", fontsize=8, frameon=True)
fig.colorbar(im, ax=ax, label="mol m$^{-3}$")
fig.tight_layout()

os.makedirs(OUT_DIR, exist_ok=True)
OUT = f"{OUT_DIR}/ic_delta_with_obs_overlay.png"
fig.savefig(OUT, dpi=140)
print(f"saved to {OUT}")

# quantitative check: |delta| at obs points vs elsewhere
j_idx = np.array([np.argmin(np.abs(lat_axis - lv)) for lv in obs_lat])
i_idx = np.array([np.argmin(np.abs(lon_axis - lv)) for lv in obs_lon360])
obs_mask = np.zeros((NY, NX), dtype=bool)
obs_mask[j_idx, i_idx] = True
wet = prior != 0
at_obs = np.abs(delta)[obs_mask & wet].mean()
away = np.abs(delta)[~obs_mask & wet].mean()
print(f"mean |delta| at obs points: {at_obs:.5f}")
print(f"mean |delta| away from obs points: {away:.5f}")
print(f"ratio: {at_obs/away:.2f}x")
