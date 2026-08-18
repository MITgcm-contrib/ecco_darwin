#!/usr/bin/env python3
"""Compute Darwin's own SOCAT pCO2 cost function J, summed over all 12
months (matching BLING's real full-year data.profiles: PCO term active
for every month, mult_profiles(m,13)=1.0). Same per-month formula as
compute_darwin_cost.py, using the darwinPCO2Snap01..12 diagnostic
snapshots (one per month, at each month's own shared SOCAT observation
timestamp -- see run_darwin_1yr/data.diagnostics for the exact
timesteps and their derivation).

Usage: compute_darwin_cost_1yr.py <cycle>
Prints total J to stdout and writes it (plus the per-month breakdown) to
darwin_cost.<cycle>.txt in DARWIN_RUN_DIR.

Requires MITgcmutils (ships with MITgcm under utils/python/MITgcmutils --
add that directory to PYTHONPATH before running this script).

Configure via environment variables:
  DARWIN_RUN_DIR      -- this pilot's run_darwin_1yr/
  DARWIN_INPUT_AD_DIR -- global_oce_biogeo_darwin/input_ad/ (for socat_pco2_clim_month*.nc)
"""
import os
import sys
import numpy as np
import netCDF4 as nc
from MITgcmutils import rdmds

DARWIN_RUN = os.environ.get("DARWIN_RUN_DIR", "../run_darwin_1yr")
DARWIN_INPUT_AD = os.environ.get("DARWIN_INPUT_AD_DIR", "../../global_oce_biogeo_darwin/input_ad")

# month -> darwinPCO2SnapNN iteration (see run_darwin_1yr/data.diagnostics)
MONTH_TIMESTEP = {
    1: 1392, 2: 4272, 3: 7056, 4: 9936, 5: 12912, 6: 15792,
    7: 18672, 8: 21552, 9: 24432, 10: 27312, 11: 30192, 12: 33072,
}


def main(cycle):
    XC = rdmds(f"{DARWIN_RUN}/XC")
    YC = rdmds(f"{DARWIN_RUN}/YC")
    lon_axis = XC[0, :]
    lat_axis = YC[:, 0]

    total_J = 0.0
    per_month = {}
    for m, ts in MONTH_TIMESTEP.items():
        ds = nc.Dataset(f"{DARWIN_INPUT_AD}/socat_pco2_clim_month{m:02d}.nc")
        obs_lon = ds.variables["prof_lon"][:]
        obs_lat = ds.variables["prof_lat"][:]
        obs_lon360 = np.where(obs_lon < 0, obs_lon + 360, obs_lon)
        obs_pco2 = ds.variables["prof_PCO"][:, 0]
        obs_weight = ds.variables["prof_PCOweight"][:, 0]

        j_idx = np.array([np.argmin(np.abs(lat_axis - lv)) for lv in obs_lat])
        i_idx = np.array([np.argmin(np.abs(lon_axis - lv)) for lv in obs_lon360])
        assert np.all(np.abs(lat_axis[j_idx] - obs_lat) < 1e-6)
        assert np.all(np.abs(lon_axis[i_idx] - obs_lon360) < 1e-6)

        field = rdmds(f"{DARWIN_RUN}/darwinPCO2Snap{m:02d}", ts)
        model_vals = field[j_idx, i_idx]
        J_month = float(np.sum(obs_weight * (model_vals - obs_pco2) ** 2))
        per_month[m] = J_month
        total_J += J_month

    with open(f"{DARWIN_RUN}/darwin_cost.{cycle:010d}.txt", "w") as fh:
        fh.write(f"{total_J!r}\n")
        for m, J_month in per_month.items():
            fh.write(f"month{m:02d}: {J_month!r}\n")

    print(f"cycle {cycle}: Darwin total J = {total_J:.6e}  "
          f"(per-month: {[f'{per_month[m]:.3e}' for m in sorted(per_month)]})")
    return total_J


if __name__ == "__main__":
    main(int(sys.argv[1]))
