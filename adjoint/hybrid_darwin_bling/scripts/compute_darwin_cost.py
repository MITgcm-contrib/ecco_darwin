#!/usr/bin/env python3
"""Compute Darwin's own SOCAT-month01 pCO2 cost function J at the current
state, replicating pkg/profiles' cost formula exactly (verified against
cost_profiles.F/cost_final.F -- see ../../darwin_bling_comparison/README.md):
    J = sum_k  weight_k * mask_k * (model_k - obs_k)**2
Every SOCAT month01 observation lands exactly on a model grid-cell
center, so this collapses to a direct grid lookup, no interpolation.

Usage: compute_darwin_cost.py <cycle>
Prints J to stdout and writes it to darwin_cost.<cycle>.txt in DARWIN_RUN_DIR.

Requires MITgcmutils (ships with MITgcm under utils/python/MITgcmutils --
add that directory to PYTHONPATH before running this script).

Configure via environment variables:
  DARWIN_RUN_DIR      -- this pilot's run_darwin/ (has XC/YC, darwinPCO2Snap)
  DARWIN_INPUT_AD_DIR -- global_oce_biogeo_darwin/input_ad/ (for socat_pco2_clim_month01.nc)
"""
import os
import sys
import numpy as np
import netCDF4 as nc
from MITgcmutils import rdmds

DARWIN_RUN = os.environ.get("DARWIN_RUN_DIR", "../run_darwin")
DARWIN_INPUT_AD = os.environ.get("DARWIN_INPUT_AD_DIR", "../../global_oce_biogeo_darwin/input_ad")

TRACERS = {
    'PO4': 0.02402,
    'DIC': 8.57,
    'ALK': 12.57,
    'O2':  2.659,
    'NO3': 0.3054,
    'FeT': 1.068e-5,
}


def main(cycle):
    XC = rdmds(f"{DARWIN_RUN}/XC")
    YC = rdmds(f"{DARWIN_RUN}/YC")
    lon_axis = XC[0, :]
    lat_axis = YC[:, 0]

    ds = nc.Dataset(f"{DARWIN_INPUT_AD}/socat_pco2_clim_month01.nc")
    obs_lon = ds.variables["prof_lon"][:]
    obs_lat = ds.variables["prof_lat"][:]
    obs_lon360 = np.where(obs_lon < 0, obs_lon + 360, obs_lon)
    obs_pco2 = ds.variables["prof_PCO"][:, 0]
    obs_weight = ds.variables["prof_PCOweight"][:, 0]

    j_idx = np.array([np.argmin(np.abs(lat_axis - lv)) for lv in obs_lat])
    i_idx = np.array([np.argmin(np.abs(lon_axis - lv)) for lv in obs_lon360])
    assert np.all(np.abs(lat_axis[j_idx] - obs_lat) < 1e-6)
    assert np.all(np.abs(lon_axis[i_idx] - obs_lon360) < 1e-6)

    field = rdmds(f"{DARWIN_RUN}/darwinPCO2Snap", 1392)
    model_vals = field[j_idx, i_idx]
    J = float(np.sum(obs_weight * (model_vals - obs_pco2) ** 2))

    with open(f"{DARWIN_RUN}/darwin_cost.{cycle:010d}.txt", "w") as fh:
        fh.write(f"{J!r}\n")
    print(f"cycle {cycle}: Darwin J = {J:.6e}")
    return J


if __name__ == "__main__":
    main(int(sys.argv[1]))
