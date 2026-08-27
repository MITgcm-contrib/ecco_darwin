#!/usr/bin/env python3
"""Compute Darwin's own SOCAT-month01 pCO2 cost function J at the current
state, replicating pkg/profiles' cost formula exactly:
    J = sum_k  weight_k * mask_k * (model_k - obs_k)**2
Every SOCAT month01 observation lands exactly on a model grid-cell
center, so this collapses to a direct grid lookup, no interpolation.
Identical logic to hybrid_darwin_bling/scripts/compute_darwin_cost.py,
just pointed at this project's own run_darwin/.

Usage: compute_darwin_cost_dic.py <cycle>
Prints J to stdout and writes it to darwin_cost.<cycle>.txt.
"""
import sys
import numpy as np
import netCDF4 as nc

sys.path.insert(0, "/Users/carrolld/Documents/research/ECCO/offline/scripts")
from MITgcmutils import rdmds

DARWIN_RUN = "/Users/carrolld/Documents/research/adjoint/MITgcm/verification/hybrid_darwin_dic/run_darwin"
DARWIN_INPUT_AD = "/Users/carrolld/Documents/research/adjoint/MITgcm/verification/global_oce_biogeo_darwin/input_ad"


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
