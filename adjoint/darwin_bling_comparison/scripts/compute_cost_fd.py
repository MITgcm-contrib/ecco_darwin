#!/usr/bin/env python3
"""Compute dJ/dtracer via finite differences of Darwin's own SOCAT-month01
pCO2 cost function, replicating pkg/profiles' cost formula:
    J = sum_k  weight_k * mask_k * (model_k - obs_k)**2
(mult_profiles=1 for the PCO term, 0 for everything else, per data.profiles).

Spatial matching: every SOCAT month01 observation lands EXACTLY on a model
grid-cell center (verified separately -- 2759/2759 exact lon+lat matches),
so pkg/profiles' bilinear interpolation collapses to a direct grid lookup
here (weight 1 at that cell, 0 at the other 3 corners) -- no interpolation
machinery needed. All matched cells are ocean (data.profiles: "ocean points
only"), so mask_k = 1 throughout.

Requires full pCO2-field snapshots (darwinPCO2Snap at iteration 1392, the
single timestep matching SOCAT month01's shared observation time) saved by
run_sensitivity_sweep_fields.sh into
pco2_fields/<TRACER>_run_{plus,minus}_darwinPCO2Snap.*

Requires MITgcmutils (ships with MITgcm under utils/python/MITgcmutils --
add that directory to PYTHONPATH before running this script).

Configure via environment variables:
  DARWIN_RUN_DIR      -- Darwin experiment dir with run/, pco2_fields/
  DARWIN_INPUT_AD_DIR -- Darwin's input_ad/ (for socat_pco2_clim_month01.nc)
"""
import os
import numpy as np
import netCDF4 as nc
from MITgcmutils import rdmds

DARWIN = os.environ.get("DARWIN_RUN_DIR", ".")
INPUT_AD = os.environ.get("DARWIN_INPUT_AD_DIR", f"{DARWIN}/input_ad")

TRACERS = {
    'PO4': 0.02402,
    'DIC': 8.57,
    'ALK': 12.57,
    'O2':  2.659,
    'NO3': 0.3054,
    'FeT': 1.068e-5,
}

XC = rdmds(f'{DARWIN}/run/XC')
YC = rdmds(f'{DARWIN}/run/YC')
lon_axis = XC[0, :]
lat_axis = YC[:, 0]

ds = nc.Dataset(f'{INPUT_AD}/socat_pco2_clim_month01.nc')
obs_lon = ds.variables['prof_lon'][:]
obs_lat = ds.variables['prof_lat'][:]
obs_lon360 = np.where(obs_lon < 0, obs_lon + 360, obs_lon)
obs_pco2 = ds.variables['prof_PCO'][:, 0]
obs_weight = ds.variables['prof_PCOweight'][:, 0]

# map each obs to its exact (j,i) grid index
j_idx = np.array([np.argmin(np.abs(lat_axis - lv)) for lv in obs_lat])
i_idx = np.array([np.argmin(np.abs(lon_axis - lv)) for lv in obs_lon360])
# sanity: confirm exact match (within float tolerance)
assert np.all(np.abs(lat_axis[j_idx] - obs_lat) < 1e-6)
assert np.all(np.abs(lon_axis[i_idx] - obs_lon360) < 1e-6)


def cost(field):
    model_vals = field[j_idx, i_idx]
    return float(np.sum(obs_weight * (model_vals - obs_pco2) ** 2))


if __name__ == '__main__':
    print(f"{'tracer':6s} {'eps':>10s} {'J_plus':>16s} {'J_minus':>16s} {'dJ/dtracer':>16s}")
    results = {}
    for name, eps in TRACERS.items():
        try:
            f_plus = rdmds(f'{DARWIN}/pco2_fields/{name}_run_plus_darwinPCO2Snap', 1392)
            f_minus = rdmds(f'{DARWIN}/pco2_fields/{name}_run_minus_darwinPCO2Snap', 1392)
        except Exception as e:
            print(f"{name:6s}  -- field not available yet ({e})")
            continue
        J_plus = cost(f_plus)
        J_minus = cost(f_minus)
        fd = (J_plus - J_minus) / (2 * eps)
        results[name] = fd
        print(f"{name:6s} {eps:10.4g} {J_plus:16.6e} {J_minus:16.6e} {fd:16.6e}")

    import json
    with open(f'{DARWIN}/cost_fd_results.json', 'w') as fh:
        json.dump(results, fh, indent=2)
    print(f"\nSaved to {DARWIN}/cost_fd_results.json")
