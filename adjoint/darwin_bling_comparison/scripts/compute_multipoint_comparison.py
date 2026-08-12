#!/usr/bin/env python3
"""38-point x 6-tracer Darwin forward-FD vs BLING adjoint comparison,
mirroring Kay's GCHP validation methodology (see readme.txt): relative
(+/-5%) IC perturbation at each point, centred FD normalized by the
ABSOLUTE perturbation actually applied (eta * c0_m, per point -- not a
shared constant), a true unperturbed baseline for the
nonlinearity/asymmetry diagnostic, and a final slope/correlation
regression of adjoint vs FD.

Darwin side: local dJ/dtracer at each point computed from that point's own
COINCIDENT SOCAT observation (every one of the 38 points was chosen from
the SOCAT obs grid, so pkg/profiles' bilinear interpolation collapses to
an exact single-cell lookup -- no interpolation machinery needed) using
the full darwinPCO2Snap field from the matched-IC baseline/plus/minus
runs (run_multi_baseline, run_multi_<TRACER>_{plus,minus} from
setup_and_launch_multipoint.sh). BLING side: adxx_ptr* (already computed
once, whole-domain, in BLING's own pilot adjoint run) converted from
control-space to physical units via sqrt(ctrl weight), per
ctrl_map_genarr3d.F's preconditioning.

Requires MITgcmutils (ships with MITgcm under utils/python/MITgcmutils --
add that directory to PYTHONPATH before running this script).

Configure via environment variables:
  DARWIN_RUN_DIR      -- Darwin experiment dir with run/, run_multi_*/
  DARWIN_INPUT_AD_DIR -- Darwin's input_ad/ (for socat_pco2_clim_month01.nc)
  BLING_DIR           -- BLING_SOCAT experiment dir (pilot_adxx_fields.npz,
                          input_ad/wt_*.bin)
"""
import json
import os
import numpy as np
import netCDF4 as nc
from MITgcmutils import rdmds

DARWIN = os.environ.get("DARWIN_RUN_DIR", ".")
INPUT_AD = os.environ.get("DARWIN_INPUT_AD_DIR", f"{DARWIN}/input_ad")
BLING = os.environ.get("BLING_DIR", "../global_oce_biogeo_bling_SOCAT")
NZ, NY, NX = 15, 64, 128

TRACERS = ['DIC', 'ALK', 'O2', 'NO3', 'PO4', 'FeT']
ICFILE = {'DIC': 'dic_darwin_init.bin', 'ALK': 'alk_darwin_init.bin', 'O2': 'o2_darwin_init.bin',
          'NO3': 'no3_darwin_init.bin', 'PO4': 'po4_darwin_init.bin', 'FeT': 'fet_darwin_init.bin'}
ADXX_KEY = {'DIC': 'ptr1', 'ALK': 'ptr2', 'O2': 'ptr3', 'NO3': 'ptr4', 'PO4': 'ptr5', 'FeT': 'ptr6'}
WT_FILE = {'DIC': 'wt_DIC.bin', 'ALK': 'wt_ALK.bin', 'O2': 'wt_O2.bin',
           'NO3': 'wt_NO3.bin', 'PO4': 'wt_PO4.bin', 'FeT': 'wt_FE.bin'}

# --- grid + SOCAT obs (for the coincident-obs lookup) ---
XC = rdmds(f'{DARWIN}/run/XC'); YC = rdmds(f'{DARWIN}/run/YC')
lon_axis = XC[0, :]; lat_axis = YC[:, 0]

ds = nc.Dataset(f'{INPUT_AD}/socat_pco2_clim_month01.nc')
obs_lon = ds.variables['prof_lon'][:]; obs_lat = ds.variables['prof_lat'][:]
obs_lon360 = np.where(obs_lon < 0, obs_lon + 360, obs_lon)
obs_pco2_all = ds.variables['prof_PCO'][:, 0]
obs_weight_all = ds.variables['prof_PCOweight'][:, 0]
obs_j = np.array([np.argmin(np.abs(lat_axis - lv)) for lv in obs_lat])
obs_i = np.array([np.argmin(np.abs(lon_axis - lv)) for lv in obs_lon360])
obs_lookup = {(int(j), int(i)): (float(obs_pco2_all[k]), float(obs_weight_all[k]))
              for k, (j, i) in enumerate(zip(obs_j, obs_i))}

# --- BLING adjoint fields + ctrl weights (already computed, whole-domain) ---
adxx = np.load(f'{BLING}/pilot_adxx_fields.npz')
wt = {t: np.fromfile(f'{BLING}/input_ad/{WT_FILE[t]}', dtype='>f4').reshape(NZ, NY, NX) for t in TRACERS}

# --- Darwin baseline field (shared across all tracers) ---
baseline_pco2 = rdmds(f'{DARWIN}/run_multi_baseline/darwinPCO2Snap', 1392)

rows = []
for t in TRACERS:
    plus_pco2 = rdmds(f'{DARWIN}/run_multi_{t}_plus/darwinPCO2Snap', 1392)
    minus_pco2 = rdmds(f'{DARWIN}/run_multi_{t}_minus/darwinPCO2Snap', 1392)
    with open(f'{DARWIN}/run_multi_{t}_plus/perturb_log.json') as fh:
        plog = json.load(fh)
    eta = plog['eta']

    for pt in plog['points']:
        j, i, base = pt['j'], pt['i'], pt['base']
        key = (j, i)
        if key not in obs_lookup:
            continue  # shouldn't happen -- points were chosen from the obs set
        obs_pco2, obs_w = obs_lookup[key]

        model0 = float(baseline_pco2[j, i])
        model_p = float(plus_pco2[j, i])
        model_m = float(minus_pco2[j, i])
        delta = eta * base  # absolute perturbation actually applied at this point

        dpco2_dtracer = (model_p - model_m) / (2 * delta)
        dJ_darwin_native = 2 * obs_w * (model0 - obs_pco2) * dpco2_dtracer  # per Darwin's native units (mmol/m3 etc)
        dJ_darwin = dJ_darwin_native * 1000.0  # -> per mol/m3, matching BLING's ctrl-space units

        # nonlinearity/asymmetry diagnostic (Kay's readme.txt, section 5)
        Jp = obs_w * (model_p - obs_pco2) ** 2
        Jm = obs_w * (model_m - obs_pco2) ** 2
        J0 = obs_w * (model0 - obs_pco2) ** 2
        asym = (Jp - J0) - (J0 - Jm)
        denom = abs(Jp - Jm)
        asym_rel = abs(asym) / denom if denom > 0 else np.nan

        wval = float(wt[t][0, j, i])
        sqrtw = np.sqrt(wval) if wval > 0 else np.nan
        a = float(adxx[ADXX_KEY[t]][0, j, i])
        dJ_bling = a * sqrtw

        rows.append(dict(tracer=t, j=j, i=i, dJ_darwin=dJ_darwin, dJ_bling=dJ_bling,
                          asym_rel=asym_rel, model0=model0, obs_pco2=obs_pco2))

if __name__ == '__main__':
    import csv
    OUT = f'{DARWIN}/multipoint_comparison_results.csv'
    with open(OUT, 'w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)
    print(f"wrote {len(rows)} rows to {OUT}")

    for t in TRACERS:
        sub = [r for r in rows if r['tracer'] == t]
        dd = np.array([r['dJ_darwin'] for r in sub])
        bb = np.array([r['dJ_bling'] for r in sub])
        finite = np.isfinite(dd) & np.isfinite(bb) & (dd != 0)
        if finite.sum() >= 2:
            slope, intercept = np.polyfit(dd[finite], bb[finite], 1)
            corr = np.corrcoef(dd[finite], bb[finite])[0, 1]
        else:
            slope, corr = np.nan, np.nan
        med_asym = np.nanmedian([r['asym_rel'] for r in sub])
        print(f"{t:5s} n={len(sub):3d}  slope={slope:8.3f}  corr={corr:6.3f}  median|asym|/|Jp-Jm|={med_asym:.4f}")
