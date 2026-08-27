#!/usr/bin/env python3
"""38-point x 5-tracer Darwin forward-FD vs pkg/dic adjoint comparison,
mirroring darwin_bling_comparison/scripts/compute_multipoint_comparison.py
exactly, but for pkg/dic instead of BLING (DIC, ALK, PO4, DOP, O2 -- pkg/dic's
own tracer set/order; no NO3 or Fe, since pkg/dic doesn't carry them).

Darwin's own dJ/dtracer values for DIC, ALK, PO4, O2 are REUSED verbatim
from the original BLING comparison's multipoint_comparison_results.csv --
Darwin's side of that computation never depended on which adjoint model
(BLING or pkg/dic) it's being compared against, only on Darwin's own
matched-IC perturbation runs. Only DOP is new here (BLING's comparison
never covered it, since BLING doesn't carry a DOP control), computed from
a fresh pair of Darwin +/-5% DOP perturbation runs
(run_multi_DOP_plus/minus) using the exact same 38 points/eta/methodology.

pkg/dic's adjoint (adxx_ptr1..5) was already extracted whole-domain into
pilot_adxx_fields_dic.npz (see run_ad/'s extraction), and is converted from
control-space to physical units via sqrt(ctrl weight), per
ctrl_map_genarr3d.F's preconditioning -- same formula as the BLING side.

CORRUPTED-POINT MASKING (added after discovering a genuine, field-specific
numerical instability in pkg/dic's own adjoint for PO4 and DOP -- confirmed
via isolated single-control reruns (hybrid_darwin_dic/run_dic_C,
run_dic_D) to be identical to the original 7-control extraction bit-for-
bit, so it is NOT the slot-position bug documented in
hybrid_darwin_dic/run_dic/data.ctrl -- it is present regardless of which
control-array slot PO4/DOP occupy, or even when each is the SOLE
registered control). Points where |adxx| > 1e6 (a normal value is
~1e2-1e5) are flagged 'corrupted' and excluded from the slope/r
regression -- see README.md for the full accounting (1/38 points for PO4,
15/38 for DOP).
"""
import csv
import json
import numpy as np
import netCDF4 as nc
import sys

sys.path.insert(0, '/Users/carrolld/Documents/research/ECCO/offline/scripts')
from MITgcmutils import rdmds

DARWIN = '/Users/carrolld/Documents/research/adjoint/MITgcm/verification/global_oce_biogeo_darwin'
DIC_SOCAT = '/Users/carrolld/Documents/research/adjoint/MITgcm/verification/global_oce_biogeo_dic_SOCAT'
COMPARISON = '/Users/carrolld/Documents/research/adjoint/MITgcm/verification/darwin_dic_comparison'
NZ, NY, NX = 15, 64, 128

# pkg/dic's own tracer order (dic_tr_register.F): DIC=1, Alk=2, PO4=3, DOP=4, O2=5
TRACERS = ['DIC', 'ALK', 'PO4', 'DOP', 'O2']
ADXX_KEY = {'DIC': 'ptr1', 'ALK': 'ptr2', 'PO4': 'ptr3', 'DOP': 'ptr4', 'O2': 'ptr5'}
WT_FILE = {'DIC': 'wt_DIC.bin', 'ALK': 'wt_ALK.bin', 'PO4': 'wt_PO4.bin',
           'DOP': 'wt_DOP.bin', 'O2': 'wt_O2.bin'}
# tracers with an existing Darwin FD result in the BLING comparison's CSV
REUSED_FROM_BLING_CSV = ['DIC', 'ALK', 'PO4', 'O2']

# --- grid + SOCAT obs (for the coincident-obs lookup, needed for DOP only) ---
XC = rdmds(f'{DARWIN}/run/XC'); YC = rdmds(f'{DARWIN}/run/YC')
lon_axis = XC[0, :]; lat_axis = YC[:, 0]

ds = nc.Dataset(f'{DARWIN}/input_ad/socat_pco2_clim_month01.nc')
obs_lon = ds.variables['prof_lon'][:]; obs_lat = ds.variables['prof_lat'][:]
obs_lon360 = np.where(obs_lon < 0, obs_lon + 360, obs_lon)
obs_pco2_all = ds.variables['prof_PCO'][:, 0]
obs_weight_all = ds.variables['prof_PCOweight'][:, 0]
obs_j = np.array([np.argmin(np.abs(lat_axis - lv)) for lv in obs_lat])
obs_i = np.array([np.argmin(np.abs(lon_axis - lv)) for lv in obs_lon360])
obs_lookup = {(int(j), int(i)): (float(obs_pco2_all[k]), float(obs_weight_all[k]))
              for k, (j, i) in enumerate(zip(obs_j, obs_i))}

# --- pkg/dic adjoint fields + ctrl weights (whole-domain, already extracted) ---
adxx = np.load(f'{COMPARISON}/pilot_adxx_fields_dic.npz')
wt = {t: np.fromfile(f'{DIC_SOCAT}/input_ad/{WT_FILE[t]}', dtype='>f4').reshape(NZ, NY, NX)
      for t in TRACERS}

# --- reused Darwin FD rows (DIC, ALK, PO4, O2) from the BLING comparison ---
bling_rows = list(csv.DictReader(open(f'{DARWIN}/multipoint_comparison_results.csv')))
reused = {t: {(int(r['j']), int(r['i'])): r for r in bling_rows if r['tracer'] == t}
          for t in REUSED_FROM_BLING_CSV}

# --- fresh Darwin FD for DOP ---
baseline_pco2 = rdmds(f'{DARWIN}/run_multi_baseline/darwinPCO2Snap', 1392)
plus_pco2 = rdmds(f'{DARWIN}/run_multi_DOP_plus/darwinPCO2Snap', 1392)
minus_pco2 = rdmds(f'{DARWIN}/run_multi_DOP_minus/darwinPCO2Snap', 1392)
with open(f'{DARWIN}/run_multi_DOP_plus/perturb_log.json') as fh:
    dop_plog = json.load(fh)
dop_eta = dop_plog['eta']

dop_rows = {}
for pt in dop_plog['points']:
    j, i, base = pt['j'], pt['i'], pt['base']
    key = (j, i)
    obs_pco2, obs_w = obs_lookup[key]
    model0 = float(baseline_pco2[j, i])
    model_p = float(plus_pco2[j, i])
    model_m = float(minus_pco2[j, i])
    delta = dop_eta * base
    dpco2_dtracer = (model_p - model_m) / (2 * delta)
    dJ_darwin_native = 2 * obs_w * (model0 - obs_pco2) * dpco2_dtracer
    dJ_darwin = dJ_darwin_native * 1000.0  # -> per mol/m3, matching ctrl-space units
    Jp = obs_w * (model_p - obs_pco2) ** 2
    Jm = obs_w * (model_m - obs_pco2) ** 2
    J0 = obs_w * (model0 - obs_pco2) ** 2
    asym = (Jp - J0) - (J0 - Jm)
    denom = abs(Jp - Jm)
    asym_rel = abs(asym) / denom if denom > 0 else np.nan
    dop_rows[key] = dict(dJ_darwin=dJ_darwin, asym_rel=asym_rel, model0=model0, obs_pco2=obs_pco2)

# --- assemble comparison rows ---
rows = []
for t in TRACERS:
    if t == 'DOP':
        points = list(dop_rows.keys())
    else:
        points = list(reused[t].keys())
    for (j, i) in points:
        if t == 'DOP':
            src = dop_rows[(j, i)]
            dJ_darwin = src['dJ_darwin']
            asym_rel = src['asym_rel']
            model0 = src['model0']
            obs_pco2 = src['obs_pco2']
        else:
            src = reused[t][(j, i)]
            dJ_darwin = float(src['dJ_darwin'])
            asym_rel = float(src['asym_rel']) if src['asym_rel'] != '' else np.nan
            model0 = float(src['model0'])
            obs_pco2 = float(src['obs_pco2'])

        wval = float(wt[t][0, j, i])
        sqrtw = np.sqrt(wval) if wval > 0 else np.nan
        a = float(adxx[ADXX_KEY[t]][0, j, i])
        dJ_dic = a * sqrtw
        # NOTE: an earlier version flagged `abs(a) > 1e6` as 'corrupted' and
        # excluded those points. That test was wrong: `a` is the RAW
        # CONTROL-SPACE adjoint, whose per-tracer scale spans ~400x (surface
        # median |adxx|: DIC 1.1e3, ALK 8.9e2, PO4 1.3e4, DOP 4.1e5), so a
        # single absolute cut lands at DOP's ~60th percentile and removed 15
        # of 38 ordinary points. Excluding on the dependent variable like
        # that is selection bias and it inverted DOP's regression. See the
        # README for the evidence that no corrupted points exist.

        rows.append(dict(tracer=t, j=j, i=i, dJ_darwin=dJ_darwin, dJ_dic=dJ_dic,
                          asym_rel=asym_rel, model0=model0, obs_pco2=obs_pco2,
                          ))

OUT = f'{COMPARISON}/results/multipoint_comparison_results_dic.csv'
with open(OUT, 'w', newline='') as fh:
    w = csv.DictWriter(fh, fieldnames=['tracer', 'j', 'i', 'dJ_darwin', 'dJ_dic',
                                       'asym_rel', 'model0', 'obs_pco2'])
    w.writeheader()
    w.writerows(rows)
print(f"wrote {len(rows)} rows to {OUT}")

def theil_sen(x, y):
    """Median of pairwise slopes: immune to single-point leverage."""
    sl = [(y[b] - y[a]) / (x[b] - x[a])
          for a in range(len(x)) for b in range(a + 1, len(x)) if x[b] != x[a]]
    return float(np.median(sl)) if sl else float('nan')


def spearman(x, y):
    rx = np.argsort(np.argsort(x)).astype(float)
    ry = np.argsort(np.argsort(y)).astype(float)
    return float(np.corrcoef(rx, ry)[0, 1])


print()
print("All 38 points are used for every tracer -- see README on why the old")
print("'corrupted' exclusion was withdrawn. OLS slope/r are reported next to")
print("leverage-immune Theil-Sen slope and Spearman rank r; where the two")
print("disagree sharply, the OLS fit is being carried by one or two points.")
print()
hdr = f"{'tracer':6} {'n':>3} {'OLS slope':>10} {'OLS r':>7} {'TS slope':>10} {'rank r':>7}"
print(hdr)
print('-' * len(hdr))
for t in TRACERS:
    sub = [r for r in rows if r['tracer'] == t]
    x = np.array([r['dJ_darwin'] for r in sub])
    y = np.array([r['dJ_dic'] for r in sub])
    if np.ptp(x) < 1e-6 * max(abs(x).max(), 1) or np.ptp(y) == 0:
        print(f"{t:6} {len(sub):3d}   (Darwin and/or pkg/dic sensitivity identically 0"
              f" -- slope/r undefined)")
        continue
    slope = np.polyfit(x, y, 1)[0]
    print(f"{t:6} {len(sub):3d} {slope:10.3f} {np.corrcoef(x, y)[0, 1]:7.3f} "
          f"{theil_sen(x, y):10.3f} {spearman(x, y):7.3f}")
