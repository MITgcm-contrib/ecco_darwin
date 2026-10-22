"""
FORTE-campaign-framework counterpart to tools/rank_observation_value.py: the SAME
Fisher-information/degrees-of-freedom-for-signal (DFS) ranking method -- see that
script's docstring for the full method, effective-sample-size correction, and
caveats, not repeated here -- applied to the Green's functions and observations
built for the climatological FORTE campaign run instead of the 2007/2008/2009
analog-year nature run. docs/osse/ and runs/osse_green/ are completely untouched;
this reads only docs/osse_forte/*.csv and runs/forte_green/green_functions_*.npz.

WHAT'S GENUINELY DIFFERENT:
  - 25 controls, not 19: discharge, clb_T, 9 marine (clb_*) boundary species, 8
    upstream (cub_*) controls (cub_TOC, cub_NO3, cub_NH4, cub_PO4, cub_dSi,
    cub_DIC, cub_ALK, cub_DIA), PLUS -- new, 2026-10-20 -- 6 REACH-LOCAL controls
    (rate_aerdeg_mouth/mid/head, lateral_toc_mouth/mid/head) testing sensitivity to
    spatially-distributed in-channel unknowns rather than just the two boundaries,
    since continuous spatial observations (airborne maps, multi-station grab
    sampling) also carry information about what's happening INSIDE the domain --
    see build_green_perturbations_forte.py's docstring for the full design. No
    special-casing needed for ANY of these in build_H_sigma: unlike mooring_offshore
    (which reads a marine boundary FORCING file directly and needs
    VAR_TO_MARINE_CONTROL), every cub_*/reach-local control's sensitivity reaches an
    interior observation purely through the model's own simulated transport, i.e.
    the ordinary `{model_var}__{control}` Green's-function lookup already in the
    interior-point loop below.
  - No per-year loop: tools/build_osse_observations_forte.py's rows all carry a
    single fixed year label ("forte", a string, not a real calendar year) and a
    single BASE_DAY=365 -- there is only one climatological pass, so every
    YEAR_START_DAY[year] dict lookup in the original collapses to one constant.
  - _marine_forcing_value reads the CLIMATOLOGICAL marine CSV (365 rows, one
    representative year) by DAY-OF-YEAR (`day_idx % 365`), not an absolute index
    into a 6935-row multi-year series -- same reasoning as build_osse_observations_
    forte.py's offshore-mooring fix.
  - Green's-function arrays (runs/forte_green/green_functions_<site>.npz) only have
    real (non-zero) content for day-of-year 0-273 (tools/build_green_functions_
    forte.py) -- doy is clipped to CAMPAIGN_END_DAY=273, not 364, though this barely
    matters in practice since days beyond 273 are never actually sampled by any
    channel in the first place.
  - The MARINE_NATIVE_DAYS~30-day effective-sample-size cap (B) is carried over
    unchanged: the climatology is built FROM real ~monthly-native ECCO-Darwin data
    (tools/build_forte_campaign_forcing.py averages across 19 years at each day-of-
    year, which are interpolated from the same recurring monthly boundaries), so
    the same physical smoothing-scale reasoning applies to the climatological
    series as to the raw interannual one.

Usage:  python tools/rank_observation_value_forte.py
Writes: docs/osse_forte/observation_value_ranking.csv,
        docs/osse_forte/observation_value_ranking_combined.csv
"""
import csv
import warnings
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent.parent
RUNS = ROOT / "runs" / "forte_green"
OSSE = ROOT / "docs" / "osse_forte"
FORC = ROOT / "forcing"
SITES = ["colville", "kuparuk", "sagavanirktok"]
MARINE_SPECIES = ["S", "DIC", "ALK", "NO3", "NH4", "PO4", "dSi", "O2", "TOC"]
UPSTREAM_CONST_SPECIES = ["NO3", "NH4", "PO4", "dSi", "DIC", "ALK", "DIA"]
REACHES = ["mouth", "mid", "head"]
LOCAL_CONTROLS = ([f"rate_aerdeg_{r}" for r in REACHES]
                  + [f"lateral_toc_{r}" for r in REACHES])
CONTROLS = (["discharge", "clb_T"] + [f"clb_{sp}" for sp in MARINE_SPECIES]
            + ["cub_TOC"] + [f"cub_{sp}" for sp in UPSTREAM_CONST_SPECIES]
            + LOCAL_CONTROLS)
N_CTRL = len(CONTROLS)
DELXI_KM = 0.1
BASE_DAY = 365          # single climatological "year" -- matches build_osse_
                        # observations_forte.py's BASE_DAY
CAMPAIGN_END_DAY = 273  # matches build_green_functions_forte.py
PRIOR_SIGMA_FRAC = 0.20

# observation variable -> (model variable in the Green's-function arrays, scale)
VAR_MAP = {
    "DOC": ("TOC", 1.0), "CDOM_proxy": ("TOC", 1.0),
    "NH4": ("NH4", 1.0), "NO3": ("NO3", 1.0), "PO4": ("PO4", 1.0), "dSi": ("dSi", 1.0),
    "DIC": ("DIC", 1.0), "ALK": ("ALK", 1.0),
    "Chl_a_mg_m3": ("DIA", 12.011 / 50.0), "Chl_a_proxy_mg_m3": ("DIA", 12.011 / 50.0),
    "salinity": ("S", 1.0), "temperature": ("T", 1.0),
    "O2": ("O2", 1.0), "SPM": ("SPM", 1.0), "depth_anomaly": ("depth", 1.0),
}
EXCLUDED_VARS = {"current_speed", "arrival_time_at_mouth_days"}
# which CONTROL a given model var's boundary corresponds to, for the mooring_offshore
# direct-boundary-forcing special case
VAR_TO_MARINE_CONTROL = {"S": "clb_S", "T": "clb_T", "O2": "clb_O2", "DIC": "clb_DIC",
                          "ALK": "clb_ALK", "NO3": "clb_NO3", "NH4": "clb_NH4",
                          "PO4": "clb_PO4", "dSi": "clb_dSi", "TOC": "clb_TOC"}
MARINE_FILE_KEY = {"S": "s", "T": "t", "O2": "o2", "DIC": "dic", "ALK": "alk",
                   "NO3": "no3", "NH4": "nh4", "PO4": "po4", "dSi": "dsi", "TOC": "toc"}


def _read_csv(site, channel):
    p = OSSE / f"{site}_{channel}.csv"
    if not p.exists():
        return []
    with open(p) as f:
        return list(csv.DictReader(f))


def _load_green(site):
    p = RUNS / f"green_functions_{site}.npz"
    if not p.exists():
        return None
    return np.load(p)


def _marine_forcing_value(site, model_var, day_idx):
    """Baseline (unperturbed) value of the marine boundary CLIMATOLOGICAL forcing
    (365 rows, one representative year) for `model_var`, indexed by DAY-OF-YEAR
    (`day_idx % 365`) -- used only for the mooring_offshore special case. A
    climatology file has no absolute-day concept (repeatYear cycles the same 365
    rows every model year), unlike the original script's multi-year series."""
    key = MARINE_FILE_KEY[model_var]
    if model_var == "T":
        fp = FORC / f"{site}_t_marine_climatology.csv"
    else:
        fp = FORC / f"{site}_{key}_marine_climatology.csv"
    if not hasattr(_marine_forcing_value, "_cache"):
        _marine_forcing_value._cache = {}
    cache_key = (site, model_var)
    if cache_key not in _marine_forcing_value._cache:
        _marine_forcing_value._cache[cache_key] = np.genfromtxt(
            open(fp, encoding="utf-8-sig"), delimiter=",")
    vals = _marine_forcing_value._cache[cache_key]
    return vals[np.asarray(day_idx) % vals.size]


MARINE_NATIVE_DAYS = 30   # ECCO-Darwin v5 marine boundary: "monthly native resolution
                          # interpolated to daily" -- carried over unchanged, see
                          # module docstring for why it still applies to the
                          # climatological series built from that same data.


def _effective_sample_weight(abs_day, cell, is_offshore, varname):
    """Per-row weight in [0,1] correcting for duplicate/highly-correlated samples --
    see tools/rank_observation_value.py's docstring for the full derivation of caps
    (A) (mechanical, daily Green's-function resolution) and (B) (physical, ~monthly
    marine-forcing update rate, offshore only). Identical logic here; only the
    absolute-day values being binned differ (a single climatological pass instead
    of a multi-year series)."""
    n = len(abs_day)
    weight = np.ones(n)
    day_int = np.round(abs_day).astype(int)

    keys_a = list(zip(cell.tolist(), day_int.tolist()))
    counts_a = {}
    for k in keys_a:
        counts_a[k] = counts_a.get(k, 0) + 1
    weight_a = np.array([1.0 / counts_a[k] for k in keys_a])

    weight_b = np.ones(n)
    if is_offshore.any():
        bin_b = (day_int[is_offshore] // MARINE_NATIVE_DAYS)
        counts_b = {}
        for b in bin_b:
            counts_b[b] = counts_b.get(b, 0) + 1
        weight_b[is_offshore] = np.array([1.0 / counts_b[b] for b in bin_b])

    return weight_a * weight_b


def build_H_sigma(rows, green, site, M):
    """(H [n x N_CTRL], sigma [n], weight [n]) for one (site, observation-channel)
    row set, skipping rows whose variable has no mapping/sensitivity (see
    EXCLUDED_VARS). `weight` is the effective-sample-size correction from
    _effective_sample_weight."""
    H_rows, sig_rows, weight_rows = [], [], []
    by_var = {}
    for r in rows:
        if r["variable"] in EXCLUDED_VARS or r["variable"] not in VAR_MAP:
            continue
        if r["sigma"] in ("", "None", None):
            continue
        by_var.setdefault(r["variable"], []).append(r)

    for varname, group in by_var.items():
        model_var, scale = VAR_MAP[varname]
        t_days = np.array([float(r["t_days"]) for r in group])
        doy = np.clip(np.round(t_days - BASE_DAY).astype(int), 0, CAMPAIGN_END_DAY)
        abs_day = BASE_DAY + doy
        sigma = np.abs(np.array([float(r["sigma"]) for r in group]))
        channels = [r["channel"] for r in group]

        offshore_mask = np.array([c == "mooring_offshore" for c in channels])
        H = np.zeros((len(group), N_CTRL))
        cell = np.full(len(group), -1)   # -1 = offshore pseudo-cell (no spatial index)

        if offshore_mask.any():
            ctrl = VAR_TO_MARINE_CONTROL.get(model_var)
            if ctrl is not None:
                ci = CONTROLS.index(ctrl)
                vals = _marine_forcing_value(site, model_var, doy[offshore_mask])
                H[offshore_mask, ci] = scale * vals

        interior_mask = ~offshore_mask
        if interior_mask.any():
            dist = np.array([float(r["distance_km"]) if r["distance_km"] not in
                              ("", "None", None) else 0.0 for r in group])
            cell_interior = np.clip(np.round(dist / DELXI_KM).astype(int), 0, M - 1)
            cell[interior_mask] = cell_interior
            for ci, control in enumerate(CONTROLS):
                key = f"{model_var}__{control}"
                if key not in green:
                    continue
                arr = green[key]   # (365, M)
                idx_doy = doy[interior_mask]
                idx_cell = cell_interior
                H[np.where(interior_mask)[0], ci] = scale * arr[idx_doy, idx_cell]

        weight = _effective_sample_weight(abs_day, cell, offshore_mask, varname)

        H_rows.append(H)
        sig_rows.append(sigma)
        weight_rows.append(weight)

    if not H_rows:
        return np.zeros((0, N_CTRL)), np.zeros(0), np.zeros(0)
    return np.vstack(H_rows), np.concatenate(sig_rows), np.concatenate(weight_rows)


def fisher_info(H, sigma, weight=None):
    """F = sum_i weight_i * (h_i h_i^T) / sigma_i^2."""
    if H.shape[0] == 0:
        return np.zeros((N_CTRL, N_CTRL))
    W = H / sigma[:, None]
    if weight is not None:
        W = W * np.sqrt(weight)[:, None]
    # See tools/rank_observation_value.py's fisher_info docstring: Apple's
    # Accelerate BLAS backend emits spurious RuntimeWarnings on this matmul even
    # when every input/output is fully finite -- suppressed with an isfinite
    # assertion immediately after, so a genuine problem is still caught loudly.
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=RuntimeWarning)
        F = W.T @ W
    assert np.isfinite(F).all(), "fisher_info produced a non-finite entry -- not the known BLAS false-positive, investigate"
    return F


def dfs_and_reduction(F, n=None):
    """n defaults to N_CTRL (F is the full 19x19 Fisher matrix). Pass a smaller n with
    F already sliced to F[:n, :n] to get the DFS as it would have been under a SMALLER
    control set -- used to reproduce, apples-to-apples, what the original 11-control
    (discharge + marine-only) analysis would have scored on this same data, see
    N_MARINE_CTRL below."""
    if n is None:
        n = N_CTRL
    P0 = np.diag([PRIOR_SIGMA_FRAC ** 2] * n)
    P0_inv = np.diag([1.0 / PRIOR_SIGMA_FRAC ** 2] * n)
    P_post = np.linalg.inv(P0_inv + F)
    dfs = n - np.trace(P_post @ P0_inv)
    reduction = np.diag(P0) - np.diag(P_post)
    return float(dfs), reduction


N_MARINE_CTRL = 2 + len(MARINE_SPECIES)   # = 11: discharge + clb_T + 9 marine species --
                                            # exactly CONTROLS[:N_MARINE_CTRL], the full
                                            # control set before upstream (cub_*) controls
                                            # were added. Slicing F to its first 11x11
                                            # block and recomputing DFS reproduces what
                                            # the original FORTE analysis (marine+discharge
                                            # sensitivity only) would have scored on this
                                            # same campaign-window data -- isolating the
                                            # TIMING effect from the upstream-control effect.
N_BOUNDARY_CTRL = N_MARINE_CTRL + 1 + len(UPSTREAM_CONST_SPECIES)   # = 19: all boundary
                                            # controls (marine+discharge+upstream), before
                                            # the 6 reach-local controls were added. Lets
                                            # the ranking isolate the REACH-LOCAL effect the
                                            # same way N_MARINE_CTRL isolates the upstream
                                            # effect -- CONTROLS[:N_BOUNDARY_CTRL] is exactly
                                            # the boundary-only control set.


CHANNEL_SOURCES = [
    ("grab", "grab", None),
    ("moorings", "mooring_inshore", "mooring_inshore"),
    ("moorings", "mooring_offshore", "mooring_offshore"),
    ("drifters", "drifter", None),
    ("airborne", "airborne", None),
    ("satellite", "swot", "swot"),
    ("satellite", "pace", "pace"),
]


def main():
    out_rows = []
    combined_by_channel = {}
    combined_by_channel_uncorrected = {}
    for site in SITES:
        green = _load_green(site)
        if green is None:
            print(f"{site}: SKIPPED -- no Green's functions (run "
                  f"tools/build_green_functions_forte.py first)")
            continue
        M = next(iter(green.values())).shape[1]
        print(f"\n=== {site} ===")
        F_all = np.zeros((N_CTRL, N_CTRL))
        for file_key, channel_label, channel_filter in CHANNEL_SOURCES:
            rows = _read_csv(site, file_key)
            if channel_filter:
                rows = [r for r in rows if r["channel"] == channel_filter]
            H, sigma, weight = build_H_sigma(rows, green, site, M)
            F = fisher_info(H, sigma, weight)
            F_uncorrected = fisher_info(H, sigma, weight=None)
            F_all += F
            dfs, reduction = dfs_and_reduction(F)
            dfs_uncorrected, _ = dfs_and_reduction(F_uncorrected)
            n_eff = float(weight.sum()) if weight.size else 0.0
            top = np.argsort(reduction)[::-1][:3]
            top_str = ", ".join(f"{CONTROLS[i]}({reduction[i]:.4f})" for i in top if reduction[i] > 1e-6)
            print(f"  {channel_label:18s} n_obs_raw={H.shape[0]:7d}  n_eff={n_eff:9.1f}  "
                  f"DFS={dfs:6.3f} (uncorrected {dfs_uncorrected:6.3f})  "
                  f"top-constrained: {top_str or '(none)'}")
            out_rows.append(dict(site=site, channel=channel_label, n_obs_raw=H.shape[0],
                                   n_obs_effective=round(n_eff, 2), dfs=round(dfs, 4),
                                   dfs_uncorrected=round(dfs_uncorrected, 4),
                                   **{f"var_reduction__{c}": round(float(v), 6)
                                      for c, v in zip(CONTROLS, reduction)}))
            combined_by_channel.setdefault(channel_label, np.zeros((N_CTRL, N_CTRL)))
            combined_by_channel[channel_label] += F
            combined_by_channel_uncorrected.setdefault(channel_label, np.zeros((N_CTRL, N_CTRL)))
            combined_by_channel_uncorrected[channel_label] += F_uncorrected
        dfs_all, _ = dfs_and_reduction(F_all)
        print(f"  {'ALL CHANNELS':18s} {'':14s}  DFS={dfs_all:6.3f}")
        out_rows.append(dict(site=site, channel="ALL_COMBINED", n_obs_raw=-1,
                               n_obs_effective=-1, dfs=round(dfs_all, 4)))

    print("\n=== combined across all 3 rivers (DFS is additive across independent river blocks) ===")
    combined_rows = []
    for channel_label, F in combined_by_channel.items():
        dfs, reduction = dfs_and_reduction(F)
        dfs_u, _ = dfs_and_reduction(combined_by_channel_uncorrected[channel_label])
        dfs_marine_only, _ = dfs_and_reduction(F[:N_MARINE_CTRL, :N_MARINE_CTRL], n=N_MARINE_CTRL)
        dfs_boundary_only, _ = dfs_and_reduction(F[:N_BOUNDARY_CTRL, :N_BOUNDARY_CTRL], n=N_BOUNDARY_CTRL)
        print(f"  {channel_label:18s} DFS={dfs:6.3f} (uncorrected {dfs_u:6.3f}, "
              f"marine+discharge-only {dfs_marine_only:6.3f}/{N_MARINE_CTRL}, "
              f"boundary-only {dfs_boundary_only:6.3f}/{N_BOUNDARY_CTRL})")
        combined_rows.append(dict(channel=channel_label, dfs=round(dfs, 4),
                                    dfs_uncorrected=round(dfs_u, 4),
                                    dfs_marine_discharge_only_11ctrl=round(dfs_marine_only, 4),
                                    dfs_boundary_only_19ctrl=round(dfs_boundary_only, 4),
                                    **{f"var_reduction__{c}": round(float(v), 6)
                                       for c, v in zip(CONTROLS, reduction)}))

    OSSE.mkdir(parents=True, exist_ok=True)
    out_path = OSSE / "observation_value_ranking.csv"
    fields = sorted({k for r in out_rows for k in r.keys()})
    with open(out_path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        for r in out_rows:
            w.writerow(r)
    print(f"\nwrote {out_path.relative_to(ROOT)}")

    combined_path = OSSE / "observation_value_ranking_combined.csv"
    combined_fields = sorted({k for r in combined_rows for k in r.keys()})
    with open(combined_path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=combined_fields)
        w.writeheader()
        for r in combined_rows:
            w.writerow(r)
    print(f"wrote {combined_path.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
