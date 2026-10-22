"""
Step 3 (the actual ranking engine) of the Green's-function/observation-value
framework. Combines:
  - the Green's-function sensitivities from tools/build_green_functions.py
    (G_c(doy, cell, var) = dX/d(fractional control), per river)
  - the synthetic observations from tools/build_osse_observations.py
    (docs/osse/<site>_<channel>.csv, with each row's own recorded instrument sigma)
into a Fisher-information / degrees-of-freedom-for-signal (DFS) ranking of which
observation TYPE most reduces uncertainty in the model's boundary/forcing inputs --
the actual research question this whole framework was built to answer.

METHOD (classic linear/Gaussian information filter, in CONTROL space -- 11 controls
per river, never more, so this is cheap regardless of how many observations there
are, unlike working in observation space which would need to invert an
n_observations x n_observations matrix):
  F_channel = sum_i  w_i h_i h_i^T / sigma_i^2     (Fisher information this channel
                                                     adds, n_controls x n_controls)
  P_post    = (P0^-1 + F_channel)^-1               (posterior covariance)
  DFS       = n_controls - trace(P_post @ P0^-1)   (Rodgers 2000 degrees of freedom
                                                     for signal -- 0 = no information,
                                                     n_controls = fully determined)
  var_reduction_per_control = diag(P0) - diag(P_post)

EFFECTIVE-SAMPLE-SIZE WEIGHT w_i (_effective_sample_weight): naively treating every
raw observation row as independent evidence badly overstates the information content
of continuously-sampled channels (mooring_inshore's 6-min cadence, drifters,
airborne's 3 flights/day) -- this was caught by direct inspection of a first version
of this script's results (mooring_inshore was outranking channels that sample far
more DISTINCT model states, driven purely by raw row count). w_i in [0,1] corrects
this with two caps, both grounded in documented facts about this project rather than
an arbitrary smoothing choice: (A) a MECHANICAL cap, applying to every channel -- the
Green's-function arrays are resolved at daily granularity only, so same-day repeat
samples of the same (variable, cell) are mathematically identical observations, not
independent ones; (B) a PHYSICAL cap, mooring_offshore only -- it reads the marine
boundary FORCING file directly, whose real native update rate is ~monthly (ECCO-
Darwin v5, documented in sites/<river>_interannual_full.py), so daily-interpolated
samples within the same ~30-day window are not new information. See
_effective_sample_weight's docstring for the full derivation.

PRIOR: P0 = diag([0.20^2] * 11) for every control -- a single, documented, equal 20%
relative-uncertainty assumption on every boundary input. Not tailored per-control
(would need literature-specific uncertainty estimates this project doesn't have yet)
-- stated plainly as the one free parameter of this whole ranking, since changing it
rescales every DFS number (though not their RELATIVE ordering, which is what
actually answers "which channel is more valuable").

h_i (one row of H per observation): for every channel except mooring_offshore, this
is G_c(doy, cell, mapped_model_var) read directly from the Green's-function arrays,
scaled by the same constant build_osse_observations.py used to derive that
observation from the model state (e.g. the Chl-a:C conversion). mooring_offshore is
a SPECIAL CASE: it observes the marine BOUNDARY FORCING directly, not a model-
propagated field, so its sensitivity to its own control is just that control's
baseline forcing VALUE at that moment (not a Green's function at all) and exactly
zero to every other control -- see module docstring in build_osse_observations.py
for why this channel is conceptually different from every other one.

NOT YET RANKED (documented gaps, not oversights):
  - drifters' "arrival_time_at_mouth_days": a nonlinear, trajectory-integrated
    quantity: its sensitivity would need re-advecting the SAME release through each
    perturbed run's velocity field (the Green's-function FIELD framework above
    doesn't apply directly to a path integral). A real, likely-valuable product
    (flagged as such in docs/osse/osse_sampling_report.pdf) -- just not implemented
    in this first version.
  - moorings' "current_speed": velocity is not one of the state variables
    build_green_functions.py differences (depth/S/T/O2/DIC/... are; velocity isn't,
    since none of the OTHER channels need it) -- would need adding it there first.

CAVEAT INHERITED FROM build_green_perturbations.py: the Green's function is built
from ONE reference year (2006) and reused for 2007/2008/2009 observations by day-of-
year, assuming the sensitivity SHAPE is stationary across different flow regimes.
This is the single biggest approximation in the whole chain -- read before trusting
small DFS differences between channels.

Usage:  python tools/rank_observation_value.py
Writes: docs/osse/observation_value_ranking.csv, prints a ranked summary table.
"""
import csv
import warnings
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent.parent
RUNS = ROOT / "runs" / "osse_green"
OSSE = ROOT / "docs" / "osse"
FORC = ROOT / "forcing"
SITES = ["colville", "kuparuk", "sagavanirktok"]
MARINE_SPECIES = ["S", "DIC", "ALK", "NO3", "NH4", "PO4", "dSi", "O2", "TOC"]
CONTROLS = ["discharge", "clb_T"] + [f"clb_{sp}" for sp in MARINE_SPECIES]
N_CTRL = len(CONTROLS)
DELXI_KM = 0.1
YEAR_START_DAY = {2006: 365, 2007: 730, 2008: 1095, 2009: 1460}
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
    """Baseline (unperturbed) value of the marine boundary forcing for `model_var`
    at absolute day `day_idx` -- used only for the mooring_offshore special case."""
    key = MARINE_FILE_KEY[model_var]
    if model_var == "T":
        fp = FORC / f"{site}_t_marine_interannual_2005-2023.csv"
    else:
        fp = FORC / f"{site}_{key}_marine_interannual_2005-2023.csv"
    if not hasattr(_marine_forcing_value, "_cache"):
        _marine_forcing_value._cache = {}
    cache_key = (site, model_var)
    if cache_key not in _marine_forcing_value._cache:
        _marine_forcing_value._cache[cache_key] = np.genfromtxt(
            open(fp, encoding="utf-8-sig"), delimiter=",")
    vals = _marine_forcing_value._cache[cache_key]
    day_idx = np.clip(day_idx, 0, vals.size - 1)
    return vals[day_idx]


MARINE_NATIVE_DAYS = 30   # ECCO-Darwin v5 marine boundary: "monthly native resolution
                          # interpolated to daily" -- see sites/<river>_interannual_
                          # full.py's own docstring. Daily-interpolated values between
                          # two real monthly data points carry no new independent
                          # information, so mooring_offshore (which reads that
                          # forcing file directly) is capped to >=1 effective sample
                          # per this many days.


def _effective_sample_weight(abs_day, cell, is_offshore, varname):
    """Per-row weight in [0,1] correcting for duplicate/highly-correlated samples,
    so that N raw samples sharing the same (varname, cell, bin) contribute exactly
    ONE effective sample's worth of Fisher information, not N. Two independent caps:

    (A) MECHANICAL, applies to every row (interior AND offshore): the Green's-
    function arrays this script reads (tools/build_green_functions.py) are resolved
    at DAILY granularity only -- two observations of the same (varname, cell) on the
    same calendar day get the LITERALLY IDENTICAL H row (verified: build_green_
    functions.py's g_by_doy has exactly one row per day-of-year). Sub-daily sampling
    (mooring's 6-min cadence, drifters, airborne's 3 flights/day) therefore cannot
    add information this linear framework can detect beyond one sample/day -- this
    is a fact about how H was built, not a physical assumption, and bins by the
    integer absolute day (NOT day-of-year, so different years are never conflated).

    (B) PHYSICAL, offshore-only: mooring_offshore reads the marine boundary FORCING
    file directly rather than a model-propagated field (see build_osse_observations.
    py's docstring), and that forcing's real native update rate is ~monthly (see
    MARINE_NATIVE_DAYS above) -- daily-interpolated values within the same ~30-day
    window are not independent new information, so offshore rows are additionally
    capped to >=1 effective sample per MARINE_NATIVE_DAYS, on top of cap (A).

    Grab/satellite rows are essentially unaffected (their own real sampling cadence
    already exceeds these windows), verified in practice: this function returns 1.0
    for the overwhelming majority of their rows."""
    n = len(abs_day)
    weight = np.ones(n)
    day_int = np.round(abs_day).astype(int)

    # Cap (A): group by (is_offshore-aware cell key, day)
    keys_a = list(zip(cell.tolist(), day_int.tolist()))
    counts_a = {}
    for k in keys_a:
        counts_a[k] = counts_a.get(k, 0) + 1
    weight_a = np.array([1.0 / counts_a[k] for k in keys_a])

    # Cap (B): offshore rows only, 30-day bins (re-weight the ALREADY-capped-to-1/day
    # offshore group further, so the two caps compose correctly: cap A alone already
    # reduces N raw daily-duplicated offshore rows to 1 effective/day; cap B then
    # further reduces up to MARINE_NATIVE_DAYS of those daily effective samples to 1).
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
    _effective_sample_weight -- see that function's docstring."""
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
        years = np.array([int(r["year"]) for r in group])
        t_days = np.array([float(r["t_days"]) for r in group])
        base_days = np.array([YEAR_START_DAY[y] for y in years])
        doy = np.clip(np.round(t_days - base_days).astype(int), 0, 364)
        abs_day = base_days + doy
        sigma = np.abs(np.array([float(r["sigma"]) for r in group]))
        channels = [r["channel"] for r in group]

        offshore_mask = np.array([c == "mooring_offshore" for c in channels])
        H = np.zeros((len(group), N_CTRL))
        cell = np.full(len(group), -1)   # -1 = offshore pseudo-cell (no spatial index)

        if offshore_mask.any():
            ctrl = VAR_TO_MARINE_CONTROL.get(model_var)
            if ctrl is not None:
                ci = CONTROLS.index(ctrl)
                abs_days_off = base_days[offshore_mask] + doy[offshore_mask]
                vals = _marine_forcing_value(site, model_var, abs_days_off)
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
    """F = sum_i weight_i * (h_i h_i^T) / sigma_i^2 -- see _effective_sample_weight
    for what `weight` corrects for. weight=None (or all-ones) reduces to the plain,
    uncorrected Fisher information."""
    if H.shape[0] == 0:
        return np.zeros((N_CTRL, N_CTRL))
    W = H / sigma[:, None]
    if weight is not None:
        W = W * np.sqrt(weight)[:, None]   # sum_i (sqrt(w_i) W_i)(sqrt(w_i) W_i)^T
                                             # = sum_i w_i W_i W_i^T, as intended
    # Apple's Accelerate BLAS backend (this numpy build's blas/lapack -- confirmed via
    # np.show_config()) emits spurious "divide by zero"/"overflow"/"invalid value"
    # RuntimeWarnings on this matmul even when every input AND the actual output are
    # fully finite -- verified directly (H, sigma, W all finite, no zero sigma
    # entries, W.T@W itself a sane finite value) before concluding this is a known
    # platform BLAS false-positive, not a real numerical problem. Suppressed ONLY
    # around this exact call, with an explicit isfinite assertion immediately after
    # so a genuine future non-finite result is still caught loudly, not hidden.
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=RuntimeWarning)
        F = W.T @ W
    assert np.isfinite(F).all(), "fisher_info produced a non-finite entry -- not the known BLAS false-positive, investigate"
    return F


def dfs_and_reduction(F):
    P0 = np.diag([PRIOR_SIGMA_FRAC ** 2] * N_CTRL)
    P0_inv = np.diag([1.0 / PRIOR_SIGMA_FRAC ** 2] * N_CTRL)
    P_post = np.linalg.inv(P0_inv + F)
    dfs = N_CTRL - np.trace(P_post @ P0_inv)
    reduction = np.diag(P0) - np.diag(P_post)
    return float(dfs), reduction


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
                  f"tools/build_green_functions.py first)")
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
            F_uncorrected = fisher_info(H, sigma, weight=None)   # for before/after
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
        print(f"  {channel_label:18s} DFS={dfs:6.3f} (uncorrected {dfs_u:6.3f})")
        combined_rows.append(dict(channel=channel_label, dfs=round(dfs, 4),
                                    dfs_uncorrected=round(dfs_u, 4),
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
