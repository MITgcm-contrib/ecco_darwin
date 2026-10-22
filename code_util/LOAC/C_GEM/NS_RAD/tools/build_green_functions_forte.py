"""
FORTE-campaign-framework counterpart to tools/build_green_functions.py: differences
the REUSED runs/forte_campaign/<site>/output.nc baseline against each perturbed run
from tools/build_green_perturbations_forte.py (runs/forte_green/pert_<control>/
<site>/output.nc) into a per-control, per-variable sensitivity field
G_c(doy, cell) = dX/d(fractional control), saved to
runs/forte_green/green_functions_<site>.npz. See tools/build_green_functions.py's
docstring for the full method (finite-difference, not a true adjoint); only what
differs is called out below. runs/osse_green/ is completely untouched.

WHAT'S GENUINELY DIFFERENT:
  - Baseline is runs/forte_campaign/<site>/output.nc (the already-built
    climatological run), NOT runs/forte_green/baseline/<site>/ -- there is no
    separate baseline run for this framework, by design (see
    build_green_perturbations_forte.py's docstring).
  - 25 controls, not 19: discharge, clb_T, 9 marine (clb_*) boundary species, 8
    upstream (cub_*) controls, PLUS -- new, 2026-10-20 -- 6 REACH-LOCAL controls
    (rate_aerdeg_mouth/mid/head, lateral_toc_mouth/mid/head) testing sensitivity to
    spatially-distributed in-channel unknowns rather than just the two boundaries --
    see build_green_perturbations_forte.py's docstring for the full design.
  - The 3 lateral_toc_* controls diff against a DIFFERENT, separate baseline
    (runs/forte_green/lateral_base/<site>/output.nc, not runs/forte_campaign/) --
    see LATERAL_CONTROLS below and build_green_perturbations_forte.py's
    "methodological note" for why (a uniform background lateral load must already
    be present in the reference run for a local +10% BOOST to be a valid finite
    difference). Every other control still diffs against runs/forte_campaign/.
  - The reference window is day-of-year 0 through CAMPAIGN_END_DAY=273 (273 days),
    NOT a full post-warmup year (0-364) -- this run was never simulated past
    day-of-year 273 (tools/run_forte_campaign.sh's MAXT_DAYS=638 = 365+273), since
    nothing in this framework's observations samples beyond the campaign season.
    The output arrays are still shape (365, M) for drop-in compatibility with the
    ranking script's indexing -- days 274-364 are simply left at zero and never
    queried (every PHASES/MOORING_WINDOW day-of-year in build_osse_observations_
    forte.py is <= 273).

Usage:  python tools/build_green_functions_forte.py
Writes: runs/forte_green/green_functions_<site>.npz (one per river)
"""
from pathlib import Path

import numpy as np
from netCDF4 import Dataset

ROOT = Path(__file__).resolve().parent.parent
BASELINE_RUNS = ROOT / "runs" / "forte_campaign"
RUNS = ROOT / "runs" / "forte_green"
LATERAL_BASELINE_RUNS = RUNS / "lateral_base"
SITES = ["colville", "kuparuk", "sagavanirktok"]
MARINE_SPECIES = ["S", "DIC", "ALK", "NO3", "NH4", "PO4", "dSi", "O2", "TOC"]
UPSTREAM_CONST_SPECIES = ["NO3", "NH4", "PO4", "dSi", "DIC", "ALK", "DIA"]
REACHES = ["mouth", "mid", "head"]
LOCAL_CONTROLS = ([f"rate_aerdeg_{r}" for r in REACHES]
                  + [f"lateral_toc_{r}" for r in REACHES])
LATERAL_CONTROLS = {f"lateral_toc_{r}" for r in REACHES}  # diff against the SEPARATE
                                                           # lateral baseline, not
                                                           # runs/forte_campaign/
CONTROLS = (["discharge", "clb_T"] + [f"clb_{sp}" for sp in MARINE_SPECIES]
            + ["cub_TOC"] + [f"cub_{sp}" for sp in UPSTREAM_CONST_SPECIES]
            + LOCAL_CONTROLS)
VARS = ["S", "T", "O2", "DIC", "ALK", "NO3", "NH4", "PO4", "dSi", "TOC", "DIA", "SPM", "depth"]
PERT_FRAC = 0.10
WARMUP_DAYS = 365
CAMPAIGN_END_DAY = 273   # matches tools/run_forte_campaign.sh's MAXT_DAYS=638


def _load(path, var):
    d = Dataset(path)
    try:
        if var not in d.variables:
            return None, None
        t = np.asarray(d.variables["time"][:]) / 86400.0
        F = np.asarray(d.variables[var][:])[:, 1:]
        return t, F
    finally:
        # Explicit close -- under memory pressure, relying on GC to release
        # netCDF4's C-level buffers (rather than Python's own array memory)
        # was leaving multi-hundred-MB files open longer than necessary across
        # 247 (var, control) reads per site.
        d.close()


def _reference_window(t, F):
    """(doy[0..273], F_sliced) for the campaign-season window of the single
    climatological run -- NOT a full year (see module docstring)."""
    mask = (t >= WARMUP_DAYS) & (t <= WARMUP_DAYS + CAMPAIGN_END_DAY)
    tt, Fy = t[mask], F[mask]
    doy = np.round(tt - WARMUP_DAYS).astype(int)
    doy = np.clip(doy, 0, CAMPAIGN_END_DAY)
    return doy, Fy


def build_site(site):
    base_path = BASELINE_RUNS / site / "output.nc"
    if not base_path.exists():
        print(f"{site}: SKIPPED -- {base_path} not found (run "
              f"tools/run_forte_campaign.sh first)")
        return
    lateral_base_path = LATERAL_BASELINE_RUNS / site / "output.nc"
    out = {}
    for var in VARS:
        t, F = _load(base_path, var)
        if F is None:
            continue
        doy, base_y = _reference_window(t, F)

        # lateral_toc_* controls diff against the SEPARATE lateral-baseline run
        # (see module docstring / LATERAL_CONTROLS above), not runs/forte_campaign/.
        base_y_lateral = None
        if lateral_base_path.exists():
            _, F_lat = _load(lateral_base_path, var)
            if F_lat is not None:
                _, base_y_lateral = _reference_window(t, F_lat)

        for control in CONTROLS:
            pert_path = RUNS / f"pert_{control}" / site / "output.nc"
            if not pert_path.exists():
                continue
            _, Fp = _load(pert_path, var)
            if Fp is None:
                continue
            _, pert_y = _reference_window(t, Fp)

            if control in LATERAL_CONTROLS:
                if base_y_lateral is None:
                    continue
                ref_y = base_y_lateral
            else:
                ref_y = base_y

            n = min(len(ref_y), len(pert_y))
            delta = pert_y[:n] - ref_y[:n]
            g = delta / PERT_FRAC   # dX / d(fractional control), one-sided
            # 365-row array (not 274) for drop-in indexing compatibility with the
            # ranking script -- days beyond CAMPAIGN_END_DAY stay zero, never queried.
            g_by_doy = np.zeros((365, g.shape[1]))
            counts = np.zeros(365)
            for i, d_ in enumerate(doy[:n]):
                g_by_doy[d_] += g[i]
                counts[d_] += 1
            counts[counts == 0] = 1
            g_by_doy /= counts[:, None]
            out[f"{var}__{control}"] = g_by_doy
    npz_path = RUNS / f"green_functions_{site}.npz"
    RUNS.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(npz_path, **out)
    print(f"  wrote {npz_path.relative_to(ROOT)} ({len(out)} (var,control) sensitivity fields)")


def main():
    for site in SITES:
        build_site(site)


if __name__ == "__main__":
    main()
