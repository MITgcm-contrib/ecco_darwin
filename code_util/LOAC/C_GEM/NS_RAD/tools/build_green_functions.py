"""
Step 2 of the Green's-function/observation-value-ranking engine: differences each
baseline/perturbed run pair from tools/build_green_perturbations.py into a per-
control, per-variable sensitivity field G_c(doy, cell) = dX/d(fractional control),
saved to runs/osse_green/green_functions_<site>.npz.

G_c is a FINITE-DIFFERENCE sensitivity (ΔX / PERT_FRAC, one-sided), not a true
adjoint Green's function -- see build_green_perturbations.py's docstring for why, and
for the "stationary in day-of-year across analog years" approximation this script's
output and tools/rank_observation_value.py both rely on (G_c is built from ONE
reference year, 2006, then reused for 2007/2008/2009 observations by day-of-year).

VARIABLES: S, T, O2, DIC, ALK, NO3, NH4, PO4, dSi, TOC, DIA, SPM, depth -- every
state variable any of the 5 synthetic observation channels actually samples (see
tools/build_osse_observations.py). CDOM_proxy/Chl_a_proxy observations reuse
G_TOC/G_DIA directly (they are linear transforms of those two fields, not separate
model state) -- handled in the ranking script, not duplicated here.

Usage:  python tools/build_green_functions.py
Writes: runs/osse_green/green_functions_<site>.npz (one per river)
"""
from pathlib import Path

import numpy as np
from netCDF4 import Dataset

ROOT = Path(__file__).resolve().parent.parent
RUNS = ROOT / "runs" / "osse_green"
SITES = ["colville", "kuparuk", "sagavanirktok"]
MARINE_SPECIES = ["S", "DIC", "ALK", "NO3", "NH4", "PO4", "dSi", "O2", "TOC"]
CONTROLS = ["discharge", "clb_T"] + [f"clb_{sp}" for sp in MARINE_SPECIES]
VARS = ["S", "T", "O2", "DIC", "ALK", "NO3", "NH4", "PO4", "dSi", "TOC", "DIA", "SPM", "depth"]
PERT_FRAC = 0.10
WARMUP_DAYS = 365


def _load(path, var):
    d = Dataset(path)
    if var not in d.variables:
        return None, None
    t = np.asarray(d.variables["time"][:]) / 86400.0
    F = np.asarray(d.variables[var][:])[:, 1:]
    return t, F


def _reference_year(t, F):
    """(doy[0..364], F_sliced) for the single post-warmup reference year (2006)."""
    mask = (t >= WARMUP_DAYS) & (t < WARMUP_DAYS + 365)
    tt, Fy = t[mask], F[mask]
    doy = np.round(tt - WARMUP_DAYS).astype(int)
    doy = np.clip(doy, 0, 364)
    return doy, Fy


def build_site(site):
    base_dir = RUNS / "baseline" / site / "output.nc"
    if not base_dir.exists():
        print(f"{site}: SKIPPED -- {base_dir} not found (run "
              f"tools/build_green_perturbations.py first)")
        return
    out = {}
    for var in VARS:
        t, F = _load(base_dir, var)
        if F is None:
            continue
        doy, base_y = _reference_year(t, F)
        for control in CONTROLS:
            pert_path = RUNS / f"pert_{control}" / site / "output.nc"
            if not pert_path.exists():
                continue
            _, Fp = _load(pert_path, var)
            if Fp is None:
                continue
            _, pert_y = _reference_year(t, Fp)
            n = min(len(base_y), len(pert_y))
            delta = pert_y[:n] - base_y[:n]
            g = delta / PERT_FRAC   # dX / d(fractional control), one-sided
            # Collapse to exactly one row per day-of-year (daily cadence already
            # gives ~1 row/day; guard against an off-by-one duplicate row).
            g_by_doy = np.zeros((365, g.shape[1]))
            counts = np.zeros(365)
            for i, d_ in enumerate(doy[:n]):
                g_by_doy[d_] += g[i]
                counts[d_] += 1
            counts[counts == 0] = 1
            g_by_doy /= counts[:, None]
            out[f"{var}__{control}"] = g_by_doy
    npz_path = RUNS / f"green_functions_{site}.npz"
    np.savez_compressed(npz_path, **out)
    print(f"  wrote {npz_path.relative_to(ROOT)} ({len(out)} (var,control) sensitivity fields)")


def main():
    for site in SITES:
        build_site(site)


if __name__ == "__main__":
    main()
