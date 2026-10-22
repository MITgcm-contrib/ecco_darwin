"""
FORTE-campaign-framework counterpart to tools/build_osse_observations.py: samples
the CLIMATOLOGICAL "FORTE campaign year" run (tools/run_forte_campaign.sh,
runs/forte_campaign/<site>/output.nc) instead of the 2007/2008/2009 analog-year
nature run (runs/osse_nature/) that script uses -- see sites/<river>_forte_
campaign.py's docstring for why a climatology, not an analog year, and CLAUDE.md
for the full reasoning trail. runs/osse_nature/ and its own docs/osse/*.csv are
completely untouched by this file's existence; this is a separate, parallel
pipeline, matching this project's "never modify a validated original, add a
parallel opt-in variant" convention used throughout.

WHAT'S GENUINELY DIFFERENT HERE, NOT JUST A PATH SWAP:
  - ONE climatological pass, not three analog years -- there is only a single
    "year" to sample (the campaign window sits in the run's second climatological
    year, absolute day 365-638), so every STUDY_YEARS loop in the original
    collapses to a single pass. Sample COUNTS are correspondingly ~3x smaller than
    the original's per-channel counts (one campaign season instead of three) --
    an honest consequence of this being a single representative year, not a
    reduction in realism.
  - The offshore mooring's marine-boundary read points at the CLIMATOLOGICAL marine
    CSVs (tools/build_forte_campaign_forcing.py's *_marine_climatology.csv, 365
    rows, one representative year) rather than the 19-year *_marine_interannual_
    2005-2023.csv series -- indexed by DAY-OF-YEAR directly (`d % 365`), not by an
    absolute day into a multi-year file, since a climatology file has no absolute-
    day concept at all (repeatYear cycles the same 365 rows every model year).
  - Satellite sampling's per-"year" window is the actual simulated span (post-
    warmup day 365 through day 638, i.e. 273 days = through day-of-year 273, the
    mooring-recovery end of the campaign season) rather than a full 365-day year
    -- this run was never simulated past day-of-year 273, unlike the original's
    full-year nature run.

CAMPAIGN PHASES USED (day-of-year, non-leap; Phase 0 excluded) -- identical to
the original script, see tools/build_osse_observations.py's docstring for the
full per-channel instrument/variable/error-model documentation (grab sampling,
moorings, drifters, airborne optical, notional satellite): that reasoning is
NOT repeated here, only what changes with the single-climatological-year design
is called out above.

OUTPUT: one CSV per (site, channel) under docs/osse_forte/<site>_<channel>.csv,
same long format as the original (a fixed "forte" year label instead of a real
calendar year).

Usage:  python tools/build_osse_observations_forte.py
"""
import csv
import math
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent.parent
NATURE = ROOT / "runs" / "forte_campaign"
FORC = ROOT / "forcing"
OUT = ROOT / "docs" / "osse_forte"
SITES = ["colville", "kuparuk", "sagavanirktok"]

YEAR_LABEL = "forte"
BASE_DAY = 365     # WARMUP_DAYS -- the campaign window sits in the run's 2nd
                   # climatological year, absolute day 365 + doy
CAMPAIGN_END_DAY = 273   # last day-of-year this run was ever simulated to (see
                         # tools/run_forte_campaign.sh's MAXT_DAYS=638 = 365+273)

# Campaign phase windows, day-of-year (non-leap). Phase 0 (under-ice, late April) is
# excluded -- the source document itself marks it infeasible/TBD.
PHASES = {
    "phase1a": (152, 156),   # Jun 1-5,  pre-breakup grab sampling
    "phase1b": (157, 171),   # Jun 6-20, post-breakup grab sampling + ice drone (not synthesized)
    "phase2":  (182, 193),   # Jul 1-12, AVIRIS-3 + moorings deployed
    "phase3":  (201, 219),   # Jul 20-Aug 7, AVIRIS-3+drone + offshore moorings
}
MOORING_WINDOW = (182, 273)   # Jul 1 - Sep 30: deployed at Phase 2, out before ice-up

C_TO_CHL = 50.0   # C:Chl-a mass ratio, commonly cited marine average (NOT 1:1)
RNG = np.random.default_rng(20270101)   # seeded for reproducibility; date is arbitrary

_CACHE = {}


def load(site, var):
    """Memoised (t_days, F[:, 1:]) from runs/forte_campaign/<site>/output.nc."""
    key = (site, var)
    if key in _CACHE:
        return _CACHE[key]
    ncp = NATURE / site / "output.nc"
    if not ncp.exists():
        _CACHE[key] = (None, None)
        return _CACHE[key]
    from netCDF4 import Dataset
    d = Dataset(ncp)
    if var not in d.variables:
        _CACHE[key] = (None, None)
        return _CACHE[key]
    t = np.asarray(d.variables["time"][:]) / 86400.0
    F = np.asarray(d.variables[var][:])
    _CACHE[key] = (t, F[:, 1:])
    return _CACHE[key]


def _delxi_km(site):
    import importlib, sys, os, warnings
    p = str(ROOT / "code")
    if p not in sys.path:
        sys.path.insert(0, p)
    for m in ("config", "sites"):
        for mod in list(sys.modules):
            if mod == m or mod.startswith(m + "."):
                del sys.modules[mod]
    os.environ["CGEM_SITE"] = site
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        import config as _cfg
    return _cfg.DELXI / 1000.0


def _stations(M):
    """5 along-channel station indices (0-based, cell 0=mouth .. M-1=upstream): the
    mouth- and upstream-adjacent cells 2 in from each boundary-forced cell (avoiding
    the openbound relaxation zone, same convention as the publication-figures
    filtering mass balance), plus 3 evenly-spaced interior stations."""
    lo, hi = 2, M - 3
    return sorted(set(int(round(x)) for x in np.linspace(lo, hi, 5)))


def _doy_mask(t, lo, hi):
    return (t >= BASE_DAY + lo) & (t <= BASE_DAY + hi)


def _add_rows(rows, site, channel, phase, t_days, distance_km, variable, value, sigma, **extra):
    noisy = float(value + RNG.normal(0.0, sigma)) if sigma else float(value)
    row = dict(site=site, channel=channel, year=YEAR_LABEL, phase=phase,
               t_days=round(float(t_days), 4),
               distance_km=None if distance_km is None else round(float(distance_km), 3),
               variable=variable, value=noisy, sigma=sigma)
    row.update(extra)
    rows.append(row)


# =====================================================================================
# 1. Grab sampling
# =====================================================================================
GRAB_VARS = [  # (nc_var, out_name, sigma_frac_or_abs, is_fraction, unit_note)
    ("TOC", "DOC", 0.05, True, "mmol C m-3"),
    ("TOC", "CDOM_proxy", 0.15, True, "mmol C m-3 (TOC proxy)"),
    ("NH4", "NH4", 0.05, True, "mmol N m-3"),
    ("NO3", "NO3", 0.05, True, "mmol N m-3"),
    ("PO4", "PO4", 0.05, True, "mmol P m-3"),
    ("dSi", "dSi", 0.05, True, "mmol Si m-3"),
    ("DIC", "DIC", 0.02, True, "mmol C m-3"),
    ("ALK", "ALK", 0.02, True, "mmol eq m-3"),
]


def sample_grab(site):
    t, _ = load(site, "S")
    if t is None:
        return []
    M = load(site, "S")[1].shape[1]
    delxi = _delxi_km(site)
    stations = _stations(M)
    rows = []
    for phase, (lo, hi) in PHASES.items():
        sample_days = sorted({lo, (lo + hi) // 2, hi})
        for sd in sample_days:
            target = BASE_DAY + sd
            idx = int(np.argmin(np.abs(t - target)))
            for st in stations:
                for nc_var, out_name, sig, is_frac, _ in GRAB_VARS:
                    _, F = load(site, nc_var)
                    if F is None:
                        continue
                    val = float(F[idx, st])
                    sigma = abs(val) * sig if is_frac else sig
                    _add_rows(rows, site, "grab", phase, t[idx], st * delxi,
                               out_name, val, sigma, station=int(st))
                _, Fdia = load(site, "DIA")
                if Fdia is not None:
                    chl = float(Fdia[idx, st]) * 12.011 / C_TO_CHL
                    _add_rows(rows, site, "grab", phase, t[idx], st * delxi,
                               "Chl_a_mg_m3", chl, abs(chl) * 0.10, station=int(st))
    return rows


# =====================================================================================
# 2. Moorings (inshore = interior model cell; offshore = marine boundary forcing)
# =====================================================================================
MOORING_VARS = [("S", "salinity", 0.1, False), ("T", "temperature", 0.2, False),
                 ("O2", "O2", 0.05, True)]
_MARINE_FILE_KEY = {"S": "s", "T": "t", "O2": "o2"}


def sample_moorings(site):
    t, _ = load(site, "S")
    if t is None:
        return []
    delxi = _delxi_km(site)
    inshore_idx = 2
    rows = []
    mask = _doy_mask(t, *MOORING_WINDOW)
    tt = t[mask]
    for nc_var, out_name, sig, is_frac in MOORING_VARS:
        _, F = load(site, nc_var)
        if F is None:
            continue
        vals = F[mask, inshore_idx]
        for ti, v in zip(tt, vals):
            sigma = abs(v) * sig if is_frac else sig
            _add_rows(rows, site, "mooring_inshore", "mooring", ti,
                       inshore_idx * delxi, out_name, v, sigma)
    _, Vel = load(site, "velocity")
    if Vel is not None:
        speed = np.abs(Vel[mask, inshore_idx])
        for ti, v in zip(tt, speed):
            _add_rows(rows, site, "mooring_inshore", "mooring", ti,
                       inshore_idx * delxi, "current_speed", v, 0.02)

    # Offshore: the marine boundary FORCING itself, not an interior cell --
    # validates the assumed boundary input (see module docstring). Reads the
    # CLIMATOLOGICAL marine CSV (365 rows, one representative year) by DAY-OF-YEAR
    # directly -- there is no absolute-day concept in a repeating single-year file.
    for nc_var, out_name, sig, is_frac in MOORING_VARS:
        key = _MARINE_FILE_KEY[nc_var]
        fp = FORC / f"{site}_{key}_marine_climatology.csv"
        if not fp.exists():
            continue
        vals = np.genfromtxt(open(fp, encoding="utf-8-sig"), delimiter=",")
        lo, hi = MOORING_WINDOW
        for d in range(lo, hi + 1):
            row_i = d % 365
            v = float(vals[row_i])
            sigma = abs(v) * sig if is_frac else sig
            _add_rows(rows, site, "mooring_offshore", "mooring", float(BASE_DAY + d),
                       None, out_name, v, sigma)
    return rows


# =====================================================================================
# 3. microSWIFT drifters (Lagrangian)
# =====================================================================================
DRIFTER_RELEASE = {"phase1a": 10, "phase1b": 10, "phase2": 5, "phase3": 5}  # of 30 total
DRIFTER_MAX_DAYS = 10.0


def _advect_one(t, Vel, S, T, delxi, M, release_day_idx, x0_km):
    """Explicit advection at the nature run's native (6-min) save cadence -- a
    deliberate simplification of the model's own 30s DELTI integration, adequate for
    an illustrative synthetic tracer (see module docstring). Returns a list of
    (t_days, x_km, S_val, T_val) samples and the arrival time at the mouth (or None
    if still in-domain after DRIFTER_MAX_DAYS)."""
    samples = []
    x = x0_km
    i = release_day_idx
    t0 = t[i]
    arrival = None
    while i < len(t) - 1 and (t[i] - t0) <= DRIFTER_MAX_DAYS:
        cell = x / delxi
        lo_c = int(np.clip(np.floor(cell), 0, M - 1))
        hi_c = int(np.clip(lo_c + 1, 0, M - 1))
        frac = np.clip(cell - lo_c, 0.0, 1.0)
        u = Vel[i, lo_c] * (1 - frac) + Vel[i, hi_c] * frac
        s_val = S[i, lo_c] * (1 - frac) + S[i, hi_c] * frac
        t_val = T[i, lo_c] * (1 - frac) + T[i, hi_c] * frac
        samples.append((t[i], x, s_val, t_val))
        dt_s = (t[i + 1] - t[i]) * 86400.0
        x = x + u * dt_s / 1000.0   # u [m/s] -> km
        if x <= 0.0:
            arrival = t[i + 1] - t0
            break
        if x >= (M - 1) * delxi:
            x = (M - 1) * delxi
        i += 1
    return samples, arrival


def sample_drifters(site):
    t, Vel = load(site, "velocity")
    _, S = load(site, "S")
    _, T = load(site, "T")
    if t is None or Vel is None:
        return []
    M = Vel.shape[1]
    delxi = _delxi_km(site)
    rows = []
    did = 0
    for phase, n in DRIFTER_RELEASE.items():
        lo, hi = PHASES[phase]
        release_day = BASE_DAY + lo
        idx0 = int(np.argmin(np.abs(t - release_day)))
        release_x = np.linspace(0.6 * (M - 1) * delxi, (M - 1) * delxi, n)
        for x0 in release_x:
            did += 1
            samples, arrival = _advect_one(t, Vel, S, T, delxi, M, idx0, x0)
            for ti, xi, sv, tv in samples:
                _add_rows(rows, site, "drifter", phase, ti, xi, "salinity",
                          sv, 0.2, drifter_id=did)
                _add_rows(rows, site, "drifter", phase, ti, xi, "temperature",
                          tv, 0.3, drifter_id=did)
            if arrival is not None:
                _add_rows(rows, site, "drifter", phase, t[idx0] + arrival, 0.0,
                          "arrival_time_at_mouth_days", arrival, None, drifter_id=did)
    return rows


# =====================================================================================
# 4. Airborne optical (AVIRIS-3 + drone), Phase 2/3 only
# =====================================================================================
AIRBORNE_PHASES = ["phase2", "phase3"]
AIRBORNE_CADENCE_H = 8.0


def sample_airborne(site):
    t, _ = load(site, "S")
    if t is None:
        return []
    delxi = _delxi_km(site)
    rows = []
    for phase in AIRBORNE_PHASES:
        lo, hi = PHASES[phase]
        n_flights = int((hi - lo) * 24.0 / AIRBORNE_CADENCE_H)
        flight_days = BASE_DAY + lo + np.arange(n_flights) * (AIRBORNE_CADENCE_H / 24.0)
        for fd in flight_days:
            idx = int(np.argmin(np.abs(t - fd)))
            for nc_var, out_name, sig in (("TOC", "CDOM_proxy", 0.15),
                                           ("SPM", "SPM", 0.15)):
                _, F = load(site, nc_var)
                if F is None:
                    continue
                vals = F[idx, :]
                for ci, v in enumerate(vals):
                    _add_rows(rows, site, "airborne", phase, t[idx], ci * delxi,
                               out_name, float(v), abs(float(v)) * sig)
            _, Fdia = load(site, "DIA")
            if Fdia is not None:
                for ci, v in enumerate(Fdia[idx, :]):
                    chl = float(v) * 12.011 / C_TO_CHL
                    _add_rows(rows, site, "airborne", phase, t[idx], ci * delxi,
                               "Chl_a_proxy_mg_m3", chl, abs(chl) * 0.20)
    return rows


# =====================================================================================
# 5. Satellite, notional (SWOT depth-anomaly + PACE optical), physically gated only
# =====================================================================================
SWOT_MIN_WIDTH_M = 100.0
SWOT_REPEAT_DAYS = 21
PACE_CLOUD_DROPOUT_PROB = 0.5


def sample_satellite(site):
    t, Dep = load(site, "depth")
    _, W = load(site, "width")
    _, Ice = load(site, "ice_frac")
    if t is None:
        return []
    delxi = _delxi_km(site)
    rows = []

    # Only one climatological pass -- the window is whatever this run actually
    # simulated past warmup (through day-of-year CAMPAIGN_END_DAY), not a full year.
    mask_window = (t >= BASE_DAY) & (t <= BASE_DAY + CAMPAIGN_END_DAY)
    window_idx = np.where(mask_window)[0]

    # SWOT: width-masked, ice-masked, ~21-day repeat
    w0 = W[0]
    swot_cells = np.where(w0 >= SWOT_MIN_WIDTH_M)[0]
    for i in window_idx[::int(SWOT_REPEAT_DAYS * 1440 / 6)]:  # ~21 days @ 6-min steps
        if Ice[i, swot_cells].mean() >= 0.5:
            continue
        for ci in swot_cells:
            if Ice[i, ci] >= 0.5:
                continue
            width = w0[ci]
            sigma_m = max(0.03, 0.10 * math.sqrt(SWOT_MIN_WIDTH_M / width))
            _add_rows(rows, site, "swot", "satellite", t[i], ci * delxi,
                       "depth_anomaly", float(Dep[i, ci]), sigma_m)

    # PACE: outermost 1-2 cells, near-daily, ice + random cloud dropout
    pace_cells = [0, 1] if Dep.shape[1] > 1 else [0]
    daily_idx = window_idx[::240]  # 240 * 6min = 1 day
    for i in daily_idx:
        if RNG.random() < PACE_CLOUD_DROPOUT_PROB:
            continue
        for ci in pace_cells:
            if Ice[i, ci] >= 0.5:
                continue
            _, Ftoc = load(site, "TOC")
            _, Fdia = load(site, "DIA")
            if Ftoc is not None:
                v = float(Ftoc[i, ci])
                _add_rows(rows, site, "pace", "satellite", t[i], ci * delxi,
                           "CDOM_proxy", v, abs(v) * 0.25)
            if Fdia is not None:
                chl = float(Fdia[i, ci]) * 12.011 / C_TO_CHL
                _add_rows(rows, site, "pace", "satellite", t[i], ci * delxi,
                           "Chl_a_proxy_mg_m3", chl, abs(chl) * 0.25)
    return rows


# =====================================================================================
def _write_csv(rows, path):
    if not rows:
        return
    OUT.mkdir(parents=True, exist_ok=True)
    fields = sorted({k for r in rows for k in r.keys()})
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        for r in rows:
            w.writerow(r)
    print(f"  wrote {path.relative_to(ROOT)} ({len(rows)} rows)")


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    channels = [("grab", sample_grab), ("moorings", sample_moorings),
                ("drifters", sample_drifters), ("airborne", sample_airborne),
                ("satellite", sample_satellite)]
    for site in SITES:
        if not (NATURE / site / "output.nc").exists():
            print(f"{site}: SKIPPED -- runs/forte_campaign/{site}/output.nc not found "
                  f"(run tools/run_forte_campaign.sh first)")
            continue
        print(f"=== {site} ===")
        for name, fn in channels:
            rows = fn(site)
            _write_csv(rows, OUT / f"{site}_{name}.csv")


if __name__ == "__main__":
    main()
