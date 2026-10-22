"""
Synthetic observations for an OSSE (Observing System Simulation Experiment) built
around the FORTE Expedition field-campaign plan (Google Doc "Document_Input_
fromTeam.docx" / "FORTE Expedition Sampling Plan - Field Campaign Summary"), sampling
the 6-min-cadence "nature run" tools/run_osse_nature.sh produces
(runs/osse_nature/<site>/output.nc) as if it were the real rivers. This is the
observation-generation half of the observation-value-ranking framework discussed
with the project owner (representer/Green's-function sensitivity + degrees-of-
freedom-for-signal ranking is the OTHER half, deferred to a follow-up tool that
consumes this one's output); this script's job is ONLY to produce realistic,
error-bearing synthetic observation sets for each of the campaign's real instrument
types, nothing more.

WHY AN ANALOG YEAR, NOT THE REAL CAMPAIGN YEAR. The real FORTE campaign is 2027
(Phase 0 "maybe 2028"); NS-RAD's forcing does not extend past 2023, and the model
has no restart/checkpoint capability (confirmed by reading main.py), so there is no
cheap way to reach a later year without simulating every day before it at the same
6-min cadence. tools/run_osse_nature.sh instead runs 2007/2008/2009 -- the cleanest
consecutive low/average/high flow-regime triplet in the record for all three rivers
simultaneously (checked directly against the discharge forcing: Colville mean Q
198/330/438 m3/s, Kuparuk 4.7/19.0/66.6, Sagavanirktok 31/77/172) -- and maps the
campaign's real day-of-year phase windows onto each of those three years. This
answers "how would a campaign shaped like this perform across different flow years,"
not literally "what would 2027 look like."

CAMPAIGN PHASES USED (day-of-year, non-leap; Phase 0 excluded -- the source document
itself marks it "TBD (consider maybe for 2028)" and infeasible due to negligible
under-ice flow):
    Phase 1a  Jun 1-5    (doy 152-156)  helicopter grab sampling, pre-breakup
    Phase 1b  Jun 6-20   (doy 157-171)  helicopter+boat grab sampling, post-breakup
    Phase 2   Jul 1-12   (doy 182-193)  boats+R/V Ukpik, AVIRIS-3, moorings deployed
    Phase 3   Jul 20-Aug 7 (doy 201-219) boats+R/V Ukpik+OceanXplorer, AVIRIS-3+drone,
                                         offshore moorings

FIVE OBSERVATION CHANNELS (the document's instrument list regrouped by sampling
geometry, which is what actually matters for an observation operator):

1. GRAB SAMPLING (`sample_grab`) -- Phases 1a/1b/2/3, 5 along-channel stations/river
   (mouth- and upstream-adjacent at k=2 cells in from the boundary-forced cell, plus
   3 evenly-spaced interior stations), 3 sampling days per phase (first/mid/last).
   Model-variable mapping, and what is deliberately NOT synthesized:
     DOC          -> TOC (+-5%)
     CDOM         -> TOC proxy again, wider error (+-15%) -- NS-RAD has no CDOM state;
                     this is an explicit, documented proxy, not a new variable
     NH4/NO3/PO4/dSi -> direct state vars (+-5%)
     DIC          -> direct (+-2%, lab-analytical precision)
     ALK          -> direct (+-2%) -- not explicit in the source doc's chemistry
                     list, but any DIC-measuring campaign also gets ALK for carbonate
                     system closure; a documented inference, not confabulation
     Chl-a/HPLC   -> DIA via a C:Chl = 50:1 (by weight) conversion (+-10%), a commonly
                     cited average marine ratio -- NOT a 1:1 identity
     POC/PON      -> SKIPPED. No NS-RAD state variable is a defensible proxy (DIA is
                     autotrophic carbon only, conflates with detrital POC if reused)
     flow cytometry, productivity incubations -> SKIPPED. Process-rate/community-
                     composition measurements with no analogous model state.
     SPM/turbidity -> NOT sampled here (reserved for the airborne channel below,
                     since the source doc lists turbidity-relevant work under
                     AVIRIS/drone, not the grab-sample chemistry list)

2. MOORINGS (`sample_moorings`) -- continuous, deployed Jul 1 - Sep 30 (doy 182-273)
   each study year (moorings go in at Phase 2, come out before ice-up -- deployed
   "for long-term measurements" per the source doc, not year-round under ice).
     INSHORE: an interior model cell (mouth-adjacent, k=2) -- continuous S, T, O2,
       current speed (|velocity|), native 6-min cadence, sonde-grade noise (T +-0.2
       degC, S +-0.1 psu, O2 +-5%, current +-0.02 m/s). A genuine interior-model
       check.
     OFFSHORE: NOT an interior model cell -- samples the marine BOUNDARY FORCING
       file itself (T/S/O2 clb series, native daily cadence) at the same window.
       This validates the ASSUMED boundary input, not interior model skill -- a
       conceptually different role from the inshore mooring, kept as a separate
       "channel" in the output for that reason.

3. DRIFTERS (`sample_drifters`) -- microSWIFT-style Lagrangian tracers, 30 total per
   study year (the source doc's "qty = 30 drifters in 2027"), released in batches at
   the start of each phase (10/10/5/5 for 1a/1b/2/3 -- front-loaded onto the freshet,
   the most dynamically interesting transport period), advected by the model's own
   native 6-min velocity field (simple explicit advection at the save cadence, a
   deliberate simplification of the model's actual 30s DELTI integration -- adequate
   for a synthetic/illustrative tracer, not claimed as a numerically converged
   trajectory), tracked until exit at the mouth or 10 days elapsed. Samples T/S along
   the trajectory (+-0.3 degC/+-0.2 psu, a mobile sensor's worse precision than a
   fixed mooring) -- NOT biogeochemistry (real microSWIFTs don't carry wet chemistry).
   Each drifter's arrival time at the mouth is also recorded as a derived scalar --
   likely a strong, cheap constraint on bulk advective/dispersion timescales, flagged
   here for the follow-up sensitivity-ranking tool to weigh accordingly.

4. AIRBORNE OPTICAL (`sample_airborne`) -- AVIRIS-3 (King Air B-200) + drone, Phase 2
   and Phase 3 ONLY (the source doc does not fly either in Phase 1; Phase 1b's drone
   work is ice-elevation structure-from-motion, a different product not synthesized
   here). Full longitudinal snapshot (every grid cell) every 8 hours throughout each
   phase window -- optical-only variables (TOC/CDOM proxy +-15%, DIA/Chl proxy
   +-20%, SPM +-15%), no T/S/nutrients (realistic: these are reflectance-based
   retrievals). An 8-hour cadence does not divide evenly into any of this coast's
   dominant tidal periods (M2 12.42 h, S2 12.00 h, O1 25.82 h, K1 23.93 h -- see
   forcing/tidal_constituents.json), so sampled tidal phase keeps advancing rather
   than re-locking; Phase 2 (12 days) and Phase 3 (19 days) are both comparable to or
   longer than the ~14.8-day M2/S2 spring-neap beat, so good spring-neap diversity
   falls out of the real campaign windows on their own -- no special-case date
   picking was needed for that (an explicit design simplification versus the
   original 5-day-burst sketch this script supersedes).

5. SATELLITE, notional (`sample_satellite`) -- SWOT (water depth as a stand-in for
   its WSE-anomaly product) and PACE (surface TOC/DIA optical proxies), kept as a
   "what if this were coordinated" channel per the project owner's explicit choice,
   since the source document only loosely implies satellite coordination, not a
   committed data source. GATED BY PHYSICAL LIMITS ONLY, not real launch dates
   (SWOT Dec 2022, PACE Feb 2024) -- the whole analog-year substitution above
   already breaks literal calendar correspondence (2007-2009 stands in for
   2027-2028), so re-imposing literal launch-date gating on top would make the
   satellite channel spuriously empty for a date-mismatch reason instead of the
   physical reasons (ice, cloud, channel width) this OSSE is meant to probe.
     SWOT: every cell with width >= 100 m (KaRIn's realistic usable-width floor;
       checked directly against this project's own geometry -- Kuparuk's upstream
       prismatic width is 58 m, essentially unresolvable, while Colville's mouth
       ~1052 m and Sagavanirktok's ~225 m are plausible), ice_frac < 0.5, ~21-day
       repeat cadence across all 3 study years. Width-dependent noise (10 cm at
       100 m, improving as sqrt(width) up to a 3 cm floor, reflecting reach-
       averaging) -- a simple heuristic, not a mission error budget.
     PACE: outermost 1-2 cells only (the mouth/plume, the only reach optically
       visible from space in a narrow shaded channel), near-daily cadence,
       ice_frac < 0.5 AND a 50% random cloud dropout (a documented simplifying
       stand-in for real cloud climatology, not an actual cloud product; seeded for
       reproducibility), TOC/DIA proxies at +-25% (coarser than the airborne
       channel's resolution/revisit).

OUTPUT: one CSV per (site, channel) under docs/osse/<site>_<channel>.csv, long format
(year, phase, t_days, distance_km or NaN, variable, value, sigma, plus channel-
specific extra columns e.g. drifter_id/is_offshore). Reads runs/osse_nature/<site>/
output.nc (built by tools/run_osse_nature.sh -- run that first) and the marine-
boundary forcing CSVs (for the offshore mooring). No model run of its own.

Usage:  python tools/build_osse_observations.py
"""
import csv
import math
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent.parent
NATURE = ROOT / "runs" / "osse_nature"
FORC = ROOT / "forcing"
OUT = ROOT / "docs" / "osse"
SITES = ["colville", "kuparuk", "sagavanirktok"]

STUDY_YEARS = [2007, 2008, 2009]          # low / average / high flow analogs
YEAR_START_DAY = {2006: 365, 2007: 730, 2008: 1095, 2009: 1460}  # see WARMUP_DAYS=365

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
    """Memoised (t_days, F[:, 1:]) from runs/osse_nature/<site>/output.nc."""
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


def _doy_mask(t, lo, hi, year):
    base = YEAR_START_DAY[year]
    return (t >= base + lo) & (t <= base + hi)


def _add_rows(rows, site, channel, year, phase, t_days, distance_km, variable, value, sigma, **extra):
    noisy = float(value + RNG.normal(0.0, sigma)) if sigma else float(value)
    row = dict(site=site, channel=channel, year=year, phase=phase, t_days=round(float(t_days), 4),
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
    for year in STUDY_YEARS:
        for phase, (lo, hi) in PHASES.items():
            if phase == "phase1b":
                pass  # grab sampling continues in 1b; the ice-focused drone work in
                      # 1b is a separate, not-synthesized product (see module docstring)
            sample_days = sorted({lo, (lo + hi) // 2, hi})
            for sd in sample_days:
                target = YEAR_START_DAY[year] + sd
                idx = int(np.argmin(np.abs(t - target)))
                for st in stations:
                    for nc_var, out_name, sig, is_frac, _ in GRAB_VARS:
                        _, F = load(site, nc_var)
                        if F is None:
                            continue
                        val = float(F[idx, st])
                        sigma = abs(val) * sig if is_frac else sig
                        _add_rows(rows, site, "grab", year, phase, t[idx], st * delxi,
                                   out_name, val, sigma, station=int(st))
                    _, Fdia = load(site, "DIA")
                    if Fdia is not None:
                        chl = float(Fdia[idx, st]) * 12.011 / C_TO_CHL
                        _add_rows(rows, site, "grab", year, phase, t[idx], st * delxi,
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
    for year in STUDY_YEARS:
        mask = _doy_mask(t, *MOORING_WINDOW, year)
        tt = t[mask]
        for nc_var, out_name, sig, is_frac in MOORING_VARS:
            _, F = load(site, nc_var)
            if F is None:
                continue
            vals = F[mask, inshore_idx]
            for ti, v in zip(tt, vals):
                sigma = abs(v) * sig if is_frac else sig
                _add_rows(rows, site, "mooring_inshore", year, "mooring", ti,
                           inshore_idx * delxi, out_name, v, sigma)
        _, Vel = load(site, "velocity")
        if Vel is not None:
            speed = np.abs(Vel[mask, inshore_idx])
            for ti, v in zip(tt, speed):
                _add_rows(rows, site, "mooring_inshore", year, "mooring", ti,
                           inshore_idx * delxi, "current_speed", v, 0.02)

        # Offshore: the marine boundary FORCING itself, not an interior cell --
        # validates the assumed boundary input (see module docstring).
        for nc_var, out_name, sig, is_frac in MOORING_VARS:
            key = _MARINE_FILE_KEY[nc_var]
            fp = FORC / f"{site}_{key}_marine_interannual_2005-2023.csv"
            if not fp.exists():
                continue
            vals = np.genfromtxt(open(fp, encoding="utf-8-sig"), delimiter=",")
            base_day = YEAR_START_DAY[year]
            lo, hi = MOORING_WINDOW
            for d in range(lo, hi + 1):
                row_i = base_day + d   # forcing CSVs are 0-indexed daily rows from 2005-01-01
                if row_i >= vals.size:
                    continue
                v = float(vals[row_i])
                sigma = abs(v) * sig if is_frac else sig
                _add_rows(rows, site, "mooring_offshore", year, "mooring", float(row_i),
                           None, out_name, v, sigma)
    return rows


# =====================================================================================
# 3. microSWIFT drifters (Lagrangian)
# =====================================================================================
DRIFTER_RELEASE = {"phase1a": 10, "phase1b": 10, "phase2": 5, "phase3": 5}  # of 30/year
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
    for year in STUDY_YEARS:
        for phase, n in DRIFTER_RELEASE.items():
            lo, hi = PHASES[phase]
            release_day = YEAR_START_DAY[year] + lo
            idx0 = int(np.argmin(np.abs(t - release_day)))
            release_x = np.linspace(0.6 * (M - 1) * delxi, (M - 1) * delxi, n)
            for x0 in release_x:
                did += 1
                samples, arrival = _advect_one(t, Vel, S, T, delxi, M, idx0, x0)
                for ti, xi, sv, tv in samples:
                    _add_rows(rows, site, "drifter", year, phase, ti, xi, "salinity",
                              sv, 0.2, drifter_id=did)
                    _add_rows(rows, site, "drifter", year, phase, ti, xi, "temperature",
                              tv, 0.3, drifter_id=did)
                if arrival is not None:
                    _add_rows(rows, site, "drifter", year, phase, t[idx0] + arrival, 0.0,
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
    for year in STUDY_YEARS:
        for phase in AIRBORNE_PHASES:
            lo, hi = PHASES[phase]
            base = YEAR_START_DAY[year]
            n_flights = int((hi - lo) * 24.0 / AIRBORNE_CADENCE_H)
            flight_days = base + lo + np.arange(n_flights) * (AIRBORNE_CADENCE_H / 24.0)
            for fd in flight_days:
                idx = int(np.argmin(np.abs(t - fd)))
                for nc_var, out_name, sig in (("TOC", "CDOM_proxy", 0.15),
                                               ("SPM", "SPM", 0.15)):
                    _, F = load(site, nc_var)
                    if F is None:
                        continue
                    vals = F[idx, :]
                    for ci, v in enumerate(vals):
                        _add_rows(rows, site, "airborne", year, phase, t[idx], ci * delxi,
                                   out_name, float(v), abs(float(v)) * sig)
                _, Fdia = load(site, "DIA")
                if Fdia is not None:
                    for ci, v in enumerate(Fdia[idx, :]):
                        chl = float(v) * 12.011 / C_TO_CHL
                        _add_rows(rows, site, "airborne", year, phase, t[idx], ci * delxi,
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

    # SWOT: width-masked, ice-masked, ~21-day repeat
    w0 = W[0]
    swot_cells = np.where(w0 >= SWOT_MIN_WIDTH_M)[0]
    for year in STUDY_YEARS:
        mask_year = (t >= YEAR_START_DAY[year]) & (t < YEAR_START_DAY[year] + 365)
        year_idx = np.where(mask_year)[0]
        for i in year_idx[::int(SWOT_REPEAT_DAYS * 1440 / 6)]:  # ~21 days @ 6-min steps
            if Ice[i, swot_cells].mean() >= 0.5:
                continue
            for ci in swot_cells:
                if Ice[i, ci] >= 0.5:
                    continue
                width = w0[ci]
                sigma_m = max(0.03, 0.10 * math.sqrt(SWOT_MIN_WIDTH_M / width))
                _add_rows(rows, site, "swot", year, "satellite", t[i], ci * delxi,
                           "depth_anomaly", float(Dep[i, ci]), sigma_m)

    # PACE: outermost 1-2 cells, near-daily, ice + random cloud dropout
    pace_cells = [0, 1] if Dep.shape[1] > 1 else [0]
    for year in STUDY_YEARS:
        mask_year = (t >= YEAR_START_DAY[year]) & (t < YEAR_START_DAY[year] + 365)
        year_idx = np.where(mask_year)[0]
        daily_idx = year_idx[::240]  # 240 * 6min = 1 day
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
                    _add_rows(rows, site, "pace", year, "satellite", t[i], ci * delxi,
                               "CDOM_proxy", v, abs(v) * 0.25)
                if Fdia is not None:
                    chl = float(Fdia[i, ci]) * 12.011 / C_TO_CHL
                    _add_rows(rows, site, "pace", year, "satellite", t[i], ci * delxi,
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
            print(f"{site}: SKIPPED -- runs/osse_nature/{site}/output.nc not found "
                  f"(run tools/run_osse_nature.sh first)")
            continue
        print(f"=== {site} ===")
        for name, fn in channels:
            rows = fn(site)
            _write_csv(rows, OUT / f"{site}_{name}.csv")


if __name__ == "__main__":
    main()
