"""
Step 1 of the Green's-function/observation-value-ranking engine (the second half of
the FORTE OSSE framework -- tools/build_osse_observations.py built the first half,
the synthetic observations themselves).

Generates one BASELINE and one PERTURBED 2-year run per (river, control), differenced
in tools/build_green_functions.py to approximate the model's LINEAR sensitivity to
each boundary/forcing input -- a finite-difference stand-in for a true adjoint
Green's function, since NS-RAD has none (it is a ported legacy Fortran model, not
auto-differentiable).

CONTROLS (11 per river): discharge, sea temperature (SEATEMP_FILE), and the marine
(`clb`) boundary for every genuinely time-varying solute in the full-forcing variant
(S, DIC, ALK, NO3, NH4, PO4, dSi, O2, TOC -- see sites/<river>_interannual_full.py's
BOUNDARY_FORCING dict). The corresponding UPSTREAM (`cub`) boundaries are NOT
perturbed here -- every solute except TOC holds its upstream boundary constant even
in the full-forcing variant (see CLAUDE.md -> "Interannual forcing"), so there is no
time-varying input there to be uncertain about in the first place.

PERTURBATION: a single ONE-SIDED +10% relative step applied to the WHOLE forcing
series (not a true delta-impulse) -- this answers "if this boundary input were
systematically biased 10% high, how would the model state respond," which is the
right question for ranking which OBSERVATIONS best catch a persistent forcing
misspecification (the actual research question), not a transient-event response.
One-sided (not central-difference) to keep the run count at 36 instead of 69 --
adequate for a ranking exercise, not claimed as a high-precision derivative.

WHY A SHORT 2-YEAR RUN, NOT THE FULL 3-ANALOG-YEAR NATURE RUN. NS-RAD has no
restart/checkpoint capability (confirmed when building tools/run_osse_nature.sh), so
reaching 2007/2008/2009 requires simulating every day before them regardless of
cadence. A full-length, 6-min-cadence PERTURBATION campaign (36 runs instead of 3)
would take many hours. Instead, each perturbation run is MAXT_DAYS=730 (the model's
TOTAL run length, NOT additive with WARMUP_DAYS -- confirmed in config.py: MAXT is
the full simulated span, WARMUP is just the within-that-span point biogeochemistry
turns on) at the fast DAILY cadence (TS=2880, same as tools/run_interannual_full.sh)
-- cheap (~minutes, not hours) because Green's functions only need to capture the
SEASONAL SHAPE of the model's sensitivity, not tidal-scale dynamics. The first 365
days (2005) are discarded as warmup; the second 365 (2006) is the single reference
year, because it is the first full post-warmup year reachable at all -- i.e. the
cheapest possible choice given the no-restart constraint.

APPROXIMATION THIS INTRODUCES (stated here once, relied on throughout the rest of
this engine): the resulting sensitivity is assumed STATIONARY IN DAY-OF-YEAR across
different calendar years -- i.e., "how a +10% DIC bias in May propagates" is assumed
similar in 2006, 2007, 2008, and 2009, even though the actual flow regime differs
across those three analog years (that is the whole reason they were chosen). This is
a real limitation of a single-reference-year sensitivity, not a free simplification,
and is the single most important caveat for interpreting the ranking results tool
(tools/rank_observation_value.py) produces. A fully rigorous version would compute a
SEPARATE reference run per analog year -- not done here due to the same no-restart
cost problem that made 2007/2008/2009 themselves expensive to reach in the first
place.

Generated scratch artifacts (regenerate with this script; not meant to be hand-
edited or committed):
  forcing/_green_scratch/<site>_<control>_pert.csv   (perturbed forcing copies)
  code/sites/_green_<site>_<control>.py              (site variants overriding
                                                        exactly one forcing file)
  runs/osse_green/baseline/<site>/output.nc
  runs/osse_green/pert_<control>/<site>/output.nc

Usage:  python tools/build_green_perturbations.py          # all 3 rivers x 11 controls
        python tools/build_green_perturbations.py colville # just one river
"""
import os
import subprocess
import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent.parent
CODE = ROOT / "code"
FORC = ROOT / "forcing"
SCRATCH = FORC / "_green_scratch"
SITES_DIR = CODE / "sites"
RUNS = ROOT / "runs" / "osse_green"

SITES = ["colville", "kuparuk", "sagavanirktok"]
MARINE_SPECIES = ["S", "DIC", "ALK", "NO3", "NH4", "PO4", "dSi", "O2", "TOC"]
PERT_FRAC = 0.10
WARMUP_DAYS = 365
MAXT_DAYS = 730   # TOTAL run length (MAXT_DAYS is not additive with WARMUP_DAYS --
                  # confirmed in config.py: MAXT = total seconds simulated, WARMUP is
                  # just the within-that-range threshold after which biogeochemistry
                  # turns on). 730 = 365 warmup + 365 actual post-warmup study year
                  # (2006) -- getting this wrong the first time (MAXT_DAYS=365, i.e.
                  # a run that is ENTIRELY warmup) is exactly why the first version
                  # of this script produced all-zero Green's functions.
TS = 2880   # daily, matches tools/run_interannual_full.sh
MAX_PARALLEL = 8
NS_RAD_PYTHON = "/opt/homebrew/anaconda3/envs/ns_rad/bin/python3"


def _site_boundary_forcing(site):
    """(DISCHARGE_FILE, SEATEMP_FILE, BOUNDARY_FORCING) for <site>_interannual_full,
    read via the real config import (not re-typed here) so this script can never
    drift from the site files."""
    for m in list(sys.modules):
        if m == "config" or m.startswith("sites") or m == "sites":
            del sys.modules[m]
    sys.path.insert(0, str(CODE))
    os.environ["CGEM_SITE"] = f"{site}_interannual_full"
    import config as cfg
    return cfg.DISCHARGE_FILE, cfg.SEATEMP_FILE, cfg.BOUNDARY_FORCING


def _perturb_csv(src_name, dst_path):
    vals = np.genfromtxt(open(FORC / src_name, encoding="utf-8-sig"), delimiter=",")
    vals_pert = vals * (1.0 + PERT_FRAC)
    np.savetxt(dst_path, vals_pert, delimiter=",", fmt="%.6f")


def controls_for(site):
    """[(control_name, override_kind, override_value), ...] for one river."""
    discharge_file, seatemp_file, boundary = _site_boundary_forcing(site)
    out = [("discharge", "DISCHARGE_FILE", discharge_file),
           ("clb_T", "SEATEMP_FILE", seatemp_file)]
    for sp in MARINE_SPECIES:
        if sp in boundary and "clb" in boundary[sp]:
            out.append((f"clb_{sp}", ("BOUNDARY_FORCING", sp, "clb"), boundary[sp]["clb"]))
    return out


def write_perturbed_forcing_and_site(site, control_name, override_kind, src_fname):
    SCRATCH.mkdir(parents=True, exist_ok=True)
    pert_fname = f"_green_scratch/{site}_{control_name}_pert.csv"
    _perturb_csv(src_fname, FORC / pert_fname)

    site_mod_name = f"_green_{site}_{control_name}"
    site_path = SITES_DIR / f"{site_mod_name}.py"
    if override_kind in ("DISCHARGE_FILE", "SEATEMP_FILE"):
        body = (
            f"# AUTO-GENERATED by tools/build_green_perturbations.py -- do not hand-edit;\n"
            f"# regenerate by rerunning that script. +{PERT_FRAC*100:.0f}% {override_kind} "
            f"perturbation for the Green's-function sensitivity engine.\n"
            f"from .{site}_interannual_full import *  # noqa: F401,F403\n"
            f"{override_kind} = {pert_fname!r}\n"
        )
    else:
        _, species, end = override_kind
        body = (
            f"# AUTO-GENERATED by tools/build_green_perturbations.py -- do not hand-edit;\n"
            f"# regenerate by rerunning that script. +{PERT_FRAC*100:.0f}% {species} "
            f"{end} perturbation for the Green's-function sensitivity engine.\n"
            f"from .{site}_interannual_full import *  # noqa: F401,F403\n"
            f"BOUNDARY_FORCING = {{k: dict(v) for k, v in BOUNDARY_FORCING.items()}}\n"
            f"BOUNDARY_FORCING[{species!r}][{end!r}] = {pert_fname!r}\n"
        )
    site_path.write_text(body)
    return site_mod_name


def run_one(site_variant, outdir):
    outdir.mkdir(parents=True, exist_ok=True)
    env = dict(os.environ)
    env.update(CGEM_SITE=site_variant, CGEM_MAXT_DAYS=str(MAXT_DAYS),
               CGEM_WARMUP_DAYS=str(WARMUP_DAYS), CGEM_TS=str(TS),
               PYTHONPATH=str(CODE), PYTHONWARNINGS="default",
               PATH=f"/opt/homebrew/anaconda3/envs/ns_rad/bin:{os.environ.get('PATH','')}")
    log = open(outdir / "run.log", "w")
    return subprocess.Popen([NS_RAD_PYTHON, str(CODE / "main.py")], cwd=outdir,
                             env=env, stdout=log, stderr=subprocess.STDOUT)


def main():
    sites = sys.argv[1:] or SITES
    jobs = []   # (label, site_variant, outdir, proc_or_None, logfile)
    for site in sites:
        jobs.append((f"{site} baseline", f"{site}_interannual_full", RUNS / "baseline" / site))
        for control_name, override_kind, src_fname in controls_for(site):
            site_variant = write_perturbed_forcing_and_site(site, control_name, override_kind, src_fname)
            jobs.append((f"{site} {control_name}", site_variant, RUNS / f"pert_{control_name}" / site))

    print(f"{len(jobs)} runs to launch ({MAXT_DAYS}+{WARMUP_DAYS} days @ daily cadence each)")
    running = []
    pending = list(jobs)
    t_start = time.time()
    while pending or running:
        while pending and len(running) < MAX_PARALLEL:
            label, variant, outdir = pending.pop(0)
            print(f"  -> launching {label}")
            proc = run_one(variant, outdir)
            running.append((label, proc))
        time.sleep(5)
        still = []
        for label, proc in running:
            if proc.poll() is None:
                still.append((label, proc))
            else:
                status = "OK" if proc.returncode == 0 else f"FAILED (rc={proc.returncode})"
                print(f"  <- {label}: {status} ({time.time()-t_start:.0f}s elapsed)")
        running = still
    print("all Green's-function perturbation runs finished")


if __name__ == "__main__":
    main()
