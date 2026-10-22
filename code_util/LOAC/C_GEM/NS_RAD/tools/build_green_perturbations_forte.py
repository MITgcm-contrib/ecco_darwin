"""
FORTE-campaign-framework counterpart to tools/build_green_perturbations.py: the
SAME one-sided +10% finite-difference perturbation design, but built on top of the
CLIMATOLOGICAL FORTE campaign run (sites/<river>_forte_campaign.py,
runs/forte_campaign/<site>/output.nc) instead of the *_interannual_full analog-year
baseline -- see tools/build_green_perturbations.py's docstring for the full
method/caveat reasoning (one-sided vs. central difference, the day-of-year-
stationarity assumption), not repeated here; only what differs is called out below.
runs/osse_green/ (the original campaign) is completely untouched.

25 CONTROLS per river (up from 19): discharge, the 9 marine (`clb`) boundaries, 8
UPSTREAM (`cub`) boundaries (TOC file-based; NO3, NH4, PO4, dSi, DIC, ALK, DIA as
plain CONSTANT overrides -- these have no time-varying source in the model yet, so
the perturbation tests "if the true upstream value differs 10% from the constant we
currently assume, how much does that matter," exactly what the campaign's own
upstream grab sampling will eventually let someone answer with real data; S and O2
excluded since nothing in the campaign plan samples them AT the upstream boundary),
and -- NEW, added 2026-10-20 at the project owner's explicit request -- 6 REACH-
LOCAL controls testing sensitivity to spatially-distributed (not just boundary)
unknowns, since continuous spatial observations (airborne/drone maps, grab sampling
at 5 stations along each estuary) also carry information about what's happening
INSIDE the domain:
  - `rate_aerdeg_mouth/mid/head`: a local +10% perturbation to the aerobic-degradation
    rate constant (kox), confined to one reach third (config.KOX_MULT_REACH, a new
    per-cell multiplier array consumed in biogeo_module.py's kernel). Baseline =
    the existing runs/forte_campaign/ (kox_mult is all-ones there by default).
  - `lateral_toc_mouth/mid/head`: a local +10% BOOST to a shared, uniform, domain-wide
    background lateral TOC load (config.LATERAL_BOOST_REACH). This control does
    NOT reuse the forte_campaign baseline -- see "methodological note" below.
Reach boundaries are coarse thirds of each river's own grid (b1=M//3, b2=2*(M//3);
mouth=(1,b1), mid=(b1+1,b2), head=(b2+1,M) -- cell 1 is the seaward/mouth end, cell M
the upstream/headwater end, same convention as distance_km/DELXI_KM indexing
elsewhere in this pipeline) -- per the project owner's explicit choice of
coarse-reach over per-station or per-grid-cell resolution.

METHODOLOGICAL NOTE on lateral_toc_*: a naive on/off version of this control
(perturbing LATERAL_INFLOW from its real value of 0 directly to something, confined
to one reach) would NOT be a valid finite difference -- dividing a 0-to-on jump by
PERT_FRAC is not "dX/d(10% nudge)", and would inflate this control's apparent
Fisher-information value relative to every other control, which IS a genuine ±10%
nudge around an always-on baseline. Fixed by introducing a dedicated, UNPERTURBED
lateral-baseline run per site (a uniform, domain-wide background lateral TOC load,
magnitude/concentration reusing sites/idealized.py's own precedent: ~15% of mean
discharge, 300.0 mmol C/m3) and perturbing it with a local +10% BOOST confined to
one reach -- a clean, dimensionless nudge fully consistent with every other control
and with rank_observation_value_forte.py's PRIOR_SIGMA_FRAC=0.20. This costs one
extra run per site (written to runs/forte_green/lateral_base/<site>/, NOT a
"control" itself -- not present in CONTROLS downstream).

WHAT'S GENUINELY DIFFERENT FROM tools/build_green_perturbations.py:
  - BASELINE IS REUSED, NOT RE-RUN. runs/forte_campaign/<site>/output.nc (built by
    tools/run_forte_campaign.sh) IS ALREADY the unperturbed climatological run this
    script needs as a baseline -- so unlike the original (which had to run its own
    3 baseline jobs), this script generates and runs ONLY the PERTURBED jobs,
    reading the baseline directly from the existing run.
  - MAXT_DAYS=638, not 730 -- matches tools/run_forte_campaign.sh (365 warmup + 273
    days of the campaign season, not a full post-warmup year), since there is no
    reason to perturb-and-simulate days this framework's observations never sample
    in the first place.
  - Site variants inherit from sites/<river>_forte_campaign.py (climatological
    forcing), not <river>_interannual_full.py (19-year real forcing) -- so the
    FILE-based perturbed forcing copies are built from the 365-row CLIMATOLOGY csvs
    (tools/build_forte_campaign_forcing.py), not the 6935-row interannual series.
  - The 7 upstream CONSTANT controls need no forcing file at all -- they override
    the site's BOUNDARIES dict entry directly (a plain Python tuple of (clb, cub)
    constants; see sites/colville.py), a genuinely different mechanism from every
    other control here, which perturbs a csv file.

Generated scratch artifacts (regenerate with this script; not meant to be hand-
edited or committed), using a DISTINCT prefix from the original campaign's scratch
files so the two can never collide:
  forcing/_fortegreen_scratch/<site>_<control>_pert.csv   (file-based controls only)
  code/sites/_fortegreen_<site>_<control>.py
  code/sites/_fortegreen_<site>_lateral_base.py           (shared baseline, 1/site)
  runs/forte_green/pert_<control>/<site>/output.nc
  runs/forte_green/lateral_base/<site>/output.nc          (shared baseline, 1/site)

Re-running this script SKIPS any (site, control) whose run already completed
(checked via its run.log's last line) -- safe to re-run after adding new controls
without redoing the original ones.

Usage:  python tools/build_green_perturbations_forte.py          # all 3 rivers x 25 controls
        python tools/build_green_perturbations_forte.py colville # just one river
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
SCRATCH = FORC / "_fortegreen_scratch"
SITES_DIR = CODE / "sites"
RUNS = ROOT / "runs" / "forte_green"

SITES = ["colville", "kuparuk", "sagavanirktok"]
MARINE_SPECIES = ["S", "DIC", "ALK", "NO3", "NH4", "PO4", "dSi", "O2", "TOC"]
# Upstream (cub) species perturbed as plain BOUNDARIES constants (no forcing file
# exists for these yet) -- matches exactly what the campaign's grab sampling
# measures at the upstream station. TOC is handled separately below (it DOES
# already have a cub forcing file, like the marine species).
UPSTREAM_CONST_SPECIES = ["NO3", "NH4", "PO4", "dSi", "DIC", "ALK", "DIA"]
# Reach-local controls' shared lateral-baseline magnitude -- reuses sites/idealized.py's
# own precedent ratio/concentration exactly, for consistency with that already-cited
# choice (LATERAL_INFLOW = 8.0 m3/s against idealized's own ~53 m3/s mean, i.e. ~15%).
LATERAL_BACKGROUND_FRAC = 0.15
LATERAL_TOC_CONC = 300.0
PERT_FRAC = 0.10
WARMUP_DAYS = 365
MAXT_DAYS = 638   # matches tools/run_forte_campaign.sh -- 365 warmup + 273 days of
                  # the campaign season (NOT a full post-warmup year; this
                  # framework's observations never sample past day-of-year 273)
TS = 12           # 6-min, matches tools/run_forte_campaign.sh (NOT the original
                  # Green's-function campaign's daily TS=2880 -- that one's baseline
                  # was a cheap NEW daily-cadence run built just for sensitivity;
                  # this one's baseline is REUSED from the already-built 6-min
                  # runs/forte_campaign/, so the perturbed runs must match its
                  # cadence to be differenced against it)
MAX_PARALLEL = 1  # lowered from 8, then 4, then 2 -- all three were OOM-killed under
                  # sustained host memory pressure from other running apps (not
                  # specific to this pipeline); fully serial is slowest but most
                  # resilient
NS_RAD_PYTHON = "/opt/homebrew/anaconda3/envs/ns_rad/bin/python3"


def _site_boundary_forcing(site):
    """(DISCHARGE_FILE, SEATEMP_FILE, BOUNDARY_FORCING, M) for <site>_forte_campaign,
    read via the real config import (not re-typed here) so this script can never
    drift from the site files."""
    for m in list(sys.modules):
        if m == "config" or m.startswith("sites") or m == "sites":
            del sys.modules[m]
    sys.path.insert(0, str(CODE))
    os.environ["CGEM_SITE"] = f"{site}_forte_campaign"
    import config as cfg
    return cfg.DISCHARGE_FILE, cfg.SEATEMP_FILE, cfg.BOUNDARY_FORCING, cfg.M


def _perturb_csv(src_name, dst_path):
    vals = np.genfromtxt(open(FORC / src_name, encoding="utf-8-sig"), delimiter=",")
    vals_pert = vals * (1.0 + PERT_FRAC)
    np.savetxt(dst_path, vals_pert, delimiter=",", fmt="%.6f")


def _reach_bounds(M):
    """{'mouth'|'mid'|'head': (i_start, i_end)} 1-indexed grid-cell range, coarse
    thirds of this river's own grid -- the project owner's explicit choice of
    spatial resolution for the reach-local controls (see module docstring). Cell 1
    is the seaward/mouth (clb) end, cell M the upstream/headwater (cub) end (same
    convention as distance_km/DELXI_KM indexing throughout the ranking pipeline) --
    'mouth' is therefore the LOWEST-index third, 'head' the HIGHEST, not the other
    way around."""
    b1, b2 = M // 3, 2 * (M // 3)
    return {"mouth": (1, b1), "mid": (b1 + 1, b2), "head": (b2 + 1, M)}


def _mean_discharge(site):
    vals = np.genfromtxt(open(FORC / f"{site}_river_discharge_climatology_m3sec.csv",
                               encoding="utf-8-sig"), delimiter=",")
    return float(vals.mean())


def write_lateral_base_site(site):
    """Write (or overwrite) the shared, UNPERTURBED lateral-baseline site variant
    that every lateral_toc_mouth/mid/head control inherits from -- a uniform, domain-
    wide background lateral TOC load, not present in runs/forte_campaign/. Exists
    only so the reach-local +10% boost is a genuine nudge around an already-nonzero
    value (see module docstring's methodological note). Not itself a "control" --
    not added to CONTROLS downstream."""
    site_mod_name = f"_fortegreen_{site}_lateral_base"
    site_path = SITES_DIR / f"{site_mod_name}.py"
    inflow = LATERAL_BACKGROUND_FRAC * _mean_discharge(site)
    body = (
        f"# AUTO-GENERATED by tools/build_green_perturbations_forte.py -- do not\n"
        f"# hand-edit; regenerate by rerunning that script. Shared, UNPERTURBED\n"
        f"# baseline for the lateral_toc_mouth/mid/head Green's-function controls: a\n"
        f"# uniform, domain-wide background lateral TOC load (no other FORTE control\n"
        f"# uses this site variant -- see that script's module docstring).\n"
        f"from .{site}_forte_campaign import *  # noqa: F401,F403\n"
        f"LATERAL_INFLOW = {inflow!r}\n"
        f"LATERAL_CONC = {{'TOC': {LATERAL_TOC_CONC!r}}}\n"
    )
    site_path.write_text(body)
    return site_mod_name


def controls_for(site):
    """[(control_name, override_kind, override_value), ...] for one river."""
    discharge_file, seatemp_file, boundary, M = _site_boundary_forcing(site)
    out = [("discharge", "DISCHARGE_FILE", discharge_file),
           ("clb_T", "SEATEMP_FILE", seatemp_file)]
    for sp in MARINE_SPECIES:
        if sp in boundary and "clb" in boundary[sp]:
            out.append((f"clb_{sp}", ("BOUNDARY_FORCING", sp, "clb"), boundary[sp]["clb"]))
    # Upstream (cub) controls -- matches exactly what the campaign's grab sampling
    # measures at the upstream station (see module docstring).
    if "TOC" in boundary and "cub" in boundary["TOC"]:
        out.append(("cub_TOC", ("BOUNDARY_FORCING", "TOC", "cub"), boundary["TOC"]["cub"]))
    for sp in UPSTREAM_CONST_SPECIES:
        out.append((f"cub_{sp}", ("BOUNDARIES_CONST", sp), None))
    # Reach-local controls -- local in-stream rate and local lateral loading, at
    # coarse upper/mid/lower-third resolution (see module docstring addendum).
    for reach_name, (i0, i1) in _reach_bounds(M).items():
        out.append((f"rate_aerdeg_{reach_name}", ("KOX_MULT_REACH", i0, i1), None))
        out.append((f"lateral_toc_{reach_name}", ("LATERAL_BOOST_REACH", i0, i1), None))
    return out


def write_perturbed_forcing_and_site(site, control_name, override_kind, src_fname):
    site_mod_name = f"_fortegreen_{site}_{control_name}"
    site_path = SITES_DIR / f"{site_mod_name}.py"

    if isinstance(override_kind, tuple) and override_kind[0] == "KOX_MULT_REACH":
        _, i0, i1 = override_kind
        mult = 1.0 + PERT_FRAC
        body = (
            f"# AUTO-GENERATED by tools/build_green_perturbations_forte.py -- do not\n"
            f"# hand-edit; regenerate by rerunning that script. Local +{PERT_FRAC*100:.0f}%\n"
            f"# aerobic-degradation-rate perturbation confined to grid cells {i0}-{i1}\n"
            f"# (one reach third) for the FORTE Green's-function engine.\n"
            f"from .{site}_forte_campaign import *  # noqa: F401,F403\n"
            f"KOX_MULT_REACH = ({i0}, {i1}, {mult!r})\n"
        )
        site_path.write_text(body)
        return site_mod_name

    if isinstance(override_kind, tuple) and override_kind[0] == "LATERAL_BOOST_REACH":
        _, i0, i1 = override_kind
        lateral_base_mod = f"_fortegreen_{site}_lateral_base"
        body = (
            f"# AUTO-GENERATED by tools/build_green_perturbations_forte.py -- do not\n"
            f"# hand-edit; regenerate by rerunning that script. Local +{PERT_FRAC*100:.0f}%\n"
            f"# lateral-TOC-loading boost confined to grid cells {i0}-{i1} (one reach\n"
            f"# third), on top of the uniform background in {lateral_base_mod}.py\n"
            f"# (written by write_lateral_base_site -- must exist before this runs).\n"
            f"from .{lateral_base_mod} import *  # noqa: F401,F403\n"
            f"LATERAL_BOOST_REACH = ({i0}, {i1}, {PERT_FRAC!r})\n"
        )
        site_path.write_text(body)
        return site_mod_name

    if isinstance(override_kind, tuple) and override_kind[0] == "BOUNDARIES_CONST":
        # Plain constant override -- no forcing file involved at all.
        species = override_kind[1]
        body = (
            f"# AUTO-GENERATED by tools/build_green_perturbations_forte.py -- do not\n"
            f"# hand-edit; regenerate by rerunning that script. +{PERT_FRAC*100:.0f}% "
            f"upstream (cub) {species} perturbation for the FORTE Green's-function\n"
            f"# engine -- a plain BOUNDARIES constant override, not a forcing file\n"
            f"# (this species' cub has no time-varying source in the model yet).\n"
            f"from .{site}_forte_campaign import *  # noqa: F401,F403\n"
            f"BOUNDARIES = dict(BOUNDARIES)\n"
            f"BOUNDARIES[{species!r}] = (BOUNDARIES[{species!r}][0], "
            f"BOUNDARIES[{species!r}][1] * {1.0 + PERT_FRAC})\n"
        )
        site_path.write_text(body)
        return site_mod_name

    SCRATCH.mkdir(parents=True, exist_ok=True)
    pert_fname = f"_fortegreen_scratch/{site}_{control_name}_pert.csv"
    _perturb_csv(src_fname, FORC / pert_fname)

    if override_kind in ("DISCHARGE_FILE", "SEATEMP_FILE"):
        body = (
            f"# AUTO-GENERATED by tools/build_green_perturbations_forte.py -- do not\n"
            f"# hand-edit; regenerate by rerunning that script. +{PERT_FRAC*100:.0f}% "
            f"{override_kind} perturbation for the FORTE Green's-function engine.\n"
            f"from .{site}_forte_campaign import *  # noqa: F401,F403\n"
            f"{override_kind} = {pert_fname!r}\n"
        )
    else:
        _, species, end = override_kind
        body = (
            f"# AUTO-GENERATED by tools/build_green_perturbations_forte.py -- do not\n"
            f"# hand-edit; regenerate by rerunning that script. +{PERT_FRAC*100:.0f}% "
            f"{species} {end} perturbation for the FORTE Green's-function engine.\n"
            f"from .{site}_forte_campaign import *  # noqa: F401,F403\n"
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


def _already_done(outdir):
    """True if this (site, control)'s run already completed successfully --
    checked via its run.log's last line, so re-running this script after adding
    new controls does not redo the ones already finished."""
    log = outdir / "run.log"
    if not log.exists():
        return False
    try:
        last = log.read_text().strip().splitlines()[-1]
    except IndexError:
        return False
    return f"t: {MAXT_DAYS:.2f} days" in last


def main():
    sites = sys.argv[1:] or SITES
    jobs = []   # (label, site_variant, outdir) -- NO baseline jobs: reused directly
                # from runs/forte_campaign/<site>/output.nc (already built)
    skipped = 0
    for site in sites:
        # Lateral baseline -- a prerequisite for the lateral_toc_* controls below,
        # not a "control" itself (not added to CONTROLS in the downstream scripts).
        lateral_base_variant = write_lateral_base_site(site)
        lateral_base_outdir = RUNS / "lateral_base" / site
        if _already_done(lateral_base_outdir):
            skipped += 1
        else:
            jobs.append((f"{site} lateral_base", lateral_base_variant, lateral_base_outdir))

        for control_name, override_kind, src_fname in controls_for(site):
            outdir = RUNS / f"pert_{control_name}" / site
            if _already_done(outdir):
                skipped += 1
                continue
            site_variant = write_perturbed_forcing_and_site(site, control_name, override_kind, src_fname)
            jobs.append((f"{site} {control_name}", site_variant, outdir))

    if skipped:
        print(f"{skipped} (site, control) runs already completed -- skipped")
    print(f"{len(jobs)} perturbed runs to launch ({MAXT_DAYS} days @ 6-min cadence "
          f"each); baseline reused from runs/forte_campaign/")
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
    print("all FORTE Green's-function perturbation runs finished")


if __name__ == "__main__":
    main()
