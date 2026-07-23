#!/usr/bin/env python3
"""
Verification harness for the IDEALIZED North Slope river (sites/idealized.py).

A reusable regression/verification fixture: it runs the model on a known analytic
input and asserts a set of physical invariants on the output, so a future change
that quietly breaks transport, the carbonate solve, the ice model, or the new
time-varying-boundary machinery is caught.

Modes
-----
  --check-only   Structural checks only, NO model run (instant, no numba needed):
                 forcings exist and are 365-long, config resolves, every
                 BOUNDARY_FORCING species/end is valid. Good first line of defence.

  (default)      --check-only, then a SHORT run (5 days) asserting the run completes,
                 output.nc is well-formed, no NaN/negative concentrations, pH in
                 range, winter ice grows, and the boundary cells match the prescribed
                 forcing. ~1-2 min including numba's first-call compile.

  --full         A long run (default 220 days: through freeze-up, the freshet
                 break-up, the DOC/alkalinity pulse and the summer salinity minimum)
                 asserting the seasonal-cycle physics as well. ~10-15 min.

Override the duration with --days N / --warmup D. Exits non-zero if any check fails.

    python tools/verify_idealized.py --check-only
    python tools/verify_idealized.py
    python tools/verify_idealized.py --full
"""
import argparse
import os
import subprocess
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.normpath(os.path.join(HERE, ".."))
CODE = os.path.join(ROOT, "code")
FORCINGS = os.path.join(ROOT, "forcing")
# The --full run IS the canonical fixture run (runs/idealized/), which the diagnostics
# PDF and movies read. The short --quick run goes to a scratch subdir so it never
# clobbers a good full run.
RUNDIR_FULL = os.path.join(ROOT, "runs", "idealized")
RUNDIR_QUICK = os.path.join(ROOT, "runs", "idealized", ".quick_check")
N_DAYS_FORCING = 365
DAY = 86400.0

# Species whose concentration must never go negative (all transported solutes,
# including the Arctic-extension tracers RDOC/CH4/N2O).
NONNEG = ["S", "DIA", "dSi", "NO3", "NH4", "PO4", "O2", "TOC", "RDOC", "CH4", "N2O",
          "SPM", "DIC", "ALK"]


class Checks:
    """Collect (name, ok, detail) results and report."""

    def __init__(self):
        self.rows = []

    def add(self, name, ok, detail=""):
        self.rows.append((name, bool(ok), detail))
        mark = "PASS" if ok else "FAIL"
        print(f"  [{mark}] {name}" + (f"  -- {detail}" if detail else ""))
        return ok

    @property
    def ok(self):
        return all(ok for _, ok, _ in self.rows)


# ---------------------------------------------------------------------------
# structural checks (no model run)
# ---------------------------------------------------------------------------
def load_site():
    """Import sites.idealized without pulling in config/numba."""
    sys.path.insert(0, CODE)
    import importlib
    return importlib.import_module("sites.idealized")


def check_structure(c):
    print("structural checks (no model run):")
    site = load_site()

    # every forcing the site names must exist and be 365 long
    files = [site.DISCHARGE_FILE, site.WATERTEMP_FILE, site.SEATEMP_FILE,
             site.WIND_FILE, site.SOLAR_FILE, site.AIRTEMP_FILE,
             site.RELHUM_FILE, site.PCO2_FILE]
    for ends in site.BOUNDARY_FORCING.values():
        files.extend(ends.values())
    for fn in files:
        p = os.path.join(FORCINGS, fn)
        if not os.path.exists(p):
            c.add(f"forcing exists: {fn}", False, "MISSING -- run build_idealized_forcings.py")
            continue
        arr = np.genfromtxt(open(p, encoding="utf-8-sig"), delimiter=",", dtype=float)
        c.add(f"forcing 365-long & finite: {fn}",
              arr.shape == (N_DAYS_FORCING,) and np.all(np.isfinite(arr)),
              f"shape {arr.shape}")

    # BOUNDARY_FORCING references valid species and ends
    valid_sp = set(site.BOUNDARIES)
    for sp, ends in site.BOUNDARY_FORCING.items():
        c.add(f"BOUNDARY_FORCING species valid: {sp}", sp in valid_sp)
        c.add(f"BOUNDARY_FORCING ends valid: {sp}",
              all(e in ("clb", "cub") for e in ends))

    c.add("IS_IDEALIZED flag set", getattr(site, "IS_IDEALIZED", False))


# ---------------------------------------------------------------------------
# run the model
# ---------------------------------------------------------------------------
def run_model(days, warmup, ts, outdir):
    os.makedirs(outdir, exist_ok=True)
    env = dict(os.environ,
               CGEM_SITE="idealized",
               CGEM_MAXT_DAYS=str(days),
               CGEM_WARMUP_DAYS=str(warmup),
               CGEM_TS=str(ts),
               CGEM_OUTPUT="nc",
               PYTHONPATH=CODE,
               PYTHONWARNINGS="ignore")
    rel = os.path.relpath(outdir, ROOT)
    print(f"\nrunning idealized model: {days} d (warmup {warmup} d, TS={ts}) in {rel}/ ...")
    log = os.path.join(outdir, "run.log")
    with open(log, "w") as f:
        r = subprocess.run([sys.executable, os.path.join(CODE, "main.py")],
                           cwd=outdir, env=env, stdout=f, stderr=subprocess.STDOUT)
    return r.returncode, os.path.join(outdir, "output.nc"), log


def _daily_mean(x, spd):
    """Bin a per-timestep series into daily means (drops a ragged final day)."""
    n = (len(x) // spd) * spd
    return x[:n].reshape(-1, spd).mean(axis=1) if n >= spd else x


def _interp_forcing(fn, times_s):
    arr = np.genfromtxt(open(os.path.join(FORCINGS, fn), encoding="utf-8-sig"),
                        delimiter=",", dtype=float)
    axis = np.linspace(0, N_DAYS_FORCING * DAY, N_DAYS_FORCING)
    tt = np.where(times_s > 365 * DAY, times_s - 365 * DAY, times_s)  # repeatYear
    return np.interp(tt, axis, arr)


def check_output(c, ncpath, days, warmup, full):
    from netCDF4 import Dataset
    ds = Dataset(ncpath)
    t = ds["time"][:].astype(float)
    M = ds.dimensions["x"].size - 1

    def field(name):
        return np.ma.filled(ds[name][:], np.nan)

    print("\noutput checks:")
    # 1. structure
    expect = {"S", "TOC", "ALK", "DIC", "NO3", "pH", "T", "FCO2",
              "ice_thickness", "velocity", "depth", "width"}
    c.add("expected variables present", expect <= set(ds.variables),
          f"missing {expect - set(ds.variables) or 'none'}")
    c.add("time steps written", t.size > 0, f"{t.size} steps")

    # 2. no NaN / negatives in interior (cells 1..M)
    for name in NONNEG:
        a = field(name)[:, 1:]
        c.add(f"{name}: finite", np.all(np.isfinite(a)),
              f"{100*np.isnan(a).mean():.2f}% nan")
        c.add(f"{name}: non-negative", np.nanmin(a) >= -1e-6,
              f"min {np.nanmin(a):.3g}")

    # 3. pH physical
    ph = field("pH")[:, 1:]
    c.add("pH within guard range [2,12]",
          np.nanmin(ph) >= 2 and np.nanmax(ph) <= 12,
          f"[{np.nanmin(ph):.2f}, {np.nanmax(ph):.2f}]")

    # 4. time-varying boundary propagation. The upstream outflow cell (M) is set
    #    almost directly from cub by openbound, so it should closely track the cub
    #    series. The downstream mouth cell (1) is tidally advected around clb, so it
    #    is only a loose proxy -- tested by daily-mean CORRELATION, not tight match.
    #    Neither can be verified on a short winter window where the inputs barely
    #    move, so the tracking tests run only in --full; quick mode just asserts the
    #    boundary cells stay inside the species' seasonal envelope (+ a tidal band).
    site = load_site()
    # Species that also receive distributed lateral loading: their upstream cell is
    # boundary + lateral + reactions, NOT a clean mirror of the cub forcing, so it is
    # verified by correlation (it still responds to the forcing) rather than tight match.
    lateral_species = set(getattr(site, "LATERAL_CONC", {}))
    steps_per_day = max(1, int(round(DAY / (np.median(np.diff(t)) if t.size > 1 else DAY))))
    for sp, ends in site.BOUNDARY_FORCING.items():
        for end, fn in ends.items():
            season = np.genfromtxt(open(os.path.join(FORCINGS, fn), encoding="utf-8-sig"),
                                   delimiter=",", dtype=float)
            lo, hi, rng = float(season.min()), float(season.max()), float(np.ptp(season))
            cell = M if end == "cub" else 1
            model = field(sp)[:, cell]
            target = _interp_forcing(fn, t)
            if end == "cub":
                # The upstream outflow cell is set almost directly from cub, so its
                # value is bounded by the imposed series (+ a small numeric band).
                band = max(2.0, 0.06 * rng)
                c.add(f"boundary {sp}.cub within seasonal envelope",
                      np.nanmin(model) >= lo - band and np.nanmax(model) <= hi + band,
                      f"model [{np.nanmin(model):.1f},{np.nanmax(model):.1f}] vs input [{lo:.1f},{hi:.1f}]")
                if full and rng > 1e-6 and np.nanstd(target) > 1e-6:
                    if sp in lateral_species:
                        # lateral loading + reactions offset the cell -> check it RESPONDS
                        md = _daily_mean(model, steps_per_day)
                        td = _daily_mean(target, steps_per_day)
                        cc = float(np.corrcoef(md, td)[0, 1]) if md.size > 2 else 0.0
                        c.add(f"boundary {sp}.cub correlates with forcing (lateral-loaded)",
                              cc > 0.5, f"r={cc:.2f}")
                    else:
                        dev = float(np.nanmax(np.abs(model - target)))
                        c.add(f"boundary {sp}.cub tracks forcing", dev <= max(3.0, 0.08 * rng),
                              f"max dev {dev:.1f} over input range {rng:.0f}")
            else:
                # The mouth (clb) cell is NOT bounded by the marine boundary: a strong
                # freshet flushes the salt intrusion out to sea and the cell goes fresh
                # (physically correct). So no envelope check -- verify instead that the
                # cell RESPONDS to the marine forcing, via daily-mean correlation.
                if full and rng > 1e-6 and np.nanstd(target) > 1e-6:
                    md = _daily_mean(model, steps_per_day)
                    td = _daily_mean(target, steps_per_day)
                    cc = float(np.corrcoef(md, td)[0, 1]) if md.size > 2 else 0.0
                    c.add(f"boundary {sp}.clb correlates with forcing (daily mean)",
                          cc > 0.5, f"r={cc:.2f}")

    # 5. ice: with a run that reaches deep winter, ice must grow; and if it spans the
    #    freshet, hydraulic break-up must clear it.
    ice = field("ice_thickness")[:, 1:]
    day = t / DAY
    wint = (day >= 1) & (day <= min(days, 40))
    if wint.any():
        c.add("ice grows in winter", np.nanmax(ice[wint]) > 0.01,
              f"max winter ice {np.nanmax(ice[wint]):.2f} m")
    if full and day.max() >= 190:
        post = ice[(day >= 185) & (day <= 205)]
        if post.size:
            c.add("ice clears after freshet break-up",
                  np.nanmean(post) < 0.02,
                  f"mean ice day 185-205 = {np.nanmean(post):.3f} m")

    # 6. full-cycle signals actually appear
    if full:
        toc_M = field("TOC")[:, M]
        if day.max() >= 165:
            c.add("freshet DOC pulse reaches upstream cell",
                  np.nanmax(toc_M[day >= 120]) > 400,
                  f"peak TOC[M] = {np.nanmax(toc_M):.0f} (baseflow ~250)")
        s1 = field("S")[:, 1]
        if day.max() >= 210:
            c.add("summer sea-ice freshening at marine cell",
                  np.nanmin(s1[day >= 180]) < 27,
                  f"min S[1] = {np.nanmin(s1):.1f} (winter ~30)")
        fco2 = field("FCO2")[:, 1:]
        c.add("FCO2 finite", np.all(np.isfinite(fco2[~np.isnan(fco2)])) if fco2.size else True)

    # 7. Arctic biogeochemistry extension (docs/arctic_biogeochemistry.md)
    if full:
        summer = day >= 170          # open-water window
        # (1) CDOM photomineralisation switches on in sunlit open water
        if "photo" in ds.variables:
            ph = field("photo")[:, 1:]
            c.add("photomineralisation active in summer",
                  np.nanmax(ph[summer]) > 0 if summer.any() else False,
                  f"max photo = {np.nanmax(ph):.2e} mmol m^-3 s^-1")
        # (2a) CH4: supersaturated river outgasses (ch4_ex > 0) once ice clears
        if "ch4_ex" in ds.variables:
            ce = field("ch4_ex")[:, 1:]
            c.add("CH4 outgasses in open water (ch4_ex > 0)",
                  np.nanmax(ce[summer]) > 0 if summer.any() else False,
                  f"max ch4_ex = {np.nanmax(ce):.2e}")
        # (2b) N2O finite and produced
        if "N2O" in ds.variables:
            n = field("N2O")[:, 1:]
            c.add("N2O finite & non-negative", np.all(np.isfinite(n)) and np.nanmin(n) >= -1e-6,
                  f"[{np.nanmin(n):.3g}, {np.nanmax(n):.3g}]")
        # (3) benthic SOD -> DIC efflux is always on (runs under ice too)
        if "sod" in ds.variables:
            sd = field("sod")[:, 1:]
            c.add("benthic SOD efflux positive", np.nanmax(sd) > 0,
                  f"max sod = {np.nanmax(sd):.2e} mmol m^-3 s^-1")
        # (4) distributed lateral loading raises interior CH4 above the river boundary
        if "CH4" in ds.variables:
            ch4 = field("CH4")
            cub_ch4 = 0.6         # idealized river-boundary CH4 [mmol m^-3]
            c.add("lateral loading elevates interior CH4",
                  np.nanmax(ch4[:, 1:M]) > cub_ch4 * 1.05,
                  f"max interior CH4 = {np.nanmax(ch4[:, 1:M]):.3g} (river bnd {cub_ch4})")
        # RDOC time-varying boundary (new-species BOUNDARY_FORCING test)
        if "RDOC" in ds.variables:
            rM = field("RDOC")[:, M]
            if day.max() >= 165:
                c.add("RDOC freshet pulse reaches upstream cell",
                      np.nanmax(rM[day >= 120]) > 600,
                      f"peak RDOC[M] = {np.nanmax(rM):.0f} (baseflow ~250)")

    ds.close()


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--check-only", action="store_true",
                    help="structural checks only, no model run")
    ap.add_argument("--full", action="store_true",
                    help="long seasonal run (default 220 d) with cycle-physics checks")
    ap.add_argument("--days", type=int, default=None, help="override run length [days]")
    ap.add_argument("--warmup", type=int, default=None, help="override warmup [days]")
    ap.add_argument("--ts", type=int, default=48, help="output cadence (timesteps)")
    args = ap.parse_args()

    c = Checks()
    check_structure(c)

    if not args.check_only:
        days = args.days if args.days is not None else (220 if args.full else 5)
        warmup = args.warmup if args.warmup is not None else (30 if args.full else 1)
        # --full writes the canonical fixture run (read by the PDF/movies); --quick
        # uses a scratch subdir so it never clobbers a good full run.
        outdir = RUNDIR_FULL if args.full else RUNDIR_QUICK
        rc, nc, log = run_model(days, warmup, args.ts, outdir)
        if not c.add("model run exit 0", rc == 0, f"see {log}"):
            print(f"\nmodel failed; tail of {log}:")
            print("".join(open(log).readlines()[-15:]))
        elif not c.add("output.nc written", os.path.exists(nc)):
            pass
        else:
            check_output(c, nc, days, warmup, args.full)

    print("\n" + ("ALL CHECKS PASSED" if c.ok else "SOME CHECKS FAILED"))
    sys.exit(0 if c.ok else 1)


if __name__ == "__main__":
    main()
