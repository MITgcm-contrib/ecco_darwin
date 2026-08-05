# Verify your downscaled regional set-up

---
## Preliminary information
**Requirement**: 
> Before proceeding to the following instructions, you will need to complete steps in the README, STEP 1, STEP 2, STEP 3 and STEP 4.\
> The anaconda python environment will be needed to run the code.

STEP4 gets the model compiled and submitted. This step answers the two questions you
actually want answered before committing supercomputer time to a production integration:

1. **Did it run sanely?** No NaN, the advective CFL condition respected, temperature and
   salinity inside physical bounds, and no runaway kinetic energy.
2. **Did the model actually use what you built?** It is entirely possible for a run to
   complete cleanly while silently ignoring your initial conditions, or while carrying no
   flow at all through one of the open boundaries. Both failure modes are quiet: the model
   does not complain, and the output looks superficially reasonable.

Both checks are performed by ``verify_run.py`` in the ``utils`` folder, and both are
self-contained — neither requires re-running the parent model.

<u>Note:</u> this is deliberately a **short** experiment. The point is to fail fast and
cheaply. A one-day integration is usually enough to expose an unstable configuration or a
dead boundary.

---
## I. Configure and run a short verification integration

Start from the run directory you prepared in STEP4, Section IV.b, and make three temporary
changes to the ``data`` namelist.

> - Integrate for a short period only. With ``deltaT=30.``, 2880 steps is one day:
```
 nTimeSteps=2880,
```
> - Ask for frequent monitor output, so there are several samples to judge drift from. This is what ``verify_run.py`` reads:
```
 monitorFreq=3600.,
 monitorSelect=2,
```
> - Dump the state the model actually starts from, so the initial condition can be checked against your pickup files:
```
 dumpInitAndLast=.TRUE.,
```

<u>Note:</u> ``monitorSelect=2`` includes the dynamic-state statistics (``dynstat_*``) used by
the health check. ``dumpInitAndLast=.TRUE.`` makes MITgcm write ``T.0000000000.data``,
``S.0000000000.data``, ``U.``, ``V.`` and ``Eta.0000000000.data`` into the run directory —
these are the files Section III compares against ``forcings/pickups/``.

Then submit as in STEP4:
```
qsub batch_file
```

When it finishes, copy the run directory's ``STDOUT.*`` file and the
``*.0000000000.data`` dumps back to the machine running your anaconda environment, or run
``verify_run.py`` directly on the supercomputer if the environment is installed there.

---
## II. Run the verification

```
conda activate downscaling
cd ecco_darwin/regions/downscaling/utils/
python3 verify_run.py -d /path/to/config_dir -n name_of_the_region\
                      -bnd EWNS -i pickup_iteration_number\
                      -r /path/to/run -plot -v
```

**Example:**

```
python3 verify_run.py -d ecco_darwin/regions/downscaling/ -n NorthSlope \
                      -bnd ENW -i 683784 -r /path/to/run -plot -v
```

To get more information about the options required for this code run ``verify_run.py -h``. Here are additional details about the options:
> - -d: Your ``config_dir`` — the same ``-d`` you gave the STEP3 scripts.
> - -n: The name of your region.
> - -bnd: Open boundaries to check, as a string of E/W/N/S. Required unless you pass ``--health``.
> - -i: The pickup iteration the run started from — the same ``-i`` you gave ``gen_pickups.py``. Required unless you pass ``--health``.
> - -r: The model run directory (holds ``STDOUT*`` and the initial-state dumps).
> - -stdout: Explicit path to a STDOUT file, if it is not at ``<run_dir>/STDOUT*`` (**optional**).
> - --health: Run **only** the health check of Section III.1. Needs nothing but a STDOUT file.
> - --inputs: Run **only** the input-consistency check of Section III.2. Needs no model run at all, so it can be used straight after STEP3.
> - -plot: Save a time-series plot of the monitored quantities under ``<output_dir>/diagnostics/`` (**optional**, requires matplotlib).
> - -o: Where to save plots (**optional**, default ``<config_dir>``).
> - -v: verbose.

With neither ``--health`` nor ``--inputs``, both run. The script exits with status **0** when
everything passes and **1** when anything FAILs, so it can gate a production job:

```
python3 verify_run.py -d ... -n ... -bnd ENW -i 683784 -r ./run && qsub production_job
```

---
## III. Reading the output

### III.1 Run health

Parsed from the ``%MON`` lines in STDOUT:

| Check | FAIL means |
|---|---|
| NaN check | A monitored quantity became NaN — the run blew up. Reduce ``deltaT``, or check the boundary forcing with `check_obcs_transport.py`. |
| CFL | ``advcfl_*`` reached 1.0 or above. Reduce ``deltaT``. A WARN at 0.5 means you are close to the limit. |
| theta / salt range | Temperature or salinity left a generous physical range. Usually a bad initial or boundary condition, not a timestep problem. |
| KE drift | Mean kinetic energy grew more than 10x. Most often inconsistent open-boundary transport — see `check_obcs_transport.py` in STEP3, Section III.2. |
| termination | No ``Execution ended Normally`` in STDOUT. The run stopped early; check the tail of the file and the job stderr. |

A healthy short run looks like this:

```
=== RUN HEALTH (STDOUT.0000) ===
  [PASS] monitor dumps  6 dump(s) parsed
  [PASS] NaN check      no NaN in any monitored quantity
  [PASS] CFL            max advcfl_uvel_max = 0.210
  [PASS] theta range    [4.200, 29.500] degC
  [PASS] salt range     [33.100, 36.900] psu
  [PASS] KE drift       mean KE 1.000e-03 -> 1.250e-03 (1.25x)
  [PASS] termination    Execution ended Normally
```

### III.2 Input consistency

Three groups of checks, none of which need the model to have run well — or at all, in the
case of the first two:

**OBCS files.** For every boundary and variable, that the file is a whole number of
``(levels x boundary points)`` records for your grid, is free of NaN/Inf, holds physically
plausible values, and — the useful one — is actually **non-zero on the wet cells of the same
interior C-grid face the model will mask it with**:

```
  [PASS] east/UVEL    5 record(s), 120 wet cell(s), 0 (0.0%) zero for the whole run
  [FAIL] west/UVEL    5 record(s), 80 wet cell(s), 80 (100.0%) zero for the whole run
                      -- this boundary carries nothing; classic wrong-C-grid-face symptom
```

A boundary reported as 100% zero carries no flow. That is the failure mode described in the
C-grid staggering note at the end of STEP3, Section III, and it is invisible in a completed
run until you notice the region is not being ventilated. Anything above 50% is reported as a
warning and is worth inspecting with `plot_boundary_velocity_profile.py`.

**Pickup files.** That each is finite and its wet-cell values are physically plausible.

**Initial state ingested.** If ``dumpInitAndLast=.TRUE.`` was set, the state the model
started from is compared against the pickup files ``gen_pickups.py`` wrote:

```
  [PASS] initial T    matches pickup_THETA (rel. diff 3.2e-08)
```

A FAIL here means the model did **not** start from the file named in ``data`` — check the
``hydrogThetaFile``/``uVelInitFile``/``pSurfInitFile`` entries (STEP4, Section II.a) and that
``nIter0=0`` with no ``pickupSuff`` line. The comparison is made at single precision, because
the model writes 32-bit dumps of a 64-bit input; a relative difference around ``1e-8`` is
expected and passes.

<u>Note:</u> if you see ``no <var>.0000000000.data dumps in the run directory``, you did not
set ``dumpInitAndLast=.TRUE.`` in Section I. The rest of the checks still run.

---
## IV. Once it passes

Revert the three temporary ``data`` changes from Section I — restore your production
``nTimeSteps``/``endTime``, and turn ``dumpInitAndLast`` back off (leaving ``monitorFreq`` at
a sensible value for a long run is good practice) — and submit the production integration.

It is worth re-running ``verify_run.py --health`` against the production STDOUT once that run
completes too; the drift and CFL checks are just as meaningful over a long integration, and
cost nothing.

---

**CONGRATULATIONS!!** Your downscaled regional set-up is built, run, and verified.
