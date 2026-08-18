# hybrid_darwin_bling

A hybrid Darwin(forward)/BLING(adjoint) **surrogate-gradient** optimization
of ocean carbon-chemistry initial conditions (DIC, alkalinity) against
SOCAT surface pCO2. Each iteration runs Darwin's own forward model to
compute the *true* objective (Darwin's own pCO2 vs. SOCAT), and separately
runs BLING's forward+adjoint to compute an *approximate* gradient -- then
feeds (Darwin's cost, BLING's gradient) to the same offline M1QN3 driver
used in `../global_oce_biogeo_bling_SOCAT/`, as if they were a consistent
pair.

This builds directly on `../darwin_bling_comparison/`'s finding that
BLING's adjoint sensitivity closely matches Darwin's own forward-FD
sensitivity specifically for DIC and alkalinity (the two tracers whose
effect on pCO2 goes through the carbonate-chemistry equations both models
share, bit-identical) -- and does *not* match well for the
biology-mediated tracers (NO3, O2, PO4, FeT). That's why this hybrid setup
is deliberately restricted to DIC+ALK only.

## Why a surrogate gradient at all

Darwin's full ecosystem (10 explicit, trait-based competing plankton
types across 7 functional groups, non-smooth Liebig-limitation `min()`
terms, threshold-based grazing/mortality switches) would be difficult or
unreliable to differentiate with TAF -- the automatic-differentiation
tool used to build BLING's adjoint. Rather than attempt that, this
experiment uses BLING's already-working, much-simpler adjoint as a stand-
in gradient for optimizing Darwin's own initial conditions, exploiting the
close agreement already validated for DIC/ALK specifically.

## Get the code

```sh
git clone https://github.com/MITgcm-contrib/ecco_darwin.git
git clone https://github.com/MITgcm/MITgcm.git
cd MITgcm
git clone https://github.com/darwinproject/darwin3.git
```
You need both `../global_oce_biogeo_darwin` and
`../global_oce_biogeo_bling_SOCAT` built and runnable first (their own
`mitgcmuv` / `mitgcmuv_ad` binaries, respectively), including
`global_oce_biogeo_darwin/scripts/match_ics_to_bling.py` already applied
so Darwin's DIC/ALK/O2/NO3/PO4/Fe ICs match BLING's own.

## Mechanism

Each cycle:
1. **BLING forward+adjoint** (`mitgcmuv_ad`) at the current control guess
   -- writes `ecco_cost_MIT_CE_000.opt<cycle>` (BLING's own cost *and*
   the packed DIC/ALK adjoint gradient, in one binary file).
2. **Apply that same control guess to Darwin's own physical IC**
   (`scripts/apply_control_to_darwin_ic*.py`) -- reproduces MITgcm's
   `ctrl_map_genarr3d.F` preconditioning transform externally:
   `physical = prior + xx_gen / sqrt(ctrl_weight)`, using the same prior
   IC and weight files (`wt_DIC.bin`/`wt_ALK.bin`) BLING's own adjoint run
   used internally.
3. **Darwin forward run** (`mitgcmuv`) on that IC.
4. **Compute Darwin's own SOCAT cost** (`scripts/compute_darwin_cost*.py`)
   -- same cost formula as BLING's own `pkg/profiles`
   (`J = sum(weight * (model-obs)^2)`), verified against
   `pkg/profiles/cost_profiles.F` / `pkg/cost/cost_final.F` in
   `../darwin_bling_comparison/README.md`.
5. **Binary-patch `ecco_cost_MIT_CE_000.opt<cycle>`'s cost field with
   Darwin's own cost** (`scripts/patch_ecco_cost*.py`) -- BLING's packed
   adjoint-gradient payload is left untouched. The file's fixed byte
   layout (Fortran unformatted sequential, big-endian: `fc` is an 8-byte
   `real*8` at a fixed offset, gradient payload immediately after) was
   verified byte-for-byte against a real 10-cycle SOCAT campaign before
   trusting this patch -- see the scripts' docstrings for the exact
   layout.
6. **`optim.x`** reads (Darwin's cost, BLING's gradient) as if
   consistent and proposes the next control guess via M1QN3.

This is a deliberate mismatch: the gradient was computed from BLING's
cost function, not Darwin's. It works to the extent that BLING's adjoint
is a good local surrogate for Darwin's own sensitivity -- which is
exactly what was validated for DIC/ALK in `../darwin_bling_comparison/`.

## Two configurations

| | `run_bling` / `run_darwin` | `run_bling_1yr` / `run_darwin_1yr` |
|---|---|---|
| Duration | 1 month (`nTimeSteps=2880`) | 1 year (`nTimeSteps=34560`) |
| SOCAT cost scope | month01 only (`mult_profiles` zeroed for months 2-12) | all 12 months (unmodified `data.profiles`) |
| `data.optim` | `numiter=5` | `numiter=1` |
| Driver | `scripts/run_hybrid_loop.sh` | `scripts/run_hybrid_loop_1yr.sh` |
| Cost/IC scripts | `*_to_darwin_ic.py`, `compute_darwin_cost.py`, `patch_ecco_cost.py` | `*_to_darwin_ic_1yr.py`, `compute_darwin_cost_1yr.py`, `patch_ecco_cost_1yr.py` |

Both restrict the active control vector to DIC+ALK only via `data.ctrl`
(`run_bling*/data.ctrl`): theta, salt, O2, NO3, PO4, and Fe are
deliberately left **unregistered** (not just zero-weighted) --
`ctrl_map_ini_genarr.F` only maps fields whose `xx_genarr3d_file` slot is
non-blank, so omitting them freezes those fields at their prior IC for
the whole run rather than letting the (poor-surrogate) adjoint gradient
move them.

The full-year run's `data.diagnostics` (`run_darwin_1yr/data.diagnostics`)
registers 12 separate one-off `darwinPCO2Snap01`..`12` diagnostic
snapshots, one per SOCAT month, each timed to that month's own shared
observation timestamp (read directly from each
`socat_pco2_clim_month*.nc`'s `prof_YYYYMMDD`/`prof_HHMMSS` -- all 12
land exactly on an integer model timestep; see the file's header comment
for the full derivation).

## Running

```sh
cd run_bling      # or run_bling_1yr
./scripts/run_hybrid_loop.sh      # or run_hybrid_loop_1yr.sh, from hybrid_darwin_bling/
```
The driver script `cd`s into `run_bling*/` itself and locates
`run_darwin*/` and `scripts/` relative to its own location, so run it
from `hybrid_darwin_bling/` (not from inside `run_bling/`).

**Disk**: the full-year config's BLING adjoint tapes alone are
~9.75GB/cycle (pure TAF checkpoint-recompute scratch, unrelated to
`pChkptFreq`/`chkptFreq`, which are set to 0 here since each cycle
cold-starts independently and doesn't need ordinary restart pickups).
`run_hybrid_loop_1yr.sh` deletes these after each cycle's cost has been
extracted and injected. Budget accordingly if you extend `numiter`.

**Timing** (on a 14-core, single-threaded-per-run machine): the 1-month
config takes ~4-4.5h/cycle (dominated by Darwin's forward run); the
full-year config takes ~2.4 days/cycle (~7.3h BLING adjoint + ~49h Darwin
forward, roughly linear in `nTimeSteps` from the 1-month benchmark).
**Each M1QN3 "iteration" in `data.optim` typically costs *two* cycles**,
not one -- the first cycle evaluates the step, the second is needed for
M1QN3's line search / cap-recognition to complete (confirmed against
`m1qn3_output.txt`'s own `iter`/`simul` counters in both runs here).

## Results

**1-month pilot** (`numiter=5`, 6 cycles 0-5): Darwin's own cost dropped
**84.5%** (40,909 -> 6,321), monotonically every cycle, no rejected
line-search steps. See `figures/cost_trajectory.png`.

**Full-year run** (`numiter=1`, 2 cycles 0-1): Darwin's own cost dropped
**10.5%** (340,153 -> 304,597) in the single completed iteration --
smaller than the pilot's relative drop, as expected: one iteration
against the full 12-month record is a much broader-scope problem than
five iterations against a single restricted month. BLING's own cost
(212,373 -> 164,280) exactly matches the reference value from the
original real 10-cycle SOCAT campaign at cycle 0, confirming the
full-year configuration here is set up correctly. See
`figures/cost_trajectory_1yr.png`. Full numeric trajectory for both runs:
`results/cost_trajectory_results.csv`.

**Spatial pattern of the IC adjustment** (`figures/ic_comparison.png`,
`figures/ic_delta_with_obs_overlay.png`): DIC increased and ALK decreased
almost everywhere (complementary effects on pCO2 in carbonate chemistry,
consistent with Darwin's baseline pCO2 sitting below the SOCAT
observations pre-optimization) -- but the adjustment is **~17x larger**
at grid cells with an actual SOCAT month-01 observation than elsewhere.
Only ~34% of ocean grid cells (2,759 of 8,192) have any SOCAT observation
in a given month -- real ship-track coverage, not a full grid -- so the
cost function only has direct gradient information there; grid cells far
from any track pick up only weak indirect adjustment via advective spread
over the run. This is expected behavior for a sparse-observation adjoint
optimization, not an artifact, and is a concrete argument for why a
full-year (12-month, denser cumulative track coverage) campaign would
likely fill in this patchiness more than any single month alone.

## Known limitations / open items

- Only 1 (full-year) or 5 (1-month) M1QN3 iterations have been run --
  neither trajectory has been carried to convergence. Extending `numiter`
  costs roughly another ~4.5h (1-month) or ~2.4 days (full-year) *per
  additional cycle*.
- The surrogate-gradient mismatch is only validated as small for DIC/ALK
  at the 38 points studied in `../darwin_bling_comparison/`; it has not
  been re-validated at the specific (evolving) states visited during this
  optimization itself.
- `data.ctrl`'s restriction to DIC+ALK means theta, salt, and the other 4
  BLING tracers are frozen throughout -- this is a deliberate scope choice
  (see "Why a surrogate gradient at all"), not a limitation to fix.
