# global_oce_biogeo_bling_SOCAT

Adjoint-based (TAF) optimization of ocean biogeochemical initial conditions
against a SOCAT surface-pCO2 climatology, built on MITgcm's
`global_oce_biogeo_bling` verification experiment (BLING biogeochemistry
coupled to the global ocean circulation model).

An M1QN3 offline quasi-Newton optimizer adjusts 8 three-dimensional initial
condition fields -- temperature, salinity, and the 6 core BLING tracers
(DIC, alkalinity, O2, NO3, PO4, Fe) -- to minimize the mismatch between
modeled and observed monthly-climatological surface pCO2, using the
MITgcm-generated adjoint (via TAF) to compute the exact gradient of the
cost function with respect to all ~445,000 control variables in a single
adjoint model integration.

See `../darwin_bling_comparison/README.md` for a companion pilot study that
compares this model's adjoint sensitivities against forward finite
differences from the Darwin ecosystem model, run on this same grid/physics.

## Get the code

```sh
git clone https://github.com/MITgcm-contrib/ecco_darwin.git
git clone https://github.com/MITgcm/MITgcm.git
cd MITgcm
git clone https://github.com/darwinproject/darwin3.git   # BLING lives in darwin3's pkg/bling
```
(or use an existing MITgcm/darwin3 checkout -- any recent checkout with pkg/bling and pkg/profiles works.)

## Why this isn't a typical MITgcm verification experiment

Standard MITgcm verification experiments are short (seconds-to-minutes),
deterministic, and checked automatically via `testreport` against a
reference `results/output.txt`. This one is not: a single forward+adjoint
model evaluation over the full 360-day model year takes **~7.3 hours on
one CPU core** (see "Performance" below), and a full optimization campaign
is dozens of such evaluations run sequentially over several days. It is
included here, at the user's request, as a faithful, documented
reproduction of the actual production campaign rather than a cut-down fast
check -- `results/output.txt` is a reference trajectory to sanity-check a
rerun against, not a pass/fail regression test `testreport` can run in CI.

## Directory layout

```
code_ad/        MITgcm source-code customizations (headers, package list, checkpointing config)
input_ad/       Runtime namelists (binary/NetCDF input files are NOT included here -- see below)
scripts/        Driver script, SOCAT preprocessing, and figure-generation scripts
results/        Reference cost-function trajectory and gradient-check numbers
figures/        Reference diagnostic figures from the campaign this was extracted from
```

**This repo intentionally excludes large binary/NetCDF files** (forcing,
initial conditions, SOCAT climatology, control weights) to keep the repo
small -- only the text namelists (`data*`) are checked in. To reproduce a
runnable `input_ad/`:

- Generic forcing files (bathymetry, wind stress, SST/SSS climatology):
  symlink from `MITgcm/verification/tutorial_global_oce_biogeo/input/`,
  following the convention used by the sibling `global_oce_biogeo_bling`
  experiment.
- BLING tracer initial conditions (`dic_init.bin`, `alk_init.bin`, etc.),
  `ones_32b.bin`, `sample_prof.nc`: come from the standard
  `global_oce_biogeo_bling` verification experiment's own `input_ad/`.
- `mah_flux_smooth.bin`: same source.
- `socat_pco2_clim_month01..12.nc`: regenerate with
  `scripts/build_socat_clim.py` from the raw SOCAT source (see that
  script's docstring for provenance/verification details) -- the raw
  source file (`SOCATv2026_tracks_gridded_monthly.nc`, ~4GB, from
  [socat.info](https://socat.info)) is not included here either.
- `wt_*.bin` control weight files: regenerate with
  `scripts/make_ctrl_weights.py` (see "Known limitations" below).

## The model configuration

- **Grid**: 128x64, 2.8125 deg, 15 levels, 2x2 tile decomposition
  (`sNx=64, sNy=32`), single MPI process (`nPx=nPy=1`) -- see "Performance."
- **Calendar**: `startDate_1=20070101`, 360-day climatological year
  (12 x 30-day months, `nTimeSteps=34560` at `deltaT=900s`). This is a
  repeating climatology, not real inter-annual 2007 conditions.
- **Packages**: GM-Redi, ptracers, gchem+BLING, profiles, ECCO/cost/ctrl/
  autodiff/grdchk (offline adjoint-optimization stack), diagnostics.
- **Control vector** (`input_ad/data.ctrl`): 8 active 3D initial-condition
  fields, `mult_genarr3d(1:8)=1.0`:

  | control | field |
  |---|---|
  | `xx_theta` | Temperature IC |
  | `xx_salt` | Salinity IC |
  | `xx_ptr1` | DIC IC |
  | `xx_ptr2` | Alkalinity IC |
  | `xx_ptr3` | O2 IC |
  | `xx_ptr4` | NO3 IC |
  | `xx_ptr5` | PO4 IC |
  | `xx_ptr6` | Fe IC |

  445,384 total control variables (8 fields x wet ocean grid points x 15
  levels x 4 tiles).
- **Cost function** (`input_ad/data.profiles`): `J = J_data + J_prior`.
  `J_data = sum(weight * (model - obs)^2)` over 12 monthly SOCAT surface
  pCO2 climatology files (`socat_pco2_clim_month01..12.nc`), where `weight
  = 1/sigma^2` is precomputed into each file's `prof_PCOweight`.
  `J_prior = sum(control perturbation^2 * background_weight)` across the 8
  control fields, currently a uniform background weight (see "Known
  limitations").
- **Checkpointing** (`code_ad/tamc.h`): 3-level TAF checkpointing,
  `nchklev_1=2, nchklev_2=60, nchklev_3=300` (product 36,000 >= 34,560
  timesteps).

## Building

This experiment needs two separate binaries: the MITgcm adjoint model
(`mitgcmuv_ad`, built via TAF) and the offline M1QN3 driver (`optim.x`,
built from `optim_m1qn3/src/`, a sibling of this MITgcm checkout -- adjust
the path below to wherever that lives in your checkout).

```sh
mkdir build_ad && cd build_ad
../../../tools/genmake2 -mods=../code_ad -of=<your optfile>
make depend
make adall            # builds mitgcmuv_ad (forward + TAF-generated adjoint)
```

```sh
cd ../../optim_m1qn3/src
make                   # builds optim.x
```

Set up a run directory (e.g. `run/`) with symlinks to `mitgcmuv_ad`,
`optim.x`, everything in `input_ad/` (plus the externally-obtained binary
inputs above), and `scripts/run_optim_loop.sh`.

**Materialize `data.optim` as a real file, not a symlink.**
`run_optim_loop.sh` edits it in place every cycle (`sed -i`), and on
macOS `sed -i` refuses to edit through a symlink ("in-place editing only
works for regular files") -- it'll fail every cycle after the first
(which happens to look fine only because `optimcycle=0` is already
correct at cycle 0 by coincidence). After symlinking everything else, run:
```sh
cp -L data.optim data.optim.tmp && rm data.optim && mv data.optim.tmp data.optim
```

## Step 1: validate the adjoint (gradient check)

Before spending days on the optimization campaign, validate the
TAF-generated adjoint gradient against a finite-difference approximation
using `pkg/grdchk`, on a cheap 1-month test window:

1. In your run directory's `data`, temporarily set `nTimeSteps=2880` (1
   month).
2. In `data.pkg`, set `useGrdchk=.TRUE.`.
3. `input_ad/data.grdchk` is already configured for a 2-test-point check
   of `xx_ptr2` (alkalinity) -- edit `grdchkvarname` to check a different
   field (e.g. `xx_ptr1` for DIC).
4. Run `./mitgcmuv_ad` directly (not through `run_optim_loop.sh` -- grdchk
   is a single self-contained run, not an optimization loop).
5. Check the RMS relative error between the adjoint and finite-difference
   gradients in the log -- expect ~1e-6 (machine precision). See
   `results/output.txt` section 1 for the reference numbers (DIC and ALK,
   both ~1.3-1.5e-6).
6. Restore `nTimeSteps=34560` and `useGrdchk=.FALSE.` before proceeding to
   the real campaign -- **do not** run the production optimization with
   grdchk enabled.

## Step 2: run the M1QN3 optimization campaign

```sh
cd run/
nohup ./run_optim_loop.sh > optim_loop_master.log 2>&1 &
```

Each pass of the loop is one full forward+adjoint model evaluation
(`mitgcmuv_ad`, ~7.3hr) followed by one `optim.x` call (seconds) that reads
the fresh cost/gradient, does one M1QN3 quasi-Newton step, and writes the
next control-vector guess. `input_ad/data.optim` ships with
`optimcycle=0, numiter=10, nfunc=100` -- a proven-safe fresh-start budget
(10 iterations, generous line-search headroom) matching the reference
campaign in `results/output.txt`.

Monitor progress via:
```sh
grep "global fc " output_optim_*.txt   # cost function per cycle
tail -f optim_loop_master.log
```

### Continuing past an iteration-cap stop

When M1QN3 hits its `numiter` cap, it reports `omode=4` and the loop exits
-- this is a budget exhaustion, not necessarily gradient convergence
(check "realized relative precision on g" in `m1qn3_output.txt` against
your `epsg` target). To continue with a larger budget:

1. Bump `numiter` (and `nfunc` with generous headroom, e.g. 3-5x
   `numiter`) in `data.optim`.
2. Set `optimcycle` to the next cycle number.
3. **Add `coldStart=.TRUE.` to the `&M1QN3` namelist block.** This is
   easy to miss: `optimcycle==0` is the only case that auto-forces a cold
   start (`optim_readparms.F`). For any other `optimcycle`, the code
   defaults to a *warm* restart, which requires an `OPWARM.opt<cycle>`
   pickup file -- but M1QN3 only writes that file when it's continuing
   normally (`reverse==1`), not when it terminates via an iteration/
   function-eval cap. Without `coldStart=.TRUE.`, `optim.x` will stop
   abnormally looking for a pickup file that was never written.
4. After the first cycle of the new "leg" completes, **remove
   `coldStart=.TRUE.`** (or set it `.FALSE.`) before the next cycle --
   otherwise every subsequent cycle keeps discarding M1QN3's accumulated
   Hessian approximation instead of building on it.
5. If you're resuming a stopped/killed loop, **archive or truncate
   `m1qn3_output.txt` first.** It's opened in append-only mode across
   restarts, so a leftover `"m1qn3: output mode"` line from a previous
   leg will cause `run_optim_loop.sh` to falsely report convergence after
   the very first cycle of the new leg (the script only checks lines
   appended since the start of each cycle, but that only helps once the
   file no longer has a stale match sitting in already-checked history --
   truncating on restart is the simplest safe move). The script keeps the
   file *existing* but empty rather than deleting it, since M1QN3 opens it
   with `status='old'`, which requires the file to be present.

## Performance

A full evaluation currently runs **single-core**: `code_ad/SIZE.h` has
`nPx=1, nPy=1` (one MPI process), `nSx=2, nSy=2` (4 tiles processed
serially within that process), and `eedata` has `nTx=1, nTy=1` (no
threading). On a 14-core test machine this leaves ~13 cores idle for the
entire ~7.3hr/cycle campaign.

Checkpointing I/O and compiler flags are not the bottleneck: the build
already uses `-O3 -ftree-vectorize -funroll-loops`, and the ~16GB of
checkpoint tape + diagnostic I/O per cycle is negligible next to 7+ hours
of compute at any plausible disk throughput. Compute (the forward
integration, TAF checkpoint-recompute sweeps, and the adjoint sweep) is
what actually costs the time.

`code_ad/SIZE.h_mpi` stages half of an MPI conversion (`nSx=nSy=1,
nPx=nPy=2`, same tile boundaries, decomposed across processes instead of
within one) -- up to a ~4x speedup on a multi-core machine is plausible,
but this build has **no MPI support compiled in at all** (no
`ALLOW_USE_MPI`, plain `gfortran` not `mpif90`). Converting requires a
from-scratch `genmake2 -mpi` reconfigure + rebuild of both the forward and
adjoint binaries, and **re-running the grdchk validation** (Step 1) before
trusting the result -- TAF-generated halo-exchange/tape-indexing code for
a real multi-process decomposition is a plausible source of new bugs that
grdchk is specifically there to catch. Treat this as a separate,
validated side project, not a change to make mid-campaign.

## Known limitations

**Control-vector weighting -- fixed in this experiment, not yet in the
campaign it was extracted from.** The production campaign this
experiment documents ran with all 8 fields' `xx_genarr3d_weight` pointed
at a uniform `ones_32b.bin` rather than field-specific weight files.
Decoding the actual optimized control vectors from that campaign showed
raw adjustments spanning **7 orders of magnitude** across fields (e.g.
~1e-11 K for temperature vs. ~1e-5 for Fe), and M1QN3's descent direction
drifting to 75-78 degrees off the true gradient by iteration 8-10 -- a
classic poorly-conditioned, unscaled multi-field optimization symptom.
Practically, temperature, salinity, and O2 were likely frozen at their
first-guess values in all but name in that campaign, even though
`mult_genarr3d` marked all 8 fields "active": the cost function decreased
correctly, but the "8 active control fields" framing overstated what was
actually moving.

`input_ad/data.ctrl` in **this** experiment now points each field at its
own physically-scaled weight file instead (generated by
`scripts/make_ctrl_weights.py`):

| field | sigma | basis |
|---|---|---|
| theta | 1.0 degC | standard ECCO-style theta IC prior |
| salt | 0.2 psu | standard ECCO-style salt IC prior |
| DIC | 0.02 mol/m^3 | ~1-2% of the field's own ~2 mol/m^3 mean |
| ALK | 0.02 mol/m^3 | same basis as DIC |
| O2 | 0.015 mol/m^3 | field max ~0.48 mol/m^3 |
| NO3 | 0.003 mol/m^3 | field max ~0.049 mol/m^3 |
| PO4 | 0.0002 mol/m^3 | field max ~0.0036 mol/m^3, Redfield-consistent with NO3 |
| Fe | 5e-7 mol/m^3 (0.5 nM) | field max ~4.3 nM |

These are defensible, physically-grounded starting values (cross-checked
against this experiment's own initial-condition field statistics), **not
a rigorously calibrated prior covariance** -- if you have a better-informed
estimate of actual IC uncertainty for your application, use it instead.
Since this changes the optimization problem itself (not just a bugfix to
identical results), **the reference trajectory in `results/output.txt`
was generated under the old uniform weighting and will not match a rerun
under these new weights** -- expect theta/salt/O2 to now move
meaningfully cycle-to-cycle where they previously didn't, and the overall
cost-decrease trajectory to differ. Re-run the gradient check (Step 1)
after this change before trusting the campaign, same as after any other
control-vector change.

**SOCAT climatology mid-month dates drift for the second half of the
year.** The preprocessing in `scripts/build_socat_clim.py` originally
computed each month's nominal date by stepping 30 days at a time from
January 15th rather than using the true calendar date, drifting 2-5 days
by December. The version shipped here fixes this and has been **verified
against the raw source data** (`SOCATv2026_tracks_gridded_monthly.nc`):
running it reproduces the exact observation count and exact per-tile
maximum recorded from the original campaign, for all 12 months (see the
table in the script's docstring and `results/output.txt`). The
already-generated `socat_pco2_clim_month*.nc` files used by the reference
campaign (not included in this repo -- see "Directory layout") still
reflect the original, drifted dates -- regenerate them with the
corrected script if you want the dates fixed. Note this only reliably
fixes `prof_date`; `prof_PCOweight` uses a formula (`1/sigma**2` from
`fco2_std_weighted`, with a floor) that matches by recollection, not by
independently reproducing the shipped files' byte values, so it may
differ slightly for very-low-observation-count cells -- see the script's
docstring.

## Reproducing the figures

`scripts/make_diagnostic_figures.py` (gradient check validation, monthly
cost breakdown, adjoint sensitivity maps) and
`scripts/make_optim_progress_figures.py` (cost function convergence,
optimized control-adjustment maps) regenerate `figures/*.png`. Both read
paths from environment variables (`BLING_RUN_OPT_DIR`, `BATHY_BIN`,
`FIGURES_OUT_DIR`) defaulting to the layout used for the reference
figures shipped here -- point them at your own run directory and pull
fresh numbers from your own logs rather than reusing the reference
campaign's hardcoded arrays (see the note at the top of each script).
