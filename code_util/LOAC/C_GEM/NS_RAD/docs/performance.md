# Speeding up NS-RAD — what was done and how it was verified

Result: **~9× faster overall, with bit-identical output.** On the 2-day verification
case, runtime fell from 37.2 s to 3.97 s (9.4×); per-site runtime for the shipped
2-year configuration dropped from **192 min to about 26 min**.

| Phase | Before | After | Speedup |
|---|---|---|---|
| Year 1 — hydrodynamics + transport | 87 steps/s | 1047 steps/s | **12.0×** |
| Year 2 — adds biogeochemistry + sediment | 63 steps/s | 370 steps/s | **5.9×** |
| Per site (840,960 timesteps) | 192 min | 26 min | **7.5×** |

Measured in two stages: the hydrodynamic kernels (2a) gave 6.1×, the density stack
(2b) a further 1.53×. All timings taken with four full simulations running
concurrently, so they are conservative.

Everything below is reproducible: `python -m cProfile -o prof.out main.py` with
`CGEM_MAXT_DAYS` / `CGEM_WARMUP_DAYS` set to isolate a phase.

## 1. Measure first

The starting estimate was guesswork and it was wrong in both directions — a naive
projection said 3.2 h, the first minutes of a real run looked like 40 h, and the
truth was 2.4 h. Timing a real run settled it:

```
CGEM_MAXT_DAYS=1 CGEM_WARMUP_DAYS=5   # never reaches WARMUP -> phase 1 only
CGEM_MAXT_DAYS=1 CGEM_WARMUP_DAYS=0   # WARMUP already passed -> phase 2 every step
```

Then `cProfile` on 1152 timesteps of phase 2 (total 29.5 s):

| tottime | cumtime | ncalls | function |
|---:|---:|---:|---|
| 6.53 | 6.90 | 71,375 | `tridag_module.coeff_a` |
| 4.76 | 4.76 | 71,375 | `tridag_module.tridag` |
| 4.48 | 4.48 | 71,375 | `uphyd_module.update` |
| 2.26 | 3.61 | 142,750 | `tridag_module.conv` |
| 1.77 | 1.77 | 783,360 | `density._seck` |
| 1.56 | 8.25 | 1,152 | `biogeo_module.biogeo` |
| 1.00 | 1.00 | 14,648,102 | `builtins.abs` |
| 0.99 | 0.99 | 783,360 | `density._dens0` |
| 0.67 | 2.99 | 470,016 | `fun_module.K1_CO2` |
| 0.47 | 3.23 | 783,360 | `density.dens` |
| 0.39 | 1.89 | 313,344 | `fun_module.K2_CO2` |
| 0.39 | 3.62 | 783,360 | `fun_module.p_bar` |

Two things fall out of the call counts, and both matter more than the times:

**The hydrodynamic solver is 61% of runtime** (6.53+4.76+4.48+2.26 = 18.0 s of 29.5).
It is called 71,375 times for 1,152 timesteps — **62 solver iterations per timestep**.
`hyd` loops `coeff_a → tridag → conv → update` until `rsum == 2.0`.

**The carbonate constants are recomputed redundantly.** `density.dens` runs 783,360
times = 680 per timestep = **5× per grid point**, because `biogeo` called
`K1_CO2` three times and `K2_CO2` twice per cell, each one going through
`p_bar → dens`. One call each suffices.

### An aside worth recording

Measuring the solver directly under both geometries:

```
new geometry (D ≈ 1–2.3 m):  44.2 iterations/timestep
old geometry (D = 15 m):     21.5 iterations/timestep
```

The observation-based depth correction **doubled the cost of the model**. Shallower
channels converge more slowly. That is physics, not a defect — but it means the
geometry fix and the need for this optimization are causally linked.

## 2. What was changed

### 2a. Hydrodynamic kernels → numba (`tridag_module.py`, `uphyd_module.py`)

`coeff_a`, `tridag`, `conv`, `update` and `new_uh` are now `@njit(cache=True)`
kernels. These are plain numeric loops over ~136 points with no Python objects, so
they are ideal JIT targets.

Two implementation details that are not optional:

- **The arrays had to become numpy.** `variables.py` now allocates the hydrodynamic
  state (`Y, E, ZZ, D, Dold, Dold2, DEPTH, Chezy, H, U, B, TH, TU, Z, C`) as
  `np.float64` instead of Python lists. Numba's support for reflected lists is slow
  and deprecated. Item assignment (`arr[i] = x`) behaves identically, so surrounding
  code is untouched — but **these must never be rebound to lists**, or the
  mutate-in-place contract that the whole model depends on breaks silently.
- **The arrays are passed as arguments, not read as globals.** Numba treats module
  globals as compile-time constants, so a kernel closing over `D` would bake in the
  contents at first compile and silently stop seeing updates. Each public function
  keeps its original signature and forwards the arrays explicitly:

  ```python
  def coeff_a(t, Qr):
      _coeff_a(Z, C, D, U, B, H, TH, Chezy, DEPTH,
               M1, M2, M3, DELTI, DELXI, G, RS, Qr)
  ```

`C` became a 2-D array; `C[j][1]` and `C[j, 1]` index it equivalently, so call sites
did not change.

### 2b. Density stack → numba (`density.py`)

`_dens0`, `_seck` and `dens` — the UNESCO-1983 density routines — are `@njit(cache=True)`.
Pure scalar arithmetic, no Python objects. `dens()` runs twice per grid point per
timestep through `fun_module.p_bar`, and profiling put the stack at ~18% of runtime
once the hydrodynamic kernels were compiled. Worth **1.53×** on its own, bit-identical.

The default arguments changed from `P=0` to `P=0.0` so numba compiles one float
specialisation rather than separate int and float ones. `0.1*0` and `0.1*0.0` are both
`0.0`, so nothing changes numerically. The unused routines (`svan`, `sigma`, `drhodt`,
`alpha`, `drhods`, `beta`) are left as plain Python and can still call the compiled
`dens()`.

### 2c. Hoist the carbonate constants (`biogeo_module.py`)

```python
k1 = K1_CO2(t, i, water_temp)
k2 = K2_CO2(t, i, water_temp)
```

computed once per grid point and reused in the `pH` call, the CO2* expression and
the flux. Both depend only on `(i, water_temp)`, which are fixed within the
iteration, so the values and the order in which they are used are unchanged. This
removes 3 of every 5 `dens()` calls.

## 3. Verification — bit-identity, not "looks right"

Before touching anything, a reference run was captured:

```bash
CGEM_MAXT_DAYS=2 CGEM_WARMUP_DAYS=1 CGEM_SITE=kuparuk python main.py
md5 *.dat | md5     # -> f2489efced2621a0f60e1bd72b9713c2
```

`WARMUP=1` of a 2-day run means the window exercises **both** phases — hydrodynamics
and transport throughout, biogeochemistry and sediment on day 2 — so all changed
code paths are covered. After the optimization the same command produced the same
hash across all 22 `.dat` files:

```
reference : f2489efced2621a0f60e1bd72b9713c2
optimized : f2489efced2621a0f60e1bd72b9713c2
```

This is the property that makes the work safe to adopt: results from before and after
are directly comparable, because they are the same bits. Any future change here
should be held to the same bar — re-capture the reference hash first, then diff.

**Two gotchas when re-running this harness today:**

1. **`CGEM_OUTPUT=dat`.** The default output format is now NetCDF (`output.nc`), so
   `md5 *.dat` finds nothing. Force the legacy `.dat` files for the hash:
   `CGEM_MAXT_DAYS=2 CGEM_WARMUP_DAYS=1 CGEM_SITE=kuparuk CGEM_OUTPUT=dat python main.py`.
   (NetCDF is a poor hash target anyway — it can embed nondeterministic metadata.)
2. **Clear `__pycache__` before capturing the reference.** The jitted kernels use
   `@njit(cache=True)`, which persists compiled code across runs in `__pycache__/*.nbc`.
   A cache left over from an earlier session — a different numba/llvm build, or code that
   has since changed — can be loaded on the *first* run and differ from a fresh compile
   by a ULP on some site's data (this bit the Tier-1 verification below: kuparuk hashed
   `02a4a142…` on a stale cache and `28cb7965…` on every clean run). `rm -rf __pycache__`
   before capturing so the reference reflects a fresh compile, not a stale artifact.

The most robust form of the check does not pin an absolute hash at all (the hash tracks
config, geometry, chemistry and boundary forcing, all of which change under active
development). Instead, A/B the two code versions **at the same commit of everything
else**: revert only the change under test to its pre-jit Python form in a copy, run both
copy and current with a freshly-cleared cache, and diff the hashes. That isolates the
change from any concurrent edits to `main.py`/`config.py`/forcings.

Coverage limit, stated honestly: the check covers days 0–2, not the freshet peak
(~4147 m³/s on Colville) or the year-1→year-2 transition. The changes are structural
rather than numerical, so bit-identity should hold generally, but it has not been
demonstrated across those regimes.

## 4. What was deliberately NOT done

- **`fastmath=True` — tried, measured, reverted.** Enabled on every jitted kernel and
  benchmarked with warm caches, three runs each:

  | | runs | median |
  |---|---|---|
  | fastmath | 3.72 / 3.79 / 3.79 s | **3.79 s** |
  | no fastmath | 3.81 / 3.86 s | **3.83 s** |

  **~1%, inside the noise floor** — and it changed the results (`fd765dbf…` vs
  `f2489efc…`). The reason is structural: fastmath only affects compiled code, and
  after the kernels above are jitted the remaining time is in pure-Python `biogeo`
  (24.7%), the transport scheme and per-cell call overhead, none of which it can
  reach. So it buys ~1% and costs the bit-identity check that has caught real
  mistakes.

  There is also a model-specific hazard. `hyd` loops on `while rsum != 2.0` — an exact
  float comparison — and `conv` decides convergence on `abs(diff) >= 1e-10`.
  Reassociation can perturb `diff` near that threshold, changing iteration counts and
  therefore results; in the worst case it changes whether convergence is reached at
  all, and that loop has no iteration cap. fastmath additionally asserts no-NaN and
  no-inf, which is an assumption rather than a guarantee — the model divides by
  `DEPTH` and `B`, and bottom-fast ice or zero flow could violate it.

  Re-enable with:
  `sed -i '' 's/@njit(cache=True)/@njit(cache=True, fastmath=True)/g' density.py tridag_module.py uphyd_module.py`
  — but re-measure before adopting, and accept that bit-identity is then gone.
- **Relaxing `TOL` (1e-10) or the `while rsum != 2.0` convergence test.** Fewer solver
  iterations is the largest theoretical win, since the loop runs 44× per timestep.
  But it changes results, so it is a modelling decision, not an optimization.
  (Note the exact float comparison `rsum != 2.0` remains a hang risk — see CLAUDE.md.)
- **Reducing the per-timestep `print`.** ~840k lines per run. Measured as negligible
  against the numeric cost, and it is the only progress signal a long run has.
- **Reducing output volume.** `TS = 12` writes every 900 s → 70,080 rows per file,
  ~4.6 GB per site. Raising `TS` would cut this substantially but changes what the
  run produces.

## 5. What is left

**Both items originally listed here are now done** (biogeo, tvd/disp_sch — see §7), and
so is a further pass on the Python orchestration layer. Current status and the next
targets are in §7 and §8.

## 6. Environment

`numba 0.66.0` + `llvmlite 0.48.0`, installed into the miniforge Python 3.13
environment (cp313 wheel, requires numpy < 2.5; numpy 2.4.3 present). `numba` was
already a documented dependency of this model — `fun_module.pH` has always been
`@njit` — it simply was not installed. First compile costs ~0.1 s per kernel;
`cache=True` persists compiled code across runs.

## 7. Third pass — the Python orchestration layer (Tier 1 + Tier 2)

By this point every per-cell numeric loop was already jitted — the §5 targets (biogeo,
`tvd`/`disp_sch`) had been done, as had `heat`, `ice` and `pH`. A fresh `cProfile`
(discounting one-time compile) then showed the cost had moved off the arithmetic and
onto **Python orchestration and a few un-jitted stragglers**:

| tottime (s) | function | nature |
|---:|---|---|
| 0.75 | `transport.transport` | pure-Python driver + pure-Python `avg` loop |
| 0.31 | `sed_module.sed` | not jitted |
| 0.22 | `tridag` (jitted) | ~64 hydro iters/step → dispatch-bound |
| 0.18 | `fun.river_dispersion` | not jitted, per-cell, every step |
| 0.16 | `fun.piston_velocity` | not jitted, per-cell |
| 0.07 | `fun.Sc` | not jitted |
| 0.04 | `posix.stat` | `forcing()` `os.path.exists`, 8×/step |

**Tier 1 (bit-identical, low-risk):**
- `transport_module` — the per-species average accumulation was a pure-Python
  `for i in range(1,M+1): avg[i]+=c[i]` loop (~1.5 B iterations/run); now `avg[1:] +=
  c[1:]`.
- `sed_module` — `sed` is now an `@njit _sed_loop` (verbatim, `wISS` inlined), same
  no-cache/read-globals convention as biogeo.
- `fun_module` — `piston_velocity` and `river_dispersion` now wrap `@njit` kernels
  (`_piston_velocity_loop`, `_river_dispersion_loop`, with `_sc`/`_d_o2` helpers); the
  Python originals are retained.
- `main` — the 8 forcing paths are resolved once before the time loop instead of
  `os.path.exists`-ing them every timestep.

**Tier 2 (bit-identical, moderate):**
- `hyd_module` — the convergence loop (`coeff_a → tridag → conv → conv → update`, ~64
  passes/step, each previously ~5 Python↔numba boundary crossings) is fused into a
  single `@njit _hyd_iterate` driver that calls the already-jitted cores natively,
  removing ~320 dispatches/step. The iteration count and final `rsum` are returned so
  the once-per-run cap warning stays in Python.

**Results** (kuparuk, warm cache, one process):

| run | original | Tier 1 | Tier 1+2 |
|---|---:|---:|---:|
| 20-day phase-2 (biogeo every step) | 37.5 s | 18.1 s (2.07×) | 16.2 s (**2.32×**) |
| 10-day phase-1 (hydro/transport only) | — | 5.94 s | 4.75 s (**1.25×** from Tier 2) |

**Verified bit-identical** on all four sites by the A/B method in §3 — reverting only
these changes to their Python form in a copy and diffing the two hashes under a fresh
cache, all else (including a concurrent `main.py`/config change adding config-driven
forcing filenames and `BOUNDARY_FORCING`) held equal. Every site matched
(orig == modified), file-by-file: transport, sed, fun and hyd each genuinely differ from
their Python form, `main.py` identical in both. This pass also uncovered the §3
stale-cache trap, which had produced a misleading one-off "difference" during Tier 1.

## 8. What is left now

1. ~~**Fuse the transport per-species loop.**~~ **Retired — not worth doing on its own.**
   Two things overtook it: the `tvd` `cold = co` aliasing defect is already fixed
   (`schemes_module.py:41`, `cold = co.copy()`, bit-identical), and Tier 1's
   avg-vectorization already removed transport's dominant cost. A fresh warm-cache
   profile puts `transport` at ~1.6% tottime (was the #1 item at 0.75 s); the whole
   transport family (transport + `openbound` + the `tvd`/`disp_sch` wrapper dispatch) is
   ~5.5%, of which a fusion recovers maybe ~4% (~1.04×). It also folds `openbound`'s
   `v[s]["clb"]/["cub"]` access into a kernel, which collides with the active
   `BOUNDARY_FORCING` work. **Decision (this session): fold it into (2) instead of doing
   it standalone.**
2. **A single jitted `step()` driver.** With every component now a jitted kernel, the
   whole per-timestep sequence (hyd → transport → heat → ice → biogeo → sed) could run
   inside one `@njit` step, leaving `main.py` to do only forcing interpolation and I/O.
   This subsumes the transport fusion and also removes the `hyd` (~11%) and `biogeo`
   (~10%) *wrapper* tottimes — the real remaining lever. Largest ceiling, largest risk.
   **Do this after the `BOUNDARY_FORCING`/boundary-refresh machinery has stabilized**,
   since the fused step must call that code; re-profile and snapshot first, verify by the
   §3 A/B method.
3. **I/O write batching** (independent, safe now). `file_module.write` + `close` are ~6%
   of warm runtime (row-by-row appends). Buffering rows and flushing in chunks is
   bit-identical and touches only `file_module` — no conflict with the boundary work.
4. **Reduce the ~64 hydro iterations/step** (warm-start from the previous step, relax
   `TOL`). The single largest remaining *compute* lever, but it **changes results** —
   validate physically, not by MD5.
