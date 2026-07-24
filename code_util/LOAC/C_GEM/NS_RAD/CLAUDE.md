# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this is

**NS-RAD — the North Slope River-Aquatic-Delta Model.** A C-GEM-based 1-D reactive-transport model of the four
main Alaskan North Slope rivers (Colville, Kuparuk, Sagavanirktok, Canning), with carbonate chemistry, air–sea
CO2/CH4/N2O exchange, a transported temperature field + surface heat budget, prognostic river ice, and an
optional Arctic biogeochemistry extension. Built on **C-GEM** (Volta et al. 2014), the underlying estuarine
reactive-transport engine; most modules still carry `translated from <name>.c` docstrings and mirror the
original C file layout one-for-one. "NS-RAD" is the project/model name; "C-GEM" refers to that engine lineage.

Vendored from `MITgcm-contrib/ecco_darwin` at `code_util/LOAC/C_GEM/North Slope`. These are plain files, not a
git checkout — there is no upstream remote to pull from.

## Running

The model is configured for **four North Slope rivers** — Colville, Kuparuk, Sagavanirktok, Canning. Each
runs as its own process from its own output directory:

```bash
tools/run_sites.sh                    # all four, in parallel -> runs/definitive/<site>/
tools/run_sites.sh kuparuk            # one site
SERIAL=1 tools/run_sites.sh           # one at a time
```

`tools/build_all.sh` is the one-command wrapper around the whole product: run the four rivers, then rebuild
every figure and stitch `docs/ns_rad_report.pdf`. Use it, not `make_report.py`, to regenerate the report — the
six section PDFs are throwaway intermediates it creates, stitches and deletes. A failing figure does not abort
the rest; it prints a PASS/FAIL summary and exits non-zero.

```bash
tools/build_all.sh                    # run 4 rivers + all figures + report + movies (~15 min)
tools/build_all.sh --figures-only     # rebuild figures from existing runs/
tools/build_all.sh --with-idealized   # also run + verify + figure the idealized fixture
tools/build_all.sh --with-regression  # also build runs/regression_bnd (+~20 min; see below)
```

**The report needs a SECOND set of runs to be complete.** `make_validation_pdf.py` reads
`runs/regression_bnd/`, which holds Kuparuk and Sagavanirktok rerun with the pure air-temperature regression
as the upstream T boundary instead of their own observed-blended series — otherwise their modelled temperature
would be scored against the very USGS record that feeds their boundary, which is near-tautological. Colville
and Canning already use the regression, so they are not rebuilt. Without these runs the report silently builds
**5 of 6 parts** and the entire model-vs-observation section is missing. Produce them once with
`tools/run_regression_bnd.sh` (~20 min; it verifies the boundary override actually took effect) or pass
`--with-regression`. They only go stale when the temperature forcing changes, which is why they are not part
of the default build.

Or directly, which is what the script does per site:

```bash
mkdir -p runs/definitive/kuparuk && cd runs/definitive/kuparuk
CGEM_SITE=kuparuk PYTHONPATH=../../../code python ../../../code/main.py
```

`runs/` is organized by **run-type, then river**: `runs/definitive/<site>/` are these operational runs (all
analysis tools read these); `runs/regression_bnd/<site>/` is the independent-boundary check consumed only by
`make_validation_pdf.py`; `runs/idealized/` is the verification fixture (see the Idealized verification site
section). Three archives hold the states two recent changes were measured against, each isolating **one**
variable: `runs/singlechannel_archive/<site>/` is pre-multichannel, `runs/multichannel_preicefix/<site>/` is
post-multichannel but pre-ice-fix, and `runs/definitive/` is current. So
`tools/compare_adoption.py` diffs the first two (width change alone, both with the old ice) and
`tools/compare_icefix.py` diffs the last two (ice change alone, both with the new width) — do not compare
`singlechannel_archive` directly against `definitive`, which would conflate them.
`runs/experiment_multichannel/` holds the multichannel sensitivity sweep and the distributary diagnostic. The
old flat `runs/<site>`, `runs_baseline/` and `runs_regression_bnd/` layout is retired.

Sites **must** run in separate directories: `main.py` deletes every `*.dat` in the cwd at startup. They also
must be separate processes — `variables.py` allocates all state at import time from `config.M`, so two sites
cannot coexist in one interpreter without a deep refactor. Since the runs are independent and slow, parallel
processes are the right answer anyway.

Smoke test without editing `config.py` (it used to require hand-editing `MAXT`, which was easy to leave in):

```bash
CGEM_MAXT_DAYS=2 CGEM_WARMUP_DAYS=1 tools/run_sites.sh
```

No build, no linter config, no dependency manifest, and no unit tests — **`tools/verify_idealized.py` is the
closest thing to a test suite** (see the Idealized verification site section) and is what to run after a change
that could affect physics. Requires `numpy`, `numba`, and `netCDF4` (output is NetCDF by default);
`matplotlib` for the figure tools, `ffmpeg` for movies, `python-docx` only to rebuild the Word guide.

Everything tunable at run time is an environment variable, read in `config.py` at import:

| var | default | effect |
|---|---|---|
| `CGEM_SITE` | `colville` | which `code/sites/<name>.py` to load |
| `CGEM_MAXT_DAYS` / `CGEM_WARMUP_DAYS` | 730 / 365 | run and warmup length — the smoke-test knobs |
| `CGEM_OUTPUT` | `nc` | `nc`, `dat`, or `both` (see Output, and the bit-identity caveat below) |
| `CGEM_TS` | 12 | output save cadence (saves every `TS*DELTI` seconds) |
| `CGEM_ICE` | `on` | `off` reverts every coupling to the legacy `previousdays` gate |
| `CGEM_MULTICHANNEL` | `on` | `off` reverts to legacy single-channel width, bit-identically |
| `CGEM_N_CHAN_UP` | derived | override the prismatic thread count, for sensitivity sweeps |
| `CGEM_WATERTEMP_FILE` | per-site | override the upstream T boundary; builds `runs/regression_bnd` |
| `CGEM_DISTANCE` | per-site | saline grid points; swept and rejected as a salinity lever (see Geometry) |
| `CGEM_FONT_SCALE` | 1.0 | figure tools only (`tools/nsrad_style.py`) |

`numba` **is installed** (`~/miniforge3`, Python 3.13, numba 0.66.0 / llvmlite 0.48.0, numpy 2.4.3) and is now
load-bearing: essentially every per-cell hot loop is `@njit` — the hydro kernels (`tridag_module`,
`uphyd_module`, and the fused `hyd_module._hyd_iterate`), `density`, `biogeo`, `schemes` (`tvd`/`disp_sch`),
`sed`, `heat`, `ice`, `pH`, and the scalar carbonate/rate/piston/dispersion helpers in `fun_module`. `main.py`
fails at import if numba is missing.

`main.py` deletes every `*.dat` in the **current working directory** at startup and writes outputs there, so
always run it from the directory you want results in. (The upstream one-variable `plot outputs.py` scratch
script has been deleted — the `tools/make_*.py` figure generators replaced it and read `output.nc`.)

At the shipped settings (`DELTI=75`, `MAXT=2 years`) this is ~840k timesteps over 136 grid points. Runtime
has come down in three optimization passes (192 min → hydro/density jit → per-cell jit → Python-orchestration
jit), each held to bit-identical output; see `docs/performance.md` for the full history, profiles, numbers and
the verification harness. `docs/performance_metrics.json` holds the latest measured per-site wall times.

**Bit-identity harness caveat (learned the hard way — see `performance.md` §3):** capture the reference with
`CGEM_OUTPUT=dat` (default output is now NetCDF) **and after `rm -rf __pycache__`**. The jitted kernels are
`@njit(cache=True)`, so a stale cross-session `.nbc` cache can be loaded on the first run and differ from a
fresh compile by a ULP on some site's data — a false "difference" that is a cache artifact, not a code change.

Two consequences of that work to be aware of when editing:

- **The hydrodynamic arrays in `variables.py` are numpy, not lists** (`Y, E, ZZ, D, Dold, Dold2, DEPTH,
  Chezy, H, U, B, TH, TU, Z, C`). Item assignment behaves the same, but never *rebind* them to lists — the
  jitted kernels require typed arrays, and the mutate-in-place contract would break silently.
- **The jitted kernels take arrays as arguments, not globals.** Numba freezes module globals as
  compile-time constants, so a kernel closing over `D` would stop seeing updates. Keep the
  thin-wrapper pattern in `tridag_module` / `uphyd_module` when adding more.

Note what a smoke test does and does not exercise: biogeochemistry and sediment are skipped entirely until
`t > WARMUP`, so a run shorter than `WARMUP` covers only hydrodynamics, transport, heat and ice.

## Idealized verification site

A fifth site, `idealized`, is a **synthetic verification fixture, not a river** —
`CGEM_SITE=idealized`. It exercises the whole model on a known, analytic, **time-varying** input over a
seasonal cycle, so a change that quietly breaks transport, the carbonate solve, the ice model or the boundary
machinery is caught by an invariant check. `config.IS_IDEALIZED` warns on every idealized run so its output is
never read as a result. Full write-up in `docs/idealized_verification.md`.

It is also the first user of a **general, file-driven time-varying boundary mechanism**. The four real rivers
hold every solute boundary constant (only temperature is seasonal, special-cased in `main.py`); a site can now
declare `BOUNDARY_FORCING = {species: {"clb"|"cub": filename}}` and `main.py` refreshes those boundaries every
timestep from a 365-day series — the same way temperature already is, exploiting the fact that
`schemes_module.openbound` reads `v[s]["clb"]/["cub"]` fresh each call. Likewise the previously-hardcoded
meteorology filenames (`WIND_FILE`, `SOLAR_FILE`, `AIRTEMP_FILE`, `RELHUM_FILE`, `PCO2_FILE`, `SEATEMP_FILE`)
are now per-site overridable. **The four real rivers are byte-for-byte unaffected:** they set none of these
attributes, so `config.py`'s `getattr` defaults reproduce the old hardcoded met files and an empty
`BOUNDARY_FORCING`. Defaults and the contract live in `sites/_baseline.py`.

```bash
python tools/build_idealized_forcings.py            # (re)generate the 13 analytic forcing/boundary CSVs
python tools/verify_idealized.py --check-only       # structural checks, no run (no numba needed)
python tools/verify_idealized.py                    # + short 5-day run, physical invariants (~1-2 min)
python tools/verify_idealized.py --full             # + 220-day seasonal run: freshet ice break-up,
                                                    #   DOC pulse, summer freshening, boundary tracking
```

Everything about the fixture is defined by one `PARAMS` block in `tools/build_idealized_forcings.py`, so it is
reproducible (the builder rewrites the CSVs byte-for-byte). To change what the fixture tests, edit `PARAMS` and
rebuild; to add a real river with a seasonal boundary, drop in a CSV and add a `BOUNDARY_FORCING` entry.

## Multi-site configuration

`config.py` selects a river at import time from `CGEM_SITE` (default `colville`) and pulls its parameters from
`code/sites/<name>.py`. Only what genuinely differs between rivers lives there — channel
geometry, tidal amplitude, the discharge filename, and the `BOUNDARIES` table of 16 species × (downstream,
upstream) concentrations. Everything else — Chezy, sediment, the whole reaction network, numerics — stays
shared in `config.py`, so a change to the biogeochemistry applies to all four sites at once.

Sites inherit from `sites/_baseline.py` and override. `init_module` no longer hardcodes the 12
`initialize_substance` calls; it loops over `config.BOUNDARIES`.

### Geometry — now observation-based

Populated from **SWORD** (width) and **USGS channel surveys** (depth). All sites declare
`GEOMETRY_IS_PLACEHOLDER = False`. Rebuild discharge with `tools/fetch_discharge.py`; the SWORD delta-width
recovery is `tools/extract_sword.py` → `docs/sword_widths.json` (see "delta-mouth width" below).

`B_lb` is the **seaward-boundary (delta-mouth) width = the raw SWORD v17b distributary SUM**; `B_ub` is the
per-channel upstream prismatic width; the flare converges `B_lb → B_ub` over `L_FLARE`.

| Site | B_lb (delta sum) | B_ub | L_FLARE | DEPTH | depth from | surveys |
|---|---|---|---|---|---|---|
| Colville | 1550 m (1310+240) | 423 m | 5.5 km | 2.25 m | `0.360·Q^0.297` | 208 |
| Kuparuk | 516 m (335+182) | 58 m | 7.0 km | 1.34 m | `0.291·Q^0.309` | 323 |
| Sagavanirktok | 553 m (320+233) | 102 m | 7.0 km | 0.98 m | `0.280·Q^0.259` | 189 |
| Canning | 797 m (562+236) | 64 m | 7.0 km | 1.11 m | `0.224·Q^0.325` | 28 |
| *shipped* | *1215* | *859* | — | *15.00 m* | — | — |

**Depth replaced 15 m with ~1–2.3 m.** The shipped 15 m was cited to the Sagavanirktok gauge, but that
gauge's own 189 ADCP surveys span 0.5–1.9 m (median 0.88 m) — the citation did not hold. Depths are now
at-a-station hydraulic geometry `D = c·Q^f` evaluated at each river's **2022 open-water mean discharge**
(the ice gate is open days 140–303, so that is the period the geometry should represent). Since
`FCO2 = RCO2/DEPTH`, this is a first-order change to the carbon flux — a ~7× increase in flux per unit
volume relative to the old value. **Caveat:** the Colville and Sagavanirktok gauges are far upstream of their
deltas, so these depths likely *understate* the estuarine channel.

**Width required a braiding correction.** Raw SWORD width made Colville and the Sagavanirktok appear to
*widen* upstream, which would have inverted the estuary funnel. They are not widening — they are braided:
Colville is single-channel at the mouth but 100% multi-channel (3–4 braids) upstream, and SWORD reports the
*sum* across braids while C-GEM models one conveyance channel. Widths are therefore
`SWORD width / n_chan_mod`, lake nodes excluded (`lakeflag == 0`), restricted to the model domain. SWORD
*nodes* (~200 m spacing) were used rather than reaches (~10 km), because node spacing matches `DELXI` and
gives ~136 points across the domain — one per model grid point. This per-channel width is what `B_ub`
(upstream prismatic) uses.

**The delta-mouth `B_lb` uses the distributary SUM, not the per-channel width — because a delta is not a
braid.** Braids split and rejoin, so summing over-counts one conveyance channel (and the sum grows upstream,
inverting the funnel — hence the per-channel correction above). But at the **delta**, the distributaries
*diverge to separate sea outlets and never rejoin*, conveying the total discharge to the coast in **parallel**
— so the estuarine/salt-exchange width is genuinely the SUM. Using the per-channel width there (with total Q)
over-estimated the mouth velocity ~`n_chan`-fold and flushed the salt intrusion to ~0. `tools/extract_sword.py`
recovers the sum from a fresh **SWORD v17b** extraction: it sums one width per distinct `reach_id` in a ~3 km
seaward band at each delta (SWORD traces the *main* channel narrow, and for Colville/Kuparuk stops 30–50 km
inland — the deltas below are unnamed `NODATA` distributaries, which the reach-based sum still captures). All
four deltas are wide two-channel systems (table above); the earlier one-off extraction had traced only the
main channel and under-represented the mouth. **Impact:** shifts velocity (mouth 0.05–0.28 m/s, well-validated
vs USGS channel measurements — see `docs/ns_rad_validation.pdf`) and the carbon flux, but **does not fix
the summer salinity** — that is bounded by the forcing, not the width: the real harmonic tide is genuinely
microtidal (~0.3 m, not `AMPL`, which is an unused fallback), so summer high flow flushes the mouth regardless
of channel width. The intrusion is otherwise a *winter* phenomenon (~25 PSU under ice, low flow). The
saltwater-intrusion driver this coast actually has — wind-driven storm surge — is now forced on all four rivers
and still does not close the gap (see "Wind-driven storm surge" below), confirming the mismatch is event-timing /
year-specific, not geometric. **Three levers have now been tested and all fail**: mouth width (SWORD sums),
storm surge, and the `distance` / Chezy-dispersion gradient — a `CGEM_DISTANCE` sweep from 1 to 68 moves the
Kuparuk summer mouth salinity by <2 PSU (median stays 0 against observed 19), so `distance = 1` is kept. The
salinity is forcing-limited (microtidal + summer flushing + seaward-lagoon sampling), not tunable in the
channel. `SWORD_MOUTH_SUM` in each `sites/<name>.py` holds the mouth width; the network is mapped on Sentinel-2
imagery in `docs/ns_rad_river_networks.pdf`.

**`B_ub` is now the SWORD per-channel median in the prismatic reach** (>7 km, main stem, braided total /
`n_chan_mod`), recovered by `tools/extract_sword.py`: Colville 423, Kuparuk 58, Sagavanirktok 102, Canning 64 m.
This replaced two borrowed/mismatched values — Sagavanirktok's `B_ub` had been *borrowed from the Canning*
(51 m) because SWORD showed its per-channel width widening upstream over a 5–20 km data gap with a flare fit of
R² = −0.65; the fresh v17b per-channel median (102 m) is a direct observation, so only the flare *length*
`L_FLARE` is still borrowed. Canning's `B_ub` dropped 132→64 m to match its per-channel scatter. The
`docs/ns_rad_geometry.pdf` width page overlays the full per-channel node scatter + median/IQR against the
flare, with the delta-sum `B_lb` marked at the mouth — the prismatic `B_ub` now sits inside the SWORD scatter
for all four.

## Width and dispersion no longer use the Savenije estuary formulation

Two coupled changes, both driven by the SWORD data. Either can be reverted via `config.py`.

### `WIDTH_MODEL = "flare"` (was: whole-domain exponential)

Width converges exponentially over the first `L_FLARE` metres, then is **prismatic** at `B_ub`. Implemented
in `init_module.width_at()`; set `WIDTH_MODEL = "expo"` to recover the original law.

The original `B(i) = B_lb·exp(-i·DELXI/LC)` fits these rivers badly — R² 0.02–0.23. The measured profiles
show a sharp drop over 4–7 km then a quasi-prismatic channel (Kuparuk: 98 → 120 → 62 → flat ~55 for 30 km).
Formal comparison on 1 km binned log-width, AIC (lower better):

| river | exponential | constant | flare+prismatic | best |
|---|---|---|---|---|
| colville | R² 0.06 / −28.3 | −29.0 | **R² 0.39 / −35.4** | flare |
| kuparuk | R² 0.23 / −74.2 | −68.9 | **R² 0.43 / −80.7** | flare |
| canning | R² 0.15 / −65.6 | −63.0 | **R² 0.35 / −70.9** | flare |
| sagavanirktok | R² 0.46 / −26.7 | −23.2 | (9 bins, gap) | untrustworthy |

Note that for Colville a *constant* width also beat the exponential — the fitted convergence was adding
nothing. The flare now converges from the SWORD delta-mouth sum to the per-channel upstream width: Colville
1550→423 m over 5.5 km; Kuparuk 516→58 over 7 km; Canning 797→132 over 7 km; Sagavanirktok 553→51 over 7 km
(upstream shape borrowed). The AIC fits above predate the delta-sum `B_lb`; they justified the flare *shape*,
which is unchanged — only the seaward endpoint was raised to the distributary sum.

### `DISPERSION_MODEL = "seo"` (was: Van der Burgh / Savenije)

`fun_module.river_dispersion()` computes Seo & Cheong (1998) per timestep from **local** width, depth and
velocity: `K = 5.915 (W/H)^0.620 (U/u*)^1.428 H u*`, with `u* = sqrt(g)·U/Chezy` taken from the model's own
friction closure rather than a channel slope — SWORD's slopes at these delta reaches are unusable (zeros,
and values as absurd as 2.5 m/m). Since `U/u* = Chezy/sqrt(g)` identically, the velocity ratio is set by
friction alone.

**Why it had to change.** Under observation-based geometry the Savenije form collapses — `beta` explodes and
dispersion is clamped to zero after a few grid points, leaving transport purely advective:

| geometry | beta | dispersion non-zero over |
|---|---|---|
| shipped (B=1215/859, D=15 m) | 0.5 | all 137 pts (27.4 km) |
| colville (D=2.25 m) | 8.2 | 16 pts (3.2 km) |
| kuparuk (D=1.34 m) | 69.5 | 3 pts (0.6 km) |
| sagavanirktok (D=0.98 m) | 138.8 | 2 pts (0.4 km) |

After the change: **350–638 m²/s at all 136 points, for every site.**

**Seo & Cheong rather than Fischer (1979).** Fischer assumes W/H ≈ 10–100; Colville runs W/H = 477, where
Fischer returns 21661 m²/s. That exceeds `DISP_MAX` and would make `disp_sch`'s Crank–Nicolson right-hand
side oscillatory: the scheme builds `r[i] = 1 − (c1+c2)·DELTI/(8·DELXI²)`, so `r ≥ 0` needs
`K ≤ 4·DELXI²/DELTI` = **2133 m²/s** at the shipped `DELTI=75`, `DELXI=200`. `DISP_MAX` enforces this and
reports once per run if it binds. **It does bind** — the Colville run log carries the warning (an earlier
note here claimed it did not, which was wrong). It is rare and confined to peak-freshet timesteps at low
water: an offline reconstruction puts it at a couple of percent of cell-timesteps at most, and the
multichannel geometry (below) makes it *less* frequent, not more. The warning fires once per run, so the log
alone tells you it happened, not how often — check with the reconstruction in `tools/compare_multichannel.py`
if it matters.

**Root cause of both.** C-GEM's width law *and* its dispersion closure both derive from Savenije's
tide-dominated alluvial-estuary theory. This config has `AMPL = 0.2 m` and `distance = 1` — negligible tides,
essentially no saline intrusion. These are rivers with short delta flares, not funnel estuaries, so the
framework was being applied outside its regime; the width misfit and the collapsing dispersion were the same
problem surfacing twice.

**Verified:** all four sites run with finite output and no negative concentrations. **Not verified:**
the effect on `FCO2` — that needs a full multi-year run.

`EL = 27175` is unchanged and common to all four. It is a modelling choice (how far upstream to simulate)
rather than an observable, kept uniform for comparability.

### `B` is the TOTAL conveyance width; dispersion uses a separate per-thread width

`CGEM_MULTICHANNEL` — **on by default**, `=off` reverts to the legacy single-channel geometry
bit-identically. Full write-up and the sensitivity analysis: `docs/multichannel_test.md`.

**The defect this fixes.** One `B(x)` was doing two physically different jobs *and changing definition
between them*: `B_lb` is the SWORD delta distributary **sum** (total across parallel paths), `B_ub` is the
SWORD **per-channel** median (one of *n* braids) — with the total discharge carried through both. In a
braided reach that divides total Q by one braid's area, over-estimating velocity by ~`n_chan` and
under-estimating residence time and water-surface area by the same factor. Symptom: a spurious 3× velocity
jet just inside the Colville mouth (0.156 → 0.467 m/s over 5.4 km), created purely by the definition change.

**The fix.** `B` is the total conveyance/surface width everywhere — prismatic end = the observed
`B_UB_TOTAL` — and `variables.B_thread = B / n_chan` feeds *only* the Seo & Cheong closure, which is a
within-channel shear process. `n_chan` blends geometrically from `N_CHAN_LB` to the derived prismatic thread
count over `L_FLARE`. Two properties hold by construction and are asserted in the verification: upstream
`B_thread` reduces **exactly** to the site's `B_ub`, and at the mouth to `B_lb / N_CHAN_LB`.

**Read the braided total; do not reconstruct it.** SWORD's raw width is *already* summed across braids, so
sites declare the observed `B_UB_TOTAL` (`tools/extract_sword.py` emits it) and `config` derives the thread
count. `B_ub × median(n_chan_mod)` is **not** the braided total — `B_ub` is a median of per-node *ratios*,
and the median of a ratio is not the ratio of medians (Colville: 1052 observed vs 1269 reconstructed).

| river | `B_ub` per-thread | `B_UB_TOTAL` observed | IQR | threads | flare is |
|---|---|---|---|---|---|
| Colville | 423 | **1052** | 792–1354 | 2.49 | ~half definitional |
| Kuparuk | 58 | 60 | 42–78 | 1.03 | real |
| Sagavanirktok | 102 | 224.5 | 146–318 | 2.20 | partly definitional |
| Canning | 64 | 63 | 53–78 | 1.00 (clamped) | real |

**What it changes.** `ZZ = B·depth`, so scaling `B` scales the cross-section and leaves `DEPTH` untouched. The
largest effect is on the **water-surface area**, so the **basin-integrated** flux roughly doubles on the
braided rivers. Measured, year 2, `tools/compare_adoption.py`:

| river | surface area | per-area gC m⁻² yr⁻¹ | basin tC yr⁻¹ |
|---|---|---|---|
| Colville | 2.13× | 93.8 → 96.8 (1.07×) | 1321 → 2900 (**2.19×**) |
| Sagavanirktok | 1.79× | 91.2 → 76.5 (**0.84×**) | 339 → 531 (1.57×) |
| Kuparuk *(control)* | 1.02× | 176.4 → 175.7 (1.00×) | 429 → 437 (1.02×) |
| Canning *(control)* | 1.00× | 105.7 → 105.8 (1.00×) | 337 → 337 (1.00×) |

**The per-area flux is NOT invariant** — a tempting but wrong inference from "`DEPTH` is unchanged". Velocity
falls with the widened cross-section, and the flow-driven piston velocity goes as **√U**
(`fun_module._piston_velocity_loop`: `kflow = sqrt(|U|·D_O2/DEPTH)`), so gas exchange per unit area falls too;
the 1.5–1.8× longer transit then shifts the along-channel DIC/ALK balance. Decomposition: Sagavanirktok
0.88× = `vp` 0.961 × chemistry 0.917; Colville 1.07× = `vp` 0.928 × chemistry 1.153. Per-volume reaction
rates barely move (Sag `denit` 0.999×, `aer_deg` 0.998×) — it is the residence time that changes.

So the braided rivers' **pre-adoption per-area fluxes were biased high**: the too-narrow cross-section made
the water flow ~`n_chan` too fast and therefore exchange gas too vigorously. Sagavanirktok is where that shows
cleanly (−16%); on Colville a compensating chemistry change masked it. Kuparuk and Canning, being genuinely
single-thread, are unchanged to 3 significant figures — which is the control that confirms the mechanism is
the velocity change and nothing else.

**Do not read the salinity as constrained.** Mouth salinity improves — the model can now reach marine values
(27.1 PSU vs a baseline that capped at 11.6) where four earlier levers could not — but it does **not** close
the observed 8–32 PSU misfit (summer median 1.40 PSU), and it is highly sensitive to the one uncertain
number: over the SWORD IQR alone the summer median spans 0.02–3.23 PSU. The *capability* to intrude is
robust; the *typical* intrusion is not.

**Verify a run used the geometry you intended.** `config` prints it into every run log — `grep multichannel
runs/definitive/<site>/run.log`. This exists because **zsh does not word-split unquoted expansions**, so
`env $EXTRA python main.py` with a multi-assignment `$EXTRA` passes one malformed assignment, the flag
evaluates false, and the run produces output *byte-identical to the old baseline* — which reads as "the
change did nothing" rather than "the flag never arrived". That invalidated a whole sweep before it was
caught. Write launches with literal, separate `VAR=value` assignments.

**`BOUNDARIES` — the riverine carbonate boundary (`cub` DIC/ALK/pH) now comes from delta-proximal USGS
discrete samples for all four rivers, and all four outgas.** Built and reproduced by
`tools/usgs_carbonate_boundary.py` from the USGS Water Quality Portal. The method is the load-bearing part:
**DIC is back-solved to reproduce the observed PAIRED-sample `pCO2`** — `pCO2` is computed per grab from
co-located (same-activity) ALK+pH+T, the median is taken, and DIC is solved from the well-constrained median
ALK and that median `pCO2`. The earlier boundaries solved DIC from median ALK + median pH taken *separately*,
which **understated the flux-relevant `pCO2`** (pH and ALK covary; the median pH is biased by high-pH low-flow
grabs while the CO₂-rich grabs drive the air–sea flux) and left three rivers spuriously near equilibrium.

| river | source (open-water paired ALK+pH) | *n* | ALK | median `pCO2` | DIC |
|---|---|---|---|---|---|
| Colville | USGS 15880000 (nr Nuiqsut, delta head) | 74 | 1319 | 658 µatm | 1337 |
| Kuparuk | USGS 15896000 (nr Deadhorse, tidewater) | 48 | 1099 | 796 µatm | 1141 |
| Sagavanirktok | Sag mainstem, 3 stations | 16 | 1609 | 666 µatm | 1628 |
| Canning | eastern-ANWR regional pool, 20 stations | 28 | 1797 | 613 µatm | 1804 |

**Kuparuk's carbonate was upgraded off the Arctic-LTER *headwater* proxy** (EDI `knb-lter-arc.10303`, ~163 km
upstream, soft-water ALK 274) **onto its own near-tidewater gauge** (delta ALK 1099, ~4× higher from downstream
mineral weathering — exactly the "lower bound for the delta" caveat the LTER note had flagged). Kuparuk's
*nutrients* (NO3/NH4/PO4/TOC) still come from LTER; Canning has no discrete carbonate of its own, so its region
is pooled and only its carbonate is sourced (all other species stay placeholder — consistent with it being the
weakest site). **Result (year-2 open-water FCO₂, gC m⁻² yr⁻¹): every river flips to a realistic net CO₂
source, 100 % of open-water cells outgassing** — Colville −3.9→**+94.6**, Kuparuk −9.3→**+177.5**,
Sagavanirktok −51.7→**+91.8**, Canning +425→**+106.6** (its placeholder had been over-outgassing). The prior
near-equilibrium/uptake was a boundary-chemistry artifact (DIC ≈ ALK), not a model deficiency; this and the
`_pbar_rho` pressure fix (Known defects) together make the carbon flux physical.

**Those four numbers predate the multi-channel adoption** (they are the single-channel geometry, retained here
because they are what isolates the carbonate-boundary effect). On the current geometry the per-area fluxes are
**Colville +96.8, Kuparuk +175.7, Sagavanirktok +76.5, Canning +105.8** — the two single-thread rivers
unchanged, the two braided ones shifted by the velocity/gas-transfer effect described under *`B` is the TOTAL
conveyance width*. Basin-integrated totals are **2900 / 437 / 531 / 337 tC yr⁻¹**.

**Canning is the weakest-constrained site**: observed geometry, but *reconstructed* discharge, *borrowed*
temperature, and only 28 surveys behind its depth.

## Forcings

`main.py` used to read from hardcoded absolute paths on the original author's machine
(`/Users/rsavelli/Documents/CMS_LOAC/...`). It now resolves `../forcing/` relative to its own
file location, so it works from any output directory, and `forcing()` raises a clear error on a missing file.
The `pCO2_barrow_2022.csv` / `pCO2_Barrow_2022.csv` case bug is fixed.

Note: `code/` used to carry stale, wrong-length copies of `watertemp.csv`, `windspeed.csv`,
`solarradiation.csv` and a dead `daily_average_weather.csv` that would raise `ValueError: fp and xp are not of
the same length` in `numpy.interp` if a path ever resolved to them. **These have been removed** — the only
forcing CSVs now live in `forcing/`, all 365 values. If that error reappears, a wrong-length file
has been reintroduced into the code dir; delete it (the real forcings are in `forcing/`).

**Meteorological forcing is shared across all four sites** — one regional record each for wind, solar,
water temperature and pCO2 (Barrow). Wind/solar/air come from NDBC PRDA2, Prudhoe Bay (raw file
`prda2h2022.txt`). This is defensible along this coastline and isolates discharge and geometry as the main
differences between sites. It does mean per-river ice timing is identical, which is a real simplification.

### Tides — per-river harmonic reconstruction

The mouth sea-surface elevation is a sum of harmonic constituents from the nearest NOAA CO-OPS station
(`tools/build_tides.py` → `forcing/tidal_constituents.json`; loaded in `config`, summed in
`fun_module.Tide`): `eta(t) = Σ Aᵢ cos(speedᵢ·(t/3600) − Gᵢ)`. Stations: Colville → Cape Halkett, Kuparuk →
Prudhoe Bay, Sagavanirktok → Cross Island, Canning → Point Thomson. This replaced the shipped single sinusoid
(`AMPL`/`pfun`, now a fallback only). **The coast is genuinely microtidal** — M2 amplitude ~0.06–0.07 m,
short-period range ~0.3–0.4 m — so the change is one of *structure* (spring-neap, diurnal inequality, per-river
phase), not amplitude; the old 0.2 m was a fair amplitude guess. **Seasonal SA/SSA constituents are excluded**
(they are the annual sea-level cycle, not tides, and would conflate with the discharge-driven seasonal signal).
The phase reference is the station's GMT equilibrium argument, so absolute calendar timing is not pinned in
this idealised climatology run — amplitudes, frequencies and relative phases are real.

### Wind-driven storm surge — the real saltwater-intrusion driver (all four rivers)

On this microtidal coast the ~0.3 m astronomical tide is **not** the dominant saltwater-intrusion mechanism —
wind-driven storm surge is (Beaufort fall storms reach 1–3 m; 2022 was moderate, peak daily-mean +0.66 m).
`fun_module.Tide` therefore adds a `SURGE` offset to the harmonic elevation, refreshed each timestep in
`main.py` from a per-site daily series (`config.SURGE_FILE`; 0 when unset). `tools/build_surge.py` builds it as
the **observed daily-mean water level** (the tide averages out) from NOAA CO-OPS. **Only Prudhoe Bay (9497645)
has continuous 2022 water level on this coast**, so **all four rivers** use `surge_prudhoe_2022_m.csv` — it is
Kuparuk's own tidal station and the **regional proxy** for the other three (a single storm surge affects the
whole Beaufort coast similarly, like the shared meteorology).

**Result of the surge (Kuparuk example):** it is the correct mechanism and moves salinity the right way — the clean
open-water surge on day 279 (surge +0.56 m, low flow) intrudes to 1.5 PSU, and the summer-mouth maximum roughly
doubles (2.7 → 5.7 PSU) — but it does **not** reproduce the observed 8–32 PSU grabs, because 2022's *largest*
surges (+0.66 m) all fall *after fall freeze-up* (ice-covered mouth → no intrusion) and the open-water surges
coincide with residual summer flow. So the salinity gap is **event-timing / year-specific**, not a model
deficiency: width, tides and even the real surge are all insufficient for 2022, and the grabs (all years, and
possibly sampled in the seaward lagoon) reflect bigger open-water surge events. The surge is kept because it is
real, station-matched, physically correct data; the effect on the 2022 climatology is just small.

### Water temperature is a TRANSPORTED field, not a scalar

`T` is species 13 in `variables.v`, with `env=1`, so the existing TVD advection and
dispersion carry it exactly as they carry salinity. It replaces the single `water_temp`
scalar that was previously applied uniformly to all 136 grid points. Every
temperature-dependent rate — `Fhet`, `Fnit`, `O2sat`, `K0/K1/K2/KB`, `p_bar` — is now
evaluated at the **local** value `Ti = c_T[i]`.

Boundary values are seasonal, so unlike salinity's fixed `clb`/`cub` they are refreshed
every timestep in `main.py` (`openbound` re-reads them per call):

- **downstream = sea temperature** — `watertemp.csv`, the PRDA2 buoy record. The file
  that was wrong as a *river* forcing is correct as a *marine boundary condition*.
- **upstream = river temperature** — the modelled `WATERTEMP_FILE`.

The values in `sites/_baseline.py` `BOUNDARIES["T"]` seed the initial profile only.

This is step (a) of `docs/ice_model_plan.md` §4b. Step (b), the **surface heat budget**,
is in (`heat_module.py`), and step (c), the **prognostic ice feedback**, is now in
(`ice_module.py`) — see both below.

### Surface heat budget (`heat_module.py`)

Applied to `v['T']` each timestep in `main.py`, **after** transport, so the atmosphere
warms/cools the field that advection has just moved:

    rho_w cp_w H dT/dt = Q_sw + Q_lw + Q_sens + Q_lat

This is what makes interior temperature physical rather than a linear blend of the two
boundaries — verified to lift mid-channel summer temperature above both end-members
under net input of ~200–300 W/m². Toggle with `config.HEAT_BUDGET`.

Humidity is now **observed**: `main.py` reads `relhum_2022_frac.csv` (Deadhorse Airport,
colocated with PRDA2; built by `tools/build_humidity.py`) and passes it per timestep as
`heat_budget`'s `rel_hum` argument. `config.REL_HUMIDITY = 0.85` survives only as that
argument's default — the comment block above it in `config.py` still describes the old
assumed-constant regime and is stale.

Three things it does NOT do well, in decreasing severity:

1. **Heat lost below freezing becomes ice** (with `ICE_MODEL` on). The per-step freeze
   energy is booked to `heat_module.ice_energy_deficit` and `ice_module` converts it to
   ice thickness via `rho_ice*L_fusion`. With `ICE_MODEL` off it is discarded and T is
   simply clamped at 0 °C (the old imitation).
2. **No cloud cover** — longwave uses clear-sky emissivity (Brutsaert), biasing cooling
   high.
3. **Daily-mean shortwave** — no diurnal cycle; fine midsummer, poor in shoulder seasons.

With `ICE_MODEL` on the budget runs **year-round** and **skips ice-covered cells** (the
slab insulates them; `ice_module` holds the under-ice water at freezing). With `ICE_MODEL`
off it falls back to being skipped whenever the crude `previousdays` gate is closed.

**No bit-identity here.** Unlike the optimization work, this changes results by design,
so the MD5 harness does not apply — validation has to be physical (see §4b).

### Prognostic ice model (`ice_module.py`)

`ice_step()` runs each timestep in `main.py`, **after** `heat_budget` (whose per-step
freeze energy it consumes) and **before** `biogeo` (which reads `ice_frac`). One prognostic
per cell, `variables.ice_thickness` [m], with `ice_frac` diagnosed from it. Enabled by
default; toggle with `CGEM_ICE=off` (`config.ICE_MODEL`), which reverts every coupling to
the legacy `previousdays` gate. Four mechanisms:

- **Freeze-up** — the freeze-clamp energy `heat_budget` books becomes new ice.
- **Conductive (Stefan) growth** — under an existing cover the base is at freezing and the
  slab conducts heat out to a surface at `min(T_air, 0)`; this builds the winter thickness
  the frazil term alone cannot.
- **Bottom-fast grounding** — the local depth limits how much *new* ice can form (once the
  slab reaches the bed there is no water left to freeze), but a falling water level never
  *removes* ice: grounded ice rests on the bed at unchanged thickness, so a grounded slab
  can be thicker than the instantaneous water depth. See the *Bottom-fast ice was a
  one-way ratchet* note below — the earlier unconditional cap was a real defect.
- **Surface melt** — warm air (sensible) + absorbed shortwave (ice albedo) thin the slab
  in spring.
- **Hydraulic break-up** — when discharge exceeds `BREAKUP_Q_FACTOR × annual-mean
  discharge` (computed once in `main.py`), the freshet surge clears the cover
  mechanically. Referenced to the mean, not winter baseflow (the ice-affected winter
  percentile is ~0); freshets run 6–30× the mean so it fires only then, and no
  air-temperature gate is needed. This is the North-Slope-specific mechanism; a purely
  thermal model breaks up weeks too late (see `docs/ice_model_plan.md` §0).

**Couplings** (all gated on `ICE_MODEL`, all no-ops when `ice_frac == 0`):
- `heat_module` skips ice-covered cells and holds the water at `T_FREEZE` (insulation).
- `biogeo_module` attenuates under-ice PAR through the slab (`k_ice_PAR`, Beer–Lambert) and
  scales O2/CO2 gas exchange by `(1 − ice_frac)` — sealed cells neither outgas nor ventilate.
- `transport_module`, `biogeo_module`, `sed_module` **conserve state year-round** instead of
  zeroing it under ice (see the retired-gate note below); erosion is scaled by `(1 − ice_frac)`.

Ice thickness and fraction are written each save step (`file_module.icewrite` →
`ice_thickness`/`ice_frac` in the NetCDF/.dat output). **Deferred:** the ice draft does not
yet reduce the hydrodynamic cross-section (`docs/ice_model_plan.md` Tier 3).

#### Bottom-fast ice was a one-way ratchet — FIXED

`ice_module` used to truncate thickness to the local depth unconditionally, every timestep:
`if h > DEPTH[i]: h = DEPTH[i]`. But `DEPTH` is the **live** water depth, and on this coast
it swings **~0.8 m per day** on tide + storm surge — comparable to the *entire* water column
on the three shallow rivers. So the cap chopped the cover roughly twice a day, and the
clipped ice was **destroyed**: when the water rose again it could only return by slow
conduction. The cover therefore ratcheted *thinner* through midwinter instead of growing.

**98% of every midwinter thinning event sat exactly at the cap** — the attribution is
unambiguous, and the behaviour was identical in the pre-multichannel runs (7.70% vs 7.54% of
steps thinning), so it long predated that work.

The fix: bottom-fast ice **grounds**, it does not melt. Depth limits how much *new* ice can
form; existing thickness is never reduced by a falling water level
(`h = DEPTH[i] if DEPTH[i] > h_prev else h_prev`). Year-2 deep winter, days 10–110
(`tools/compare_icefix.py`, both sides on the same multi-channel geometry):

| river | depth | thinning steps | reversals | domain-mean ice | peak ice | ice destroyed |
|---|---|---|---|---|---|---|
| Colville | 2.11 m | 1.13% → **0** | 8/100 → **0** | 1.29→1.64 ⇒ 1.29→**2.07** | 1.96 → 2.12 | 0.54 → **0** |
| Kuparuk | 1.20 m | 3.41% → **0** | 25/100 → **0** | 0.89→0.80 ⇒ 1.17→**1.39** | 1.20 → 1.40 | 1.53 → **0** |
| Canning | 0.97 m | 4.95% → **0** | 29/100 → **0** | 0.74→0.60 ⇒ 1.02→**1.17** | 1.02 → 1.19 | 1.89 → **0** |
| Sagavanirktok | 0.85 m | 7.54% → **0** | 35/100 → **0** | 0.69→0.53 ⇒ 0.94→**1.06** | 0.94 → 1.07 | 2.08 → **0** |

Ice now **grows** through winter on all four; before, three of the four were *losing* ice
through midwinter. The severity scaled with shallowness, as the mechanism predicts.

**The carbon results are completely unaffected — verified, not assumed.** Per-area and
basin-integrated `FCO2` are identical to four significant figures on every river (ratio
1.000), because gas exchange is gated on the *binary* `ice_frac`, which does not move: the
slab never approached `ICE_FORM_THRESH`, so `ice_frac` still flips exactly twice per cell per
year and the open-water fraction is unchanged to 0.01%. (Canning has one cell at 16.2 km that
freezes, briefly reopens in a late-October warm spell and refreezes — 274 flips, physical,
and identical before and after.) What *was* wrong is the slab thickness itself, which feeds
under-ice PAR through Beer–Lambert. `runs/multichannel_preicefix/` holds the pre-fix runs.

### Provenance of the temperature forcing — do not use `watertemp.csv` as a river input

`watertemp.csv` is **sea-water** temperature from the PRDA2 buoy (r = 0.9922 against `prda2h2022.txt`,
min −2.00 °C — the freezing point of *seawater*, which fresh river water cannot reach). It is retained only
as provenance and is **no longer read**. Using it as a river forcing did two things:

- biased every temperature-dependent rate cold by 6–8 °C;
- drove the ice gate off *coastal sea ice*, which clears weeks after river breakup — the gate opened day 173
  while Colville and Kuparuk peak day 153, discarding 45% / 40% of annual water volume **including the entire
  freshet**.

The model now reads `river_watertemp_2022_degC.csv`, built by `tools/build_river_temp.py` from

    T_river = max(0, 1.432 * (T_air_10day + 4.0))

fitted to USGS parameter 00010. Key properties, all verified:

| | old (sea) | new (river) |
|---|---|---|
| ice gate open | day 173–271 (99 d) | day 140–303 (164 d) |
| Colville annual volume discarded | 45.1% | 0.2% |
| Kuparuk annual volume discarded | 39.9% | 0.2% |
| mean T over open water | 1.20 °C | 8.02 °C |
| `Fhet` respiration | 0.691 | 1.378 (**1.99×**) |
| `Fnit` nitrification | 0.556 | 1.666 (**3.00×**) |
| `O2sat` | 441.8 | 370.0 mmol/m³ (−16.3%) |
| sub-zero river days | 267 | 0 |

**Caveats.** The relation is *regional*, fitted on Kuparuk and applied to all four rivers — Kuparuk is the
only gauge whose temperature air temperature explains (R² = 0.706, RMSE 2.45 °C). The Sagavanirktok is
decoupled (constrained fits to its own data give **negative** R²; it is mountain-fed, groundwater- and
aufeis-influenced), and Colville and Canning have no 00010 record at all. Validation against observations:
Kuparuk RMSE 2.45 °C / bias +0.05; Sagavanirktok RMSE 4.05 °C / bias **+2.61**. So this trades a −6 to −8 °C
cold bias for a smaller warm one — a clear improvement in seasonal timing, but not a good thermal model.
The zero-crossing is *constrained* (slope-only fit through T0 = −4 °C) because an unconstrained fit
extrapolated to 297 open-water days on the Sagavanirktok. See `docs/ice_model_plan.md` §0.

`airtemp_2022_degC.csv` is also produced and is currently unused — it is the input a real ice model needs.

### Discharge

Rebuild all four with `python tools/fetch_discharge.py` (USGS NWIS, parameterCd 00060, statCd 00003,
cfs→m³/s via `0.3048**3`). That script reproduces the original Colville file **byte-for-byte**, which is the
check that one recipe is being applied to all four sites.

| Site | Gauge | 2022 mean | Caveat |
|---|---|---|---|
| Colville | USGS 15875000 (at Umiat) | 238.8 m³/s | Gauged far upstream of the delta — **understates** mouth discharge |
| Kuparuk | USGS 15896000 (nr Deadhorse) | 63.8 m³/s | Near tidewater; best-constrained of the four |
| Sagavanirktok | USGS 15908000 (nr Pump Sta 3) | 47.6 m³/s | ~130 km inland, under half the basin — **understates** badly |
| Canning | *none* | 45.4 m³/s | **RECONSTRUCTED** — see below |

**Canning has no 2022 observations.** USGS 15955000 ran only 2008-06-23 to 2012-09-30. The forcing is
Hulahula River (USGS 15980000, nearest gauge active in 2022) scaled by 2.971 — the ratio of means over their
731-day common record, log-space r = 0.87. Ratio-of-means rather than median-of-daily-ratios (~3.7) because
it preserves annual water volume, which the transport and carbonate budgets are sensitive to. `config.py`
warns at import. Do not present Canning at the same confidence as the other three.

Roughly 250 of 365 days each year carry the USGS `A:e` (approved, estimated) flag — the ice-affected winter
months are modelled, not measured. The ice gate suppresses biogeochemistry across most of that period anyway.
In 2022 water temperature is below zero **269 days**, so only ~96 days do active carbon work.

## Architecture

**Global mutable state.** `variables.py` allocates every array at import time, sized from `config.M`. Modules
`from variables import ...` and mutate those lists **in place** — there is no state object, and nothing is
passed by return value. A module that rebinds rather than mutates (`x = [...]` instead of `x[i] = ...`)
silently detaches from the shared state. `config.py` is the single source of all parameters; `M` is derived
from `EL / DELXI` and forced even, with `M1..M3` as the offset variants the schemes index by.

**Per-timestep call graph**, driven by the loop in `main.py`:

```
main() ── exfread() × 9+ ──> forcings interpolated to time t
      ├─ recompute dispersion[] from Qr (Seo & Cheong)
      ├─ hyd(t, Qr) ──> new_bc → [coeff_a → tridag → conv → update]* → new_uh
      ├─ transport(t, previousdays) ──> per species: openbound → tvd → disp_sch
      ├─ heat_budget(t, ...) ──> surface heat flux onto v['T'] (skips ice cells)
      ├─ ice_step(t, ...) ──> grow/melt ice_thickness, diagnose ice_frac
      └─ if t > WARMUP:
           ├─ biogeo(t, Uw_sal, Uw_tid, v['T']['c'], pCO2, I0, previousdays)
           └─ sed(t, previousdays)
```

Nine fixed `exfread` calls over eight distinct files (discharge, wind ×2, solar, air temp, relative humidity,
pCO2, river temp, sea temp), plus one per `BOUNDARY_FORCING` entry and one for `SURGE_FILE` when set.
`windspeed.csv` is read twice, into `Uw_sal` and `Uw_tid`, so the saline/tidal wind split in
`fun_module.piston_velocity` is fed identical values.

`exfread` used to re-open and `genfromtxt`-parse the whole CSV on **every** call — one parse per call ×
~840k timesteps. It now caches parsed series in `file_module._FORCING_CACHE`, keyed by filename; the
interpolation still runs per timestep, so results are bit-for-bit unchanged. The rolling ice-sum is cached
alongside. If you add a forcing that changes on disk mid-run, this cache is what will bite you.

`hyd` iterates the tridiagonal solve to convergence (`while rsum != 2.0` — an exact float comparison on the
sum of two `conv()` flags, which will spin forever if either never returns exactly 1.0).

**Staggered grid.** Odd indices hold concentrations and cross-sections; even indices hold velocities and
fluxes. This is why `tvd` runs `for j in range(1, M2+1, 2)` writing `co[j+1]`, then a second pass over odd `j`.
Any new scheme must respect the same parity or it will silently read velocity slots as concentrations.

**Species registry.** `variables.v` is a dict of 16 species (`DIA, dSi, NO3, NH4, PO4, O2, TOC, RDOC, CH4,
N2O, S, SPM, DIC, ALK, pH, T`), each with concentration `c`, boundary values `clb`/`cub`, running `avg`, and
flux accumulators. `RDOC/CH4/N2O` are the Arctic-extension tracers (see below); `T` is the transported
temperature field. The `env` flag gates transport: `env=1` species are advected and dispersed; `pH` is set to
`env=0` in `init_module` because it is diagnosed from S/DIC/ALK rather than transported. Adding a transported
tracer is just three edits — the `names` list, the `BOUNDARIES` table, and (if it reacts) the `biogeo` loop.
Boundary concentrations come from the active site's `BOUNDARIES` table (`sites/<name>.py`).

**Arctic biogeochemistry extension** (`docs/arctic_biogeochemistry.md`). Four process groups layered on the
shipped reaction network, **additive** and mostly gated to open water like the rest: (1) a refractory +
chromophoric DOC pool `RDOC` with slow oxidation AND **CDOM photomineralisation** (`photo`, light-driven →
DIC; `TOC` stays the labile pool); (2) **CH4 and N2O** cycling with methanotrophy, N2O yields from
nitrification/denitrification, and air-water gas exchange (`fun_module._sc_ch4/_sc_n2o/_ch4_eq/_n2o_eq`,
Wanninkhof 2014 / Wiesenburg-Guinasso 1979 / Weiss-Price 1980); (3) a **benthic** DIC/alkalinity/CH4 efflux
via a first-order SOD closure (runs under ice, so DIC/CH4 build up and vent at ice-out); (4) **distributed
lateral loading** (`lateral_module.py`) injecting tundra/thermokarst inputs along the channel. The extension
reactions are **opt-in** via `config.ARCTIC_BGC` (default `False`; the idealized fixture sets it `True`): with
it off, `RDOC/CH4/N2O` advect as inert passive tracers and the DIC/O2/ALK/NO3 updates add exactly `0.0`, so the
four real rivers' existing fields are **bit-identical** to the shipped network until their RDOC/CH4/N2O boundary
chemistry is constrained (they carry placeholder values today; lateral is also off, `LATERAL_INFLOW` unset).
All parameters + citations are in `config.py` under "ARCTIC BIOGEOCHEMISTRY EXTENSION". Verified end-to-end by
`tools/verify_idealized.py --full`; diagnostics on page 4 of `docs/idealized_verification.pdf`.

**Ice gating.** With `ICE_MODEL` on (default) the prognostic `ice_module` is the single ice authority
(see "Prognostic ice model" above): `transport`/`biogeo`/`sed` **conserve** state year-round and couple to
`ice_frac`, and `previousdays` no longer gates them. `previousdays` survives only as the **legacy** gate used
when `CGEM_ICE=off`. In that legacy path: `file_module.exfread` returns two values, the interpolated forcing
and `previousdays` — a rolling `nbday_ice`-day cumulative sum of the same series; only the water-temperature
call's second return is meaningful, so `main.py` relies on the **`P_WATERTEMP` (river temperature) call being
the last one that binds `previousdays`**. The calls after it (`P_SEATEMP`, the boundary-forcing loop,
`P_SURGE`) discard the second return into `_` for exactly this reason — note that `P_SEATEMP` *is*
`watertemp.csv`, so "last file read" and "last `previousdays`" are no longer the same call. Reordering these
`exfread` calls, or binding `previousdays` in one of the trailing ones, will silently break the legacy ice
logic. When `previousdays <= 0` (legacy only), `transport`, `biogeo`, and `sed` zero their state rather than
integrating.

**Carbonate chemistry** lives in `fun_module`: `_k0/_k1/_k2/_kb` (Millero 1995, pressure-corrected via
`density.dens`), then the **mol/kg** H+ solve `_h_solve_kg` (Follows et al. 2006), then CO2* and `FCO2` in
`biogeo_module`. **Solved consistently in mol/kg (ported from upstream C-GEM v2):** `biogeo` converts the
mmol/m³ state (DIC/ALK) to mol/kg via the local in-situ density before speciation, so it shares the mol/kg
basis of the constants, then converts CO2* and the Henry term back to mmol/m³ for the flux. Sign convention:
**`FCO2 > 0` is outgassing** (sea→air) — flipped from the original `< 0`; `biogeo` drops the leading minus on
`RCO2` and carries `-FCO2` in the DIC budget, so DIC evolution is numerically unchanged.

**The mol/kg fix is a first-order correction to `FCO2`, not cosmetic.** The pre-v2 code (which we inherited)
subtracted the atmospheric-saturation term `K0·pCO2` in **mol/kg** from `co2s` in **mmol/m³** — ~10⁶× too
small — so it was effectively negligible and `FCO2 ≈ vp·co2s`: *always outgassing, proportional to dissolved
CO₂, ignoring the air–sea gradient*. With the units reconciled the flux became the real disequilibrium
residual — so `FCO2` is now **genuinely sensitive to the boundary DIC/ALK**, and its sign is set by whether the
river carries a respiratory CO₂ excess (DIC > carbonate alkalinity). With the *delta-proximal* carbonate
boundaries now sourced for all four rivers (see `BOUNDARIES` above), the observed water `pCO2` is 580–800 µatm
(supersaturated) and **all four rivers outgas** (+77 to +176 gC m⁻² yr⁻¹ on the current multi-channel
geometry; +92 to +178 on the single-channel geometry those boundaries were commissioned against, 100 % of
open-water cells) — the
physically expected behaviour for Arctic permafrost rivers. *(Historical note: with the earlier boundaries
these three read near-equilibrium/uptake and Canning over-outgassed; that was a boundary-chemistry artifact
of DIC ≈ ALK, not the flux code. Two things had to be right for the sign to be trustworthy — the mol/kg
saturation term here, and the `_pbar_rho` pressure unit fix in Known defects.)*

**Output.** A single self-describing **`output.nc`** per run directory by default
(`config.OUTPUT_FORMAT`/`CGEM_OUTPUT = nc`), written when `(t / (TS*DELTI)) % 1 == 0`. All the analysis and
figure tools read this. The legacy path writes one tab-separated **`.dat`** per variable — appended one row
per save, time in column 0 followed by `M` values, so each file is a time × distance matrix — and is still
available via `CGEM_OUTPUT=dat` (or `both`), which is what the bit-identity harness uses. Either way the write
calls are the same: state variables through `file_module.transwrite` (filename from `v[name]["name"]`),
process rates through `Rates` (filename hardcoded at the `biogeo` call site), ice through `icewrite`, and
`close_output()` at the end — which the NetCDF path **needs** in order to flush.

## Known defects in the vendored code

Inherited from upstream, not introduced locally. Left as-is so far.

- **`fun_module.pH` — FIXED (three defects).** (1) The call site passed salinity as the H+ guess and the
  guess as salinity (args swapped); now in the right order. (2) The loop reset `hg = Hplus` on every pass, so
  it never iterated (and wasted 49× work under `@njit`); `hg` is now initialised once and fed back. (3) The
  result-changing one: after ice-out the gate leaves DIC/ALK **tiny-positive** (e.g. 1e-7, 1e-38) while
  transport re-establishes them, and the solve divided by a ~0 carbonate alkalinity to return a huge H+ —
  producing non-physical pH (values of −25, −inf littered the fields). Now guarded three ways: a 1.0 mmol/m³
  floor on DIC/ALK, the existing `cag ≤ 0` / `dummy < 0` checks, and a hard physical-range bound — any H+
  implying pH outside [2, 12] returns a neutral pH 7 (a bounded transient over a few ice-out timesteps).
  Verified: run-wide pH now sits in ~5–9.4 with zero cells outside [2, 12] (was thousands). Open-water `FCO2`
  is essentially unchanged — the garbage was confined to the ice-recovery transient, which `FCO2` gates to
  zero anyway.
- **Carbonate unit inconsistency — FIXED (mol/kg solve, ported from C-GEM v2).** The old solve mixed mmol/m³
  state with mol/kg constants (the borate term, and — more consequentially — the `K0·pCO2` saturation term in
  `FCO2`, ~10⁶× too small). It is now solved consistently in mol/kg with a density conversion; see "Carbonate
  chemistry" above. This was a first-order correction to `FCO2` — it made the flux the real air–sea
  disequilibrium residual rather than always-outgassing, so the sign then followed the boundary DIC/ALK (now
  delta-sourced; all four rivers outgas). Not the small change the earlier note anticipated.
- **`fun_module._pbar_rho` carbonate pressure — FIXED (Pa→bar unit error).** The mid-column hydrostatic
  pressure was returned as `(depth/2·rho·G)·0.1` — Pascals scaled by `0.1` — instead of the Pa→bar factor
  `1e-5`. At these rivers' ~1.3 m depth that gave `pb ≈ 657 bar` instead of `≈ 0.066 bar` (10⁴× too large),
  inflating the Millero-1995 pressure corrections on `K1/K2/KB` ~1.9× and biasing the **diagnosed pH ~0.3–0.5
  units low** (median ~6.8 vs the physical ~7.7 in freshwater cells). **`FCO2` is essentially unaffected
  (<2%):** the flux uses `co2s` from the *same* pressure-corrected constants and the H⁺ that the solve shifts
  to conserve alkalinity, so the `K1/K2` inflation cancels out of the flux — the bug corrupts only the pH
  **output** field, not the air–sea flux. Verified by driving the flux with correctly-solved DIC/ALK (98% sign
  agreement with the stored `FCO2`) and by watching `c_pH` drift under the buggy pressure in an instrumented
  run. Because `FCO2` moves slightly, output is **not** bit-identical to pre-fix runs but is physically
  equivalent for the flux. Fixed at `code/fun_module.py` (commit `d825cec`).
- **`init_module` line 33** hardcodes `50` in the Chezy ramp (`(i - 50) / (M - 50)`) while the surrounding
  conditional keys off `distance`, which is `1` in the current config. For `1 <= i < 50` the numerator is
  negative and Chezy exceeds `Chezy_lb`.
- **`transport_module`** does `names = list(v.keys())` then `for names in names:`, rebinding the list to each
  key mid-loop. Works only because the iterator was already constructed.
- **`schemes_module.tvd`** sets `cold = co`, aliasing rather than copying, so the "old" values it reads back
  are partly updated. It also shadows the `cold` imported from `variables`.
- **`fun_module.I0`** is dead. `biogeo` takes `I0` as a forcing argument from `main`, and never imports the
  function of the same name — so the synthetic diurnal light curve it computes is unused.
- **`transport_module`** rebinds `v[name]["c"] = [0.0] * (M + 1)` when ice-gated rather than zeroing in
  place. Access through the `v` dict stays correct, but any module holding a direct reference to the old
  list would detach — the general hazard called out under *Global mutable state*.

## Configuration provenance

`config.py` is North Slope; `config_OLD.py` is the prior site (160 km estuary, 3.5 m tides) and is dead code
kept for reference — nothing imports it. Diffing the two is the fastest way to see which parameters are
site-specific: geometry, `distance`, tidal amplitude, sediment erosion constants (retuned to Clark et al.
2020/2022), and the phytoplankton parameters (retuned to Le Fouest et al. 2013, which also introduced the
`Chla2CMIN`/`KE` photoacclimation terms absent from the old version).

The forcing CSVs are single-column daily series for 2022, no header, and `exfread` interpolates them onto
`linspace(0, 365*86400, 365)` — so they must be exactly 365 values long. With `repeatYear = 1` the model wraps
`t` back by one year, making a multi-year run a repeated climatology rather than a continuous series.

`code/` also carries a `.idea/` PyCharm directory from upstream, still named for the earlier
`code_python_FCO2` project.

## Where the other documentation lives

This file is the per-decision developer reference — architecture, provenance, and why each non-obvious choice
was made. It is not the only doc, and the others are not redundant with it:

| doc | what it is for |
|---|---|
| `README.md` (root) | the outward-facing project page — what NS-RAD is, quick start, layout, citation |
| `docs/README.md` | index of `docs/`, and **which files are hand-written vs generated** — read before editing anything in there |
| `docs/FEATURES.md` | catalog of everything NS-RAD adds on top of upstream C-GEM, as *was → now* |
| `docs/model_description.md` | the technical model description: governing equations, parameterizations, parameter tables |
| `docs/performance.md`, `ice_model_plan.md`, `arctic_biogeochemistry.md`, `idealized_verification.md` | the deep write-ups this file summarizes |

**Only the Markdown files in `docs/` are hand-written.** The PDFs, `.docx`, movies, `.json` and `.png` there
are generated artifacts — regenerate with `tools/build_all.sh`, never hand-edit.

**Documentation caveat.** `readme.txt` refers to `./code_python_FCO2` and `./code_python`; the actual
directories are `code` and `forcing`. It is the upstream C-GEM readme and its paths do not describe this tree.
It is still the authoritative citation list — Volta et al. 2014/2016, Regnier et al. 2002/2013, Follows et al.
2006, Laruelle et al. 2017.

Both previously-noted drifts in the sibling docs are now fixed: `docs/FEATURES.md` no longer advertises a
`CGEM_CARB_UNITS` variable (it does not exist in `config.py`), and the **v17c** SWORD references throughout
the code, tools and docs are corrected to **v17b** — the extraction reads `na_sword_v17b.nc`, so v17b was
always the right label.
