# NS-RAD run catalog

NS-RAD has no restart/checkpoint capability (every run simulates from day zero), so
each distinct question this project asks ends up as its **own run directory** under
`runs/`, built from its own site variant and read by its own set of figure/report
tools. This file catalogs every run type that currently exists: what it is, what it's
for, exactly how to (re)build it, and what reads it downstream. All of `runs/` is
git-ignored (regenerable, multi-GB) — this file is the map back to how to regenerate
any of it.

Three rivers are runnable: **Colville, Kuparuk, Sagavanirktok** (`CGEM_SITE=colville`
etc.). Canning is excluded everywhere below (`EL=0` placeholder, non-runnable) unless
named explicitly. See `CLAUDE.md` for the full architecture/provenance reference this
file assumes as background — it is not repeated here.

## Quick reference

| run directory | site variant(s) | built by | what it's for |
|---|---|---|---|
| `runs/definitive/<site>/` | `<site>` | `tools/run_sites.sh` | the operational baseline run every "current state of the model" report reads |
| `runs/regression_bnd/<site>/` | `<site>` + `CGEM_WATERTEMP_FILE` override | `tools/run_regression_bnd.sh` | independent upstream-T-boundary check for `make_validation_pdf.py` |
| `runs/idealized/` | `idealized` | `tools/verify_idealized.py` | synthetic verification fixture, not a real river |
| `runs/interannual/<site>/` | `<site>_interannual` | `tools/run_interannual.sh` | 44-year (1980-2023), discharge+DOC only multi-year |
| `runs/interannual_full/<site>/` | `<site>_interannual_full` | `tools/run_interannual_full.sh` | 19-year (2005-2023), **every** forcing category genuinely multi-year |
| `runs/definitive_full/<site>/` | `<site>_interannual_full` (6-min, 730-day slice) | `tools/run_definitive_full.sh` | real 2005-2006 at 6-min cadence, for Fig3's tidal-aliasing calibration only |
| `runs/osse_nature/<site>/` | `<site>_interannual_full` (one continuous run spanning 2005-2009) | `tools/run_osse_nature.sh` | "truth" run for the original (analog-year) OSSE |
| `runs/osse_green/pert_<control>/<site>/` | `_green_<site>_<control>` (11 controls) | `tools/build_green_perturbations.py` | boundary-only Green's-function perturbations, analog-year framework |
| `runs/forte_campaign/<site>/` | `<site>_forte_campaign` | `tools/run_forte_campaign.sh` | climatological "typical campaign year" run |
| `runs/forte_green/pert_<control>/<site>/` + `lateral_base/<site>/` | `_fortegreen_<site>_<control>` (25 controls + 3 baselines) | `tools/build_green_perturbations_forte.py` | boundary + reach-local Green's-function perturbations, FORTE framework |

## 1. Operational single-year runs

### `runs/definitive/<site>/`
The baseline run for "what does the model currently say about this river" — single
repeating 2022-pattern year, 6-minute output cadence. Every diagnostics/validation
report and the main combined PDF read this.
```bash
tools/run_sites.sh                    # all 3 rivers, parallel
tools/run_sites.sh kuparuk            # one river
CGEM_MAXT_DAYS=2 CGEM_WARMUP_DAYS=1 tools/run_sites.sh   # smoke test
```

### `runs/regression_bnd/<site>/`
Kuparuk and Sagavanirktok rerun with the pure air-temperature regression
(`river_watertemp_2022_degC.csv`) as the upstream T boundary instead of their own
observed-blended series, so `make_validation_pdf.py`'s temperature skill check isn't
tautological (scoring the model against the very data that fed its own boundary).
Colville/Canning already use the regression, so they aren't rebuilt.
```bash
tools/run_regression_bnd.sh
```

### `runs/idealized/`
Not a real river — a synthetic fixture with analytic, time-varying forcing used to
catch regressions in transport, carbonate, ice, or the boundary machinery. `config.
IS_IDEALIZED` warns on every such run so its output is never mistaken for a result.
```bash
python tools/build_idealized_forcings.py       # (re)generate the 13 analytic CSVs
python tools/verify_idealized.py --check-only  # structural only, no run
python tools/verify_idealized.py               # + 5-day physical-invariant run
python tools/verify_idealized.py --full        # + 220-day seasonal run
```

## 2. Long-term multi-year runs

### `runs/interannual/<site>/`
44 real years (1980-2023) of discharge + riverine DOC forcing (PWBM); every other
forcing category (meteorology, marine boundary, tides) still repeats one 2022-typical
year. Daily output cadence (not 6-min — a 44-year run at 6-min would be tens of GB per
site). Feeds `docs/ns_rad_interannual.pdf`.
```bash
python tools/build_interannual_forcings.py     # rebuild discharge/DOC forcing first
tools/run_interannual.sh colville kuparuk sagavanirktok
python tools/make_interannual_pdf.py
```

### `runs/interannual_full/<site>/`
19 real years (2005-2023) — the shortest common window across every independent
forcing source (wind/air-temp from NDBC PRDA2 is the binding constraint; station
didn't exist before 2005). **Every** forcing category is genuinely multi-year here:
discharge+DOC (PWBM), wind/air-temp (PRDA2, hourly), humidity (Deadhorse ISD, hourly),
solar (NOAA GML Barrow, hourly), marine boundary chemistry — S/DIC/ALK/NO3/NH4/PO4/
dSi/O2/TOC/sea-T (ECCO-Darwin v5), storm surge (NOAA CO-OPS 9497645), river
temperature (regional air->water regression), and tides (harmonic, computed
analytically — no file needed). Daily output cadence, same reasoning as above. This is
the run `docs/ns_rad_full_project_report.pdf` Part 1 characterizes, and the run the
FORTE campaign climatology (section 3 below) is built from.
```bash
python tools/build_interannual_met.py          # + build_interannual_{humidity,solar_barrow,surge,river_temp,marine}.py
tools/run_interannual_full.sh colville kuparuk sagavanirktok
python tools/make_publication_figures_full.py  # manuscript-style figure PDFs
```

### `runs/definitive_full/<site>/`
A narrow-purpose run: the `_interannual_full` site variants re-run for only the first
730 days (real 2005-2006) at the normal 6-minute cadence, instead of the full 19-year
run's daily cadence. Exists solely to give `make_publication_figures_full.py`'s Fig 3
(estuarine filtering / tidal-aliasing correction) a 6-min-resolution reference from the
*same* 2005-2023 forcing window — `runs/definitive/` is a repeated-2022-climatology
run and cannot serve that role.
```bash
tools/run_definitive_full.sh
```

## 3. FORTE campaign framework

### `runs/forte_campaign/<site>/`
A single climatological "typical campaign year": every forcing input from
`runs/interannual_full`'s 19-year record, averaged day-of-year (or hour-of-year) across
all 19 years into one smooth repeating year. Built this way rather than from one
analog historical year or the real (unreachable, 2027+) campaign dates — see
`tools/build_forte_campaign_forcing.py`'s docstring and `docs/ns_rad_full_project_
report.pdf` Part 2 for the full rationale. `MAXT_DAYS=638` (365-day warm-up + 273 days
covering the full real campaign season, day-of-year 0-273), 6-minute cadence (needed
for the Green's-function work in section 4).
```bash
python tools/build_forte_campaign_forcing.py   # climatology CSVs from runs/interannual_full
tools/run_forte_campaign.sh
```

## 4. Green's-function / OSSE sensitivity analysis

Two **parallel, independent** frameworks exist — the original analog-year one and the
newer FORTE-climatology one. Neither touches the other's files.

### Original (analog-year) framework
- **`runs/osse_nature/<site>/`** — one continuous "truth" run per site spanning real
  2005-2009, from which 3 analog study years (2007 low-flow / 2008 average / 2009
  high-flow) are sliced for sampling. `tools/run_osse_nature.sh`.
- **`docs/osse/<site>_<channel>.csv`** — synthetic sampling of that truth by each of
  FORTE's 5 planned measurement types (grab, moorings, drifters, airborne, satellite).
  `python tools/build_osse_observations.py`.
- **`runs/osse_green/pert_<control>/<site>/`** — 11 boundary-only Green's-function
  perturbation runs (discharge + 9 marine + 1 marine-T) per site, ±10% one-sided
  finite difference, baseline is its own separate run in the same tree.
  `python tools/build_green_perturbations.py` then `python tools/build_green_
  functions.py`.
- **Ranking**: `python tools/rank_observation_value.py` ->
  `docs/osse/observation_value_ranking*.csv`.
- **Report**: `python tools/make_osse_report.py` -> `docs/osse/osse_sampling_report.
  pdf`.

### FORTE (climatology) framework — current, most detailed analysis
- **`docs/osse_forte/<site>_<channel>.csv`** — same 5 channels, sampled against the
  `runs/forte_campaign` climatological run instead of an analog year.
  `python tools/build_osse_observations_forte.py`.
- **`runs/forte_green/pert_<control>/<site>/`** — **25** controls per site (up from
  11): discharge + 9 marine (`clb_*`) + 8 upstream (`cub_*`) boundary controls, PLUS 6
  **reach-local** controls new to this round — `rate_aerdeg_{mouth,mid,head}` (local
  aerobic-degradation rate) and `lateral_toc_{mouth,mid,head}` (local lateral organic-
  carbon loading), each confined to one coarse third of the channel. Most controls
  reuse `runs/forte_campaign/<site>/` directly as baseline (no separate baseline run);
  the 3 `lateral_toc_*` controls instead diff against **`runs/forte_green/
  lateral_base/<site>/`**, a dedicated uniform-background-loading baseline (see
  `tools/build_green_perturbations_forte.py`'s "methodological note" for why an
  on/off perturbation would otherwise be an invalid finite difference).
  ```bash
  python tools/build_green_perturbations_forte.py   # launches all pending (site, control) runs
  python tools/build_green_functions_forte.py       # -> runs/forte_green/green_functions_<site>.npz
  ```
  This step is the most compute-intensive in the whole project: 25 controls x 3 sites
  (75 perturbation runs) + 3 `lateral_base` baselines = 78 full 638-day runs.
  `build_green_perturbations_forte.py` skips any (site, control) whose run already
  completed, so it is always safe to re-run after adding new controls. `MAX_PARALLEL`
  inside that script currently defaults to `1` (fully serial) — lowered from `8` during
  development after repeated OOM kills under this host's memory pressure from other
  running applications; raise it if running on a host with more headroom.
- **Ranking**: `python tools/rank_observation_value_forte.py` ->
  `docs/osse_forte/observation_value_ranking*.csv`.
- **Report**: `python tools/make_forte_campaign_report.py` ->
  `docs/osse_forte/forte_campaign_report.pdf` (22 pages, FORTE-only).

### The combined report
`python tools/make_full_project_report.py` -> `docs/ns_rad_full_project_report.pdf`
(34 pages) — reads `runs/interannual_full` directly (Part 1), plus everything the
FORTE framework above produces (Parts 2-5), into one standalone, start-to-finish
document. Requires Parts 1-5's underlying runs/CSVs to already exist; builds no model
output of its own.

## 5. Historical comparison archives

Pre-date this round of work; kept for diffing specific isolated changes, not rebuilt:

- `runs/singlechannel_archive/<site>/` — pre-multichannel-geometry state.
- `runs/multichannel_preicefix/<site>/` — post-multichannel, pre-ice-grounding-fix.
- `runs/experiment_multichannel/` — the multichannel sensitivity sweep + distributary
  diagnostic.

`tools/compare_adoption.py` and `tools/compare_icefix.py` diff these against
`runs/definitive/` to isolate one variable change at a time — see CLAUDE.md for which
diffs which.
