# Idealized verification experiment

A fifth, synthetic "site" — `idealized` — that exercises the whole model on a
**known, analytic, time-varying** input over a full seasonal cycle. It is a testing
fixture, **not a river**: nothing in it is measured, and its output must never be
interpreted as a result. `config.IS_IDEALIZED` makes the model emit a warning to that
effect on every idealized run.

Its purpose is regression/verification: run it, and a change that quietly breaks
transport, the carbonate solve, the ice model, or the boundary machinery is caught by
a physical-invariant check rather than by a surprised reader months later.

## What makes it "idealized"

| aspect | idealized fixture | the four real rivers |
|---|---|---|
| geometry | round numbers: 30 km, 2 m deep, a 400→100 m flare over 8 km | observation-based (SWORD + USGS) |
| discharge | analytic baseflow + a Gaussian spring freshet | USGS gauge records |
| meteorology | analytic sinusoids (own wind/solar/air/humidity/pCO₂) | shared regional records |
| air temperature | crosses 0 °C twice, so ice freezes up and breaks up | observed 2022 |
| **boundary chemistry** | **time-varying** — DOC/ALK/DIC/NO₃/salinity follow the season | **constant** (only T is seasonal) |
| tides | single-sinusoid fallback (`AMPL`/`pfun`), not in the constituent JSON | per-river harmonic reconstruction |

Every series is a closed-form function of day-of-year, built from one `PARAMS` block
in `tools/build_idealized_forcings.py`, so the fixture is fully reproducible — re-run
the builder and the CSVs come back byte-for-byte.

## The feature under test: time-varying boundaries

The four real rivers hold every solute boundary **constant** in time; only temperature
is seasonal, and `main.py` special-cases it. The idealized site is the first user of a
**general, file-driven** mechanism that drives *any* species' downstream (`clb`) and/or
upstream (`cub`) boundary from a forcing series, refreshed every timestep exactly as
temperature is. It is declared per-site:

```python
# sites/idealized.py
BOUNDARY_FORCING = {
    "S":   {"clb": "idealized_S_clb_marine.csv"},   # marine end freshens in summer
    "TOC": {"cub": "idealized_TOC_cub_river.csv"},  # riverine DOC flushed on freshet
    "NO3": {"cub": "idealized_NO3_cub_river.csv"},  # nitrate winter-high
    "ALK": {"cub": "idealized_ALK_cub_river.csv"},  # snowmelt dilution
    "DIC": {"cub": "idealized_DIC_cub_river.csv"},
}
```

This works because `schemes_module.openbound` already reads `v[s]["clb"]`/`["cub"]`
fresh on every call. `main.py` resolves each `(species, end)` to an absolute path once,
then updates `v[species][end]` each timestep from the interpolated series. The
`BOUNDARIES` table still seeds `t = 0` and supplies the constant value for any end not
listed. The mechanism is generic: a real river can adopt it the moment a seasonal
boundary record exists — just drop in a 365-day CSV and add a `BOUNDARY_FORCING` entry.

The chemistry is tied to the hydrograph so the test is physically legible:
- **DOC (TOC)** is flushed **high** on the freshet (tundra carbon), 250 → 750 mmol m⁻³.
- **Alkalinity and DIC** are **diluted** by low-alkalinity snowmelt on the freshet.
- **Nitrate** is winter-high, drawn down through the open-water season.
- **Marine salinity** freshens ~30 → 24 PSU under summer sea-ice melt (a `clb` test).

![forcings](idealized_forcings_preview.png)

## Running it

Build the forcings once (already committed; rebuild after editing `PARAMS`):

```bash
python tools/build_idealized_forcings.py            # writes forcing/idealized_*.csv
python tools/build_idealized_forcings.py --plot     # + docs/idealized_forcings_preview.png
```

The canonical fixture run lives at `runs/idealized/` and is produced by the verify
harness below (`--full`). The diagnostics PDF and movies read it from there. To run
the site by hand in its own directory:

```bash
mkdir -p runs/idealized && cd runs/idealized
CGEM_SITE=idealized CGEM_MAXT_DAYS=220 CGEM_WARMUP_DAYS=30 \
  PYTHONPATH=../../code python ../../code/main.py
```

(`tools/run_sites.sh idealized` also works, but it writes under `runs/definitive/` with
the four real rivers — the fixture's home is `runs/idealized/`, so prefer the harness.)

## Verifying

`tools/verify_idealized.py` runs the fixture and asserts physical invariants. It exits
non-zero if any check fails, so it drops straight into CI or a pre-commit hook.

```bash
python tools/verify_idealized.py --check-only   # structural only, no run, no numba needed
python tools/verify_idealized.py                # + short 5-day run (~1-2 min w/ compile)
python tools/verify_idealized.py --full         # + 220-day seasonal run (~10-15 min)
```

Checks, by mode:

- **structural** (always): every named forcing exists, is 365 values and finite; every
  `BOUNDARY_FORCING` species/end is valid; `IS_IDEALIZED` is set.
- **quick run**: the run exits 0 and writes a well-formed `output.nc`; no NaN/Inf and no
  negative concentrations in the interior; pH inside the guard range; winter ice grows;
  each boundary cell stays inside its species' seasonal envelope (plus a tidal band for
  the mouth cell).
- **full run** (adds): the freshet **clears the ice** (hydraulic break-up); the DOC
  pulse **reaches the upstream cell** (proves the `cub` mechanism); the marine cell
  **freshens in summer** (proves the `clb` mechanism); `cub` cells track their forcing
  to tolerance and the tidal `clb` cell correlates with its forcing on daily means;
  `FCO₂` is finite.

Why tracking is full-mode only: on a short winter window the boundary inputs barely
move, so there is nothing to track — the quick run can only confirm values stay in a
plausible envelope. The seasonal signals (freshet, break-up, freshening) are what
actually prove the time-varying machinery drives the model.

## Diagnostics PDF and movies

After a `--full` run (output at `runs/idealized/`):

```bash
python tools/make_idealized_verification_pdf.py     # docs/idealized_verification.pdf (5 pages)
python tools/make_movies.py --run runs/idealized idealized   # docs/movies/idealized_*.mp4
```

The PDF's verification-results table is generated by importing `verify_idealized` and
running its checks on the same `output.nc`, so the report can never disagree with the
harness. Its five pages: (1) the analytic forcings + the full check table; (2) the
time-varying boundary overlays (the headline); (3) distance × time Hovmöllers of the
core fields and ice; (4) the **Arctic biogeochemistry extension** — the new RDOC/CH₄/N₂O
tracers, photomineralisation, and the CH₄/N₂O ice-out venting; (5) the freshet/ice/
carbonate story. The two movies are the along-channel profile animation and the
Hovmöller with a sweeping day cursor.

The fixture also exercises the **Arctic biogeochemistry extension**
(`docs/arctic_biogeochemistry.md`): the refractory-DOC (RDOC), CH₄ and N₂O tracers, a
time-varying RDOC river boundary (a new-species `BOUNDARY_FORCING` test), and
`LATERAL_INFLOW = 8 m³ s⁻¹` of distributed tundra loading — so photomineralisation,
CH₄/N₂O gas exchange, benthic efflux and lateral loading are all verified end to end.

## Files

| file | role |
|---|---|
| `code/sites/idealized.py` | the site: idealized geometry, forcings, `BOUNDARY_FORCING` |
| `tools/build_idealized_forcings.py` | generates the 13 analytic CSVs from one `PARAMS` block |
| `tools/verify_idealized.py` | runs the fixture and asserts invariants |
| `tools/make_idealized_verification_pdf.py` | 5-page summary-diagnostics PDF (incl. Arctic tracers) |
| `forcing/idealized_*.csv` | the generated forcings + boundary series |
| `runs/idealized/output.nc` | the canonical full-seasonal fixture run |
| `docs/idealized_verification.pdf` | the diagnostics report |
| `docs/movies/idealized_*.mp4` | profile + Hovmöller animations |
| `docs/idealized_forcings_preview.png` | overview of all 13 signals |

The shared plumbing that made this possible — per-site meteorology filenames and the
`BOUNDARY_FORCING` hook — lives in `sites/_baseline.py` (defaults), `config.py`
(resolution, unchanged for the real rivers via `getattr` defaults) and `main.py` (the
per-timestep refresh loop). The four real rivers are byte-for-byte unaffected: they set
none of the new attributes, so they get the old hardcoded met files and an empty
`BOUNDARY_FORCING`.
