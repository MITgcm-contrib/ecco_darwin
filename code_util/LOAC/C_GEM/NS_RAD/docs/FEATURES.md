# NS-RAD — features added since the original C-GEM code

NS-RAD started from the upstream **C-GEM North Slope** code vendored from
`MITgcm-contrib/ecco_darwin` (`code_util/LOAC/C_GEM/North Slope`) — a single-site,
tide-estuary reactive-transport model (Volta et al. 2014) that shipped as plain Python
files with hardcoded paths and one river's configuration.

This document catalogs everything added on top of that starting point. Each item notes
**was → now** where relevant and points to the detailed write-up. The authoritative
per-decision reference is [`CLAUDE.md`](../CLAUDE.md).

---

## 1. Multi-river configuration & run infrastructure

- **Four rivers, one codebase** — `config.py` selects a river at import from `CGEM_SITE`;
  only what differs (geometry, tides, discharge file, boundary chemistry) lives in
  `code/sites/<river>.py`. *Was:* a single hardcoded site.
- **Portable forcing paths** — resolved relative to the code tree, not the original
  author's absolute machine paths; a clear error on a missing file.
- **NetCDF output** (`output.nc`, compact, self-describing) as the default, with the
  legacy tab-separated `.dat` path preserved (`CGEM_OUTPUT=dat|both`).
- **Environment-variable controls** — `CGEM_MAXT_DAYS` / `CGEM_WARMUP_DAYS` (smoke tests
  without editing `config.py`), `CGEM_TS` (output cadence), `CGEM_OUTPUT`, `CGEM_ICE`,
  `CGEM_MULTICHANNEL`, `CGEM_N_CHAN_UP`, `CGEM_DISTANCE`, `CGEM_FONT_SCALE`.
- **Multi-river runner** — `tools/run_sites.sh` runs each river as its own process in its
  own output directory (required: state is global, allocated at import), in parallel.
- **Run outputs organized by type** — `runs/definitive/`, `runs/regression_bnd/`,
  `runs/idealized/`.

## 2. Channel geometry & hydrodynamics

- **Observation-based geometry** — channel **width** from SWORD v17b nodes (per-channel,
  braiding/distributary corrections) and **depth** from USGS ADCP surveys via at-a-station
  hydraulic geometry `D = c·Q^f` (~1–2 m). *Was:* a single hybrid placeholder (15 m deep).
  See CLAUDE.md → *Geometry*.
- **Flare + prismatic width law** (`WIDTH_MODEL="flare"`) — width converges over a short
  delta flare then is prismatic, which fits these rivers by AIC. *Was:* whole-domain
  Savenije exponential (poor fit, R² 0.02–0.23).
- **Multi-channel geometry** (`CGEM_MULTICHANNEL`, on by default) — `B` is the **total**
  conveyance and water-surface width (summed over parallel braids/distributaries, with the
  prismatic end read directly from SWORD's raw width), while a separate per-thread width
  `B/n_chan` drives the shear-dispersion closure. *Was:* one width field whose definition
  changed along the domain — a distributary **sum** at the mouth but **one of n** braids
  upstream — while carrying the total discharge through both, which over-estimated braided
  velocity by ~`n_chan` and under-estimated water-surface area by the same factor. Roughly
  doubles the basin-integrated carbon flux on the two braided rivers (Colville 2.19×,
  Sagavanirktok 1.57×) and leaves the two single-thread rivers unchanged; also shifts their
  per-area flux (−16% to +7%) via the √U dependence of the gas-transfer velocity.
  `=off` reverts bit-identically. See [`multichannel_test.md`](multichannel_test.md).
- **Seo & Cheong (1998) longitudinal dispersion** (`DISPERSION_MODEL="seo"`) computed per
  timestep from local width/depth/velocity, with `u*` from the model's own friction.
  *Was:* Van der Burgh / Savenije estuary form, which collapsed to zero under realistic
  shallow geometry (purely advective transport). See CLAUDE.md → *Width and dispersion*.

## 3. Observed 2022 forcings

- **Per-river discharge** — USGS NWIS gauges (byte-reproducible builder); **Canning
  reconstructed** from a donor gauge (no 2022 record).
- **River water temperature** — modelled from air temperature (regional regression) +
  USGS 00010 observations for Kuparuk/Sag. Fixes a −6 to −8 °C cold bias and a broken ice
  gate the shipped sea-water `watertemp.csv` caused.
- **Meteorology** — NDBC PRDA2 (Prudhoe Bay) wind/solar/air; Deadhorse Airport humidity;
  Barrow atmospheric pCO₂.
- **Per-river tides** — multi-constituent **harmonic reconstruction** from the nearest
  NOAA CO-OPS station (`tools/build_tides.py`). *Was:* a single shipped sinusoid.
- **Wind-driven storm surge** — a `SURGE` offset added to the marine boundary elevation
  (`tools/build_surge.py`), Prudhoe Bay as the **regional proxy for all four rivers** — the
  actual saltwater-intrusion driver on this microtidal coast.
- **Boundary chemistry** — Arctic LTER (Kuparuk headwater proxy) and WQP grabs
  (Colville, Sagavanirktok); DIC solved from observed pH + alkalinity.

## 4. Temperature, heat & ice — new physics

- **Transported temperature field** — `T` is species 13, advected/dispersed like salinity,
  with time-varying sea/river boundaries. *Was:* a single scalar applied uniformly.
- **Surface heat budget** (`heat_module.py`) — `ρ cp H dT/dt = Q_sw+Q_lw+Q_sens+Q_lat`
  warms/cools the interior; makes mid-channel temperature physical, not a boundary blend.
- **Prognostic river-ice model** (`ice_module.py`) — freeze-up from the heat-budget
  deficit, conductive (Stefan) growth, surface melt, and **hydraulic (freshet) break-up**
  (the North-Slope-specific mechanism), with **bottom-fast grounding**: depth limits how much
  new ice can form, but a falling water level never destroys ice already there. *Was:*
  thickness truncated to the live water depth every step — a one-way ratchet against a
  ~0.8 m/day tide+surge swing that thinned the cover through midwinter instead of growing it
  (98% of midwinter thinning sat at that cap). Couplings: insulates the
  heat budget, shuts O₂/CO₂ gas exchange, attenuates under-ice PAR, and **conserves** state
  year-round instead of zeroing it under ice. See [`ice_model_plan.md`](ice_model_plan.md).

## 5. Carbonate chemistry — unit-correct solve + fixes

- **mol/kg carbonate system** — DIC/ALK converted to mol/kg via in-situ density before the
  Follows et al. (2006) speciation, matching the mol/kg dissociation constants
  (`CARBONATE_UNITS`). *Was:* the legacy path mixed mmol/m³ state with mol/kg constants and
  patched only the borate term.
- **pH-solver defect fixes** — the iteration actually iterates (guess hoisted out of the
  loop), and negative/zero-H⁺ pathologies at ice-out are guarded (floors + physical bound
  → neutral fallback), eliminating the non-physical pH the vendored code produced.

## 6. Arctic biogeochemistry extension *(opt-in)*

Four process groups that represent what actually controls Arctic land-to-ocean
CO₂/CH₄/N₂O — gated behind `config.ARCTIC_BGC` (default off; real-river fields stay
bit-identical to the shipped network when off). Full write-up:
[`arctic_biogeochemistry.md`](arctic_biogeochemistry.md).

- **Refractory/chromophoric DOC (`RDOC`) + CDOM photomineralisation** — a second DOC pool
  with slow microbial oxidation *and* sunlight-driven photo-oxidation to DIC (Cory et al.
  2014); `TOC` remains the labile pool.
- **Methane and nitrous oxide** — `CH4`/`N2O` tracers with methanotrophy, N₂O yields from
  nitrification/denitrification, and air–water gas exchange (Wanninkhof 2014 Schmidt
  numbers; Wiesenburg-Guinasso / Weiss-Price solubilities).
- **Benthic efflux** — a sediment-oxygen-demand closure returning DIC / alkalinity / CH₄
  from the bed (runs under ice → builds up and vents at ice-out).
- **Distributed lateral loading** (`lateral_module.py`) — tundra/thermokarst inputs entering
  along the channel, not only at the upstream boundary (`LATERAL_INFLOW`/`LATERAL_CONC`).

The species registry grew from 12 → **16** transported tracers.

## 7. General time-varying boundary conditions

- **`BOUNDARY_FORCING`** — any species' downstream (`clb`) and/or upstream (`cub`) boundary
  can follow a 365-day forcing series, refreshed every timestep. *Was:* only temperature was
  seasonal (special-cased); all solute boundaries were constant.
- **Per-site meteorology filenames** — the shared met records are now per-site overridable,
  so a site can supply its own wind/solar/air/humidity/pCO₂/sea-T (used by the idealized
  fixture); the four real rivers are unaffected via `getattr` defaults.

## 8. Verification & testing

- **Idealized verification fixture** — a fifth, synthetic site (`idealized`) driven by
  analytic, time-varying forcing and boundaries; exercises the whole model end-to-end with
  a known input. See [`idealized_verification.md`](idealized_verification.md).
- **Invariant test harness** — `tools/verify_idealized.py` runs the fixture and asserts a
  suite of physical invariants (mass conservation, ice freeze-up/break-up, boundary
  tracking, gas-flux signs, the Arctic extension); exits non-zero on any failure. **82/82
  checks pass.**

## 9. Performance

- **Numba JIT** of the hot loops (hydrodynamic kernels, density stack, transport schemes,
  biogeochemistry, sediment, heat, ice, carbonate) — held **bit-identical** to the
  pre-optimization output at each step. ~14–15 min per river (2-year run), four in parallel
  ~15 min. See [`performance.md`](performance.md) for the history, profiles, and the
  bit-identity harness.
- **Forcing cache** — parsed forcing series are cached instead of re-parsed every timestep.

## 10. Tooling & reproducibility

- **One-command pipeline** — `tools/build_all.sh` runs all four rivers then regenerates
  every PDF and movie (PASS/FAIL/SKIP summary; `--figures-only`, `--with-idealized`).
- **Data-build scripts** — reproducible builders for every forcing: `fetch_discharge.py`,
  `build_tides.py`, `build_surge.py`, `build_river_temp.py`, `build_humidity.py`,
  `lter_boundary.py`, `build_boundary_chem.py`, `build_idealized_forcings.py`.
- **Figure/PDF/movie generators** — `make_diagnostics_pdf`, `make_validation_pdf`,
  `make_geometry_pdf`, `make_summary_figures`, `make_river_maps`, `make_schematic_pdf`,
  `make_idealized_verification_pdf`, `make_movies`, `make_report` (combined), plus a shared
  publication/presentation style (`nsrad_style.py`: embedded fonts, presentation sizing).
- **Run-instructions Word guide** — `tools/make_run_instructions_docx.py`.

## 11. Documentation

- **`CLAUDE.md`** — the full developer guide: architecture, provenance, and every
  non-obvious decision, with "was → now" reasoning throughout.
- **`README.md`** (repo landing), **`docs/README.md`** (docs index), and the topic docs:
  this file, `arctic_biogeochemistry.md`, `idealized_verification.md`, `performance.md`,
  `ice_model_plan.md`.

## 12. Fixed defects in the vendored code

Documented in CLAUDE.md → *Known defects*. Highlights: the `fun_module.pH` solver (three
defects — swapped args, non-iterating loop, ice-out blow-up), the `pCO2_barrow`/`pCO2_Barrow`
filename-case bug, the hardcoded absolute forcing paths, and several documented (deliberately
un-patched) unit inconsistencies (e.g. the boron term) left as flagged physics passes.

---

*Project rebranded to **NS-RAD** (North Slope River-Aquatic-Delta Model); "C-GEM" now refers
to the underlying estuarine reactive-transport engine (Volta et al. 2014) it is built on.*
