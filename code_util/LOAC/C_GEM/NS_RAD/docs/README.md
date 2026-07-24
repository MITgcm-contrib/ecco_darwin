# `docs/` — index

Two kinds of file live here:

- **Design docs** (Markdown) — hand-written and the source of truth. Edit these directly.
- **Everything else** (PDFs, `.docx`, movies, `.json`, `.png`) — **generated artifacts**,
  reproducible from a model run. Do not hand-edit; regenerate with `tools/build_all.sh`
  (see [Regenerating](#regenerating) below).

**None of the generated PDFs, movies or the `.docx` are tracked in git** — they are gitignored
build output (the report alone is ~32 MB) and will not exist in a fresh clone until you run
`tools/build_all.sh`. The small `.json` provenance files *are* tracked, because the site
configs and this documentation cite numbers from them.

## Design docs (hand-written)

**Model reference:** [`model_description.md`](model_description.md) — full technical description (equations, components, parameterizations, parameter tables), in the style of the MITgcm/Darwin package docs.

**Start here:** [`FEATURES.md`](FEATURES.md) — the full catalog of everything NS-RAD adds on top of the original upstream C-GEM code.


| File | What it covers |
|---|---|
| [`performance.md`](performance.md) | Optimization history — three passes (hydro/density jit → per-cell jit → Python-orchestration jit), ~2.3× overall, all bit-identical; the verification harness and its cache/A-B caveats; what is left. |
| [`ice_model_plan.md`](ice_model_plan.md) | River/estuary ice model: the sea-ice-vs-river-ice framing problem, forcing data, the tiered implementation, and the prognostic thermal coupling. |
| [`arctic_biogeochemistry.md`](arctic_biogeochemistry.md) | The Arctic biogeochemistry extension — refractory/chromophoric DOC + CDOM photomineralisation, CH4 & N2O cycling, benthic (SOD) efflux, and distributed lateral loading; opt-in via `config.ARCTIC_BGC`. |
| [`idealized_verification.md`](idealized_verification.md) | The `idealized` verification fixture — a synthetic, time-varying, analytically-known site that exercises transport, the carbonate solve, the ice model, and the `BOUNDARY_FORCING` machinery. |
| [`multichannel_test.md`](multichannel_test.md) | Multi-channel geometry (ADOPTED, on by default): why `B` is the total conveyance width while dispersion uses a separate per-thread width, the SWORD braided-total extraction, the thread-count sensitivity, and the delta-distributary diagnostic that was *not* adopted. |

## Generated reports (regenerable)

| File | What it is | Built by |
|---|---|---|
| **`ns_rad_report.pdf`** | **The deliverable — the six section PDFs below, stitched into one document.** The only report PDF kept on disk after a build; not tracked in git. | `tools/make_report.py` |
| `idealized_verification.pdf` | Idealized-fixture report (kept separate — not part of the combined report) | `tools/make_idealized_verification_pdf.py` |
| `NS-RAD_running_instructions.docx` | Step-by-step run guide | `tools/make_run_instructions_docx.py` |
| `movies/<site>_evolution.mp4`, `movies/<site>_hovmoller.mp4` | Along-channel and distance×time animations | `tools/make_movies.py` |
| `idealized_forcings_preview.png` | Preview of the idealized forcing series | `tools/build_idealized_forcings.py` |

The six section PDFs — `ns_rad_model_summary`, `ns_rad_model_schematic`, `ns_rad_geometry`,
`ns_rad_river_networks`, `ns_rad_diagnostics`, `ns_rad_validation` — are built as **throwaway
intermediates** by their `tools/make_*.py` scripts, stitched into `ns_rad_report.pdf`, then
deleted. They are gitignored; rebuild the report with `tools/build_all.sh --figures-only`
(which recreates them, stitches, and removes them), not `make_report.py` on its own.

## Data files (tool inputs / outputs, regenerable)

| File | Role | Produced by |
|---|---|---|
| `performance_metrics.json` | Per-site wall times for the summary panel | `tools/bench_timing.sh` |
| `validation_obs.json` | Observations plotted in the validation PDF | (validation data prep) |
| `usgs_velocity_obs.json` | USGS velocity observations | (geometry/validation prep) |
| `sword_widths.json` | SWORD v17 channel-width extraction | (one-off SWORD extraction) |

## Regenerating

One command rebuilds the model run and every artifact above:

```bash
tools/build_all.sh                  # run the 4 rivers, then all figure PDFs, the combined
                                    #   report, and the movies (~15 min; model run dominates)
tools/build_all.sh --figures-only   # skip the model run; rebuild figures from existing runs/
tools/build_all.sh --with-idealized # ALSO run + verify + figure the idealized fixture
SERIAL=1 tools/build_all.sh         # run the four rivers one at a time
```

Individual artifacts can be rebuilt with their `tools/make_*.py` script directly (see the
tables above). `tools/make_report.py` only stitches the six existing report PDFs together —
it does not itself run the model or make figures.

> Note: this folder is a flat mix of authored and generated files by design (kept simple
> while the model is under active development). The Markdown docs are the only files here
> meant to be edited by hand.
