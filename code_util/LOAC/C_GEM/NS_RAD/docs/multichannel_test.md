# Multi-channel geometry

**Status: ADOPTED — on by default.** `CGEM_MULTICHANNEL=off` reverts to the legacy
single-channel geometry, and that path is bit-identical to the pre-adoption code
(verified across all 36 output variables on a 3-day Colville A/B with `CGEM_OUTPUT=dat`
and cleared `__pycache__`). The definitive runs use the multi-channel geometry; the
pre-adoption runs are kept at `runs/singlechannel_archive/` and
`tools/compare_adoption.py` diffs the two.

**Adopted for the surface-area / basin-flux correction, which is robust. NOT adopted as a
fix for the salinity misfit, which it reduces but does not close and does not constrain —
see the Sensitivity section before quoting any salinity number.**

The second change tested alongside this one — resolving the delta into separate
distributary runs — was **not** adopted. It remains a diagnostic (`sites/colville_main.py`,
`colville_minor.py`); the channels cannot exchange water and each is extended over the
full domain. It earns its place by showing that lumping is not flux-neutral.

## The problem

C-GEM has one width field `B(x)`, and it is asked to do two physically different jobs:

| role | wants | used by |
|---|---|---|
| conveyance / water surface | total width summed over **parallel** threads | continuity, momentum, cross-section, surface area |
| shear dispersion | **per-thread** width (a within-channel process) | `fun_module.river_dispersion`, Seo & Cheong (1998) |

Worse, in the flare `B`'s *definition silently changes* between them:

- `B_lb` = the SWORD delta distributary **SUM** — total across parallel paths
- `B_ub` = the SWORD **per-channel** median — *one of n* parallel braids

and the model carries the **total** discharge through both. In a braided reach that
divides total Q by one braid's area, over-estimating velocity by ~n_chan and
under-estimating residence time and water-surface area by the same factor.

## The two changes

### 1. Width-role separation (`CGEM_MULTICHANNEL=on`)

`B` becomes the total conveyance/surface width everywhere — its prismatic end is the
observed braided total `B_UB_TOTAL` — and a separate per-thread width

    B_thread(x) = B(x) / n_chan(x)

feeds the dispersion closure alone. `n_chan` blends geometrically from `N_CHAN_LB`
(delta distributaries) to the derived prismatic thread count over `L_FLARE`.

Two consistency properties, both verified:
- upstream `B_thread` reduces **exactly** to the site's existing `B_ub` (423 / 58 / 102 / 64 m);
- at the mouth `B_thread` = `B_lb / N_CHAN_LB`, the mean distributary width.

`ZZ = B * depth`, so scaling `B` scales the cross-section and leaves `DEPTH` untouched.
Velocity and residence time move; the per-area gas exchange (`FCO2 = RCO2/DEPTH`) does
not. What moves is the **area** any basin-integrated flux is multiplied by.

### 2. Distributary resolution (`sites/colville_main.py`, `colville_minor.py`)

A two-outlet delta run as two separate sites, exploiting the fact that the model already
runs one process per site — the multi-**site** machinery doubling as multi-**channel**
machinery, with no change to the single-channel solver. `Q_FRACTION` gives each channel
its conveyance share. See `sites/_colville_delta.py` for the full derivation.

The partition is **derived, not assumed**. Manning with shared slope and roughness gives
`Q_i ~ W_i·D_i^(5/3)`; Colville's own at-a-station geometry gives `D = 0.360·Q^0.297`.
Eliminating D:

    Q_i^(1 - 5f/3) ~ W_i   ->   Q_i ~ W_i^1.98

so discharge partitions as nearly the **square** of width. The minor channel is 18% of
the mouth width but takes only 3.3% of the flow, and its depth then follows from that
small share:

| channel | W (m) | Q share | Q (m³/s) | D (m) | U_mouth |
|---|---|---|---|---|---|
| main | 1310 | 0.967 | 468.2 | 2.236 | 0.160 |
| minor | 240 | 0.034 | 16.3 | 0.824 | 0.082 |

A lumped 1550 m channel cannot express that spread.

## The data: read the braided total, don't reconstruct it

`tools/extract_sword.py` now also emits the **prismatic braided total** — the median RAW
(un-divided) SWORD width beyond the flare. SWORD's raw width is *already* the sum across
braids, so this is a direct observation.

| river | `B_ub` (per-channel) | raw braided total | IQR | implied threads | median `n_chan_mod` |
|---|---|---|---|---|---|
| Colville | 423 | **1052** | 792–1354 | 2.49 | 3.00 |
| Kuparuk | 58 | 60 | 42–78 | 1.03 | 1.00 |
| Sagavanirktok | 102 | 224.5 | 146–318 | 2.20 | 2.00 |
| Canning | 64 | 63 | 53–78 | 0.98 | 1.00 |

**The implied thread count is not the median `n_chan_mod`, and that matters.** `B_ub` is a
median of per-node *ratios*, and the median of a ratio is not the ratio of medians:
Colville's reconstruction `423 × 3.0 = 1269 m` overstates the observed 1052 m by 20%.
An earlier pass here used `423 × 3.5 = 1480 m`, overstating it by 40% and inflating every
downstream result. Sites therefore declare the **observed** `B_UB_TOTAL` and `config`
derives the thread count from it — the reconstruction is out of the loop entirely.

Consequence for the "flare is an artifact" claim: Colville runs 1550 → 1052 m, so roughly
**half** its flare is real convergence and half is the definition change — not "very
nearly all," as the first pass suggested. Kuparuk (60 vs 58) and Canning (63 vs 64) are
genuinely single-thread, so their flares are entirely real; they are clean controls.

## Results

Colville, year 2, `tools/compare_multichannel.py` (reads `runs/experiment_multichannel/`).
The analysis reproduces the published baseline flux to within 1% (93.8 vs +94.6
gC m⁻² yr⁻¹), so these sit on the same footing as the documented numbers.

| run | U_mouth | S summer med | S max | % open water >1 PSU | gC m⁻² yr⁻¹ | tC yr⁻¹ | surface km² |
|---|---|---|---|---|---|---|---|
| baseline (lumped) | 0.211 | **0.00** | 11.6 | 3.3 | 93.8 | 1321 | 14.15 |
| multichannel (n=2.49) | 0.210 | 1.40 | 27.1 | 42.2 | 96.8 | 2903 | 30.19 |
| main distributary | 0.243 | 3.58 | 27.1 | 48.3 | 90.5 | 2579 | 28.66 |
| minor distributary | 0.088 | 0.00 | 26.5 | 25.3 | 59.2 | 73 | 1.38 |
| delta average | 0.236 | **3.42** | 27.1 | 47.2 | 89.1 | 2653 | 30.04 |

### Sensitivity to the thread count — the load-bearing result

| n_chan_up | B_up (m) | S summer med | S max | intrusion (km) | gC m⁻² yr⁻¹ | tC yr⁻¹ |
|---|---|---|---|---|---|---|
| 1.00 | 423 | 0.00 | 11.6 | 0.07 | 93.8 | 1321 |
| 1.87 | 791 | 0.02 | 26.9 | 0.85 | 96.5 | 2263 |
| **2.49** | **1053** | **1.40** | **27.1** | **1.46** | **96.8** | **2903** |
| 3.20 | 1354 | 3.23 | 27.1 | 1.94 | 96.4 | 3619 |
| 3.50 | 1480 | 3.85 | 27.1 | 2.09 | 96.2 | 3916 |

**Read this table before trusting anything above it.** Across the *observational* IQR
(1.87–3.20) the summer median salinity spans 0.02 → 3.23 PSU. The typical intrusion is
therefore essentially undetermined by the SWORD data — it is not a tuned parameter, but
its observational spread covers nearly the whole range of outcomes.

What **is** robust across that spread:

- **Basin-integrated flux rises**, by 1.7–3.0× (2.2× at the observed value), tracking
  surface area. This is the main practical result.
- **The model becomes *able* to intrude to full marine salinity.** `S_max` is 26.9–27.1
  for every n ≥ 1.87 against the baseline's 11.6 — the baseline could not reach the marine
  boundary value at all. The *capability* is robust; the *typical* value is not.

### Adoption: all four rivers (`tools/compare_adoption.py`)

| river | surface area | per-area gC m⁻² yr⁻¹ | basin tC yr⁻¹ | threads |
|---|---|---|---|---|
| Colville | 2.13× | 93.8 → 96.8 (1.07×) | 1321 → 2900 (**2.19×**) | 2.49 |
| Sagavanirktok | 1.79× | 91.2 → 76.5 (**0.84×**) | 339 → 531 (1.57×) | 2.20 |
| Kuparuk *(control)* | 1.02× | 176.4 → 175.7 (1.00×) | 429 → 437 (1.02×) | 1.03 |
| Canning *(control)* | 1.00× | 105.7 → 105.8 (1.00×) | 337 → 337 (1.00×) | 1.00 |

**Correction to an earlier claim in this document: the per-area flux is NOT invariant.**
The reasoning that `ZZ = B·depth` leaves `DEPTH` untouched, so `FCO2 = RCO2/DEPTH` cannot
move, is incomplete — it covers the `/DEPTH` conversion but ignores that `RCO2` itself
depends on velocity. The flow-driven piston velocity goes as **√U**
(`fun_module._piston_velocity_loop`: `kflow = sqrt(|U|·D_O2/DEPTH)`), so widening the
channel to its true cross-section slows the water and reduces gas exchange per unit area;
the 1.5–1.8× longer transit then shifts the along-channel DIC/ALK balance. Decomposition:

| | transfer velocity `vp` | CO₂ excess (DIC−ALK) | net |
|---|---|---|---|
| Sagavanirktok | 0.961× | 0.917× | 0.88× |
| Colville | 0.928× | 1.153× | 1.07× |
| Kuparuk | 0.999× | 0.998× | 1.00× |

Per-volume reaction rates barely move (Sag `denit` 0.999×, `aer_deg` 0.998×) — it is the
residence time, not the kinetics, that changes.

**The braided rivers' pre-adoption per-area fluxes were therefore biased high**: their
too-narrow cross-section made the water flow ~`n_chan` too fast and exchange gas too
vigorously. Sagavanirktok shows it cleanly (−16%); on Colville a compensating chemistry
change masked it. The two single-thread rivers are unchanged to three significant figures,
which is the control confirming the mechanism is the velocity change and nothing else.

### What this does and does not do for the salinity misfit

It helps, more than any lever tried before it (mouth width, storm surge, and the
`distance`/Chezy sweep each moved summer median salinity by <2 PSU and left it at 0), and
it is the only one that lets the mouth reach marine values. But at the observed geometry
the summer median is 1.40 PSU lumped / 3.42 PSU resolved, against observed grabs of
8–32 PSU. **The misfit is reduced, not closed.**

Resolving the two distributaries adds about as much again as the width-role fix (3.42 vs
1.40 at the same total conveyance), because Seo & Cheong scales as `(W/H)^0.62` and is
therefore nonlinear in width — splitting one 1550 m channel into 1310 + 240 m is not
flux-neutral.

**A hypothesis this refuted.** The experiment was designed expecting the narrow, slow
minor distributary (U_mouth 0.088 vs 0.243 m/s) to be the one holding salt. It is the
*freshest* (summer median 0.00 vs the main channel's 3.58). Intrusion here is
dispersion-dominated, and the wide fast channel disperses salt inland harder than the
narrow slow one. The gain comes from removing the spurious near-mouth velocity jet and
from the width-nonlinearity above — not from the velocity-contrast mechanism proposed.

## Caveats

- **The distributary runs are a diagnostic, not a delta model.** Each channel is extended
  over the full 27 km domain because C-GEM needs an upstream boundary; the stream-tube
  construction keeps that reach mass- and velocity-consistent with the trunk, but the two
  channels **cannot exchange water**. The minor channel also acquires its own convergence
  jet by construction (240 → 35 m).
- **Sagavanirktok is the weakest site** — braided, with a 5–20 km SWORD coverage gap.
  Indicative only.
- **`DISP_MAX` capping is pre-existing**, not introduced here: the baseline Colville run
  carries the same once-per-run warning. CLAUDE.md's claim that it "does not currently
  bind for any site" is stale.

## Verifying a run actually used the geometry you think

`config` prints the multichannel state into every run log, and prints a distinct warning
if `CGEM_MULTICHANNEL` is set to anything that is not exactly `on`. **Check it before
reading any result:**

```bash
grep multichannel runs/experiment_multichannel/<run>/run.log
#   [multichannel] ON  n_chan: 2.00 (mouth) -> 2.49 (prismatic); B_ub total = 1053 m ...
```

This exists because of a real failure. **zsh does not word-split unquoted expansions**, so
`env CGEM_SITE=x $EXTRA python3 main.py` with `EXTRA="CGEM_MULTICHANNEL=on CGEM_N_CHAN_UP=1.87"`
passes *one* assignment, `CGEM_MULTICHANNEL="on CGEM_N_CHAN_UP=1.87"`, which is not `on` —
so the flag silently evaluates false and the run produces output **byte-identical to the
baseline**. That reads as "the change did nothing" rather than "the flag never arrived,"
and it invalidated a whole sweep before being caught (all three sweep points returned the
baseline value to the last digit, which no physical sensitivity does). Write each launch
with literal, separate `VAR=value` assignments.
