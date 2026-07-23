# River / estuary ice model for NS-RAD — implementation plan

Status: PROGNOSTIC ICE MODEL IMPLEMENTED (`ice_module.py`). The full-cover thermodynamic
model is live; only Tier 3's hydraulic cross-section feedback is deferred.

  - **Tier 0a (gate FCO2 under ice) — DONE**, and now superseded: gas exchange is scaled
    by the real ice fraction (below), so FCO2 and O2_ex fall to zero under cover by
    construction rather than by an else-branch patch. Plus the pH-conversion guard.
  - **§4b (a) transported temperature — DONE.** Advected field, not a scalar.
  - **§4b (b) surface heat budget — DONE** (`heat_module.py`). Now runs YEAR-ROUND under
    the ice model and skips ice-covered cells (insulation, §4b(c) item 2).
  - **§4b (c) ice as a latent-heat store — DONE** (`ice_module.py`): (1) latent-heat clamp
    — the per-step freeze energy `heat_module` books is converted to ice via rho_ice*L_f;
    (2) insulation — the atmosphere is decoupled under ice and the water is held at
    T_FREEZE; (3) melt buffering — surface melt (sensible + absorbed shortwave) thins the
    slab through spring before water can warm. Plus conductive (Stefan) growth to build the
    winter thickness, a bottom-fast thickness cap, and HYDRAULIC BREAKUP on the freshet
    surge (§5 decision 3 — the key North-Slope mechanism).
  - **Tier 1 (light + gas coupling) — DONE.** Under-ice PAR attenuated through the slab
    (`k_ice_PAR`); piston-velocity gas exchange scaled by `(1-ice_frac)`.
  - **Tier 2 (conserve state under ice) — DONE.** The winter zeroing in
    transport/biogeo/sed is retired under `ICE_MODEL`; state is conserved and keeps
    reacting, so respiration builds under-ice DIC that vents at breakup. Erosion is scaled
    by the open-water fraction. Legacy zeroing remains as the `CGEM_ICE=off` fallback.
  - `previousdays` is retired as the physical gate when `ICE_MODEL` is on (§5 decision 5);
    `ice_frac`/`ice_thickness` are the single authority. It still gates only the legacy path.
  - Tier 0b (river temperature forcing) — DONE earlier; see CLAUDE.md.

  - **DEFERRED — Tier 3 hydraulic coupling.** Ice draft is capped at the local depth
    (bottom-fast), but that cap does NOT yet feed back into the hydrodynamic cross-section
    `D[i]` or velocity — the flow solver is ice-unaware. This is the one invasive coupling
    (it touches the tridiagonal solve) and is held for a separate, carefully-verified pass.
    In practice near-zero winter discharge already makes under-ice advection small.

## 0. The framing problem: the model's "ice" is SEA ice, not river ice

This has to come first, because it invalidates the premise of the existing gate and reframes
everything below.

`forcing/watertemp.csv` is not river temperature. It is NDBC buoy **PRDA2 sea-water
temperature** at Prudhoe Bay, verified against `prda2h2022.txt`:

- correlation r = **0.9922**, mean absolute difference 0.093 °C, identical min/max
- it reaches **−2.00 °C**, the freezing point of *seawater*. Fresh river water freezes at 0 °C
  and cannot reach this value. 230 of 365 days sit below −1.5 °C.

So `previousdays` — the rolling sum that gates `transport`, `biogeo` and `sed` — is tracking
**landfast sea-ice conditions in the Beaufort coastal zone**, and applying them to four rivers.

### Consequence 1: the spring freshet is discarded

Coastal sea water stays sub-zero until late June; North Slope rivers break up and flood in late
May. The gate opens on day 173, but:

| Site | Peak Q | Peak day | Gate at peak | Annual volume discarded |
|---|---|---|---|---|
| Colville | 4147 m³/s | 153 | **CLOSED** | **45.1%** |
| Kuparuk | 1387 m³/s | 153 | **CLOSED** | **39.9%** |
| Sagavanirktok | 286 m³/s | 219 | open | 27.5% |
| Canning | 429 m³/s | 174 | open | 17.6% |

For Colville and Kuparuk the single largest hydrological event of the Arctic year — the freshet,
which carries the bulk of the annual terrestrial carbon and sediment flux — occurs ~20 days
*before* the model permits any transport or biogeochemistry, and is zeroed.

### Consequence 2: every temperature-dependent rate is biased cold

USGS logs actual river temperature (parameter 00010) at two of the four gauges. Comparing the
Sagavanirktok record against the sea temperature the model substitutes, over 93 overlapping days:

- river mean **8.63 °C**, sea-substitute mean **2.62 °C** — bias **+6.00 °C**, max +12.74 °C
- on day 163 the river is at 7.90 °C and flowing while the gate, reading −0.06 °C, is still shut

Propagated through the existing rate laws:

| Term | River | Sea-substitute | Error |
|---|---|---|---|
| `Fhet` respiration (Q10 = 2.75) | 1.466 | 0.798 | understated **1.84×** |
| `Fnit` nitrification (Q10 = 5) | 1.837 | 0.698 | understated **2.63×** |
| `O2sat` | 364.6 | 425.0 mmol/m³ | overstated **+16.6%** |

Carbonate constants `K0/K1/K2/KB` take the same biased temperature.

### Why river ice is a different model from sea ice

- freezes at 0 °C, not −1.8 °C; no brine, no salinity-depressed freezing point
- breakup is **mechanical and hydraulic** — driven by the snowmelt surge, ice jams and
  overtopping — not by gradual thermal decay. It is abrupt, and weeks earlier than coastal
  sea-ice clearance.
- **bottom-fast ice**: shallow delta channels freeze to the bed, so the channel is hydraulically
  absent rather than merely capped. The 2022 records show this — Kuparuk and Canning hit
  0.00 m³/s.
- **aufeis** (icing) is a river phenomenon, regionally significant on the Sagavanirktok.

None of these have sea-ice analogues, and a thermal sea-ice formulation will get breakup timing
wrong in exactly the way the current gate does.

### What this means for the forcing set

River water temperature is needed but is **not** available year-round: USGS 00010 gives 93 days
(Sagavanirktok) and 55 days (Kuparuk) in 2022, summer only, and essentially nothing for Colville
or Canning. So river temperature has to be *modelled* from air temperature and validated against
those summer windows, rather than prescribed. This is added scope that the sea-ice framing hid.

## 1. What the model does today

There is no ice in the model. The only representation is `previousdays`, a rolling
`nbday_ice`-day (10 d) cumulative sum of water temperature returned as the second value from
`file_module.exfread`. When it goes negative, three modules take a branch:

| Module | Line | Behaviour when "iced" |
|---|---|---|
| `transport_module` | 21 | rebinds every `env=1` species' `c` to zeros |
| `biogeo_module` | 110–132 | zeroes all 11 rate arrays **and** all 12 species |
| `sed_module` | 48 | zeroes SPM |

Three consequences, in increasing order of severity:

1. **State is destroyed, not conserved.** Under ice the water column is set to zero rather than
   held. There is no under-ice DIC accumulation from respiration, no nutrient build-up, and at
   thaw everything restarts from zero and re-grows from the boundary values. The freshet carbon
   pulse — a defining feature of Arctic river systems — cannot be represented.
2. **Gas exchange is not actually gated.** `FCO2` is computed at `biogeo_module.py:93`, *above*
   the `if previousdays > 0` branch, and is absent from the list the else-branch zeroes.
   `O2_ex` **is** zeroed (line 119). Oxygen exchange is ice-gated; CO2 exchange is not.
3. **The ungated flux has the wrong sign.** Because DIC has been zeroed, the water appears
   infinitely CO2-undersaturated, so `FCO2` flips positive (uptake). Single-point check at
   `DEPTH=15 m`, `vp=1e-5 m/s`, `pCO2=420 uatm`:

   | Condition | co2s | FCO2 | Direction |
   |---|---|---|---|
   | open water, DIC=1800 | 0.0030 | −1.81e-09 | outgassing |
   | ice, DIC=0 | 0.0000 | **+2.14e-11** | **uptake** |

   Per timestep this is ~1% of the open-water magnitude, but it persists 269 days against 96
   open-water days, and `plot outputs.py` masks only exact zeros, so it survives into plots and
   annual integrals.

Item 2 is a one-line fix. Items 1 and 3 are what an ice model is actually for.

## 2. Forcing data

**Available now.** `forcing/prda2h2022.txt` is NDBC buoy PRDA2 (Prudhoe Bay), 85,987
10-minute records for 2022, columns `ATMP` (air) and `WTMP` (water). Air temperature is 6.2%
missing, range −40.0 to +20.8 °C, mean −10.33 °C. This is everything a degree-day ice model
needs and it is already in the tree — it just is not read. `main.py` has a bare
`# Air Temperature` comment where the call would go.

Derived quantities for 2022 (daily means, gap-filled by interpolation):

- 256 days below 0 °C
- autumn freeze-up (first 5-day sustained sub-zero after the summer peak): **day 267 (25 Sep)**
- accumulated freezing degree-days to 31 Dec: 1036 °C·d; over a full winter: 3895 °C·d

Stefan's law `h = alpha * sqrt(FDD)` (h in cm, FDD in °C·d):

| alpha | regime | end-Dec | full winter |
|---|---|---|---|
| 1.7 | snow-covered river | 55 cm | 106 cm |
| 2.2 | moderate snow | 71 cm | 137 cm |
| 2.7 | windswept / clear | 87 cm | **169 cm** |

**Caution on this calibration:** the 1.5–2.0 m figure usually quoted for the North Slope is
*landfast sea ice*, which is the wrong target (see §0). River ice is typically thinner in the
thalweg and, in shallow channels, grows until it reaches the bed and stops. The alpha bracket
above is therefore a starting point to be calibrated against river observations or bottom-fast
behaviour, not a validated value — and picking alpha to hit 1.7 m would be fitting to the wrong
quantity.

Freeze-up is the one place the two agree: the air-temperature estimate (day 267) and the existing
sea-water proxy (gate closes day 271) fall within 4 days, which is expected since both rivers and
the coastal zone respond to the same autumn cooling. **Breakup is where the framing fails
completely** — the proxy opens on day 173 (22 Jun) whereas the rivers peak on day 153 and are
already at 7.9 °C by day 163. That is not a 1–3 week refinement; it discards 40–45% of annual
flow on the two largest rivers.

**Needed but absent.** Snow depth on ice, which controls the effective alpha and dominates light
transmission. Options, cheapest first: (a) fold into a tuned alpha and a constant snow
attenuation — adequate for tier 1; (b) a prescribed seasonal snow climatology; (c) ERA5-Land or
SNOTEL, a new external dependency. Recommend (a) initially, revisit if under-ice light matters
to the result.

## 3. Tiered implementation

Tiers are ordered so each is independently useful and independently verifiable. Stop at any
point.

### Tier 0 — correct the existing defects (no new physics)

Two changes, both small, both removing known-wrong signal. Worth doing regardless of whether the
rest proceeds.

**0a. Gate `FCO2` consistently with `O2_ex`.** Zero `FCO2[i]` in the else-branch of
`biogeo_module`. Removes the spurious wrong-signed winter uptake. One line, no calibration.

**0b. Replace the sea-water temperature forcing.** This is the higher-value fix and is
independent of any ice model. `watertemp.csv` should be river temperature, not PRDA2 sea water.
Constraining it at ≥ 0 °C also removes the physically impossible sub-zero fresh water.

Fitted `T_river = max(0, a·T_air + b)` against the USGS 00010 records, using an n-day running
mean of air temperature to represent thermal inertia:

| River | best window | fit | R² | RMSE | vs sea substitute |
|---|---|---|---|---|---|
| Kuparuk | 10-day | `1.483·T_air + 5.52` | 0.707 | 2.45 °C | sea: RMSE 9.39, bias −8.30 → **3.5× better** |
| Sagavanirktok | 3-day | `0.352·T_air + 7.27` | 0.165 | 2.23 °C | sea: RMSE 6.89, bias −6.00 → **3.1× better** |

Read this honestly: the air-temperature model is *mediocre*, and for the Sagavanirktok it is
**poor** — R² = 0.16 means air temperature does not explain that river's summer temperature. The
Sag is mountain-fed with strong snowmelt and groundwater influence (and is the aufeis river), so
it is decoupled from local atmospheric forcing. Kuparuk behaves as expected, improving
monotonically with smoothing (R² 0.54 → 0.71 from 1-day to 10-day means); the Sag does not,
which is itself diagnostic.

The case for adopting it anyway is that it is unbiased and ~3× more accurate than a forcing that
is currently 6–8 °C cold. It is an interim correction, not a thermal model. A proper treatment
needs an equilibrium-temperature or heat-budget formulation with discharge dependence, and
Colville and Canning have no 00010 record at all to fit against — they would have to borrow a
neighbouring river's relation, with the same donor-transfer caveat as the Canning discharge.

Note 0b will *move the ice gate* as a side effect, since `previousdays` is computed from whatever
`watertemp.csv` contains. That is the intended correction — but it means 0b changes the open-water
window and the annual integral on its own, and should be run and inspected before any ice physics
is layered on top. Expect the freshet to appear.

### Tier 1 — ice as a diagnostic scalar

Introduce ice thickness and cover fraction as diagnosed quantities that modulate exchange, while
leaving the state-zeroing behaviour alone.

New forcing:
- `tools/build_airtemp.py` → `forcing/airtemp.csv`, 365 daily means from PRDA2 `ATMP`,
  same format as the existing forcings (CRLF, 2 dp, no header). Reuse the gap-fill and
  365-value contract from `fetch_discharge.py`.
- Read in `main.py` where the `# Air Temperature` comment already sits. **Must not be read
  last** — `previousdays` from the water-temperature call is what the ice gate consumes, and the
  last `exfread` wins. Insert before the `watertemp.csv` call.

New state in `variables.py` (sized `M+1`, written at odd/concentration indices to respect the
staggered grid):
- `ice_h[i]` — ice thickness [m]
- `ice_frac[i]` — areal cover fraction [0–1]
- `fdd[i]` — running freezing degree-days [°C·d]

New `ice_module.py`, called from `main.py` **before** `biogeo`:

```
ice(t, air_temp, water_temp)
  per i:
    accumulate/decay fdd from air_temp
    ice_h  = stefan_alpha * sqrt(fdd)          growth
           - melt_rate * positive degree days   decay
    ice_frac = smooth ramp in ice_h over [h_open, h_full]
```

New parameters in `config.py` (shared — ice physics is regional, not per-river):
`stefan_alpha`, `melt_rate`, `h_open`, `h_full`, `k_ice_PAR`, `albedo_ice`.

Coupling:
- `fun_module.piston_velocity`: `vp[i] *= (1.0 - ice_frac[i])`. This is the physically correct
  "lid" and makes both `O2_ex` and `FCO2` fall to zero under full cover *by construction*,
  superseding the tier-0 patch rather than duplicating it.
- `biogeo_module` light: attenuate incident `I0` by ice and snow before the existing
  Beer–Lambert water-column integral — `I0_under = I0 * (1-albedo_ice) * exp(-k_ice_PAR * ice_h)`
  applied on the `ice_frac` portion. This is what permits an under-ice spring bloom instead of a
  hard on/off.
- Output `ice_h` and `ice_frac` via `file_module.Rates` so they are plottable alongside `FCO2`.

At this tier the state-zeroing is still present, so winter DIC is still zero — but the flux is
now suppressed by cover rather than being computed from destroyed state.

### Tier 2 — conserve state under ice (the real physics change)

Replace zeroing with conservation. This is the substantive scientific change and the one that
alters the seasonal cycle.

- `transport_module`: under ice, stop rebinding `c` to zeros. Either continue transport with
  reduced cross-section (tier 3) or hold concentrations fixed.
- `biogeo_module`: under ice, keep integrating respiration/remineralisation (which continue at
  low temperature via the existing `Fhet` Arrhenius term) but with production suppressed by the
  tier-1 light attenuation, and exchange suppressed by the tier-1 piston velocity. Delete the
  else-branch state zeroing.
- `sed_module`: hold SPM rather than zeroing.

Expected consequence, and the reason to want it: DIC accumulates under ice through winter from
continued respiration with no ventilation, then vents at breakup. That produces a spring
outgassing pulse the model currently cannot generate. For an Arctic CO2 study this is likely
first-order.

Watch for: with `previousdays <= 0` no longer zeroing, any latent instability in the transport
scheme during near-zero-flow winter conditions becomes visible. Winter discharge is ~0.1–3 m³/s;
check CFL and the `hyd` convergence loop (`while rsum != 2.0`, an exact float comparison) under
those conditions before trusting long runs.

### Tier 3 — hydraulic coupling

- Ice draft reduces the flow cross-section: `D[i] -= B[i] * ice_h[i] * (rho_ice/rho_w)` where
  ice is floating.
- **Bottom-fast ice** — North Slope specific and important. Shallow delta channels freeze to the
  bed; discharge goes to zero and the channel is hydraulically absent. The 2022 records already
  show this: Kuparuk and Canning reach 0.00 m³/s in winter. Trigger when `ice_h >= DEPTH`, then
  set the local cross-section and transport to zero for those cells only — which is the *correct*
  version of the blanket zeroing the model does today, applied per grid cell rather than globally.
- This tier interacts with the geometry work in progress (SWORD widths, USGS depths). Do it after
  the geometry is real, since bottom-fast triggering depends entirely on having correct depths —
  and the shipped `DEPTH = 15 m` would never trigger it, while the gauge-observed ~1–2 m would
  trigger it readily.

### Tier 4 — optional refinements

- Solute exclusion / brine rejection: growing ice rejects salt, DIC and ALK into the residual
  water column, concentrating them. Matters for under-ice carbonate chemistry.
- Explicit breakup dynamics and ice-jam backwater.
- Snow forcing from reanalysis.
- Aufeis (icing) in the Sagavanirktok, which is regionally significant.

## 4. Validation

- **Freeze-up / breakup timing.** The USGS daily records already downloaded carry qualifier
  flags; the estimated (`A:e`) periods approximate ice-affected windows per river and give a
  per-site observational check on modelled freeze-up and breakup dates. Verify the flag semantics
  before relying on them — `A:e` means "approved, estimated", which correlates with but is not
  strictly a declaration of ice.
- **Ice thickness.** Compare peak modelled thickness against the 1.5–2.0 m regional expectation;
  Alaska DOT ice-road and lake-ice records exist for the Prudhoe area if a harder target is
  wanted.
- **Internal consistency.** Modelled freeze-up should stay near the day-267 air-temperature
  estimate and the day-271 water-temperature proxy.
- **Regression.** Tier 0 and tier 1 should leave the 96 open-water days essentially unchanged;
  any large open-water difference means a coupling was applied year-round by mistake. Keep a
  saved `FCO2.dat` from a pre-change run for diffing.

## 4b. Coupling ice to river temperature (prognostic thermal model)

Everything in §3 treats temperature as an *input*: `watertemp` is read from file and ice
is diagnosed from it. Making ice **affect** temperature requires inverting that
dependency — temperature becomes a prognostic state with a heat budget, and ice becomes
a store of latent heat that buffers it. This is a larger change than tiers 0–3 and is
described separately here.

### Why the current setup cannot express the feedback

Three structural facts, all verified in the code:

1. **`water_temp` is a scalar.** One value for all 136 grid points (`main.py:83` →
   `biogeo(..., water_temp, ...)`); there is no `water_temp[i]` anywhere. A river
   carrying a freshet down a cold channel has a strong longitudinal temperature
   gradient that the model currently cannot represent at all.
2. **It has no memory.** `T_river = max(0, 1.432·(T_air_10day + 4.0))` is instantaneous;
   the only inertia is the 10-day smoothing of air temperature. Heat capacity, and the
   latent-heat buffering that keeps a freezing river pinned at 0 °C, are absent.
3. **Ice is diagnosed downstream of temperature**, so it cannot feed back by
   construction.

### The three changes required

**(a) Temperature becomes a transported state variable. — IMPLEMENTED.** Add `T` to the species
registry in `variables.py` with `env=1`, so the existing TVD advection and dispersion
carry it exactly as they carry salinity. Consumers change from the scalar `water_temp`
to `T[i]`: 46 references in `fun_module`, 14 in `biogeo_module`.

Two things make this much cheaper than it sounds:
  - the jitted scalar helpers (`_fhet`, `_o2sat`, `_k1`, …) *already* take temperature
    as an argument, so inside the compiled loop this is a one-token change per call;
  - the species registry is already numpy, so adding a transported scalar costs nothing
    structurally.

Boundary conditions fall out neatly: the **downstream** boundary is sea temperature —
which is exactly what the deprecated `watertemp.csv` contains. The file that was wrong
as a river forcing is correct as a marine boundary condition.

*Implementation notes (done):* `T` is species 13 with `env=1`; `main.py` refreshes
`v['T']['clb']` (sea) and `v['T']['cub']` (river) each timestep since `openbound`
re-reads them per call; `piston_velocity` and the jitted `biogeo` kernel take the field
and evaluate at `Ti = c_T[i]`. `Fhet`/`Fnit` moved back inside the grid loop — they had
been hoisted out only because temperature was uniform. Verified: T develops a monotonic
mouth-to-head gradient (-1.63 → 12.0 degC) mirroring salinity (24.58 → 0) under
identical forcing.

*Known consequence:* with no heat budget yet, T is a pure mixing tracer, so the channel
interior is colder than the old uniform scalar and rates near the mouth are depressed.
This is expected at step (a) and is what step (b) fixes. Do not interpret absolute
temperatures from an (a)-only run.

**(b) A surface heat budget. — IMPLEMENTED (heat_module.py).** Replace the regression with

    rho*cp*H * dT/dt  =  Q_sw(1-albedo) + Q_lw + Q_sens + Q_lat   [+ advection, from (a)]

Data status for each term:
  - `Q_sw` — have (`solarradiation.csv`), needs an ice/snow-dependent albedo
  - `Q_sens` — have (wind 97% valid, air temp 94% valid, bulk formula)
  - `Q_lw` — derivable from air temperature with a cloud parameterisation
  - `Q_lat` — **BLOCKED: no humidity.** PRDA2's `DEWP` field is entirely empty. Options
    are (i) assume saturation over water/ice, crude but standard; (ii) source ERA5-Land,
    a new external dependency; (iii) fold it into an effective bulk coefficient tuned to
    the USGS temperature records. None is free, and latent heat is not a small term in
    an Arctic surface budget.

Stiffness is not a concern: with H ≈ 1 m, `rho*cp*H ≈ 4.2e6 J m^-2 K^-1`, so
Q ≈ 100 W m^-2 gives ~2 K/day — comfortable at `DELTI = 75 s`.

**(c) Ice as a latent-heat store — this is the actual feedback.** Four couplings, in
descending order of importance:

  1. **Latent-heat clamp.** Once T reaches freezing, further heat loss forms ice instead
     of cooling water: `dh = -Q_net·dt / (rho_ice·L_f)`, with T pinned at T_freeze. This
     is the single most important term, and the one the current `max(0, …)` clamp
     imitates without any of the physics.
  2. **Insulation.** Under ice, atmospheric exchange is replaced by conduction through
     the ice/snow slab, so under-ice water sits near 0 °C almost regardless of air
     temperature. Presently the regression drives water temperature to 0 for the whole
     winter by clamping — the right answer for the wrong reason, and it gives no
     freeze-up or breakup dynamics.
  3. **Melt buffering.** At breakup, incoming heat melts ice *before* warming water,
     delaying spring warming by days to weeks. Entirely absent now, and it directly
     affects when biogeochemistry switches on.
  4. **Albedo.** Ice and snow reflect most shortwave, so `Q_sw` into the water column
     collapses under cover and recovers sharply at breakup.

Note (1) and (3) are what make freeze-up and breakup *hysteretic* rather than a
threshold crossing — which is the physical reason a diagnostic temperature model gets
breakup timing wrong.

### What this would buy

- A **longitudinal temperature gradient**, currently absent entirely.
- **Per-river ice timing.** Today ice is regional and identical across all four sites;
  with a heat budget it would respond to each river's own depth, discharge and velocity.
- **Removal of the Sagavanirktok warm bias.** Its +2.6 °C error exists because its
  temperature is borrowed from a regression fitted on the Kuparuk. Computing temperature
  from physics removes the need to borrow — and the 93-day Sagavanirktok and 55-day
  Kuparuk USGS records stop being *fitting* targets and become *validation* targets,
  which is a much stronger position.
- **Bottom-fast ice becomes expressible** (tier 3), since ice thickness would be a real
  state rather than a degree-day estimate.

### Effort and risk, honestly

The heat budget itself is standard and small. The cost is concentrated in (a): promoting
temperature from a scalar to a transported field touches ~60 call sites across
`fun_module`, `biogeo_module`, `main` and the registry. It is mechanical, and the
bit-identity harness in `docs/performance.md` does **not** help here — this changes
results by design, so it needs physical validation instead: freeze-up/breakup dates
against the USGS ice-affected flags, and summer temperature against the 00010 records.

The honest risk is scope creep. (a) alone is a substantial refactor that changes every
temperature-dependent rate in the model; (b) has a genuine data gap; (c) introduces
hysteresis that will move the open-water window and hence the annual carbon budget.
Doing (a) first, validating that transported temperature reproduces the current results
when driven by the existing forcing at both boundaries, is the natural checkpoint before
committing to (b) and (c).

## 5. Decisions needed before implementing

1. **Do tier 0b first, on its own.** Replacing sea-water temperature with modelled river
   temperature is not an ice-model feature — it is a forcing correction that stands alone, and it
   will change the open-water window and every rate. Everything downstream should be assessed
   against a post-0b baseline, not the current one.
2. **How far up the tiers to go.** Tier 0 is unambiguous. Tier 1 is well-constrained. Tier 2
   changes the seasonal cycle and is a scientific judgement, not a bug fix.
3. **Breakup mechanism.** A purely thermal ice model will reproduce the current failure — late
   breakup — because river breakup is hydraulic. Recommend triggering breakup on discharge rise
   (freshet onset) rather than, or in addition to, ice melting through. This is the single most
   consequential design choice in the plan.
4. **Whether ice is per-site or regional.** One PRDA2 air-temperature record covers all four
   rivers, so ice *timing* would be identical across sites. But bottom-fast ice depends on local
   depth, so the ice **state** would still differ per river once tier 3 lands. Regional forcing
   with per-site response is probably the right compromise.
5. **Whether to retire `previousdays`.** Once ice is explicit, running both invites disagreement
   about when winter starts. Recommend the ice model becomes the single authority, but this
   touches the `transport` and `sed` signatures.
