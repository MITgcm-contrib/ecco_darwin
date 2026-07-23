# Arctic biogeochemistry extension

The NS-RAD (North Slope River-Aquatic-Delta Model) Arctic biogeochemistry extension:
four process groups added on top of the shipped C-GEM reaction network, motivated by
what actually controls Arctic land-to-ocean CO₂/CH₄/N₂O. Everything is **additive** and
**opt-in**: the original `TOC` pool keeps its role as labile/bioavailable DOC, and the
extension reactions run only where a site sets `ARCTIC_BGC = True` (`config.py`). With the
switch off — the default, and what the four real rivers use — `RDOC/CH4/N2O` advect as
inert passive tracers and every DIC/O₂/ALK/NO₃ update adds exactly `0.0`, so those rivers'
existing fields are **bit-identical** to the shipped network until their boundary chemistry
is constrained. Lateral loading has its own switch (`LATERAL_INFLOW`, also off by default).
Every active term is gated to open water and the ice/`react` branch like the existing
reactions. The idealized verification fixture sets both switches and tests all four groups.

Implemented in `code/biogeo_module.py` (reaction terms), `code/fun_module.py` (gas
solubility/Schmidt helpers), `code/lateral_module.py` (distributed loading),
`code/variables.py` (registry + rate arrays), `code/config.py` (parameters), and the
`sites/` boundary tables. Verified by `tools/verify_idealized.py --full`.

## New state variables

Three transported tracers join the registry (`variables.names`, all `env=1`):

| tracer | meaning | units |
|---|---|---|
| `RDOC` | refractory + chromophoric DOC | mmol C m⁻³ |
| `CH4`  | dissolved methane | mmol m⁻³ |
| `N2O`  | dissolved nitrous oxide | mmol m⁻³ |

`TOC` is retained as the **labile** DOC pool. Seven process-rate diagnostics are written
to the output alongside the existing ones: `rdoc_ox`, `photo`, `ch4_ox`, `ch4_ex`,
`n2o_prod`, `n2o_ex`, `sod`.

## 1. Refractory DOC + CDOM photomineralisation

Arctic rivers carry mostly *refractory* (old, permafrost-derived) DOC, and its
conversion to CO₂ is driven as much by **sunlight** as by microbes — photochemistry can
rival or exceed bacterial DIC production in Arctic surface waters (Cory et al. 2014,
*Science* 345:925). Two terms act on `RDOC` (both → DIC):

```
rdoc_ox = krefr · f_het(T) · RDOC/(RDOC+K_RDOC) · O2/(O2+KO2_ox)      slow microbial oxidation
photo   = PHOTO_EFF · RDOC/(RDOC+K_RDOC) · I_abs / DEPTH              CDOM photomineralisation
          I_abs = I_use − I_bottom   (PAR absorbed in the column; ice-attenuated via I_use)
```

A fraction `f_photo_lab` of the photoproduct returns as labile DOC (`TOC`), the rest goes
straight to DIC (Cory & Kling 2018). Photolysis self-gates under ice because `I_use` is
already attenuated through the slab.

Parameters (`config.py`): `krefr = 1e-5` mmol m⁻³ s⁻¹ (~60× slower than labile `kox`;
refractory turnover months–years, Hansell 2013), `K_RDOC = 100`, `PHOTO_EFF = 5e-8`
(tuned so summer surface photo-oxidation is a sizeable fraction of respiration),
`f_photo_lab = 0.2`.

## 2. Methane and nitrous oxide

Two dissolved greenhouse gases with a production/consumption term and an air–water flux.
The gas-transfer velocity is Schmidt-number rescaled off the model's existing O₂ piston
velocity, `vp_gas = vp · √(Sc_ref/Sc_gas)`, with freshwater Schmidt numbers from
Wanninkhof (2014) and equilibrium solubilities from Wiesenburg & Guinasso (1979, CH₄) and
Weiss & Price (1980, N₂O). Sign convention matches `FCO2`: **positive flux = outgassing**.

```
ch4_ox  = k_ch4ox · f_het · CH4/(CH4+K_ch4) · O2/(O2+K_ch4O2)              methanotrophy → DIC (−2 O2)
ch4_ex  = openfac · (vp_ch4/DEPTH) · (CH4 − CH4_eq(T,S,pCH4_atm))          air–water flux
n2o_prod= 0.5·(y_nit·nit + y_denit·N_denit)                               yield from nitrif + denit
n2o_ex  = openfac · (vp_n2o/DEPTH) · (N2O − N2O_eq(T,S,pN2O_atm))          air–water flux
```

Parameters: `k_ch4ox = 1.2e-6` (turnover days–weeks, Bussmann 2013), `pCH4_atm = 1.9 ppm`,
`pN2O_atm = 335 ppb` (NOAA GML 2022), `y_nit = 0.3 %` (Beaulieu et al. 2011, *PNAS*),
`y_denit = 0.5 %`. Solubilities verified against Bunsen coefficients (CH₄ ≈4.8 nmol L⁻¹,
N₂O ≈20 nmol L⁻¹ at 0 °C, present-day atmosphere).

## 3. Benthic (sediment) efflux

The shipped sediment module is SPM erosion/deposition only. In 1–2 m channels the benthos
is a large part of the metabolic budget, so a first-order **sediment-oxygen-demand**
closure adds a per-area flux (÷ depth to a volumetric source):

```
sod       = k_sod   · f_het · O2/(O2+K_sod)          / DEPTH     → DIC (aerobic benthic resp)
b_denit   = k_bdenit· f_het · NO3/(NO3+KNO3)         / DEPTH     → ALK (+1 eq/mol), NO3 loss
b_methano = k_methano·f_het · (1 − O2/(O2+K_sod))    / DEPTH     → CH4 (anaerobic)
```

Benthic terms are **not** gated by open-water fraction — sediment respiration continues
under ice, so DIC/CH₄ build up beneath the cover and vent at ice-out. `k_sod = 2.3e-4`
mmol O₂ m⁻² s⁻¹ (~20 mmol m⁻² d⁻¹, within the 5–50 river/estuary SOD range, Cai & Sayles
1996). The fuller version (a prognostic sediment-OM pool fed by settling POC) is future
work; this SOD closure is the minimal step.

## 4. Distributed lateral loading

Tundra/thermokarst/tributary inputs enter *along* the channel, not only at the upstream
boundary. `code/lateral_module.py` nudges each declared species toward the lateral
concentration in proportion to the local inflow fraction, each timestep after transport:

```
dc_i = (q_lat_i / V_i) · (c_lat − c_i) · dt ,   V_i = D_i · DELXI
```

The active site sets `LATERAL_INFLOW` [m³ s⁻¹] (spread uniformly into `variables.q_lat`
during init) and `LATERAL_CONC` (the lateral water's chemistry). This is the minimal
*solute-source* form — it adds mass but does not yet grow the downstream discharge (a
hydrodynamic change left for a later pass). **Off for the four real rivers** (unconstrained),
so their results are unchanged.

## Verification

`tools/verify_idealized.py --full` runs the idealized fixture (which enables all four
groups: RDOC/CH4/N2O boundaries, `LATERAL_INFLOW = 8 m³ s⁻¹`, a time-varying RDOC river
boundary) and asserts:

- all new tracers finite and non-negative;
- **photomineralisation** switches on in sunlit open water;
- the supersaturated river **outgasses CH₄** (`ch4_ex > 0`) once the ice clears;
- **N₂O** is finite and produced;
- **benthic SOD** efflux is positive;
- **lateral loading** raises interior CH₄ above the river boundary value;
- the **RDOC** freshet pulse reaches the upstream cell (new-species `BOUNDARY_FORCING`).

Diagnostics are in `docs/idealized_verification.pdf` (a dedicated Arctic-tracers page)
and the species registry on page 3 of `docs/ns_rad_model_schematic.pdf`.

## What is still simplified

- Photomineralisation uses a single broad-band apparent efficiency, not a wavelength-
  resolved apparent quantum yield spectrum — `PHOTO_EFF` needs site CDOM absorption data.
- The benthic term is a prescribed SOD, not a diagenetic sediment-OM pool.
- Lateral loading adds solute mass without incrementing discharge.
- CH₄/N₂O gas transfer reuses the O₂ piston velocity rescaled by Schmidt number (the
  molecular-diffusivity current term is only approximately rescaled).
- Real-river boundary values for RDOC/CH₄/N₂O are placeholder Arctic numbers pending data.
