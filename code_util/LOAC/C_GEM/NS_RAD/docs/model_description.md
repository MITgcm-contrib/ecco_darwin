# NS-RAD model description

**North Slope River-Aquatic-Delta Model** — a one-dimensional (along-channel), reactive
transport model of Arctic river–estuary–delta systems, built on the **C-GEM** estuarine
biogeochemical engine (Volta et al., 2014). NS-RAD couples channel hydrodynamics, solute
and heat transport, an aquatic carbonate system, a pelagic ecosystem/biogeochemistry
network, a prognostic river-ice model, and an optional Arctic land-to-ocean
biogeochemistry extension (dissolved CH₄/N₂O, chromophoric DOC photochemistry, benthic
exchange, and distributed lateral loading).

This document is the technical reference: the governing equations, every component and
its parameterization, and the parameter tables. For the developer/architecture guide see
[`CLAUDE.md`](../CLAUDE.md); for the catalog of what NS-RAD adds to upstream C-GEM see
[`FEATURES.md`](FEATURES.md).

## Contents

1. [Overview and governing equations](#1-overview-and-governing-equations)
2. [Model domain and grid](#2-model-domain-and-grid)
3. [Channel geometry](#3-channel-geometry)
4. [Hydrodynamics](#4-hydrodynamics)
5. [Longitudinal dispersion](#5-longitudinal-dispersion)
6. [Transport scheme](#6-transport-scheme)
7. [Water temperature and surface heat budget](#7-water-temperature-and-surface-heat-budget)
8. [Prognostic river ice](#8-prognostic-river-ice)
9. [Aquatic carbonate system](#9-aquatic-carbonate-system)
10. [Pelagic biogeochemistry](#10-pelagic-biogeochemistry)
11. [Suspended sediment](#11-suspended-sediment)
12. [Arctic biogeochemistry extension](#12-arctic-biogeochemistry-extension)
13. [Boundary conditions and forcing](#13-boundary-conditions-and-forcing)
14. [Numerical implementation](#14-numerical-implementation)
15. [Parameter tables](#15-parameter-tables)
16. [References](#16-references)

---

## 1. Overview and governing equations

NS-RAD resolves the along-channel structure of a river/estuary as a sequence of
cross-sectionally averaged cells. The prognostic variables are the water level and
velocity (hydrodynamics), and a set of transported scalars (temperature, salinity,
suspended matter, and the biogeochemical species). Each scalar $C$ obeys a
one-dimensional advection–dispersion–reaction equation,

$$\frac{\partial (A\,C)}{\partial t} \;+\; \frac{\partial (Q\,C)}{\partial x}
\;=\; \frac{\partial}{\partial x}\!\left(A\,K\,\frac{\partial C}{\partial x}\right)
\;+\; A\,\sum_j S_j(C,\dots),$$

where $x$ is distance along the channel, $A(x,t)$ the wetted cross-sectional area,
$Q=A\,U$ the volumetric discharge, $U$ the section-averaged velocity, $K$ the
longitudinal dispersion coefficient, and $\sum_j S_j$ the net of the biogeochemical,
gas-exchange, benthic, and lateral source/sink terms acting on $C$. Operator splitting is
used within each timestep: hydrodynamics → transport (advection then dispersion) → surface
heat flux → ice → reactions.

The model is **cross-sectionally averaged** and therefore does not resolve vertical
stratification or lateral structure. This is appropriate for the shallow (1–2 m),
river-dominated, microtidal North Slope systems for which it is configured, but is a
limitation in strongly stratified estuaries.

## 2. Model domain and grid

The domain spans an estuarine length $E_L$ (default 27.175 km) discretized on a uniform,
**staggered** grid with spacing $\Delta x = 200$ m, giving $M = E_L/\Delta x$ cells
(forced even; $M=136$). On the staggered grid:

- **odd indices** carry scalar concentrations and cross-sectional areas ($C$, $A$),
- **even indices** carry velocities and fluxes ($U$, $fl$).

Time integration uses a fixed step $\Delta t = 75$ s. Output is written every $T_S$ steps
(default $T_S=12$, i.e. every 900 s). The seaward end ($x=0$) is the **marine** boundary;
the landward end ($x=E_L$) is the **riverine** boundary.

## 3. Channel geometry

**Width.** NS-RAD uses a *flare + prismatic* width law (`WIDTH_MODEL="flare"`): the width
converges exponentially from the seaward width $B_{lb}$ to the prismatic upstream width
$B_{ub}$ over a short delta-flare length $L_F$, and is constant thereafter,

$$B(s) = \begin{cases}
B_{lb}\,\exp(-s/L_C), & 0 \le s < L_F, \quad L_C = -L_F / \ln(B_{ub}/B_{lb}),\\
B_{ub}, & s \ge L_F,
\end{cases}$$

with $s$ the distance from the seaward boundary. This replaces the whole-domain Savenije
exponential of standard C-GEM, which fits these rivers poorly; per-river geometry is set
from SWORD v17c node widths and USGS ADCP depth surveys. Setting `WIDTH_MODEL="expo"`
recovers the original exponential law.

**Depth and cross-section.** Depth is prescribed from at-a-station hydraulic geometry
$D = c\,Q^{f}$ evaluated at the open-water mean discharge (from the USGS surveys); the
reference cross-section is $A_0 = B\,D$, and the total wetted area $A = A_0 + B\,\eta$
evolves with the free-surface elevation $\eta$.

## 4. Hydrodynamics

The section-averaged Saint-Venant equations (continuity + momentum) are solved implicitly
each timestep for water level and velocity, with friction from a Chézy closure (Chézy
coefficient ramped from $C_z^{lb}=60$ at the mouth to $C_z^{ub}=40$ upstream). The
tridiagonal system is iterated to convergence (`hyd_module`).

**Marine boundary elevation.** The seaward water level is the sum of an astronomical tide
and a wind-driven storm surge,

$$\eta(t) = \sum_i A_i \cos\!\big(\omega_i\, t/3600 - G_i\big) \;+\; \eta_{\text{surge}}(t),$$

a multi-constituent harmonic reconstruction (amplitudes $A_i$, speeds $\omega_i$, phases
$G_i$ from the nearest NOAA CO-OPS station), plus $\eta_{\text{surge}}$, the observed
daily-mean sea-level residual. On the microtidal Beaufort coast (tidal range ~0.3 m) the
**surge is the dominant saltwater-intrusion driver**; Prudhoe Bay is used as the regional
surge proxy for all four rivers. If no constituent data is available, a single-sinusoid
fallback is used.

## 5. Longitudinal dispersion

The dispersion coefficient is computed **per timestep from local hydraulics** using the
Seo & Cheong (1998) river formula,

$$K = 5.915\left(\frac{W}{H}\right)^{0.620}\left(\frac{U}{u_*}\right)^{1.428} H\,u_*,
\qquad u_* = \sqrt{g}\,\frac{U}{C_z},$$

with $W$ the local width, $H$ the local depth, $U$ the local velocity, and shear velocity
$u_*$ taken from the model's own Chézy friction (so $U/u_* = C_z/\sqrt{g}$; no channel
slope is required, which matters because SWORD slopes at these delta reaches are
unusable). This replaces the Van der Burgh / Savenije estuary form, which collapses to zero
under realistic shallow geometry, leaving transport purely advective. $K$ is capped at
$K_{\max} = 4\,\Delta x^2/\Delta t \approx 2133$ m² s⁻¹ for numerical stability of the
dispersion scheme (reported once per run if it binds).

## 6. Transport scheme

Advection uses a **Total-Variation-Diminishing (TVD)** scheme with a flux limiter
(parity-aware on the staggered grid). Dispersion uses a **Crank–Nicolson** implicit scheme
solved with the Thomas (tridiagonal) algorithm. Boundary concentrations are imposed by an
upwind `openbound` operator that reads the marine ($c_{lb}$) or riverine ($c_{ub}$) value
depending on the sign of the boundary velocity. Every species with the transport flag
`env=1` is advected and dispersed; `pH` is diagnosed (`env=0`).

Under the prognostic ice model, transport runs **year-round**: solute state is conserved
and advected under a floating ice cover rather than zeroed each winter (the legacy gate).

## 7. Water temperature and surface heat budget

Temperature is a **transported scalar** ($T$, `env=1`), carried by the same TVD +
dispersion machinery as salinity, with time-varying boundary values (sea temperature at
the mouth, river temperature upstream). After transport, a surface heat budget warms/cools
the interior each timestep (`heat_module`),

$$\rho_w\, c_{p,w}\, H\, \frac{\partial T}{\partial t}
= Q_{sw} + Q_{lw} + Q_{sens} + Q_{lat},$$

with net shortwave $Q_{sw}$ (albedo $\alpha_w=0.06$), net longwave $Q_{lw}$ (clear-sky
Brutsaert emissivity), and turbulent sensible $Q_{sens}$ and latent $Q_{lat}$ heat from
bulk aerodynamic formulae. This makes interior temperature physically forced rather than a
mixing blend of the two boundaries. Energy removed below freezing is booked to ice
formation (§8). *Caveats:* relative humidity is observed (Deadhorse ISD); cloud cover is
not represented; shortwave is a daily mean.

## 8. Prognostic river ice

A prognostic per-cell ice thickness $h_i$ (`ice_module`) with areal fraction $f_i$
diagnosed from it. Four mechanisms:

- **Freeze-up** — heat removed below the freezing point by the surface budget is converted
  to new ice via the latent heat of fusion, $\Delta h = Q_{\text{freeze}}\Delta t /
  (\rho_{ice} L_f)$.
- **Conductive (Stefan) growth** — under an existing cover the base is at freezing and the
  slab conducts heat to a surface at $\min(T_{air},0)$: $\rho_{ice}L_f\,dh/dt = k_{ice}
  (T_{freeze}-T_{surf})/h$. Capped at the local depth (**bottom-fast**).
- **Surface melt** — warm air (sensible) plus absorbed shortwave (ice albedo $\alpha_{ice}
  =0.6$) thin the slab in spring.
- **Hydraulic (freshet) break-up** — when discharge exceeds a multiple
  ($\text{BREAKUP\_Q\_FACTOR}=3$) of the annual-mean discharge, the freshet surge clears the
  cover mechanically. This is the North-Slope-specific mechanism; a purely thermal model
  breaks up weeks too late.

**Couplings** (all no-ops where $f_i=0$): the ice slab insulates the heat budget (skips
ice-covered cells, holds the water at freezing); scales O₂/CO₂/CH₄/N₂O gas exchange by the
open-water fraction $(1-f_i)$; attenuates under-ice PAR through the slab (Beer–Lambert,
$k_{ice,PAR}$); and scales sediment erosion by $(1-f_i)$. State is **conserved** under ice,
so respiration builds up under-ice DIC/CH₄ that vents at break-up.

## 9. Aquatic carbonate system

Dissolved inorganic carbon (DIC) and total alkalinity (ALK) are prognostic; pH, aqueous
CO₂, and the air–sea CO₂ flux are diagnosed. The speciation is solved consistently in
**mol kg⁻¹**: the mmol m⁻³ state is converted via the local in-situ density $\rho$ (from a
$T,S,p$ equation of state), so it shares the mol kg⁻¹ basis of the dissociation constants
$K_0,K_1,K_2$ (carbonic acid; Millero 1995 / Weiss 1974) and $K_B$ (boric acid). Hydrogen
ion is found by the Follows et al. (2006) iteration,

$$[\mathrm{H^+}] = \tfrac12\Big[(\gamma-1)K_1 + \sqrt{(1-\gamma)^2 K_1^2 - 4K_1K_2(1-2\gamma)}\Big],
\qquad \gamma = \frac{\mathrm{DIC}}{\mathrm{Alk_C}},$$

with carbonate alkalinity $\mathrm{Alk_C} = \mathrm{ALK} - [\mathrm{B(OH)_4^-}]$, iterated
to convergence ($\le 50$ iterations). Aqueous CO₂ is
$[\mathrm{CO_2^*}] = \mathrm{DIC}\,/\,(1 + K_1/[\mathrm{H^+}] + K_1K_2/[\mathrm{H^+}]^2)$,
and pH $= -\log_{10}[\mathrm{H^+}]$. Guards return a neutral pH 7 for degenerate
(near-zero) carbonate states, a bounded transient over a few ice-out steps.

**Air–sea CO₂ flux** (sign convention: $F_{CO_2}>0$ is outgassing):

$$F_{CO_2} = (1-f_i)\,\frac{v_{p,CO_2}}{H}\big([\mathrm{CO_2^*}] - K_0\,p\mathrm{CO_2^{atm}}\big),
\qquad v_{p,CO_2} = 0.915\,v_p,$$

with $v_p$ the O₂ piston velocity (§10) and $p\mathrm{CO_2^{atm}}$ the atmospheric partial
pressure (Barrow).

## 10. Pelagic biogeochemistry

The C-GEM reaction network resolves one phytoplankton group (diatoms, `DIA`), dissolved
inorganic nutrients (`NO3`, `NH4`, `PO4`, `dSi`), dissolved oxygen (`O2`), a labile
organic-matter pool (`TOC`), and the carbonate system (`DIC`, `ALK`, diagnosed `pH`). All
temperature-dependent rates use per-cell Arrhenius factors $f_{het}=2.75^{(T-5)/10}$
(heterotrophic) and $f_{nit}=5^{(T-5)/10}$ (nitrification).

### 10.1 Primary production

Light attenuation combines background and suspended-matter terms,
$K_D = K_{D1} + K_{D2}(1000\,\mathrm{SPM}+100)$, giving a Beer–Lambert profile
$I(z)=I_{use}\exp(-K_D z)$ where $I_{use}$ is the surface PAR (equal to incident PAR in
open water, ice-attenuated under a cover). A Platt-style photosynthesis–irradiance curve
$P = P_b^{\max}\big(1-e^{-I/E_k}\big)$ is integrated analytically over the water column
(a difference of exponential-integral terms), with a photoacclimated $E_k$ from a variable
chlorophyll:carbon ratio. Nutrient limitation is multiplicative (Liebig via product form),

$$L_{nut} = \frac{dSi}{dSi+K_{dSi}}\cdot\frac{NO_3+NH_4}{NO_3+NH_4+K_N}\cdot\frac{PO_4}{PO_4+K_{PO4}}.$$

Gross and net production, and phytoplankton mortality, are

$$\mathrm{GPP} = P_b^{\max}\,\mathrm{DIA}\,L_{nut}\,\mathcal{I},\qquad
\mathrm{NPP} = \frac{\mathrm{GPP}}{H}(1-k_{exc})(1-k_{grw}) - k_{maint}\,\mathrm{DIA},\qquad
\mathrm{PD} = k_{mort}\,\mathrm{DIA},$$

where $\mathcal I$ is the depth-integrated light-limited growth. NPP is partitioned into
nitrate- and ammonium-based uptake ($\mathrm{NPP_{NO3}}$, $\mathrm{NPP_{NH4}}$) by an
ammonium-preference factor $f_{NH4}=NH_4/(10+NH_4)$; silica uptake is
$\mathrm{Si_{cons}}=r_{Si}\,\mathrm{NPP}$.

### 10.2 Organic-matter degradation and nitrogen cycling

Heterotrophic degradation of the labile pool proceeds aerobically and, at low oxygen, by
denitrification; ammonium is nitrified:

$$\begin{aligned}
\text{aerobic:}\quad & R_{ox} = k_{ox}\,f_{het}\,\frac{TOC}{TOC+K_{TOC}}\,\frac{O_2}{O_2+K_{O2}^{ox}},\\
\text{denitrification:}\quad & R_{den} = k_{den}\,f_{het}\,\frac{TOC}{TOC+K_{TOC}}\,\frac{K_{inO2}}{O_2+K_{inO2}}\,\frac{NO_3}{NO_3+K_{NO3}},\\
\text{nitrification:}\quad & R_{nit} = k_{nit}\,f_{nit}\,\frac{O_2}{O_2+K_{O2}^{nit}}\,\frac{NH_4}{NH_4+K_{NH4}}.
\end{aligned}$$

### 10.3 Gas exchange (piston velocity)

The gas-transfer velocity combines a current-shear (surface-renewal) and a wind term,

$$v_p = \underbrace{\sqrt{|U|\,D_{O_2}/H}}_{\text{shear}}
\;+\; \underbrace{\tfrac{1}{3.6\times10^5}\,0.31\,U_w^2\,(Sc/660)^{-1/2}}_{\text{wind}},$$

with $D_{O_2}(T)$ the molecular diffusivity and $Sc(T,S)$ the Schmidt number. Oxygen
exchange is $O_{2,ex} = (1-f_i)\,(v_p/H)\,(O_2^{sat}-O_2)$; CO₂ uses $0.915\,v_p$ (§9); the
Arctic gases rescale $v_p$ by their Schmidt numbers (§12).

### 10.4 Reaction stoichiometry

Applied to each cell per timestep (Redfield ratios $r_N=16/106$, $r_P=1/106$,
$r_{Si}=16/80$; alkalinity changes follow Wolf-Gladrow et al. 2007):

$$\begin{aligned}
\Delta\mathrm{DIA} &= (\mathrm{NPP}-\mathrm{PD})\,\Delta t, &
\Delta\mathrm{dSi} &= -\mathrm{Si_{cons}}\,\Delta t,\\
\Delta\mathrm{NO_3} &= \big(-\tfrac{94.4}{106}R_{den}+R_{nit}-r_N\mathrm{NPP_{NO3}}\big)\Delta t, &
\Delta\mathrm{NH_4} &= \big(r_N(R_{ox}-\mathrm{NPP_{NH4}})-R_{nit}\big)\Delta t,\\
\Delta\mathrm{PO_4} &= r_P(R_{ox}+R_{den}-\mathrm{NPP})\,\Delta t, &
\Delta\mathrm{TOC} &= (-R_{ox}-R_{den}+\mathrm{PD})\,\Delta t,\\
\Delta\mathrm{O_2} &= \big(-R_{ox}+\mathrm{NPP_{NH4}}+\tfrac{138}{106}\mathrm{NPP_{NO3}}-2R_{nit}+O_{2,ex}\big)\Delta t, &
\Delta\mathrm{DIC} &= \big(-F_{CO_2}-\mathrm{NPP}+R_{ox}+R_{den}\big)\Delta t,\\
\Delta\mathrm{ALK} &= \big(\tfrac{15}{106}R_{ox}+\tfrac{93.4}{106}R_{den}-2R_{nit}-\tfrac{15}{106}\mathrm{NPP_{NH4}}+\tfrac{17}{106}\mathrm{NPP_{NO3}}\big)\Delta t. &&
\end{aligned}$$

Net ecosystem metabolism is reported as $\mathrm{NEM}=\mathrm{NPP}-R_{ox}-R_{den}$ (minus
the Arctic-extension respiration terms when enabled).

## 11. Suspended sediment

Suspended particulate matter (`SPM`) is eroded and deposited according to the bed shear
stress $\tau_b = \rho_w g U^2 / C_z^2$: erosion where $\tau_b>\tau_{ero}$ (rate $M_{ero}$),
deposition where $\tau_b<\tau_{dep}$ with a concentration-dependent settling velocity
$w_s = w_{\max}\,\mathrm{SPM}/(\mathrm{SPM}+k_{ISS})$ (Clark et al. 2020, 2022). Erosion is
scaled by open-water fraction under ice. SPM feeds back on light attenuation (§10.1).

## 12. Arctic biogeochemistry extension

Four process groups added on top of the shipped network, motivated by what controls Arctic
land-to-ocean CO₂/CH₄/N₂O. **Opt-in** via `config.ARCTIC_BGC` (default off): when off, the
three added tracers advect as inert passive tracers and every term below is exactly zero,
so real-river carbon fields are bit-identical to the shipped C-GEM network. Detailed
write-up and citations: [`arctic_biogeochemistry.md`](arctic_biogeochemistry.md).

### 12.1 Refractory DOC and CDOM photomineralisation

A second DOC pool `RDOC` (refractory + chromophoric) is oxidized slowly by microbes and
photochemically by sunlight, both to DIC (`TOC` remains the labile pool):

$$R_{RDOC} = k_{refr}\,f_{het}\,\frac{RDOC}{RDOC+K_{RDOC}}\,\frac{O_2}{O_2+K_{O2}^{ox}},
\qquad
R_{photo} = \Phi_{photo}\,\frac{RDOC}{RDOC+K_{RDOC}}\,\frac{I_{abs}}{H},$$

where $I_{abs}=I_{use}-I_{bottom}$ is the PAR absorbed in the column and $\Phi_{photo}$ an
apparent photomineralisation efficiency (Cory et al. 2014). A fraction $f_{lab}$ of the
photoproduct returns as labile DOC, the rest to DIC. Photolysis self-gates under ice
because $I_{use}$ is already ice-attenuated.

### 12.2 Methane and nitrous oxide

`CH4` and `N2O` are transported gases with production/consumption and air–water exchange
(sign convention: flux $>0$ is outgassing). Schmidt numbers use Wanninkhof (2014) and
equilibrium solubilities Wiesenburg & Guinasso (1979, CH₄) / Weiss & Price (1980, N₂O):

$$\begin{aligned}
R_{CH_4,ox} &= k_{ch4}\,f_{het}\,\frac{CH_4}{CH_4+K_{ch4}}\,\frac{O_2}{O_2+K_{ch4}^{O2}}, &
F_{CH_4} &= (1-f_i)\,\frac{v_{p,CH_4}}{H}\big(CH_4 - CH_4^{eq}\big),\\
P_{N_2O} &= \tfrac12\big(y_{nit}R_{nit} + y_{den}N_{den}\big), &
F_{N_2O} &= (1-f_i)\,\frac{v_{p,N_2O}}{H}\big(N_2O - N_2O^{eq}\big),
\end{aligned}$$

with $N_{den}=\tfrac{94.4}{106}R_{den}$ the N reduced by denitrification, and
$v_{p,gas}=v_p\sqrt{Sc_{ref}/Sc_{gas}}$. Methanotrophy consumes 2 O₂ per CH₄ and produces
DIC.

### 12.3 Benthic (sediment) exchange

A first-order sediment-oxygen-demand closure (per-area flux ÷ depth), returning DIC,
alkalinity, nitrate loss, and CH₄ from the bed. Benthic terms run **under ice** (sediment
respiration continues), so they feed the ice-out vent:

$$\mathrm{SOD} = \frac{k_{sod}\,f_{het}}{H}\,\frac{O_2}{O_2+K_{sod}},\quad
B_{den} = \frac{k_{bden}\,f_{het}}{H}\,\frac{NO_3}{NO_3+K_{NO3}},\quad
B_{CH_4} = \frac{k_{meth}\,f_{het}}{H}\Big(1-\frac{O_2}{O_2+K_{sod}}\Big).$$

SOD adds DIC (aerobic benthic respiration); benthic denitrification adds alkalinity
(~1 eq per mol NO₃ reduced) and removes nitrate; methanogenesis sources CH₄.

### 12.4 Distributed lateral loading

Tundra/thermokarst/tributary inputs entering *along* the channel (not only at the upstream
boundary) are represented as a mixing source applied after transport,

$$\Delta C_i = \frac{q_{lat,i}}{V_i}\,(c_{lat}-C_i)\,\Delta t, \qquad V_i = A_i\,\Delta x,$$

where the total lateral inflow $Q_{lat}$ is spread uniformly over the domain and $c_{lat}$
is the lateral water's chemistry (per species). This is the minimal *solute-source* form —
it adds mass but does not yet grow the downstream discharge (a hydrodynamic change left for
future work). Off by default (`LATERAL_INFLOW=0`).

## 13. Boundary conditions and forcing

**Boundaries.** Each transported species has a marine ($c_{lb}$) and riverine ($c_{ub}$)
boundary value from the active site's `BOUNDARIES` table. Temperature boundaries are
refreshed every timestep from the sea/river temperature forcings. Any species can be made
**time-varying** at either end via `BOUNDARY_FORCING` (a per-site map of species → 365-day
CSV series), refreshed each timestep.

**Forcing (2022 climatology, repeated).** Per-river discharge (USGS gauges; Canning
reconstructed); river temperature (air-temperature regression + USGS 00010); shared
meteorology (NDBC PRDA2 wind/solar/air, Deadhorse ISD humidity, Barrow pCO₂); per-river
harmonic tides (NOAA CO-OPS) + regional storm surge (Prudhoe Bay); boundary chemistry
(Arctic LTER for Kuparuk, WQP grabs for Colville/Sagavanirktok). See
[`CLAUDE.md`](../CLAUDE.md) for the full provenance and per-river caveats.

## 14. Numerical implementation

- **Staggered grid**, odd = scalars/areas, even = velocities/fluxes (§2).
- **Hydrodynamics**: implicit, iterated tridiagonal solve to convergence.
- **Advection**: explicit TVD with a flux limiter (parity-aware).
- **Dispersion**: implicit Crank–Nicolson, Thomas tridiagonal solve; capped at $K_{\max}$.
- **Timestep** $\Delta t=75$ s; **cell size** $\Delta x=200$ m.
- **Global mutable state**: all arrays are allocated once (`variables.py`) and mutated in
  place; there is no state object.
- **Performance**: the per-cell hot loops (hydrodynamic kernels, density stack, transport
  schemes, biogeochemistry, sediment, heat, ice, carbonate) are `@njit`-compiled with
  Numba, held bit-identical to the pure-Python reference (see
  [`performance.md`](performance.md)).
- **Output**: compact NetCDF by default (`output.nc`, time × distance), with a legacy
  tab-separated `.dat` path.

## 15. Parameter tables

### 15.1 Domain and numerics

| Symbol | Parameter | Value | Units |
|---|---|---|---|
| $E_L$ | estuarine length | 27 175 | m |
| $\Delta x$ | cell size | 200 | m |
| $\Delta t$ | timestep | 75 | s |
| $M$ | grid cells | 136 | – |
| $K_{\max}$ | dispersion cap | 2133 | m² s⁻¹ |
| — | pH iterations | 50 | – |

### 15.2 Hydrodynamics and sediment

| Symbol | Parameter | Value | Units | Source |
|---|---|---|---|---|
| $C_z^{lb},C_z^{ub}$ | Chézy (mouth, upstream) | 60, 40 | m^½ s⁻¹ | Panchenko & Alabyan 2022 |
| $M_{ero}$ | erosion coefficient | 5.79×10⁻⁵ | mg m⁻² s⁻¹ | Clark et al. 2020 |
| $\tau_{ero}$ | erosion shear stress | 0.005 | N m⁻² | Clark et al. 2022 |
| $\tau_{dep}$ | deposition shear (mouth, up) | 0.4, 1.0 | N m⁻² | — |
| $w_{\max}$ | max settling velocity | 2.31×10⁻⁵ | m s⁻¹ | Clark et al. 2022 |
| $k_{ISS}$ | settling half-saturation | 0.051 | g L⁻¹ | Clark et al. 2022 |

### 15.3 Phytoplankton and light

| Symbol | Parameter | Value | Units | Source |
|---|---|---|---|---|
| $P_b^{\max}$ | max photosynthetic rate | 1.39×10⁻⁵ | s⁻¹ | Le Fouest et al. 2013 |
| $\alpha$ | photosynthetic efficiency | 2.0 | mgC (mgChl)⁻¹ (E m⁻²d⁻¹)⁻¹ | Le Fouest et al. 2013 |
| $E_K$ | photoacclimation | 8.0 | E m⁻² d⁻¹ | Le Fouest et al. 2013 |
| $\theta_{min}$ | min Chl:C | 0.0125 | g g⁻¹ | Le Fouest et al. 2013 |
| $k_{maint}$ | maintenance rate | 4.6×10⁻⁷ | s⁻¹ | — |
| $k_{mort}$ | mortality | 1.16×10⁻⁷ | s⁻¹ | Le Fouest et al. 2013 |
| $k_{exc}$ | excretion | 0.05 | – | — |
| $k_{grw}$ | growth cost | 0.29 | – | — |
| $K_{D1}$ | background attenuation | 1.3 | m⁻¹ | — |
| $K_{D2}$ | SPM attenuation | 0.06 | mg⁻¹ m⁻¹ | — |

### 15.4 Nutrient/oxygen half-saturations

| Symbol | Parameter | Value | Units |
|---|---|---|---|
| $K_{dSi}$ | silica | 1.07 | µM Si |
| $K_{PO4}$ | phosphate | 0.2 | µM P |
| $K_{NH4}$ | ammonium | 0.5 | µM N |
| $K_{NO3}$ | nitrate | 1.0 | µM N |
| $K_N$ | dissolved N | 1.13 | µM N |
| $K_{TOC}$ | organic matter | 186.25 | µM C |
| $K_{O2}^{ox}$ | O₂ (oxic degradation) | 31.0 | µM O₂ |
| $K_{O2}^{nit}$ | O₂ (nitrification) | 51.25 | µM O₂ |
| $K_{inO2}$ | O₂ (denitrif. inhibition) | 33.0 | µM O₂ |

### 15.5 Reaction rates and Redfield ratios

| Symbol | Parameter | Value | Units | Source |
|---|---|---|---|---|
| $k_{ox}$ | aerobic degradation | 6.08×10⁻⁴ | µM C s⁻¹ | — |
| $k_{den}$ | denitrification | 5.05×10⁻⁴ | µM C s⁻¹ | — |
| $k_{nit}$ | nitrification | 3.47×10⁻⁶ | µM C s⁻¹ | Le Fouest et al. 2013 |
| $r_N$ | N:C Redfield | 16/106 | mol/mol | Redfield |
| $r_P$ | P:C Redfield | 1/106 | mol/mol | Redfield |
| $r_{Si}$ | Si:C | 16/80 | mol/mol | — |

### 15.6 Arctic biogeochemistry extension

| Symbol | Parameter | Value | Units | Source |
|---|---|---|---|---|
| $k_{refr}$ | RDOC oxidation max rate | 1.0×10⁻⁵ | mmol m⁻³ s⁻¹ | Hansell 2013 |
| $K_{RDOC}$ | RDOC half-saturation | 100 | mmol m⁻³ | — |
| $\Phi_{photo}$ | photomineralisation efficiency | 5.0×10⁻⁸ | mmol m⁻² s⁻¹ / (µE m⁻² s⁻¹) | Cory et al. 2014 (tuned) |
| $f_{lab}$ | photoproduct labile fraction | 0.2 | – | Cory & Kling 2018 |
| $k_{ch4}$ | methanotrophy max rate | 1.2×10⁻⁶ | mmol m⁻³ s⁻¹ | Bussmann 2013 |
| $K_{ch4}$ | CH₄ half-saturation | 0.5 | mmol m⁻³ | — |
| $K_{ch4}^{O2}$ | O₂ half-sat (methanotrophy) | 5.0 | mmol m⁻³ | — |
| $p\mathrm{CH_4^{atm}}$ | atmospheric CH₄ | 1.9×10⁻⁶ | atm | NOAA GML 2022 |
| $y_{nit}$ | N₂O yield (nitrification) | 0.003 | mol/mol | Beaulieu et al. 2011 |
| $y_{den}$ | N₂O yield (denitrification) | 0.005 | mol/mol | — |
| $p\mathrm{N_2O^{atm}}$ | atmospheric N₂O | 3.35×10⁻⁷ | atm | NOAA GML 2022 |
| $k_{sod}$ | sediment O₂ demand max rate | 2.3×10⁻⁴ | mmol O₂ m⁻² s⁻¹ | Cai & Sayles 1996 |
| $K_{sod}$ | O₂ half-sat (SOD) | 10 | mmol m⁻³ | — |
| $k_{bden}$ | benthic denitrification | 5.0×10⁻⁵ | mmol N m⁻² s⁻¹ | — |
| $k_{meth}$ | benthic methanogenesis | 2.0×10⁻⁸ | mmol CH₄ m⁻² s⁻¹ | — |

## 16. References

- **Beaulieu, J. J., et al. (2011).** Nitrous oxide emission from denitrification in stream
  and river networks. *PNAS* 108, 214–219.
- **Bussmann, I. (2013).** Distribution of methane in the Lena Delta and Buor-Khaya Bay,
  Russia. *Biogeosciences* 10, 4641–4652.
- **Cai, W.-J., & Sayles, F. L. (1996).** Oxygen penetration depths and fluxes in marine
  sediments. *Mar. Chem.* 52, 123–131.
- **Clark, J. B., et al. (2020, 2022).** Sediment and biogeochemical dynamics of the
  Mackenzie/Arctic river systems. *JGR Biogeosciences.*
- **Cory, R. M., et al. (2014).** Sunlight controls water column processing of carbon in
  arctic fresh waters. *Science* 345, 925–928.
- **Cory, R. M., & Kling, G. W. (2018).** Interactions between sunlight and microorganisms
  influence dissolved organic matter degradation. *Limnol. Oceanogr. Lett.* 3, 102–116.
- **Follows, M. J., et al. (2006).** On the solution of the carbonate chemistry system in
  ocean biogeochemistry models. *Ocean Modelling* 12, 290–301.
- **Hansell, D. A. (2013).** Recalcitrant dissolved organic carbon fractions. *Ann. Rev.
  Mar. Sci.* 5, 421–445.
- **Le Fouest, V., et al. (2013).** The fate of riverine nutrients on Arctic shelves.
  *Biogeosciences* 10, 4785–4800.
- **Millero, F. J. (1995).** Thermodynamics of the carbon dioxide system in the oceans.
  *Geochim. Cosmochim. Acta* 59, 661–677.
- **Regnier, P., et al. (2002, 2013).** Reactive-transport modelling of estuarine
  biogeochemistry / land–ocean carbon.
- **Seo, I. W., & Cheong, T. S. (1998).** Predicting longitudinal dispersion coefficient in
  natural streams. *J. Hydraul. Eng.* 124, 25–32.
- **Volta, C., et al. (2014).** C-GEM (v 1.0): a new, cost-efficient biogeochemical model
  for estuaries. *Geosci. Model Dev.* 7, 1271–1295.
- **Wanninkhof, R. (2014).** Relationship between wind speed and gas exchange over the ocean
  revisited. *Limnol. Oceanogr. Methods* 12, 351–362.
- **Weiss, R. F., & Price, B. A. (1980).** Nitrous oxide solubility in water and seawater.
  *Mar. Chem.* 8, 347–359.
- **Wiesenburg, D. A., & Guinasso, N. L. (1979).** Equilibrium solubilities of methane in
  water and seawater. *J. Chem. Eng. Data* 24, 356–360.
- **Wolf-Gladrow, D. A., et al. (2007).** Total alkalinity: the explicit conservative
  expression and its application to biogeochemical processes. *Mar. Chem.* 106, 287–300.

---

*NS-RAD is built on the C-GEM engine (Volta et al. 2014). This document describes the
North-Slope configuration and its extensions; see `FEATURES.md` for the change log relative
to upstream C-GEM.*
