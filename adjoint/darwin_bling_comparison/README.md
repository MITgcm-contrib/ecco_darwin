# darwin_bling_comparison

A pilot study comparing **Darwin's forward finite-difference sensitivity**
against **BLING's adjoint sensitivity** of the same SOCAT surface-pCO2 cost
function, at matched grid points, on identical grid/physics/forcing. The
two ecosystem models differ substantially (Darwin: 10 explicit,
trait-based competing plankton types across 7 functional groups; BLING:
2-3 implicit biomass pools with Liebig single-limiting-nutrient growth),
so this is not a strict adjoint-correctness gradient check (see "Relation
to Kay's methodology" below) -- it's a test of how far two independently-
built ecosystem models agree on the sensitivity of ocean pCO2 to
perturbations in their shared state variables.

**Bottom line**: after fixing two real bugs (a 1000x unit mismatch and an
initial-condition provenance mismatch, both described below), Darwin and
BLING agree closely on the sensitivity of pCO2 to DIC and alkalinity
(direct carbonate-chemistry pathway, same equations in both models,
verified bit-identical) across 38 globally-distributed points: slope
1.11-1.13, correlation 0.999. They diverge substantially on the
biology-mediated tracers (NO3, O2, PO4, FeT) for reasons that are mostly
mechanistically understood (see "Results" below) but not fully resolved
for iron.

## Get the code

```sh
git clone https://github.com/MITgcm-contrib/ecco_darwin.git
git clone https://github.com/MITgcm/MITgcm.git
cd MITgcm
git clone https://github.com/darwinproject/darwin3.git
```
You'll need both `../global_oce_biogeo_darwin` and
`../global_oce_biogeo_bling_SOCAT` built and runnable first (see their
own READMEs), including `scripts/match_ics_to_bling.py` from the Darwin
side.

## Relation to Kay's methodology

This pilot's finite-difference design mirrors an internal reference
methodology (a GCHP CO2 adjoint-validation writeup, not included in this
repo) for validating an adjoint gradient against forward finite
differences:

- **Centered (two-sided) FD**: `(J+ - J-) / (2*delta)`, not one-sided --
  cancels the leading nonlinear term.
- **Perturb the IC at a single cell** (or, in the 38-point extension,
  many well-separated cells simultaneously), nothing else.
- **Compare against an actual cost-function value**, not a raw field
  value: `J = sum(weight * (model - obs)^2)`, matching pkg/profiles'
  own cost formula exactly (verified against `pkg/profiles/cost_profiles.F`
  and `pkg/cost/cost_final.F` -- no missing factor of 1/2 or 2 anywhere
  in the chain).
- **Relative (not absolute) perturbation size**: `delta = eta * c0_m`
  (eta=5%), normalized by the actual perturbation applied at each point
  -- not a shared constant -- avoiding both a common normalization bug
  and negative-concentration blowups for tracers spanning a wide dynamic
  range (NO3 alone varies ~600x across the 38 validation points).
  See `scripts/make_multipoint_perturbed_ic.py`.
  An earlier version of this pilot used fixed *absolute* eps values
  (`scripts/make_perturbed_ic.py`, `run_sensitivity_sweep_fields.sh`) at
  a single point -- fine there, but doesn't generalize across points with
  very different baseline concentrations, hence the switch to relative
  perturbations for the 38-point extension.
- **Nonlinearity/asymmetry diagnostic**: `(J+-J0) - (J0-J-)` using a true
  unperturbed baseline, computed in `compute_multipoint_comparison.py`.
  Caveat found during this pilot: because `J` is itself quadratic in the
  model-obs residual, this ratio can look large even when the underlying
  FD gradient estimate is fine, specifically at points where the
  model-obs residual happens to be small (the ratio's denominator
  shrinks, an expected artifact of squaring, not evidence of a bad FD
  estimate) -- verified by checking that high-asymmetry-ratio points
  still showed excellent Darwin/BLING agreement.

**Where this pilot's setup differs from a strict adjoint-validation
check**, and why the results can't be held to a strict slope=1/r=1 bar
the way a same-model FD-vs-adjoint check can:

- **Cross-model, not self-consistency.** A same-model check (GCHP-vs-GCHP)
  validates that the adjoint *code* is correct. This pilot compares two
  *different* models, so a mismatch can mean a real scientific difference
  between the models rather than a bug in either.
- **BLING's adjoint field is in control-preconditioned units, not
  physical units.** MITgcm's `ctrl_map_genarr3d.F` divides the raw
  control by `sqrt(ctrl weight)` before adding it to the physical field
  (`fld = fld_prior + xx_gen/sqrt(w)`), so the `adxx_*` output is
  `dJ/d(xx_gen)`, not `dJ/d(field)`. Converting requires multiplying by
  `sqrt(w)` from the field's own `wt_*.bin` weight file -- see
  `compute_multipoint_comparison.py`. This preconditioning doesn't exist
  in every adjoint system (a reference GCHP CO2-adjoint setup, for
  example, reports `SpeciesAdj_CO2` directly in physical mol/mol with no
  conversion needed), so it's easy to miss if porting this methodology
  elsewhere.

## How to reproduce

1. **Match Darwin's ICs to BLING's** (`global_oce_biogeo_darwin/scripts/match_ics_to_bling.py`)
   -- see "Bug #2" below for why this matters.
2. **Single-point pilot** (optional, superseded by step 3, kept for
   reference): `scripts/run_sensitivity_sweep_fields.sh` (6 tracers x
   +-eps at grid point j=40,i=97) -> `scripts/compute_cost_fd.py` ->
   `scripts/plot_fd_vs_adjoint.py`.
3. **38-point extension** (the real result):
   `scripts/setup_and_launch_multipoint.sh` (6 tracers x +-5% relative
   perturbation, simultaneously at 38 well-separated global points, plus
   1 unperturbed baseline = 13 concurrent MITgcm runs) ->
   `scripts/compute_multipoint_comparison.py` ->
   `scripts/plot_multipoint_comparison.py`.
4. **IC sanity check**: `scripts/plot_initial_conditions.py` (maps +
   global-mean vertical profiles, Darwin vs BLING, for all 6 shared
   tracers) -- run before/after step 1 to see the fix's effect.

The 38 validation points were chosen from SOCAT's own observation grid
cells (every SOCAT month-01 observation lands exactly on a model grid
cell center -- confirmed 2759/2759 exact matches -- so pkg/profiles'
bilinear spatial interpolation collapses to a direct single-cell lookup,
no interpolation machinery needed for the Darwin-side cost calculation),
with a minimum Chebyshev separation of 10 grid cells so that
simultaneously perturbing all 38 in one run doesn't let one point's
advective "footprint" contaminate another's local cost contribution
(BLING's own adjoint field shows sensitivity spreading several cells via
advection over a 30-day integration -- confirmed by inspecting the raw
`adxx_ptr1` field directly).

## Bugs found and fixed along the way

**Bug #1 -- 1000x unit mismatch (dominant cause of an early ~9,430x,
wrong-signed mismatch).** BLING's DIC initial-condition file is in
mol/m^3; Darwin's is in mmol/m^3. Confirmed directly: BLING's IC at a
test point was 1.99 mol/m^3, Darwin's was 2027 mmol/m^3 -- the same
physical value, 1000x apart in stored units. All unit conversions in
`compute_cost_fd.py`/`compute_multipoint_comparison.py` explicitly
account for this (`* 1000.0` comments mark every conversion point).

**Bug #2 -- initial-condition provenance mismatch (the remaining ~2x,
sign-correct-but-off discrepancy, after fixing Bug #1).** Darwin's
DIC/ALK/O2/NO3/PO4/Fe initial conditions came from a completely
different, never-cross-validated source than BLING's: Darwin's were
regridded from a prior Darwin LLC270 ecosystem run
(`global_oce_biogeo_darwin/scripts/regrid_darwin_ics.py`); BLING's ship
with the stock MITgcm BLING tutorial. At one test point these disagreed
by 7.1% in ALK and 1.9% in DIC already at t=0 (and this barely grew over
a 30-day integration: 7.07%->7.19%, 1.85%->2.21%) -- i.e. the mismatch
was a data-provenance problem, not divergent model physics. Verified
mechanistically: both models' carbonate-chemistry solvers
(`DARWIN_CALC_PCO2_APPROX` / `CALC_PCO2_APPROX`) were confirmed
line-for-line identical (same OCMIP2/Mehrbach-Millero lineage, same
inputs `t,s,DIC,PO4,SiO2,ALK` + equilibrium constants, bit-identical
output for identical input), so feeding both the same DIC/ALK state
reproduces d(pCO2)/dDIC to within ~15% -- matching the empirical
~2x-before/~1.1x-after result almost exactly.
`global_oce_biogeo_darwin/scripts/match_ics_to_bling.py` fixes this by
overwriting Darwin's IC for the 6 shared tracers with BLING's own
(unit-converted), and the full 38-point sweep was rerun after the fix.

## Results

| Tracer | slope (BLING/Darwin) | correlation | pathway |
|---|---|---|---|
| DIC | 1.11 | 0.999 | direct carbonate chemistry |
| ALK | 1.13 | 0.999 | direct carbonate chemistry |
| NO3 | 0.20 (0.60 at NO3-limited points) | 0.60 | biology-mediated |
| O2 | -- (Darwin's FD is exactly 0 at all 38 points) | -- | no direct pathway |
| PO4 | -0.12 | -0.92 | biology-mediated |
| FeT | 0.09 | 0.47 | biology-mediated |

**DIC and ALK converge cleanly** -- both go through the verified-identical
carbonate-chemistry solver, so once the IC mismatch (Bug #2) was fixed,
there was no remaining structural reason for them to disagree.

**O2 is exactly zero in Darwin's forward sensitivity, structurally, not
just "not enough time."** `DARWIN_CALC_PCO2_APPROX`'s full input list
(`t, s, DIC, PO4, SiO2, ALK` + equilibrium constants) has no O2 term at
all -- standard OCMIP2 carbonate chemistry simply doesn't take dissolved
oxygen as an input, in either model. Any real effect would have to be
indirect (O2 -> respiration/remineralization -> DIC/ALK), which Darwin's
slower multi-species community dynamics apparently didn't complete in 30
days; BLING's tiny nonzero adjoint value (0.14) suggests its simpler,
faster-responding biology let a trace of that indirect signal through.

**NO3's divergence is explained by BLING's Liebig (single-limiting-
nutrient) growth law.** `bling_bio_nitrogen.F` computes growth as
`min(NO3_lim, PO4_lim, Fe_lim)` -- only the most-limiting nutrient
controls growth at any given point. Cross-referencing against BLING's own
diagnosed limiting-nutrient field (`BLGNLIM`/`BLGPLIM`/`BLGFELIM`, added
via a one-off diagnostic run, `global_oce_biogeo_bling_SOCAT` with those
fields enabled in `data.diagnostics`) confirmed this directly: at the 24
of 38 points where NO3 actually is the binding nutrient in BLING,
Darwin/BLING correlation is 0.60; at the 14 points where Fe is binding
instead, it *inverts* to -0.49. Perturbing NO3 where it isn't locally
limiting does essentially nothing in BLING's growth law, so Darwin's
different (non-gated, community-competition-based) response shows up as
pure disagreement there.

**PO4 never binds at any of the 38 points** (`PO4_lim` is never the
`min()` of the three) -- so BLING's own PO4 sensitivity is close to a
noise floor everywhere in this sample, while Darwin's structurally
different community model still produces a real, if small, response.
Comparing a real signal against near-pure noise is a plausible
explanation for the negative correlation, more than a genuine sign
inversion in the underlying biology.

**FeT's divergence does NOT fit the same Liebig-gating story** -- unlike
NO3, correlation is essentially the same whether Fe is locally binding
(0.47) or not (0.50). FeT sensitivities are also enormous and wildly
scattered in absolute terms (10^3 to 10^7, with several points differing
in sign), consistent with iron being the least-constrained, most
model-specific nutrient cycle in ocean biogeochemistry generally
(scavenging kinetics, ligand complexation, dust-deposition sensitivity
all vary a great deal between independently-built Fe cycle
formulations). **This is the main open item**: pinning it down further
would mean comparing both models' Fe-cycle source terms directly (the
way the carbonate-chemistry solvers were compared for Bug #2), which
hasn't been done.

## Directory layout

```
scripts/    All analysis/setup scripts (see docstrings for env-var config)
figures/    Comparison plots (see below) + IC sanity-check plots
results/    Small result files: multipoint_comparison_results.csv,
            cost_fd_results.json, sensitivity_sweep_results.txt,
            pilot_adxx_fields.npz (BLING's raw adxx_* fields, extracted
            once via MITgcmutils.rdmds, so the 38-point comparison never
            needs to re-run BLING's adjoint)
```

**Figures**: `darwin_vs_bling_multipoint.png` (the main result: 38-point
map + DIC/ALK regression + all-6-tracer full context, Okabe-Ito
colorblind-safe palette), `darwin_vs_bling_fd_vs_adjoint.png` (earlier
single-point version), `initial_conditions_maps.png` /
`initial_conditions_profiles.png` (Darwin vs BLING IC sanity check).
