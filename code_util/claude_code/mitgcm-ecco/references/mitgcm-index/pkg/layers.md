# pkg/layers

Diagnostics of transport in isopycnal/temperature/salinity layers (residual overturning).

**runtime switch:** `useLAYERS`-style flag in `data.pkg` (check exact name in packages_boot.F)
**reads:** `data.layers`
**manual:** `doc/examples/examples.rst`, `doc/examples/reentrant_channel/reentrant_channel.rst`
**adjoint support files:** layers_ad_diff.list

## Namelist parameters
### LAYERS_PARM01
- `layers_G`
- `layers_taveFreq`
- `layers_diagFreq`
- `LAYER_nb`
- `layers_kref`
- `useBOLUS`
- `layers_bolus`
- `layers_name`
- `layers_bounds` — boundaries of tracer layers
- `layers_krho`

## CPP options (defaults as shipped)
- `LAYERS_UFLUX` (define, LAYERS_OPTIONS.h) — Compute isopycnal tranports in the U direction?
- `LAYERS_VFLUX` (define, LAYERS_OPTIONS.h) — Compute isopycnal tranports in the V direction?
- `LAYERS_THICKNESS` (define, LAYERS_OPTIONS.h) — Keep track of layer thicknesses?
- `LAYERS_THERMODYNAMICS` (undef, LAYERS_OPTIONS.h) — Do water mass thermodynamics?
- `LAYERS_FINEGRID_DIAPYCNAL` (undef, LAYERS_OPTIONS.h) — Use refined grid for diapycnal terms? (gives worse results)
- `LAYERS_MNC` (undef, LAYERS_OPTIONS.h) — The MNC stuff is too complicated
- `LAYERS_PRHO_REF` (define, LAYERS_OPTIONS.h) — Allow use of potential density as a layering field.
- `LAYERS_MSE` (undef, LAYERS_OPTIONS.h) — Allow use of Moist Static Energy as a coordinate (relevant in the atmosphere)

## Headers
- `LAYERS.h` — --   Header for LAYERS package. By Ryan Abernathey. --   For computing volume fluxes in isopyncal layers
- `LAYERS_OPTIONS.h` — CPP options file for LAYERS package Use this file for selecting options within package "LAYERS"
- `LAYERS_P2SHARE.h` — BOP
- `LAYERS_SIZE.h` — * Compiled-in size options for the LAYERS package * - Just as you have to define Nr in SIZE.h, you must define the number of vertical layers for isopy

## Routines (17)
`layers_calc.F`, `layers_calc_divergence.F`, `layers_check.F`, `layers_diagnostics_init.F`, `layers_diapycnal.F`, `layers_fill.F`, `layers_fluxcalc.F`, `layers_init_fixed.F`, `layers_init_varia.F`, `layers_locate.F`, `layers_mnc_init.F`, `layers_output.F`, `layers_readparms.F`, `layers_thermodynamics.F`, `layers_wsurf_tr.F`

## Called from outside the package
- `LAYERS_FILL` ← `model/src/diags_oceanic_surf_flux.F:155,193`
- `LAYERS_CALC` ← `model/src/do_the_model_io.F:236`
- `LAYERS_FILL` ← `model/src/impldiff.F:385`
- `LAYERS_CHECK` ← `model/src/packages_check.F:470`
- `LAYERS_INIT_FIXED` ← `model/src/packages_init_fixed.F:609`
- `LAYERS_INIT_VARIA` ← `model/src/packages_init_variables.F:503`
- `LAYERS_READPARMS` ← `model/src/packages_readparms.F:384`
- `LAYERS_WSURF_TR` ← `model/src/thermodynamics.F:161`
- `LAYERS_FILL` ← `pkg/diagnostics/diagnostics_fill_state.F:413,438,727,752`
- `LAYERS_FILL` ← `pkg/generic_advdiff/gad_advection.F:853,854,1087`
- `LAYERS_FILL` ← `pkg/generic_advdiff/gad_calc_rhs.F:320,369,449,498`
- `LAYERS_FILL` ← `pkg/generic_advdiff/gad_implicit_r.F:344,442`
- `LAYERS_FILL` ← `pkg/generic_advdiff/gad_som_advect.F:471,472,666`

## Verification experiments compiling it (3)
`cfc_example` `exp4` `tutorial_reentrant_channel`
