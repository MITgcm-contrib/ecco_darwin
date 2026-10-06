# pkg/bbl

Bottom boundary layer: a thin bottom layer with its own T/S exchanging with the interior and downslope (dense-overflow representation; mutually exclusive with down_slope).

**runtime switch:** `useBBL`-style flag in `data.pkg` (check exact name in packages_boot.F)
**reads:** `data.bbl`
**adjoint support files:** bbl_ad_diff.list

## README
```
Package `BBL` is a simple bottom boundary layer scheme.  The initial
motivation is to allow dense water that forms on the continental shelf around
Antarctica (High Salinity Shelf Water) in the CS510 configuration to sink to
the bottom of the model domain and to become a source of Antarctic Bottom
Water.  The bbl package aims to address the following two limitations of
package down_slope:

(i) In `pkg/down_slope`, dense water cannot flow down-slope unless there is a
step, i.e., a change of vertical level in the bathymetry.  In pkg/bbl, dense
water can flow even on a slight incline or flat bottom.

(ii) In `pkg/down_slope`, dense water is diluted as it flows into grid cells
whose thickness depends on model configuration, typically much thicker than a
bottom boundary layer.  In pkg/bbl, dense water is contained in a thin
sub-layer and hence able to preserve its tracer properties.

Specifically, the bottommost wet grid cell of thickness

    thk = hFacC(kBot) * drF(kBot),

with properties `tracer`, and density `rho` is divided in two sub-levels:

1. A bottom boundary layer with T/S tracer properties `bbl_tracer`,
density `bbl_rho`, and thickness `bbl_eta`.

2. A residual thickness `resThk = thk - bbl_eta` with tracer properties

    resTracer = ( tracer * thk - bbl_tracer * bbl_eta ) / resThk

such that the volume integral of `bbl_tracer` and `resTracer` is consistent with
```

## Namelist parameters
### BBL_PARM01
- `bbl_wvel` — default vertical entrainment velocity (m/s)
- `bbl_hvel`
- `bbl_initEta` — default initial thickness of BBL (m)
- `bbl_thetaFile`
- `bbl_saltFile`
- `bbl_etaFile`

## Headers
- `BBL.h` — bbl_wvel    :: default vertical entrainment velocity (m/s) bbl_hvvel   :: default horizontal velocity of BBL (m/s) bbl_initEta :: default initial thic
- `BBL_OPTIONS.h` — CPP options file for BBL Use this file for selecting options within package "BBL"

## Routines (17)
`bbl_calc_rho.F`, `bbl_calc_rhs.F`, `bbl_check.F`, `bbl_diagnostics_init.F`, `bbl_diagnostics_state.F`, `bbl_init_fixed.F`, `bbl_init_varia.F`, `bbl_read_pickup.F`, `bbl_readparms.F`, `bbl_tendency_apply.F`, `bbl_write_pickup.F`

## Called from outside the package
- `BBL_TENDENCY_APPLY_S` ← `model/src/apply_forcing.F:977`
- `BBL_TENDENCY_APPLY_T` ← `model/src/apply_forcing.F:745`
- `BBL_CALC_RHO` ← `model/src/do_oceanic_phys.F:750`
- `BBL_CALC_RHS` ← `model/src/do_oceanic_phys.F:1083`
- `BBL_DIAGNOSTICS_STATE` ← `model/src/do_statevars_diags.F:90`
- `BBL_TENDENCY_APPLY_S` ← `model/src/external_forcing.F:807`
- `BBL_TENDENCY_APPLY_T` ← `model/src/external_forcing.F:603`
- `BBL_CHECK` ← `model/src/packages_check.F:258`
- `BBL_INIT_FIXED` ← `model/src/packages_init_fixed.F:341`
- `BBL_INIT_VARIA` ← `model/src/packages_init_variables.F:280`
- `BBL_READPARMS` ← `model/src/packages_readparms.F:216`
- `BBL_WRITE_PICKUP` ← `model/src/packages_write_pickup.F:132`

## Verification experiments compiling it (1)
`global_oce_latlon`
