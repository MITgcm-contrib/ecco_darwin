# pkg/frazil

Frazil ice formation: removes supercooling and adjusts heat/salt (used with shelfice/seaice set-ups).

**runtime switch:** `useFRAZIL`-style flag in `data.pkg` (check exact name in packages_boot.F)
**adjoint support files:** frazil_ad_diff.list

## Headers
- `FRAZIL.h` — FrazilForcingT : frazil temperature forcing, > 0 increases theta [W/m^2]
- `FRAZIL_OPTIONS.h` — CPP options file for FRAZIL Use this file for selecting options within package "Frazil"

## Routines (5)
`frazil_calc_rhs.F`, `frazil_diagnostics_init.F`, `frazil_init_fixed.F`, `frazil_init_varia.F`, `frazil_tendency_apply.F`

## Called from outside the package
- `FRAZIL_TENDENCY_APPLY_T` ← `model/src/apply_forcing.F:697`
- `FRAZIL_CALC_RHS` ← `model/src/do_oceanic_phys.F:372`
- `FRAZIL_TENDENCY_APPLY_T` ← `model/src/external_forcing.F:555`
- `FRAZIL_INIT_FIXED` ← `model/src/packages_init_fixed.F:505`
- `FRAZIL_INIT_VARIA` ← `model/src/packages_init_variables.F:422`

## Verification experiments compiling it (1)
`global_oce_latlon`
