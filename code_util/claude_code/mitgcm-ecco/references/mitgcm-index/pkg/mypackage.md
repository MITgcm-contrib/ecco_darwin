# pkg/mypackage

TEMPLATE package: skeleton showing how to write a new package (readparms, check, init, diagnostics, pickup, tendencies). Start new packages from it.

**runtime switch:** `useMYPACKAGE`-style flag in `data.pkg` (check exact name in packages_boot.F)
**reads:** `data.mypackage`
**manual:** `doc/contributing/contributing.rst`
**adjoint support files:** mypackage_ad_diff.list

## Namelist parameters
### MYPACKAGE_PARM01
- `myPa_MNC`
- `myPa_StaV_Cgrid`
- `myPa_Tend_Cgrid`
- `myPa_applyTendT`
- `myPa_applyTendS`
- `myPa_applyTendU`
- `myPa_applyTendV`
- `myPa_doSwitch1`
- `myPa_doSwitch2`
- `myPa_index1`
- `myPa_index2`
- `myPa_param1`
- `myPa_param2`
- `myPa_string1`
- `myPa_string2`
- `myPa_Scal1File`
- `myPa_Scal2File`
- `myPa_VelUFile`
- `myPa_VelVFile`
- `myPa_Surf1File`
- `myPa_Surf2File`

## CPP options (defaults as shipped)
- `MYPACKAGE_3D_STATE` (define, MYPACKAGE_OPTIONS.h) — to reduce memory storage, disable unused array with those CPP flags :
- `MYPACKAGE_2D_STATE` (define, MYPACKAGE_OPTIONS.h)
- `MYPACKAGE_TENDENCY` (define, MYPACKAGE_OPTIONS.h)
- `MYPA_SPECIAL_COMPILE_OPTION1` (undef, MYPACKAGE_OPTIONS.h)
- `MYPA_SPECIAL_COMPILE_OPTION2` (define, MYPACKAGE_OPTIONS.h)

## Headers
- `MYPACKAGE.h` — BOP
- `MYPACKAGE_OPTIONS.h` — BOP

## Routines (14)
`mypackage_calc_rhs.F`, `mypackage_check.F`, `mypackage_diagnostics_init.F`, `mypackage_diagnostics_state.F`, `mypackage_init_fixed.F`, `mypackage_init_varia.F`, `mypackage_mnc_init.F`, `mypackage_read_pickup.F`, `mypackage_readparms.F`, `mypackage_tendency_apply.F`, `mypackage_write_pickup.F`

## Called from outside the package
- `MYPACKAGE_TENDENCY_APPLY_S` ← `model/src/apply_forcing.F:985`
- `MYPACKAGE_TENDENCY_APPLY_T` ← `model/src/apply_forcing.F:753`
- `MYPACKAGE_TENDENCY_APPLY_U` ← `model/src/apply_forcing.F:189`
- `MYPACKAGE_TENDENCY_APPLY_V` ← `model/src/apply_forcing.F:378`
- `MYPACKAGE_CALC_RHS` ← `model/src/do_oceanic_phys.F:1090`
- `MYPACKAGE_DIAGNOSTICS_STATE` ← `model/src/do_statevars_diags.F:120`
- `MYPACKAGE_TENDENCY_APPLY_S` ← `model/src/external_forcing.F:815`
- `MYPACKAGE_TENDENCY_APPLY_T` ← `model/src/external_forcing.F:611`
- `MYPACKAGE_TENDENCY_APPLY_U` ← `model/src/external_forcing.F:140`
- `MYPACKAGE_TENDENCY_APPLY_V` ← `model/src/external_forcing.F:279`
- `MYPACKAGE_CHECK` ← `model/src/packages_check.F:508`
- `MYPACKAGE_INIT_FIXED` ← `model/src/packages_init_fixed.F:646`
- `MYPACKAGE_INIT_VARIA` ← `model/src/packages_init_variables.F:550`
- `MYPACKAGE_READPARMS` ← `model/src/packages_readparms.F:423`
- `MYPACKAGE_WRITE_PICKUP` ← `model/src/packages_write_pickup.F:245`

## Verification experiments compiling it (1)
`hs94.1x64x5`
