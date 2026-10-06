# pkg/icefront

Vertical ice-front (tidewater glacier face) melt parameterisation at side walls.

**runtime switch:** `useICEFRONT`-style flag in `data.pkg` (check exact name in packages_boot.F)
**reads:** `data.icefront`
**manual:** `doc/examples/examples.rst`

## Namelist parameters
### ICEFRONT_PARM01
- `rhoIcefront`
- `ICEFRONTkappa`
- `ICEFRONTlatentHeat`
- `ICEFRONTHeatCapacity_Cp`
- `ICEFRONTthetaSurface`
- `applyIcefrontTendT`
- `applyIcefrontTendS`
- `ICEFRONTdepthFile`
- `ICEFRONTlengthFile`
### ICEFRONT_EXF_PARM02
- `SGRunOffFile`  _[ifdef ALLOW_EXF]_
- `SGRunOffperiod`  _[ifdef ALLOW_EXF]_
- `SGRunOffStartTime`  _[ifdef ALLOW_EXF]_
- `SGRunOffstartdate1`  _[ifdef ALLOW_EXF]_
- `SGRunOffstartdate2`  _[ifdef ALLOW_EXF]_
- `SGRunOffconst`  _[ifdef ALLOW_EXF]_
- `SGRunOff_inscal`  _[ifdef ALLOW_EXF]_
- `SGRunOff_remov_intercept`  _[ifdef ALLOW_EXF]_
- `SGRunOff_remov_slope`  _[ifdef ALLOW_EXF]_

## Headers
- `ICEFRONT.h` — BOP
- `ICEFRONT_OPTIONS.h` — CPP options file for ICEFRONT package. Use this file for selecting options within the ICEFRONT package.

## Routines (8)
`icefront_check.F`, `icefront_diagnostics_init.F`, `icefront_init_fixed.F`, `icefront_init_varia.F`, `icefront_readparms.F`, `icefront_tendency_apply.F`, `icefront_thermodynamics.F`

## Called from outside the package
- `ICEFRONT_TENDENCY_APPLY_S` ← `model/src/apply_forcing.F:945`
- `ICEFRONT_TENDENCY_APPLY_T` ← `model/src/apply_forcing.F:713`
- `ICEFRONT_THERMODYNAMICS` ← `model/src/do_oceanic_phys.F:539`
- `ICEFRONT_TENDENCY_APPLY_S` ← `model/src/external_forcing.F:775`
- `ICEFRONT_TENDENCY_APPLY_T` ← `model/src/external_forcing.F:571`
- `ICEFRONT_CHECK` ← `model/src/packages_check.F:373`
- `ICEFRONT_INIT_FIXED` ← `model/src/packages_init_fixed.F:496`
- `ICEFRONT_INIT_VARIA` ← `model/src/packages_init_variables.F:413`
- `ICEFRONT_READPARMS` ← `model/src/packages_readparms.F:291`

## Verification experiments compiling it (1)
`isomip`
