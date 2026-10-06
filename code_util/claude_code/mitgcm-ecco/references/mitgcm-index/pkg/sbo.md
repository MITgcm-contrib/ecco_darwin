# pkg/sbo

Solid-body ocean diagnostics: angular momentum, centre of mass (Earth rotation studies).

**runtime switch:** `useSBO`-style flag in `data.pkg` (check exact name in packages_boot.F)
**reads:** `data.sbo`
**manual:** `doc/examples/examples.rst`

## Namelist parameters
### SBO_PARM01
- `sbo_taveFreq`
- `sbo_monFreq` — SBO monitor frequency           (s)

## Headers
- `SBO.h` — Basic header for SBO
- `SBO_OPTIONS.h` — CPP options file for SBO package. Use this file for selecting options within the SBO package.

## Routines (5)
`sbo_calc.F`, `sbo_check.F`, `sbo_output.F`, `sbo_readparms.F`, `sbo_rho.F`

## Called from outside the package
- `SBO_CALC` ← `model/src/do_the_model_io.F:181`
- `SBO_OUTPUT` ← `model/src/do_the_model_io.F:182`
- `SBO_CHECK` ← `model/src/packages_check.F:451`
- `SBO_READPARMS` ← `model/src/packages_readparms.F:358`

## Verification experiments compiling it (1)
`global_ocean.90x40x15`
