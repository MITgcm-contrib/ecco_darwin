# pkg/matrix

Transport-matrix method: extracts explicit/implicit tracer transport matrices.

**pkg_depend:** +ptracers -gchem  (`+` requires, `-` excludes)
**runtime switch:** `useMATRIX`-style flag in `data.pkg` (check exact name in packages_boot.F)
**reads:** `data.matrix`

## Namelist parameters
### MATRIX_PARM01
- `expMatrixWriteTime`
- `impMatrixWriteTime`

## Headers
- `MATRIX.h` — 
- `MATRIX_OPTIONS.h` — Matrix pkg Options & Macros go here

## Routines (7)
`matrix_init_varia.F`, `matrix_output.F`, `matrix_readparms.F`, `matrix_store_tendency.F`, `matrix_write_grid.F`, `matrix_write_tendency.F`

## Called from outside the package
- `MATRIX_OUTPUT` ← `model/src/do_the_model_io.F:219`
- `MATRIX_INIT_VARIA` ← `model/src/packages_init_variables.F:363`
- `MATRIX_READPARMS` ← `model/src/packages_readparms.F:271`
- `MATRIX_STORE_TENDENCY_IMP` ← `model/src/tracers_correction_step.F:123`
- `MATRIX_STORE_TENDENCY_EXP` ← `pkg/ptracers/ptracers_integrate.F:440`

## Verification experiments compiling it (1)
`matrix_example`
