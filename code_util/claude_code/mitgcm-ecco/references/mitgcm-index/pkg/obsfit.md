# pkg/obsfit

Model-data comparison for generic (unstructured) observations in the cost function.

**pkg_depend:** +cal  (`+` requires, `-` excludes)
**runtime switch:** `useOBSFIT`-style flag in `data.pkg` (check exact name in packages_boot.F)
**reads:** `data.obsfit`
**manual:** `doc/examples/examples.rst`, `doc/ocean_state_est/ocean_state_est.rst`
**adjoint support files:** obsfit_ad_diff.list

## Namelist parameters
### OBSFIT_NML
- `obsfitDir`
- `obsfitFiles`
- `mult_obsfit`
- `obsfit_facmod`
- `obsfitDoNcOutput`
- `obsfitDoGenGrid`

## CPP options (defaults as shipped)
- `OBSFIT_USE_MDSFINDUNITS` (undef, OBSFIT_OPTIONS.h) — To use file units between 9 and 99 (seems to conflict with NF_OPEN some times, but is needed when using g77)
- `ALLOW_OBSFIT_EXCLUDE_CORNERS` (undef, OBSFIT_OPTIONS.h)

## Headers
- `OBSFIT.h` — BOP
- `OBSFIT_OPTIONS.h` — BOP
- `OBSFIT_SIZE.h` — BOP

## Routines (29)
`obsfit_active_file.F`, `obsfit_active_file_ad.F`, `obsfit_active_file_control.F`, `obsfit_active_file_g.F`, `obsfit_cost.F`, `obsfit_cost_final.F`, `obsfit_findunit.F`, `obsfit_ini_io.F`, `obsfit_init_equifiles.F`, `obsfit_init_fixed.F`, `obsfit_init_varia.F`, `obsfit_inloop.F`, `obsfit_nc_utils.F`, `obsfit_read_obs.F`, `obsfit_readparms.F`, `obsfit_sampling.F`

## Called from outside the package
- `OBSFIT_INI_IO` ← `model/src/ini_model_io.F:240`
- `OBSFIT_INIT_FIXED` ← `model/src/packages_init_fixed.F:399`
- `OBSFIT_INIT_VARIA` ← `model/src/packages_init_variables.F:583`
- `OBSFIT_READPARMS` ← `model/src/packages_readparms.F:353`
- `OBSFIT_COST` ← `model/src/the_main_loop.F:762`
- `OBSFIT_INLOOP` ← `model/src/the_main_loop.F:694,759`
- `OBSFIT_NC_CLOSE` ← `model/src/the_model_main.F:777`
- `OBSFIT_COST_FINAL` ← `pkg/cost/cost_final.F:92`

## Verification experiments compiling it (1)
`global_oce_biogeo_bling`
