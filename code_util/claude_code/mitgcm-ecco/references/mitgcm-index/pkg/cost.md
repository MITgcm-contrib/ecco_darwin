# pkg/cost

Cost-function framework for adjoint/state estimation (accumulates and writes the cost).

**in groups:** adjoint
**runtime switch:** `useCOST`-style flag in `data.pkg` (check exact name in packages_boot.F)
**reads:** `data.cost`
**manual:** `doc/autodiff/autodiff.rst`, `doc/examples/global_oce_optim/global_oce_optim.rst`, `doc/examples/tracer_adjsens/tracer_adjsens.rst`
**adjoint support files:** cost_ad_diff.list

## Namelist parameters
### COST_NML
- `mult_atl`
- `mult_test`
- `mult_tracer`
- `multTheta`
- `multSalt`
- `multUvel`
- `multVvel`
- `multEtan`
- `mult_depth`
- `mult_temp_tut`  _[ifdef ALLOW_COST_HFLUXM]_
- `mult_hflux_tut`  _[ifdef ALLOW_COST_HFLUXM]_
- `lastinterval`
- `cost_mask_file`

## CPP options (defaults as shipped)
- `ALLOW_COST_STATE_FINAL` (undef, COST_OPTIONS.h)
- `ALLOW_COST_VECTOR` (undef, COST_OPTIONS.h)
- `ALLOW_COST_ATLANTIC_HEAT` (undef, COST_OPTIONS.h) — >>> Cost function contributions
- `ALLOW_COST_ATLANTIC_HEAT_DOMASS` (undef, COST_OPTIONS.h)
- `ALLOW_COST_TEST` (undef, COST_OPTIONS.h)
- `ALLOW_COST_TSQUARED` (undef, COST_OPTIONS.h)
- `ALLOW_COST_TRACER` (undef, COST_OPTIONS.h)
- `ALLOW_COST_TEMP` (undef, COST_OPTIONS.h) — List these options here: -  User needs to provide "cost_temp.F"  code before defining following option:
- `ALLOW_COST_HFLUXM` (undef, COST_OPTIONS.h) — -  User needs to provide "cost_hflux.F" code before defining following option:
- `ALLOW_DIC_COST` (undef, COST_OPTIONS.h) — -  The following option contains some hacks (reset cost-function):
- `ALLOW_THSICE_COST_TEST` (undef, COST_OPTIONS.h)
- `ALLOW_COST_SHELFICE` (undef, COST_OPTIONS.h)
- `ALLOW_COST_STREAMICE` (undef, COST_OPTIONS.h)

## Headers
- `COST_OPTIONS.h` — BOP
- `COST_TAP_ADJ.h` — HEADER COST_TAP_ADJ
- `adcost.h` — HEADER ADCOST Header for model-data comparison; adjoint part. started: Christian Eckert eckert@mit.edu  06-Apr-2000 changed: Christian Eckert eckert@m
- `cost.h` — HEADER COST Header for model-data comparison. The individual cost function contributions are multiplied by factors mult_"var" which allow to switch of
- `g_cost.h` — HEADER G_COST Header for model-data comparison; tangent linear part. started: Christian Eckert eckert@mit.edu  06-Apr-2000 changed: Christian Eckert e

## Routines (18)
`cost_accumulate_mean.F`, `cost_atlantic_heat.F`, `cost_check.F`, `cost_copy_file.F`, `cost_dependent_init.F`, `cost_depth.F`, `cost_driver.F`, `cost_final.F`, `cost_final_restore.F`, `cost_final_store.F`, `cost_init_fixed.F`, `cost_init_varia.F`, `cost_readparms.F`, `cost_state_final.F`, `cost_test.F`, `cost_tile.F`, `cost_tracer.F`, `cost_vector.F`

## Called from outside the package
- `COST_TILE` ← `model/src/forward_step.F:1163`
- `COST_INIT_VARIA` ← `model/src/initialise_varia.F:268`
- `COST_CHECK` ← `model/src/packages_check.F:424`
- `COST_INIT_FIXED` ← `model/src/packages_init_fixed.F:360`
- `COST_READPARMS` ← `model/src/packages_readparms.F:328`
- `COST_DRIVER` ← `model/src/the_main_loop.F:769`
- `COST_FINAL` ← `model/src/the_main_loop.F:774`
- `COST_DEPENDENT_INIT` ← `model/src/the_model_main.F:651`
- `COST_FINAL_RESTORE` ← `model/src/the_model_main.F:690`
- `COST_FINAL_STORE` ← `model/src/the_model_main.F:680`
- `COST_DEPENDENT_INIT` ← `pkg/openad/the_model_main.F:181`

## Verification experiments compiling it (16)
`1D_ocean_ice_column` `bottom_ctrl_5x5` `global_oce_biogeo_bling` `global_oce_latlon` `global_ocean.90x40x15` `global_ocean.cs32x15` `halfpipe_streamice` `hs94.1x64x5` `isomip` `lab_sea` `obcs_ctrl` `offline_exf_seaice` `tutorial_dic_adjoffline` `tutorial_global_oce_biogeo` `tutorial_global_oce_optim` `tutorial_tracer_adjsens`
