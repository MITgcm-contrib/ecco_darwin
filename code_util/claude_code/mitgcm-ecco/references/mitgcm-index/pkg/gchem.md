# pkg/gchem

Geochemistry driver: interface between ptracers and BGC packages (dic, bling, cfc, darwin); separate forcing step, surface forcing, chemistry calls.

**pkg_depend:** +ptracers  (`+` requires, `-` excludes)
**runtime switch:** `useGCHEM`-style flag in `data.pkg` (check exact name in packages_boot.F)
**reads:** `data.gchem`
**manual:** `doc/examples/cfc_offline/cfc_offline.rst`, `doc/examples/global_oce_biogeo/global_oce_biogeo.rst`, `doc/outp_pkgs/outp_pkgs.rst`
**adjoint support files:** gchem_ad_diff.list

## Namelist parameters
### GCHEM_PARM01
- `useCFC` — flag to turn on/off CFC pkg
- `useDIC` — flag to turn on/off DIC pkg
- `useBLING` — flag to turn on/off BLING pkg
- `useSPOIL` — flag to turn on/off SPOIL pkg
- `useDARWIN` — flag to turn on/off darwin pkg
- `fileName1`
- `fileName2`
- `fileName3`
- `fileName4`
- `fileName5`
- `gchem_int1`
- `gchem_int2`
- `gchem_int3`
- `gchem_int4`
- `gchem_int5`
- `nsubtime` — number of chemistry timesteps per deltaTtracer (default 1)
- `gchem_rl1`
- `gchem_rl2`
- `gchem_rl3`
- `gchem_rl4`
- `gchem_rl5`
- `gchem_ForcingPeriod` — periodic forcing parameter specific for gchem (secs)
- `gchem_ForcingCycle` — periodic forcing parameter specific for gchem (secs)
- `gchem_secondsPerYear` — used for gchem_insolation only (secs)
- `tIter0`

## CPP options (defaults as shipped)
- `GCHEM_SEPARATE_FORCING` (define, GCHEM_OPTIONS.h) — o Allow separated update of Geo-Chemistry and Advect-Diff (fractional time-stepping type) for some gchem tracers
- `GCHEM_ADD2TR_TENDENCY` (undef, GCHEM_OPTIONS.h) — o Allow single update of some gchem tracers, adding Geo-Chemistry tendency to Advect-Diff tendency
- `GCHEM_ADD2TR_TENDENCY` (define, GCHEM_OPTIONS.h)
- `GCHEM_ADD2TR_TENDENCY` (define, GCHEM_OPTIONS.h)

## Headers
- `GCHEM.h` — BOP
- `GCHEM_FIELDS.h` — BOP
- `GCHEM_OPTIONS.h` — BOP
- `GCHEM_SIZE.h` — BOP

## Routines (15)
`gchem_add_tendency.F`, `gchem_calc_tendency.F`, `gchem_check.F`, `gchem_cons.F`, `gchem_diagnostics_init.F`, `gchem_fields_load.F`, `gchem_forcing_sep.F`, `gchem_init_fixed.F`, `gchem_init_vari.F`, `gchem_insolation.F`, `gchem_output.F`, `gchem_readparms.F`, `gchem_surfmean.F`, `gchem_tr_register.F`, `gchem_write_pickup.F`

## Called from outside the package
- `GCHEM_OUTPUT` ← `model/src/do_the_model_io.F:225`
- `GCHEM_CALC_TENDENCY` ← `model/src/forward_step.F:690`
- `GCHEM_FORCING_SEP` ← `model/src/forward_step.F:1081`
- `GCHEM_FIELDS_LOAD` ← `model/src/load_fields_driver.F:243`
- `GCHEM_CHECK` ← `model/src/packages_check.F:308`
- `GCHEM_INIT_FIXED` ← `model/src/packages_init_fixed.F:438`
- `GCHEM_INIT_VARI` ← `model/src/packages_init_variables.F:348`
- `GCHEM_READPARMS` ← `model/src/packages_readparms.F:256`
- `GCHEM_WRITE_PICKUP` ← `model/src/packages_write_pickup.F:164`
- `GCHEM_INSOLATION` ← `pkg/bling/bling_light.F:166`
- `GCHEM_INSOLATION` ← `pkg/dic/bio_export.F:71`
- `GCHEM_ADD_TENDENCY` ← `pkg/ptracers/ptracers_apply_forcing.F:73`

## Verification experiments compiling it (6)
`cfc_example` `global_oce_biogeo_bling` `so_box_biogeo` `tutorial_cfc_offline` `tutorial_dic_adjoffline` `tutorial_global_oce_biogeo`
