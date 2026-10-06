# pkg/ocn_compon_interf

Ocean-side interface of the coupler for coupled atmosphere-ocean runs.

**runtime switch:** `useOCN_COMPON_INTERF`-style flag in `data.pkg` (check exact name in packages_boot.F)
**reads:** `data.cpl`

## Namelist parameters
### CPL_OCN_PARAM
- `cpl_earlyExpImpCall`  _[ifdef COMPONENT_MODULE]_
- `useImportHFlx` — True => use the Imported HeatFlux from couler  _[ifdef COMPONENT_MODULE]_
- `useImportFW` — True => use the Imported Fresh Water flux fr cpl  _[ifdef COMPONENT_MODULE]_
- `useImportTau` — True => use the Imported Wind-Stress from couler  _[ifdef COMPONENT_MODULE]_
- `useImportSLP` — True => use the Imported Sea-level Pressure  _[ifdef COMPONENT_MODULE]_
- `useImportRunOff` — True => use the Imported RunOff flux from coupler  _[ifdef COMPONENT_MODULE]_
- `useImportSIce` — True => use the Imported Sea-Ice mass as ice-loading  _[ifdef COMPONENT_MODULE]_
- `useImportThSIce` — True => use the Imported thSIce state var from coupler  _[ifdef COMPONENT_MODULE]_
- `useImportSltPlm` — True => use the Imported Salt-Plume flux from coupler  _[ifdef COMPONENT_MODULE]_
- `useImportFice` — True => use the Imported Seaice fraction (DIC-only)  _[ifdef COMPONENT_MODULE]_
- `useImportCO2` — True => use the Imported atmos. CO2 from coupler  _[ifdef COMPONENT_MODULE]_
- `useImportWSpd` — True => use the Imported surf. Wind speed from coupler  _[ifdef COMPONENT_MODULE]_
- `cpl_taveFreq`  _[ifdef COMPONENT_MODULE]_
- `cpl_snapshot_mnc`  _[ifdef COMPONENT_MODULE]_
- `cpl_timeave_mnc`  _[ifdef COMPONENT_MODULE]_

## Headers
- `CPL_PARAMS.h` — Header file for Coupling component interface this version is specific to 1 component (ocean)
- `OCNCPL.h` — Variables shared between coupling layer and ocean component. These variables are used in the ocean component. Grid variables have already been mapped/
- `OCN_CPL_OPTIONS.h` — Package-specific Options & Macros go here

## Routines (21)
`cpl_diagnostics_fill.F`, `cpl_diagnostics_init.F`, `cpl_exch_configs.F`, `cpl_export_import_data.F`, `cpl_import_cplparms.F`, `cpl_ini_vars.F`, `cpl_init.F`, `cpl_init_fixed.F`, `cpl_readparms.F`, `cpl_register.F`, `cpl_write_pickup.F`, `ocn_apply_import.F`, `ocn_check_cplconfig.F`, `ocn_cpl_diags.F`, `ocn_cpl_read_pickup.F`, `ocn_export_data.F`, `ocn_export_fields.F`, `ocn_export_ocnconfig.F`, `ocn_import_atmconfig.F`, `ocn_import_fields.F`, `ocn_store_my_data.F`

## Called from outside the package
- `CPL_REGISTER` ← `eesupp/src/eeboot.F:159`
- `CPL_INIT` ← `eesupp/src/eeboot_minimal.F:171`
- `OCN_APPLY_IMPORT` ← `model/src/do_oceanic_phys.F:347`
- `OCN_EXPORT_DATA` ← `model/src/do_oceanic_phys.F:484`
- `CPL_EXPORT_IMPORT_DATA` ← `model/src/forward_step.F:587`
- `CPL_EXCH_CONFIGS` ← `model/src/initialise_fixed.F:290`
- `CPL_IMPORT_CPLPARMS` ← `model/src/initialise_fixed.F:125`
- `CPL_INIT_FIXED` ← `model/src/packages_init_fixed.F:618`
- `CPL_INI_VARS` ← `model/src/packages_init_variables.F:544`
- `CPL_READPARMS` ← `model/src/packages_readparms.F:409`
- `CPL_WRITE_PICKUP` ← `model/src/packages_write_pickup.F:231`
- `CPL_DIAGNOSTICS_FILL` ← `pkg/atm_compon_interf/cpl_export_import_data.F:97`
- `CPL_DIAGNOSTICS_INIT` ← `pkg/atm_compon_interf/cpl_init_fixed.F:26`

## Verification experiments compiling it (1)
`cpl_aim+ocn`
