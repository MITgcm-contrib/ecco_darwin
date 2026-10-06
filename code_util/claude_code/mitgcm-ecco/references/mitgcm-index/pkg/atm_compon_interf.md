# pkg/atm_compon_interf

Atmosphere-side interface of the coupler (cpl_aim+ocn style atmosphere-ocean coupling).

**runtime switch:** `useATM_COMPON_INTERF`-style flag in `data.pkg` (check exact name in packages_boot.F)
**reads:** `data.cpl`

## Namelist parameters
### CPL_ATM_PARAM
- `cpl_earlyExpImpCall`  _[ifdef COMPONENT_MODULE]_
- `cpl_oldPickup` — restart from an old pickup (= until checkpoint 59h)  _[ifdef COMPONENT_MODULE]_
- `useImportMxlD` — True => use Imported Mix.Layer Detph from coupler  _[ifdef COMPONENT_MODULE]_
- `useImportSST` — True => use the Imported SST from coupler  _[ifdef COMPONENT_MODULE]_
- `useImportSSS` — True => use the Imported SSS from coupler  _[ifdef COMPONENT_MODULE]_
- `useImportVsq` — True => use the Imported Surf. velocity^2  _[ifdef COMPONENT_MODULE]_
- `useImportThSIce` — True => use the Imported thSIce state var from coupler  _[ifdef COMPONENT_MODULE]_
- `useImportFlxCO2` — True => use the Imported air-sea CO2 flux from coupler  _[ifdef COMPONENT_MODULE]_
- `cpl_atmSendFrq` — Frequency^-1 for sending data to coupler (s)  _[ifdef COMPONENT_MODULE]_
- `maxNumberPrint` — max number of printed Export/Import messages  _[ifdef COMPONENT_MODULE]_

## Headers
- `ATMCPL.h` — Variables shared between atmos. component to coupler layer. These variables are used in the atmos component. Grid variables have already been mapped/i
- `ATM_CPL_OPTIONS.h` — Package-specific Options & Macros go here
- `CPL_PARAMS.h` — Header file for Coupling component interface this version is specific to 1 component (atmos)

## Routines (17)
`atm_apply_import.F`, `atm_check_cplconfig.F`, `atm_cpl_read_pickup.F`, `atm_export_atmconfig.F`, `atm_export_fields.F`, `atm_export_fld.F`, `atm_get_atmconfig.F`, `atm_import_fields.F`, `atm_import_ocnconfig.F`, `atm_store_aim_fields.F`, `atm_store_aim_wndstr.F`, `atm_store_dynvars.F`, `atm_store_land.F`, `atm_store_my_data.F`, `atm_store_surfflux.F`, `atm_store_thsice.F`, `cpl_output.F`

## Called from outside the package
- `ATM_STORE_MY_DATA` ← `pkg/aim_v23/aim_do_physics.F:211`
- `ATM_APPLY_IMPORT` ← `pkg/aim_v23/aim_surf_bc.F:353`
- `ATM_APPLY_IMPORT` ← `pkg/atm_phys/atm_phys_driver.F:159`
- `ATM_STORE_MY_DATA` ← `pkg/atm_phys/atm_phys_driver.F:530`

## Verification experiments compiling it (1)
`cpl_aim+ocn`
