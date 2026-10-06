# pkg/cfc

CFC-11/CFC-12 ocean tracers with air-sea gas exchange (OCMIP protocol) via gchem/ptracers.

**pkg_depend:** +gchem  (`+` requires, `-` excludes)
**runtime switch:** `useCFC`-style flag in `data.pkg` (check exact name in packages_boot.F)
**reads:** `data.cfc`
**manual:** `doc/examples/cfc_offline/cfc_offline.rst`
**adjoint support files:** cfc_ad_check_lev1_dir.h, cfc_ad_check_lev2_dir.h, cfc_ad_check_lev3_dir.h, cfc_ad_check_lev4_dir.h, cfc_ad_diff.list

## Namelist parameters
### CFC_FORCING
- `atmCFC_inpFile` — file name of Atmospheric CFC time series (ASCII file)
- `atmCFC_recSepTime` — time spacing between 2 records of atmos CFC [s]
- `atmCFC_timeOffset` — time offset for atmos CFC (cfcTime = myTime + offSet)
- `atmCFC_yNorthBnd` — Northern Lat boundary for interpolation [y-unit]
- `atmCFC_ySouthBnd` — Southern Lat boundary for interpolation [y-unit]
- `CFC_windFile` — file name of wind speeds
- `CFC_atmospFile` — file name of atmospheric pressure
- `CFC_iceFile` — file name of seaice fraction
- `CFC_forcingPeriod` — record spacing time for CFC forcing (seconds)
- `CFC_forcingCycle` — periodic-cycle freq for CFC forcing (seconds)

## Headers
- `CFC.h` — schmidt number coefficients
- `CFC_ATMOS.h` — BOP
- `CFC_SIZE.h` — BOP
- `cfc_ad_check_lev1_dir.h` — CADJ STORE AtmosCFC11   = comlev1, key = ikey_dynamics CADJ STORE AtmosCFC12   = comlev1, key = ikey_dynamics ADJ STORE Atmosp       = comlev1, key = 
- `cfc_ad_check_lev2_dir.h` — CADJ STORE AtmosCFC11   = tapelev2, key = ilev_2 CADJ STORE AtmosCFC12   = tapelev2, key = ilev_2 ADJ STORE Atmosp       = tapelev2, key = ilev_2 ADJ 
- `cfc_ad_check_lev3_dir.h` — CADJ STORE AtmosCFC11   = tapelev3, key = ilev_3 CADJ STORE AtmosCFC12   = tapelev3, key = ilev_3 ADJ STORE Atmosp       = tapelev3, key = ilev_3 ADJ 
- `cfc_ad_check_lev4_dir.h` — CADJ STORE AtmosCFC11   = tapelev4, key = ilev_4 CADJ STORE AtmosCFC12   = tapelev4, key = ilev_4 ADJ STORE Atmosp       = tapelev4, key = ilev_4 ADJ 

## Routines (10)
`cfc11_forcing.F`, `cfc11_surfforcing.F`, `cfc12_forcing.F`, `cfc12_surfforcing.F`, `cfc_atmos.F`, `cfc_check.F`, `cfc_fields_load.F`, `cfc_param.F`, `cfc_readparms.F`, `cfc_tr_register.F`

## Called from outside the package
- `CFC11_FORCING` ← `pkg/gchem/gchem_calc_tendency.F:114`
- `CFC12_FORCING` ← `pkg/gchem/gchem_calc_tendency.F:121`
- `CFC_CHECK` ← `pkg/gchem/gchem_check.F:231`
- `CFC_FIELDS_LOAD` ← `pkg/gchem/gchem_fields_load.F:47`
- `CFC_ATMOS` ← `pkg/gchem/gchem_init_fixed.F:38`
- `CFC_PARAM` ← `pkg/gchem/gchem_init_fixed.F:36`
- `CFC_READPARMS` ← `pkg/gchem/gchem_readparms.F:152`
- `CFC_TR_REGISTER` ← `pkg/gchem/gchem_tr_register.F:75`

## Verification experiments compiling it (2)
`cfc_example` `tutorial_cfc_offline`
