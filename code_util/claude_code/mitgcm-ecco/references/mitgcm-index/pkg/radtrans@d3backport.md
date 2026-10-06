# pkg/radtrans  (d3backport: ~/Documents/research/debug/darwin3)

---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+---- BOP

**vs its upstream base (merge-base, see README):** new (not in its upstream base)

**runtime switch:** `useRADTRANS`-style flag in `data.pkg` (check exact name in packages_boot.F)
**reads:** `data.radtrans`

## Namelist parameters
### RADTRANS_FORCING_PARAMS
- `RT_Edfile` — downward direct irradiance below sea surface [W/m^2] per waveband, not taking into account ice cover
- `RT_Esfile` — downward diffuse irradiance below sea surface [W/m^2] per waveband, not taking into account ice cover
- `RT_E_mask`
- `RT_E_period`
- `RT_E_RepCycle`
- `RT_E_startTime`
- `RT_E_startdate1`
- `RT_E_startdate2`
- `RT_Ed_const`
- `RT_Ed_exfremo_intercept`
- `RT_Ed_exfremo_slope`
- `RT_inscal_Ed`
- `RT_Es_const`
- `RT_Es_exfremo_intercept`
- `RT_Es_exfremo_slope`
- `RT_inscal_Es`
- `RT_icefile` — fraction of the sea surface covered by ice used to reduce incoming irradiances
- `RT_iceperiod`
- `RT_iceRepCycle`
- `RT_iceStartTime`
- `RT_icestartdate1`
- `RT_icestartdate2`
- `RT_iceconst`
- `RT_ice_exfremo_intercept`
- `RT_ice_exfremo_slope`
- `RT_icemask`
- `RT_inscal_ice`
### RADTRANS_INTERP_PARAMS
- `RT_E_lon0`  _[ifdef USE_EXF_INTERPOLATION]_
- `RT_E_lat0`  _[ifdef USE_EXF_INTERPOLATION]_
- `RT_E_nlon`  _[ifdef USE_EXF_INTERPOLATION]_
- `RT_E_nlat`  _[ifdef USE_EXF_INTERPOLATION]_
- `RT_E_lon_inc`  _[ifdef USE_EXF_INTERPOLATION]_
- `RT_E_interpMethod`  _[ifdef USE_EXF_INTERPOLATION]_
- `RT_E_lat_inc`  _[ifdef USE_EXF_INTERPOLATION]_
- `RT_ice_lon0`  _[ifdef USE_EXF_INTERPOLATION]_
- `RT_ice_lat0`  _[ifdef USE_EXF_INTERPOLATION]_
- `RT_ice_nlon`  _[ifdef USE_EXF_INTERPOLATION]_
- `RT_ice_nlat`  _[ifdef USE_EXF_INTERPOLATION]_
- `RT_ice_lon_inc`  _[ifdef USE_EXF_INTERPOLATION]_
- `RT_ice_interpMethod`  _[ifdef USE_EXF_INTERPOLATION]_
- `RT_ice_lat_inc`  _[ifdef USE_EXF_INTERPOLATION]_
### RADTRANS_PARAMS
- `RT_refract_water` — refractive index of water
- `RT_rmud_max` — cutoff for inverse cosine of solar zenith angle
- `RT_wbEdges` — waveband edges [nm]
- `RT_wbRefWLs` — reference wavelengths for wavebands [nm]
- `RT_kmax` — maximum depth index for radtrans computations
- `RT_useOASIMrmud` — flag for using cosine of solar zenith angle from oasim pkg
- `RT_useMeanCosSolz` — flag for using mean daytime cosine of solar zenith angle
- `RT_useNoonSolz` — flag for using noon solar zenith angle; if false use angle at actual time
- `RT_sfcIrrThresh` — minimum irradiance for radiative transfer computations [W/m^2]
- `RT_oasimWgt`
### RADTRANS_DEPENDENT
- `RT_wbWidths`

## CPP options (defaults as shipped)
- `RADTRANS_DIAG_SOLUTION` (undef, RADTRANS_OPTIONS.h) — fill diagnostics for radiative transfer solution parameters

## Headers
- `RADTRANS_EXF_PARAMS.h` — BOP
- `RADTRANS_FIELDS.h` — BOP
- `RADTRANS_OPTIONS.h` — BOP
- `RADTRANS_PARAMS.h` — BOP
- `RADTRANS_SIZE.h` — BOP

## Routines (13)
`radtrans_calc.F`, `radtrans_check.F`, `radtrans_declination_spencer.F`, `radtrans_diagnostics_init.F`, `radtrans_fields_load.F`, `radtrans_init_fixed.F`, `radtrans_init_varia.F`, `radtrans_monitor.F`, `radtrans_readparms.F`, `radtrans_rmud_below.F`, `radtrans_solve.F`, `radtrans_solve_tridiag.F`, `radtrans_solz_daytime.F`

## Called from outside the package
- `RADTRANS_FIELDS_LOAD` ← `model/src/load_fields_driver.F:198`
- `RADTRANS_CHECK` ← `model/src/packages_check.F:308`
- `RADTRANS_INIT_FIXED` ← `model/src/packages_init_fixed.F:416`
- `RADTRANS_INIT_VARIA` ← `model/src/packages_init_variables.F:343`
- `RADTRANS_READPARMS` ← `model/src/packages_readparms.F:252`
- `RADTRANS_CALC` ← `pkg/darwin/darwin_light_radtrans.F:177`
