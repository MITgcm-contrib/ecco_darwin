# pkg/bulk_force

Simple bulk-formula surface forcing (older alternative to exf).

**runtime switch:** `useBULK_FORCE`-style flag in `data.pkg` (check exact name in packages_boot.F)
**reads:** `data.blk`

## Namelist parameters
### BULKF_CONST
- `rhoA` — density of air [kg/m^3]
- `rhoFW` — density of fresh water [kg/m^3]
- `cpAir` — specific heat of air [J/kg/K]
- `Lvap` — latent heat of vaporization at 0.oC [J/kg]
- `Lfresh` — latent heat of melting of pure ice [J/kg]
- `Tf0kel` — Freezing temp of fresh water in Kelvin = 273.15
- `Rgas` — gas constant for dry air   [J/kg/K]
- `xkar` — von Karman constant  [-]
- `stefan` — Stefan-Boltzmann constant [W/m^2/K^4]
- `zref` — reference height [m] for transfer coefficient
- `zwd` — height [m] of near-surface wind-speed input data
- `zth` — height [m] of near-surface air-temp. & air-humid. input
- `cDrag_1`
- `cDrag_2`
- `cDrag_3`
- `cStantonS` — coefficients used to evaluate Stanton number (for sensib. Heat Flx), under Stable / Unstable stratification
- `cStantonU`
- `cDalton` — coefficient used to evaluate Dalton number (for Evap)
- `umin` — minimum wind speed used in bulk-formulae [m/s]
- `humid_fac` — dry-air - water-vapor molecular mass ratio (minus one) (used to calculate virtual temp.)
- `saltQsFac` — reduction of sat. vapor pressure over salty water
- `gamma_blk` — adiabatic lapse rate
- `atm_emissivity`
- `ocean_emissivity`
- `snow_emissivity`
- `ice_emissivity`
- `FWIND0` — ratio of near-sfc wind to lowest-level wind  _[ifdef ALLOW_FORMULA_AIM]_
- `CHS` — heat exchange coefficient over sea  _[ifdef ALLOW_FORMULA_AIM]_
- `VGUST` — wind speed for sub-grid-scale gusts  _[ifdef ALLOW_FORMULA_AIM]_
- `DTHETA` — Potential temp. gradient for stability correction  _[ifdef ALLOW_FORMULA_AIM]_
- `dTstab` — potential temp. increment for stability function derivative  _[ifdef ALLOW_FORMULA_AIM]_
- `FSTAB` — Amplitude of stability correction (fraction)  _[ifdef ALLOW_FORMULA_AIM]_
- `ocean_albedo` — ocean surface albedo [0-1]
### BULKF_PARM01
- `useFluxFormula_AIM` — set to T when using AIM flux formula rather than the default formula (LANL)
- `blk_nIter` — Number of iterations to find turbulent transfer coeff.
- `calcWindStress` — True to calculate Wind-Stress from surface wind
- `blk_taveFreq`
- `AirTempFile`
- `AirHumidityFile`
- `RainFile`
- `SolarFile`
- `LongwaveFile`
- `UWindFile`
- `VWindFile`
- `RunoffFile`
- `WSpeedFile`
- `QnetFile`
- `EmPFile`
- `CloudFile`
- `airPotTempFile`
### BULKF_PARM02
- `qnet_off`  _[ifdef CONSERV_BULKF]_
- `empmr_off`  _[ifdef CONSERV_BULKF]_
- `conservcycle`  _[ifdef CONSERV_BULKF]_

## CPP options (defaults as shipped)
- `CONSERV_BULKF` (undef, BULK_FORCE_OPTIONS.h)
- `ALLOW_FORMULA_AIM` (define, BULK_FORCE_OPTIONS.h) — allow to use of AIM surface flux formulation (S/R BULKF_FORMULA_AIM) rather than the default (S/R BULKF_FORMULA_LANL)

## Headers
- `BULKF.h` — variable for forcing using bulk formula FORCING FIELD VARIABLES - Mandatory: tair      :: air temperature (K) qair      :: specific humidity at surfac
- `BULKF_CONSERV.h` — /==========================================================\ Header for Bulk formula conservation variables \=========================================
- `BULKF_INT.h` — Intermediate variables for bulk forcing
- `BULKF_PARAMS.h` — BULK_PARAMS.h Header file for BULK_FORCE package parameters: - basic parameter ( I/O frequency, etc ...) - physical constants
- `BULK_FORCE_OPTIONS.h` — Package-specific Options & Macros go here

## Routines (11)
`bulkf_ave.F`, `bulkf_fields_load.F`, `bulkf_flux_adjust.F`, `bulkf_forcing.F`, `bulkf_formula_aim.F`, `bulkf_formula_lanl.F`, `bulkf_formula_lay.F`, `bulkf_init_varia.F`, `bulkf_output.F`, `bulkf_readparms.F`, `bulkf_sh2rh_aim.F`

## Called from outside the package
- `BULKF_FORCING` ← `model/src/forward_step.F:551`
- `BULKF_FIELDS_LOAD` ← `model/src/load_fields_driver.F:202`
- `BULKF_INIT_VARIA` ← `model/src/packages_init_variables.F:301`
- `BULKF_READPARMS` ← `model/src/packages_readparms.F:226`
- `BULKF_FORMULA_AIM` ← `pkg/thsice/thsice_get_bulkf.F:96`
- `BULKF_FORMULA_LANL` ← `pkg/thsice/thsice_get_bulkf.F:120`
- `BULKF_FORMULA_LAY` ← `pkg/thsice/thsice_get_bulkf.F:128`
- `BULKF_FLUX_ADJUST` ← `pkg/thsice/thsice_step_fwd.F:410`

## Verification experiments compiling it (1)
`global_ocean.cs32x15`
