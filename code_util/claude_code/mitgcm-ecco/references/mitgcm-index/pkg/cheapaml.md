# pkg/cheapaml

Cheap atmospheric mixed layer: simple prognostic atmospheric boundary layer above the ocean for surface fluxes.

**runtime switch:** `useCHEAPAML`-style flag in `data.pkg` (check exact name in packages_boot.F)
**reads:** `data.cheapaml`
**manual:** `doc/examples/examples.rst`

## Namelist parameters
### CHEAPAML_CONST
- `cheapaml_ntim`
- `cheapaml_mask_width`
- `cheapaml_h`
- `cheapaml_kdiff`
- `cheap_tauRelax` — main relaxation time-scale (in sec) for atm T & Q used with cheapMask (if provided) or over land
- `cheap_tauRelaxOce` — relaxation time-scale (in sec) for atm T & Q, used over ocean if cheapMask is not provided
- `cdrag_1`
- `cdrag_2`
- `cdrag_3`
- `rhoa`
- `cpair`
- `stefan`
- `gasR` — gas constant
- `xkar` — von Karman constant
- `dsolms` — Solar variation at Southern boundary
- `dsolmn` — Solar variation at Northern boundary
- `zu`
- `zt`
- `zq`
- `xphaseinit` — user input initial phase of year relative to mid winter. e.g. xphaseinit = pi implies time zero is mid summer.
- `gamma_blk` — atmospheric adiabatic lapse rate
- `humid_fac` — humidity factor for computing virtual potential temperature
- `p0` — surface pressure in mb
- `cheap_pr1` — precipitation time constant
- `cheap_pr2` — precipitation time constant
- `cheapaml_taurelax`
- `cheapaml_taurelaxocean`
### CHEAPAML_PARM01
- `periodicExternalForcing_cheap`
- `externForcingPeriod_cheap`
- `externForcingCycle_cheap`
- `AirTempFile`
- `SolarFile`
- `UWindFile`
- `VWindFile`
- `TrFile`
- `QrFile`
- `AirQFile`
- `UStressFile`
- `VStressFile`
- `WaveHFile`
- `WavePFile`
- `TracerFile`
- `TracerRfile`
- `cheapMaskFile`
- `cheap_hFile`
- `cheap_clFile`
- `cheap_dlwFile`
- `cheap_prFile`
### CHEAPAML_PARM02
- `cheapamlXperiodic` — domain (including land) is periodic in X dir
- `cheapamlYperiodic` — domain (including land) is periodic in Y dir
- `useFreshWaterFlux` — option to include evap+precip  (on  by default)
- `useFluxLimit` — use flux limiting advection    (off by default)
- `FluxFormula`
- `WaveModel`
- `useStressOption` — use stress option              (off by default)
- `useCheapTracer` — use passive tracer option      (off by default)
- `useTimeVarBLH` — use time varying BL height option (off by default)
- `useClouds` — use clouds option              (off by default)
- `useDLongWave` — use imported downward longwave  (off by default)
- `usePrecip` — use imported precipitation (off by default)
- `useRelativeWind` — use relative wind (off by default)

## CPP options (defaults as shipped)
- `INCONSISTENT_WIND_LOCATION` (undef, CHEAPAML_OPTIONS.h) — to reproduce old results, with inconsistent wind location, grid-cell center and grid-cell edges (C-grid).

## Headers
- `CHEAPAML.h` — #ifdef ALLOW_CHEAPAML
- `CHEAPAML_OPTIONS.h` — BOP

## Routines (16)
`cheapaml.F`, `cheapaml_calc_rhs.F`, `cheapaml_coare3_flux.F`, `cheapaml_copy_edges.F`, `cheapaml_diagnostics_init.F`, `cheapaml_fields_load.F`, `cheapaml_init_fixed.F`, `cheapaml_init_varia.F`, `cheapaml_lanl_flux.F`, `cheapaml_read_pickup.F`, `cheapaml_readparms.F`, `cheapaml_seaice.F`, `cheapaml_timestep.F`, `cheapaml_write_pickup.F`

## Called from outside the package
- `CHEAPAML` ← `model/src/forward_step.F:564`
- `CHEAPAML_FIELDS_LOAD` ← `model/src/load_fields_driver.F:225`
- `CHEAPAML_INIT_FIXED` ← `model/src/packages_init_fixed.F:262`
- `CHEAPAML_INIT_VARIA` ← `model/src/packages_init_variables.F:316`
- `CHEAPAML_READPARMS` ← `model/src/packages_readparms.F:236`
- `CHEAPAML_WRITE_PICKUP` ← `model/src/packages_write_pickup.F:140`

## Verification experiments compiling it (1)
`cheapAML_box`
