# verification/lab_sea

## README (first 40 lines)
```
Labrador Sea Region with Sea-Ice
=========================================

### Primary test Overview:
This example sets up a small (20x16x23) Labrador Sea experiment
coupled to a dynamic thermodynamic sea-ice model (MITgcm Documentation 8.6.2).

The domain of integration spans $`[280, 320]^\circ`$E and $`[46, 78]^\circ`$N.
Horizontal grid spacing is 2 degrees.
The 23 vertical levels and the bathymetry file

```
  bathyFile      = 'bathy.labsea1979'
```
are obtained from the the 2$`^\circ`$ ECCO configuration.

Integration is initialized from annual-mean Levitus climatology

```
 hydrogThetaFile = 'LevCli_temp.labsea1979'
 hydrogSaltFile  = 'LevCli_salt.labsea1979'
```

Surface salinity relaxation is to the monthly mean Levitus climatology

```
 saltClimFile    = 'SSS.labsea1979'
```

Forcing files are a 1979-1999 monthly climatology computed from the
NCEP reanalysis (see [`SEAICE_PARAMS.h`](https://github.com/MITgcm/MITgcm/blob/master/pkg/seaice/SEAICE_PARAMS.h) for units and signs)

```
  uwindFile      = 'u10m.labsea1979'  # 10-m zonal wind
  vwindFile      = 'v10m.labsea1979'  # 10-m meridional wind
  atempFile      = 'tair.labsea1979'  # 2-m air temperature
  aqhFile        = 'qa.labsea1979'    # 2-m specific humidity
  lwdownFile     = 'flo.labsea1979'   # downward longwave radiation
  swdownFile     = 'fsh.labsea1979'   # downward shortwave radiation
  precipFile     = 'prate.labsea1979' # precipitation
... (truncated)
```

## Build variants (code*/)
- **code**: packages.conf = `oceanic cd_code exf seaice salt_plume ptracers longstep diagnostics mnc`
  - expanded: mom_common mom_fluxform mom_vecinv generic_advdiff debug mdsio rw monitor gmredi kpp cd_code exf seaice salt_plume ptracers longstep diagnostics mnc
  - SIZE.h: grid 20x16x23; sNx=10, sNy=8, OLx=4, OLy=4, nSx=2, nSy=2, nPx=1, nPy=1, Nr=23; has SIZE.h_mpi
  - option/size headers: CPP_OPTIONS.h, DIAGNOSTICS_SIZE.h, GMREDI_OPTIONS.h, SEAICE_OPTIONS.h
- **code_ad**: packages.conf = `oceanic cd_code down_slope exf seaice salt_plume diagnostics mnc ecco adjoint`
  - expanded: mom_common mom_fluxform mom_vecinv generic_advdiff debug mdsio rw monitor gmredi kpp cd_code down_slope exf seaice salt_plume diagnostics mnc ecco autodiff cost ctrl grdchk
  - SIZE.h: grid 20x16x23; sNx=10, sNy=8, OLx=4, OLy=4, nSx=2, nSy=2, nPx=1, nPy=1, Nr=23; has SIZE.h_mpi
  - option/size headers: AUTODIFF_OPTIONS.h, CPP_OPTIONS.h, CTRL_OPTIONS.h, CTRL_SIZE.h, DIAGNOSTICS_SIZE.h, EXF_OPTIONS.h, GMREDI_OPTIONS.h, SEAICE_OPTIONS.h
  - modified/extra source: tamc.h
- **code_tap**: packages.conf = `oceanic cd_code down_slope exf seaice salt_plume diagnostics mnc ecco adjoint tapenade`
  - expanded: mom_common mom_fluxform mom_vecinv generic_advdiff debug mdsio rw monitor gmredi kpp cd_code down_slope exf seaice salt_plume diagnostics mnc ecco autodiff cost ctrl grdchk tapenade
  - SIZE.h: grid 20x16x23; sNx=10, sNy=8, OLx=4, OLy=4, nSx=2, nSy=2, nPx=1, nPy=1, Nr=23; has SIZE.h_mpi
  - option/size headers: AUTODIFF_OPTIONS.h, CPP_OPTIONS.h, CTRL_OPTIONS.h, CTRL_SIZE.h, DIAGNOSTICS_SIZE.h, EXF_OPTIONS.h, GMREDI_OPTIONS.h, SEAICE_OPTIONS.h

## Input variants (input*/)
- **input**: data.pkg on: useGMRedi, useKPP, useEXF, useCAL, useSEAICE, useDiagnostics, useMNC
  - data: deltaTmom=3600.0, deltaTtracer=3600.0, endTime=36000., startTime=3600.0, eosType='JMD95Z', implicitFreeSurface=.TRUE., usingSphericalPolarGrid=.TRUE., momStepping=.TRUE.
  - namelist files: data data.cal data.diagnostics data.exf data.exf_YearlyFields data.gmredi data.kpp data.mnc data.pkg data.seaice data_YearlyFields eedata
- **input.fd**: data.pkg on: useGMRedi, useKPP, useEXF, useCAL, useSEAICE, useDiagnostics
  - namelist files: data.pkg data.seaice eedata
- **input.hb87**: data.pkg on: useGMRedi, useKPP, useEXF, useCAL, useSEAICE, useDiagnostics
  - data: deltaTmom=3600.0, deltaTtracer=3600.0, endTime=36000., startTime=0.0, eosType='JMD95Z', implicitFreeSurface=.TRUE., usingSphericalPolarGrid=.TRUE., momStepping=.TRUE.
  - namelist files: data data.pkg data.seaice eedata
- **input.longstep**: data.pkg on: useGMRedi, useKPP, useDiagnostics, usePTRACERS
  - data: deltaTmom=3600.0, endTime=93600., startTime=21600., eosType='POLY3', implicitFreeSurface=.TRUE., usingSphericalPolarGrid=.TRUE.
  - namelist files: data data.diagnostics data.gmredi data.kpp data.longstep data.pkg data.ptracers prepare_run
- **input.natl_box**: data.pkg on: useKPP, useDiagnostics
  - data: deltaTmom=3600.0, endTime=93600., startTime=21600., eosType='POLY3', implicitFreeSurface=.TRUE., usingSphericalPolarGrid=.TRUE.
  - namelist files: data data.diagnostics data.kpp data.pkg eedata
- **input.salt_plume**: data.pkg on: useGMRedi, useKPP, useEXF, useCAL, useSEAICE, useSALT_PLUME
  - data: deltaTmom=3600.0, deltaTtracer=3600.0, endTime=36000., startTime=3600.0, eosType='JMD95Z', implicitFreeSurface=.TRUE., usingSphericalPolarGrid=.TRUE., tempAdvScheme=7, saltAdvScheme=7, momStepping=.TRUE.
  - namelist files: data data.pkg data.salt_plume data.seaice eedata
- **input_ad**: data.pkg on: useGMRedi, useKPP, useEXF, useSEAICE, useDOWN_SLOPE, useECCO, useGrdchk
  - data: deltaTmom=3600.0, deltaTtracer=3600.0, nTimeSteps=4, startTime=0.0, eosType='JMD95Z', implicitFreeSurface=.TRUE., usingSphericalPolarGrid=.TRUE., tempAdvScheme=30, saltAdvScheme=30, momStepping=.TRUE.
  - namelist files: data data.autodiff data.cal data.cost data.ctrl data.down_slope data.ecco data.err data.exf data.gmredi data.grdchk data.kpp data.optim data.pkg data.seaice eedata prepare_run
- **input_ad.noseaice**: data.pkg on: useGMRedi, useKPP, useEXF, useDiagnostics, useMNC, useECCO, useGrdchk
  - data: deltaTmom=3600.0, deltaTtracer=3600.0, nTimeSteps=12, startTime=0.0, eosType='JMD95Z', implicitFreeSurface=.TRUE., usingSphericalPolarGrid=.TRUE., tempAdvScheme=30, saltAdvScheme=30, momStepping=.TRUE.
  - namelist files: data data.diagnostics data.ecco data.mnc data.pkg
- **input_ad.noseaicedyn**: data.pkg on: useGMRedi, useKPP, useEXF, useSEAICE, useDOWN_SLOPE, useSALT_PLUME, useECCO, useGrdchk
  - data: deltaTmom=3600.0, deltaTtracer=3600.0, nTimeSteps=12, startTime=0.0, eosType='JMD95Z', implicitFreeSurface=.TRUE., usingSphericalPolarGrid=.TRUE., tempAdvScheme=30, saltAdvScheme=30, momStepping=.TRUE.
  - namelist files: data data.grdchk data.pkg data.salt_plume data.seaice
- **input_tap**: data.pkg on: useGMRedi, useKPP, useEXF, useSEAICE, useDOWN_SLOPE, useCAL, useMNC, useECCO, useGrdchk
  - namelist files: data.mnc data.pkg prepare_run
- **input_tap.noecco**: data.pkg on: useGMRedi, useKPP, useEXF, useSEAICE, useDOWN_SLOPE, useCAL, useGrdchk
  - namelist files: data.cost data.ctrl data.pkg data.seaice

## Reference results
`output.fd.txt` `output.hb87.txt` `output.longstep.txt` `output.natl_box.txt` `output.salt_plume.txt` `output.txt` `output_adm.noseaice.txt` `output_adm.noseaicedyn.txt` `output_adm.txt` `output_tap_adj.noecco.txt` `output_tap_adj.txt` `output_tlm.txt.gz`

Run: `cd verification; ./testreport -of <optfile> -t lab_sea` (add `-adm`/`-tlm`/`-tap` for adjoint variants).
