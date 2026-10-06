# verification/offline_exf_seaice

## README (first 40 lines)
```
Seaice-only verification experiment in idealized periodic channel
-----------------------------------------------------------------

1) main forward experiment (code, input)

  Re-entrant zonally periodic channel (80x42 grid points) with just level (Nr=1)
   uniform resolution (5.km, 10m), solid Southern boundary with triangular shape
   coastline ("bathy_3c.bin")

  Use seaice (dynamics & thermodynamics from pkg/thsice) with EXF (see data.pkg)
   with initial ice thickness of 0.2 m (but no snow)
   (thSIceThick_InitFile='const+20.bin', in "input/data.ice")
  Initial seaice concentration is 100 % everywhere
   (thSIceFract_InitFile='const100.bin', in "input/data.ice")
  and seaice is initially at rest.

  At runtime turn off time-stepping in 'data', PARM01, using:
    momStepping  = .FALSE.,
    saltStepping = .FALSE.,
    tempAdvection=.FALSE.,
  And just keep surface temp relaxation (tauRelax = 1 month) toward fixed SST:
   in data.exf :
  > climsstperiod      = 0.0,
  > climsstTauRelax    = 2592000.,
  >  climsstfile       = 'tocn.bin',

 Forcing:
  None of the forcing vary with time; Most of the input files have been
   generated using matlab script "input/gendata.m".
  SST relaxation field is uniform in X, parabolic function of Y with
   maximum close to Southern boundary.

  Atmospheric air temp is uniform in Y, and only vary with X (~sin(2.pi.x/Lx))
   with an amplitude of 4.K ('tair_4x.bin');
  Uses constant Relative Humidity (70%, file 'qa70_4x.bin')
  constant and uniform downward shortwave (100.W/m2, 'dsw_100.bin'),
                       downward longwave (250.W/m^2, 'dlw_250.bin'),
                       zonal wind (10.m/s, 'windx.bin'),
  no meridional wind, no precip.

... (truncated)
```

## Build variants (code*/)
- **code**: packages.conf = `gfd exf -cal seaice thsice diagnostics`
  - expanded: mom_common mom_fluxform mom_vecinv generic_advdiff debug mdsio rw monitor exf seaice thsice diagnostics
  - SIZE.h: grid 80x42x1; sNx=40, sNy=21, OLx=3, OLy=3, nSx=2, nSy=2, nPx=1, nPy=1, Nr=1; has SIZE.h_mpi
  - option/size headers: DIAGNOSTICS_SIZE.h, SEAICE_OPTIONS.h, SEAICE_SIZE.h
- **code_ad**: packages.conf = `gfd exf -cal obcs seaice thsice diagnostics adjoint`
  - expanded: mom_common mom_fluxform mom_vecinv generic_advdiff debug mdsio rw monitor exf obcs seaice thsice diagnostics autodiff cost ctrl grdchk
  - SIZE.h: grid 80x42x1; sNx=40, sNy=21, OLx=3, OLy=3, nSx=2, nSy=2, nPx=1, nPy=1, Nr=1; has SIZE.h_mpi
  - option/size headers: COST_OPTIONS.h, CTRL_SIZE.h, DIAGNOSTICS_SIZE.h, EXF_OPTIONS.h, OBCS_OPTIONS.h, SEAICE_OPTIONS.h, SEAICE_SIZE.h, THSICE_OPTIONS.h
  - modified/extra source: MDSIO_BUFF_WH.h, tamc.h

## Input variants (input*/)
- **input**: data.pkg on: useEXF, useSEAICE, useThSIce, useDiagnostics
  - data: deltaT=900.0, nTimeSteps=24, startTime=0.0, eosType='LINEAR', implicitFreeSurface=.TRUE., usingCartesianGrid=.TRUE., momStepping=.FALSE.
  - namelist files: data data.diagnostics data.exf data.ice data.pkg data.seaice eedata
- **input.dyn_ellnnfr**: data.pkg on: useEXF, useSEAICE, useDiagnostics
  - data: deltaT=1800.0, nTimeSteps=12, startTime=0.0, eosType='LINEAR', implicitFreeSurface=.TRUE., usingCartesianGrid=.TRUE., momStepping=.FALSE.
  - namelist files: data data.diagnostics data.exf data.pkg data.seaice
- **input.dyn_jfnk**: data.pkg on: useEXF, useSEAICE, useThSIce, useDiagnostics
  - data: deltaT=1800.0, nTimeSteps=12, startTime=0.0, eosType='LINEAR', implicitFreeSurface=.TRUE., usingCartesianGrid=.TRUE., momStepping=.FALSE.
  - namelist files: data data.diagnostics data.exf data.ice data.pkg data.seaice
- **input.dyn_lsr**: data.pkg on: useEXF, useSEAICE, useDiagnostics
  - data: deltaT=1800.0, nTimeSteps=12, startTime=0.0, eosType='LINEAR', implicitFreeSurface=.TRUE., usingCartesianGrid=.TRUE., momStepping=.FALSE.
  - namelist files: data data.diagnostics data.exf data.pkg data.seaice
- **input.dyn_mce**: data.pkg on: useEXF, useSEAICE, useDiagnostics
  - data: deltaT=1800.0, nTimeSteps=12, startTime=0.0, eosType='LINEAR', implicitFreeSurface=.TRUE., usingCartesianGrid=.TRUE., momStepping=.FALSE.
  - namelist files: data data.diagnostics data.exf data.pkg data.seaice
- **input.dyn_paralens**: data.pkg on: -
  - data: deltaT=1800.0, nTimeSteps=12, startTime=0.0, eosType='LINEAR', implicitFreeSurface=.TRUE., usingCartesianGrid=.TRUE., momStepping=.FALSE.
  - namelist files: data data.diagnostics data.exf data.ice data.seaice
- **input.dyn_teardrop**: data.pkg on: -
  - data: deltaT=1800.0, nTimeSteps=12, startTime=0.0, eosType='LINEAR', implicitFreeSurface=.TRUE., usingCartesianGrid=.TRUE., momStepping=.FALSE.
  - namelist files: data data.diagnostics data.exf data.ice data.seaice
- **input.thermo**: data.pkg on: useEXF, useSEAICE, useDiagnostics
  - data: deltaT=3600.0, nTimeSteps=120, startTime=0.0, eosType='LINEAR', implicitFreeSurface=.TRUE., usingCartesianGrid=.TRUE., momStepping=.FALSE.
  - namelist files: data data.diagnostics data.exf data.pkg data.seaice
- **input.thsice**: data.pkg on: useEXF, useThSIce, useDiagnostics
  - data: deltaT=3600.0, nTimeSteps=120, startTime=0.0, eosType='LINEAR', implicitFreeSurface=.TRUE., usingCartesianGrid=.TRUE., momStepping=.FALSE.
  - namelist files: data data.diagnostics data.exf data.ice data.pkg
- **input_ad**: data.pkg on: useEXF, useSEAICE, useDiagnostics, useGrdchk
  - data: deltaT=3600.0, nTimeSteps=120, startTime=0.0, eosType='LINEAR', implicitFreeSurface=.TRUE., usingCartesianGrid=.TRUE., momStepping=.FALSE.
  - namelist files: data data.autodiff data.cost data.ctrl data.diagnostics data.exf data.grdchk data.optim data.pkg data.seaice eedata prepare_run
- **input_ad.obcs**: data.pkg on: useEXF, useSEAICE, useOBCS, useDiagnostics, useGrdchk
  - data: deltaT=3600.0, nTimeSteps=12, startTime=0.0, eosType='LINEAR', implicitFreeSurface=.TRUE., usingCartesianGrid=.TRUE., momStepping=.FALSE.
  - namelist files: data data.autodiff data.ctrl data.diagnostics data.exf data.grdchk data.obcs data.pkg data.seaice
- **input_ad.thsice**: data.pkg on: useEXF, useThSIce, useGrdchk
  - data: deltaT=3600.0, nTimeSteps=60, startTime=0.0, eosType='LINEAR', implicitFreeSurface=.TRUE., usingCartesianGrid=.TRUE., momStepping=.FALSE.
  - namelist files: data data.grdchk data.ice data.pkg

## Reference results
`output.dyn_ellnnfr.txt` `output.dyn_jfnk.txt` `output.dyn_lsr.txt` `output.dyn_mce.txt` `output.dyn_paralens.txt` `output.dyn_teardrop.txt` `output.thermo.txt` `output.thsice.txt` `output.txt` `output_adm.obcs.txt` `output_adm.thsice.txt` `output_adm.txt` `output_tlm.obcs.txt.gz` `output_tlm.thsice.txt.gz` `output_tlm.txt.gz`

Run: `cd verification; ./testreport -of <optfile> -t offline_exf_seaice` (add `-adm`/`-tlm`/`-tap` for adjoint variants).
