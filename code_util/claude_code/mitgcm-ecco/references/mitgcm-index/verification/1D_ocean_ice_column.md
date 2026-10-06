# verification/1D_ocean_ice_column

## Build variants (code*/)
- **code**: packages.conf = `gfd kpp exf seaice`
  - expanded: mom_common mom_fluxform mom_vecinv generic_advdiff debug mdsio rw monitor kpp exf seaice
  - SIZE.h: grid 1x1x23; sNx=1, sNy=1, OLx=2, OLy=2, nSx=1, nSy=1, nPx=1, nPy=1, Nr=23
  - option/size headers: CPP_OPTIONS.h, DIAGNOSTICS_SIZE.h, SEAICE_OPTIONS.h
- **code_ad**: packages.conf = `gfd kpp exf seaice ecco autodiff cost ctrl grdchk`
  - expanded: mom_common mom_fluxform mom_vecinv generic_advdiff debug mdsio rw monitor kpp exf seaice ecco autodiff cost ctrl grdchk
  - SIZE.h: grid 1x1x23; sNx=1, sNy=1, OLx=2, OLy=2, nSx=1, nSy=1, nPx=1, nPy=1, Nr=23
  - option/size headers: CPP_OPTIONS.h, CTRL_OPTIONS.h, CTRL_SIZE.h, DIAGNOSTICS_SIZE.h, EXF_OPTIONS.h, SEAICE_OPTIONS.h
  - modified/extra source: tamc.h

## Input variants (input*/)
- **input**: data.pkg on: useKPP, useEXF, useCAL, useSEAICE
  - data: deltaTtracer=3600.0, nTimeSteps=10, startTime=0.0, eosType='JMD95Z', implicitFreeSurface=.TRUE., usingCartesianGrid=.TRUE., tempAdvScheme=30, saltAdvScheme=30, momStepping=.TRUE., useRealFreshWaterFlux=.TRUE.
  - namelist files: data data.cal data.err data.exf data.kpp data.pkg data.seaice eedata
- **input_ad**: data.pkg on: useKPP, useEXF, useSEAICE, useECCO, useGrdchk
  - data: deltaTtracer=3600.0, nTimeSteps=10, startTime=0.0, eosType='JMD95Z', implicitFreeSurface=.TRUE., usingCartesianGrid=.TRUE., tempAdvScheme=30, saltAdvScheme=30, momStepping=.TRUE., useRealFreshWaterFlux=.TRUE.
  - namelist files: data data.autodiff data.cal data.cost data.ctrl data.ecco data.exf data.grdchk data.kpp data.optim data.pkg data.seaice eedata prepare_run

## Reference results
`output.txt` `output_adm.txt` `output_tlm.txt.gz`

Run: `cd verification; ./testreport -of <optfile> -t 1D_ocean_ice_column` (add `-adm`/`-tlm`/`-tap` for adjoint variants).
