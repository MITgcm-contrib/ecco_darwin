# verification/wetdry_overflood@wetdry

**vs its upstream base (merge-base):** new (not in its upstream base)

## Build variants (code*/)
- **code**: packages.conf = `gfd diagnostics exf cal seaice wetdry overflood`
  - expanded: mom_common mom_fluxform mom_vecinv generic_advdiff debug mdsio rw monitor diagnostics exf cal seaice wetdry overflood
  - SIZE.h: grid 80x40x5; sNx=20, sNy=20, OLx=3, OLy=3, nSx=4, nSy=2, nPx=1, nPy=1, Nr=5
  - option/size headers: CPP_OPTIONS.h, DIAGNOSTICS_SIZE.h
  - modified/extra source: EXCH.h

## Input variants (input*/)
- **input**: data.pkg on: useWETDRY, useOVERFLOOD, useDiagnostics
  - data: nIter0=0, nonlinFreeSurf=4, eosType='LINEAR', usingCartesianGrid=.TRUE., momStepping=.TRUE.
  - namelist files: data data.diagnostics data.overflood data.pkg data.wetdry eedata
- **input.icefld**: data.pkg on: useWETDRY, useOVERFLOOD, useSEAICE, useEXF, useCAL, useDiagnostics
  - data: nIter0=0, nonlinFreeSurf=4, eosType='LINEAR', usingCartesianGrid=.TRUE., momStepping=.TRUE.
  - namelist files: data data.cal data.diagnostics data.exf data.overflood data.pkg data.seaice data.wetdry
- **input.melt**: data.pkg on: useWETDRY, useOVERFLOOD, useSEAICE, useEXF, useCAL, useDiagnostics
  - data: nIter0=0, nonlinFreeSurf=4, eosType='LINEAR', usingCartesianGrid=.TRUE., momStepping=.TRUE.
  - namelist files: data data.cal data.diagnostics data.exf data.overflood data.pkg data.seaice

## Reference results
`output.icefld.txt` `output.melt.txt` `output.txt`

Run: `cd verification; ./testreport -of <optfile> -t wetdry_overflood` (add `-adm`/`-tlm`/`-tap` for adjoint variants).
