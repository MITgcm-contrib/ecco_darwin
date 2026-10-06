# verification/wad_iceground@wad

**vs its upstream base (merge-base):** new (not in its upstream base)

## Build variants (code*/)
- **code**: packages.conf = `gfd diagnostics exf cal seaice wad`
  - expanded: mom_common mom_fluxform mom_vecinv generic_advdiff debug mdsio rw monitor diagnostics exf cal seaice wad
  - SIZE.h: grid 60x10x5; sNx=30, sNy=10, OLx=3, OLy=3, nSx=2, nSy=1, nPx=1, nPy=1, Nr=5
  - option/size headers: CPP_OPTIONS.h, DIAGNOSTICS_SIZE.h, SEAICE_OPTIONS.h
  - modified/extra source: EXCH.h

## Input variants (input*/)
- **input**: data.pkg on: useWAD, useSEAICE, useEXF, useCAL, useDiagnostics
  - data: nIter0=0, nonlinFreeSurf=4, eosType='LINEAR', usingCartesianGrid=.TRUE., momStepping=.TRUE., useRealFreshWaterFlux=.TRUE.
  - namelist files: data data.cal data.diagnostics data.exf data.pkg data.seaice data.wad eedata
- **input.dyn**: data.pkg on: -
  - data: nIter0=0, nonlinFreeSurf=4, eosType='LINEAR', usingCartesianGrid=.TRUE., momStepping=.TRUE., useRealFreshWaterFlux=.TRUE.
  - namelist files: data data.diagnostics data.exf data.seaice data.wad

## Reference results
`output.dyn.txt` `output.txt`

Run: `cd verification; ./testreport -of <optfile> -t wad_iceground` (add `-adm`/`-tlm`/`-tap` for adjoint variants).
