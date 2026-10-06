# verification/wad_mudflat@wad

**vs its upstream base (merge-base):** new (not in its upstream base)

## Build variants (code*/)
- **code**: packages.conf = `gfd obcs kpp ptracers rbcs diagnostics exf cal seaice wad sediment`
  - expanded: mom_common mom_fluxform mom_vecinv generic_advdiff debug mdsio rw monitor obcs kpp ptracers rbcs diagnostics exf cal seaice wad sediment
  - SIZE.h: grid 100x60x5; sNx=25, sNy=30, OLx=3, OLy=3, nSx=4, nSy=2, nPx=1, nPy=1, Nr=5
  - option/size headers: CPP_OPTIONS.h, DIAGNOSTICS_SIZE.h, PTRACERS_SIZE.h
  - modified/extra source: EXCH.h, balzano_get_etan.F, obcs_calc.F, obcs_wad_offshore.F

## Input variants (input*/)
- **input**: data.pkg on: useOBCS, useKPP, usePTRACERS, useRBCS, useWAD
  - data: nIter0=0, nonlinFreeSurf=4, eosType='LINEAR', usingCartesianGrid=.TRUE., momStepping=.TRUE.
  - namelist files: data data.kpp data.obcs data.pkg data.ptracers data.rbcs data.wad eedata
- **input.ice**: data.pkg on: useOBCS, useKPP, usePTRACERS, useRBCS, useEXF, useCAL, useSEAICE, useWAD
  - data: nIter0=0, nonlinFreeSurf=4, eosType='LINEAR', usingCartesianGrid=.TRUE., momStepping=.TRUE.
  - namelist files: data data.cal data.exf data.pkg data.seaice data.wadob
- **input.rain**: data.pkg on: -
  - data: nIter0=0, nonlinFreeSurf=4, eosType='LINEAR', usingCartesianGrid=.TRUE., momStepping=.TRUE., useRealFreshWaterFlux=.TRUE.
  - namelist files: data
- **input.rfwf0**: data.pkg on: -
  - data: nIter0=0, nonlinFreeSurf=4, eosType='LINEAR', usingCartesianGrid=.TRUE., momStepping=.TRUE., useRealFreshWaterFlux=.TRUE.
  - namelist files: data
- **input.sed**: data.pkg on: useOBCS, useKPP, usePTRACERS, useRBCS, useWAD, useSEDIMENT
  - namelist files: data.pkg data.ptracers data.sediment
- **input.sedw**: data.pkg on: useOBCS, useKPP, usePTRACERS, useRBCS, useWAD, useSEDIMENT
  - namelist files: data.pkg data.ptracers data.sediment

## Reference results
`output.ice.txt` `output.rain.txt` `output.rfwf0.txt` `output.sed.txt` `output.sedw.txt` `output.txt`

Run: `cd verification; ./testreport -of <optfile> -t wad_mudflat` (add `-adm`/`-tlm`/`-tap` for adjoint variants).
