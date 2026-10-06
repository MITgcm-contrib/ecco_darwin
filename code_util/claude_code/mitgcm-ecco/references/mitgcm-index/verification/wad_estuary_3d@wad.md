# verification/wad_estuary_3d@wad

**vs its upstream base (merge-base):** new (not in its upstream base)

## Build variants (code*/)
- **code**: packages.conf = `gfd obcs kpp gmredi ptracers diagnostics exf cal seaice wad`
  - expanded: mom_common mom_fluxform mom_vecinv generic_advdiff debug mdsio rw monitor obcs kpp gmredi ptracers diagnostics exf cal seaice wad
  - SIZE.h: grid 40x30x5; sNx=20, sNy=15, OLx=3, OLy=3, nSx=2, nSy=2, nPx=1, nPy=1, Nr=5
  - option/size headers: CPP_OPTIONS.h, DIAGNOSTICS_SIZE.h, PTRACERS_SIZE.h
  - modified/extra source: EXCH.h, balzano_get_etan.F, obcs_calc.F, obcs_wad_offshore.F

## Input variants (input*/)
- **input**: data.pkg on: useOBCS, useKPP, useGMRedi, useWAD
  - data: nIter0=0, nonlinFreeSurf=4, eosType='LINEAR', usingCartesianGrid=.TRUE., momStepping=.TRUE.
  - namelist files: data data.gmredi data.kpp data.obcs data.pkg data.wad eedata
- **input.carry**: data.pkg on: -
  - data: nIter0=0, nonlinFreeSurf=4, eosType='LINEAR', usingCartesianGrid=.TRUE., momStepping=.TRUE.
  - namelist files: data data.wad
- **input.ice**: data.pkg on: useOBCS, useEXF, useCAL, useSEAICE, useWAD
  - data: nIter0=0, nonlinFreeSurf=4, eosType='LINEAR', usingCartesianGrid=.TRUE., momStepping=.TRUE.
  - namelist files: data data.cal data.exf data.pkg data.seaice data.wadob
- **input.ptr**: data.pkg on: useOBCS, useKPP, useGMRedi, usePTRACERS, useDiagnostics, useWAD
  - namelist files: data.diagnostics data.pkg data.ptracers

## Reference results
`output.carry.txt` `output.ice.txt` `output.ptr.txt` `output.txt`

Run: `cd verification; ./testreport -of <optfile> -t wad_estuary_3d` (add `-adm`/`-tlm`/`-tap` for adjoint variants).
