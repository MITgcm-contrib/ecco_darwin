# verification/wad_mangrove@wad

**vs its upstream base (merge-base):** new (not in its upstream base)

## Build variants (code*/)
- **code**: packages.conf = `gfd obcs kpp ptracers diagnostics wad sediment mangrove`
  - expanded: mom_common mom_fluxform mom_vecinv generic_advdiff debug mdsio rw monitor obcs kpp ptracers diagnostics wad sediment mangrove
  - SIZE.h: grid 100x50x5; sNx=25, sNy=25, OLx=3, OLy=3, nSx=4, nSy=2, nPx=1, nPy=1, Nr=5
  - option/size headers: CPP_OPTIONS.h, DIAGNOSTICS_SIZE.h, PTRACERS_SIZE.h
  - modified/extra source: EXCH.h, MGTIDE.h, balzano_get_etan.F, mgtide_read.F, obcs_calc.F

## Input variants (input*/)
- **input**: data.pkg on: useOBCS, useKPP, useWAD, useMANGROVE
  - data: nIter0=0, nonlinFreeSurf=4, eosType='LINEAR', usingCartesianGrid=.TRUE., momStepping=.TRUE.
  - namelist files: data data.kpp data.mangrove data.mgtide data.obcs data.pkg data.wad eedata
- **input.avic**: data.pkg on: -
  - namelist files: data.mangrove
- **input.none**: data.pkg on: useOBCS, useKPP, useWAD
  - namelist files: data.pkg
- **input.patchy**: data.pkg on: -
- **input.sedw**: data.pkg on: useOBCS, useKPP, usePTRACERS, useWAD, useSEDIMENT, useMANGROVE
  - namelist files: data.mangrove data.pkg data.ptracers data.sediment
- **input.steady**: data.pkg on: -
  - data: nIter0=0, nonlinFreeSurf=4, eosType='LINEAR', usingCartesianGrid=.TRUE., momStepping=.TRUE.
  - namelist files: data data.mangrove data.mgtide data.obcs

## Reference results
`output.avic.txt` `output.none.txt` `output.patchy.txt` `output.sedw.txt` `output.steady.txt` `output.txt`

Run: `cd verification; ./testreport -of <optfile> -t wad_mangrove` (add `-adm`/`-tlm`/`-tap` for adjoint variants).
