# verification/wad_flat_xz@wad

**vs its upstream base (merge-base):** new (not in its upstream base)

## Build variants (code*/)
- **code**: packages.conf = `gfd obcs kpp wad`
  - expanded: mom_common mom_fluxform mom_vecinv generic_advdiff debug mdsio rw monitor obcs kpp wad
  - SIZE.h: grid 100x1x10; sNx=20, sNy=1, OLx=3, OLy=3, nSx=5, nSy=1, nPx=1, nPy=1, Nr=10
  - option/size headers: CPP_OPTIONS.h
  - modified/extra source: balzano_get_etan.F, obcs_calc.F

## Input variants (input*/)
- **input**: data.pkg on: useWAD
  - data: nIter0=0, nonlinFreeSurf=4, eosType='LINEAR', usingCartesianGrid=.TRUE., momStepping=.TRUE.
  - namelist files: data data.pkg data.wad eedata
- **input.tide**: data.pkg on: useOBCS, useKPP, useWAD
  - data: nIter0=0, nonlinFreeSurf=4, eosType='LINEAR', usingCartesianGrid=.TRUE., momStepping=.TRUE.
  - namelist files: data data.kpp data.obcs data.pkg

## Reference results
`output.tide.txt` `output.txt`

Run: `cd verification; ./testreport -of <optfile> -t wad_flat_xz` (add `-adm`/`-tlm`/`-tap` for adjoint variants).
