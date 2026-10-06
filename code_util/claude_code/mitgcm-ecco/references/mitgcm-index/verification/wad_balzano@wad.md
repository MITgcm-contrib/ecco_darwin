# verification/wad_balzano@wad

**vs its upstream base (merge-base):** new (not in its upstream base)

## Build variants (code*/)
- **code**: packages.conf = `gfd obcs wad diagnostics`
  - expanded: mom_common mom_fluxform mom_vecinv generic_advdiff debug mdsio rw monitor obcs wad diagnostics
  - SIZE.h: grid 140x1x1; sNx=28, sNy=1, OLx=3, OLy=3, nSx=5, nSy=1, nPx=1, nPy=1, Nr=1
  - option/size headers: CPP_OPTIONS.h
  - modified/extra source: balzano_get_etan.F, obcs_calc.F

## Input variants (input*/)
- **input**: data.pkg on: useOBCS, useWAD
  - data: deltaT=10., nTimeSteps=12960, nIter0=0, nonlinFreeSurf=4, select_rStar=2, eosType='LINEAR', implicitFreeSurface=.TRUE., usingCartesianGrid=.TRUE., saltAdvScheme=33, momStepping=.TRUE.
  - namelist files: data data.obcs data.pkg data.wad eedata
- **input.pool**: data.pkg on: -
- **input.step**: data.pkg on: -
- **input.zstar**: data.pkg on: -
  - data: deltaT=10., nTimeSteps=12960, nIter0=0, nonlinFreeSurf=4, select_rStar=0, eosType='LINEAR', implicitFreeSurface=.TRUE., usingCartesianGrid=.TRUE., saltAdvScheme=33, momStepping=.TRUE.
  - namelist files: data

## Reference results
`output.pool.txt` `output.step.txt` `output.txt` `output.zstar.txt`

Run: `cd verification; ./testreport -of <optfile> -t wad_balzano` (add `-adm`/`-tlm`/`-tap` for adjoint variants).
