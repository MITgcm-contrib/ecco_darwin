# verification/wad_thacker_1d@wad

**vs its upstream base (merge-base):** new (not in its upstream base)

## Build variants (code*/)
- **code**: packages.conf = `gfd wad diagnostics`
  - expanded: mom_common mom_fluxform mom_vecinv generic_advdiff debug mdsio rw monitor wad diagnostics
  - SIZE.h: grid 250x1x1; sNx=50, sNy=1, OLx=3, OLy=3, nSx=5, nSy=1, nPx=1, nPy=1, Nr=1
  - option/size headers: CPP_OPTIONS.h

## Input variants (input*/)
- **input**: data.pkg on: useWAD
  - data: deltaT=4.947464852, nTimeSteps=1360, nIter0=0, nonlinFreeSurf=4, select_rStar=2, eosType='LINEAR', implicitFreeSurface=.TRUE., usingCartesianGrid=.TRUE., saltAdvScheme=33, momStepping=.TRUE.
  - namelist files: data data.pkg data.wad eedata
- **input.carry**: data.pkg on: -
  - namelist files: data.wad
- **input.zstar**: data.pkg on: -
  - data: deltaT=4.947464852, nTimeSteps=1360, nIter0=0, nonlinFreeSurf=4, select_rStar=0, eosType='LINEAR', implicitFreeSurface=.TRUE., usingCartesianGrid=.TRUE., saltAdvScheme=33, momStepping=.TRUE.
  - namelist files: data

## Reference results
`output.carry.txt` `output.txt` `output.zstar.txt`

Run: `cd verification; ./testreport -of <optfile> -t wad_thacker_1d` (add `-adm`/`-tlm`/`-tap` for adjoint variants).
