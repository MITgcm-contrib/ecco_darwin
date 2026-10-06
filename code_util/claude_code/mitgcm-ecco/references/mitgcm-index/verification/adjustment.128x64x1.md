# verification/adjustment.128x64x1

## Build variants (code*/)
- **code**: packages.conf = `(inherits/none)`
  - SIZE.h: grid 128x64x1; sNx=64, sNy=32, OLx=2, OLy=2, nSx=2, nSy=2, nPx=1, nPy=1, Nr=1; has SIZE.h_mpi

## Input variants (input*/)
- **input**: data.pkg on: -
  - data: deltaT=450.0, nTimeSteps=24, nIter0=0, eosType='IDEALG', implicitFreeSurface=.TRUE., usingSphericalPolarGrid=.TRUE., buoyancyRelation='ATMOSPHERIC'
  - namelist files: data data.pkg eedata

## Reference results
`output.txt`

Run: `cd verification; ./testreport -of <optfile> -t adjustment.128x64x1` (add `-adm`/`-tlm`/`-tap` for adjoint variants).
