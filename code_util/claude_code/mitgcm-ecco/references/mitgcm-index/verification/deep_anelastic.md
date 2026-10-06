# verification/deep_anelastic

## Build variants (code*/)
- **code**: packages.conf = `gfd diagnostics`
  - expanded: mom_common mom_fluxform mom_vecinv generic_advdiff debug mdsio rw monitor diagnostics
  - SIZE.h: grid 1x160x120; sNx=1, sNy=40, OLx=2, OLy=2, nSx=1, nSy=4, nPx=1, nPy=1, Nr=120; has SIZE.h_mpi
  - option/size headers: CPP_OPTIONS.h, DIAGNOSTICS_SIZE.h

## Input variants (input*/)
- **input**: data.pkg on: useDiagnostics
  - data: deltaT=300., nTimeSteps=18, nIter0=0, eosType='LINEAR', implicitFreeSurface=.TRUE., nonHydrostatic=.TRUE., usingSphericalPolarGrid=.TRUE., tempAdvScheme=77
  - namelist files: data data.diagnostics data.pkg eedata
- **input.vecinv**: data.pkg on: -
  - data: deltaT=300., nTimeSteps=18, nIter0=0, eosType='LINEAR', implicitFreeSurface=.TRUE., nonHydrostatic=.TRUE., usingSphericalPolarGrid=.TRUE., tempAdvScheme=77
  - namelist files: data data.diagnostics

## Reference results
`output.txt` `output.vecinv.txt`

Run: `cd verification; ./testreport -of <optfile> -t deep_anelastic` (add `-adm`/`-tlm`/`-tap` for adjoint variants).
