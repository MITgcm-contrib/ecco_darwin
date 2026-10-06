# verification/short_surf_wave

## Build variants (code*/)
- **code**: packages.conf = `gfd diagnostics`
  - expanded: mom_common mom_fluxform mom_vecinv generic_advdiff debug mdsio rw monitor diagnostics
  - SIZE.h: grid 52x1x50; sNx=13, sNy=1, OLx=2, OLy=2, nSx=4, nSy=1, nPx=1, nPy=1, Nr=50; has SIZE.h_mpi
  - option/size headers: CPP_OPTIONS.h, DIAGNOSTICS_SIZE.h

## Input variants (input*/)
- **input**: data.pkg on: useDiagnostics
  - data: deltaT=5.e-3, nTimeSteps=11, nIter0=1, eosType='LINEAR', implicitFreeSurface=.TRUE., nonHydrostatic=.TRUE., usingCartesianGrid=.TRUE.
  - namelist files: data data.1it data.diagnostics data.pkg eedata

## Reference results
`output.txt`

Run: `cd verification; ./testreport -of <optfile> -t short_surf_wave` (add `-adm`/`-tlm`/`-tap` for adjoint variants).
